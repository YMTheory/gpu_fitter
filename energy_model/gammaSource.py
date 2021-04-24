import numpy as np
from ROOT import TH1D, TF1
import uproot as up
import elecLoader as eloader
import random
import time
from numba import cuda, float32
from numba.cuda.random import create_xoroshiro128p_states, xoroshiro128p_uniform_float32
#import matplotlib.pyplot as plt



def loadPrmBeta(filename):
    maxN = 100
    secBetaArr, secAntiBetaArr = np.zeros((5000, maxN)), np.zeros((5000, maxN))
    nLine = 0
    with open(filename) as f:
        for lines in f.readlines():
            line = lines.strip("\n")
            data = line.split(" ")
            counta = 0
            countb = 0
            for i in data:
                if "a" in i:
                    tmp = list(i)
                    tmp.pop()
                    j = ''.join(tmp)
                    secAntiBetaArr[nLine, counta] = float(j)
                    counta+=1
                if "b" in i:
                    tmp = list(i)
                    tmp.pop()
                    j = ''.join(tmp)
                    secBetaArr[nLine, countb] = float(j)
                    countb+=1
            nLine += 1
    return secBetaArr, secAntiBetaArr

def loadSimTruth(filename):
    evt = up.open(filename)["evt"]
    totpe = evt["totalPE"].array()
    return totpe


def GammaPePred(secBeta, secAntiBeta):
    # two-times sampling alg here
    sc = 1.00
    calcpe = []
    ctime = []
    sample_id = np.random.randint(0, 5000, size=5000)
    for i in sample_id:
        st = time.time()
        tmppe = 0
        for j in secBeta[i]:
            if j == 0:
                break
            tmppe += random.gauss(eloader.getNPE(j), eloader.getSPE(j)*sc)
        for j in secAntiBeta[i]:
            if j == 0:
                break
            tmppe += random.gauss(eloader.getNPE(j), eloader.getSPE(j)*sc) + 2*random.gauss(660.8, 27.07*sc)
        calcpe.append(tmppe)
        et = time.time()
        ctime.append(et-st)
    return calcpe, ctime


def preSampling(secBeta, secAntiBeta):
    st = time.time()
    Nevents = 5000
    mu_arr, sigma_arr = [], []
    sample_id = np.random.randint(0, 5000, size=Nevents)
    for ii in sample_id:
        tmpmu, tmpsigma = 0, 0
        for j in secBeta[ii, :]:
            if j == 0:
                break
            tmpmu += eloader.getNPE(j)
            tmpsigma += eloader.getSPE(j)**2
        for j in secAntiBeta[ii, :]:
            if j == 0:
                break
            tmpmu += ( eloader.getNPE(j) + 660.8*2 )
            tmpsigma += ( eloader.getSPE(j)**2 + 2*27.07**2 )
        tmpsigma = np.sqrt(tmpsigma)
        mu_arr.append(tmpmu)
        sigma_arr.append(tmpsigma)
    et = time.time()
    mu_arr = np.array(mu_arr)
    sigma_arr = np.array(sigma_arr)
    print("Time consumed during pre-sampling: %s s" %(et-st))
    return mu_arr, sigma_arr


@cuda.jit
def preSampling_inGPU(secBeta, secAntiBeta, indexarr, npearr, spearr, muarr, sigmaarr):
    # shared memory usage :
    sNPE = cuda.shared.array(shape=(1, 845), dtype=float32)
    sSPE = cuda.shared.array(shape=(1, 845), dtype=float32)

    x, y = cuda.grid(2)

    # threading positioning
    tx = cuda.threadIdx.x
    ty = cuda.blockIdx.x
    bw = cuda.blockDim.x
    pos = tx + ty * bw

    if pos < muarr.size:
        ii = indexarr[pos]
        tmpmu, tmpsigma = 0, 0
        for j in secBeta[ii, :]:
            if j == 0:
                break
            tmpmu += eloader.getNPE(j)
            tmpsigma += eloader.getSPE(j)**2
        for j in secAntiBeta[ii, :]:
            if j == 0:
                break
            tmpmu += ( eloader.getNPE(j) + 660.8*2 )
            tmpsigma += ( eloader.getSPE(j)**2 + 2*27.07**2 )
        tmpsigma = np.sqrt(tmpsigma)
        muarr[pos] = tmpmu
        sigmaarr[pos] = tmpsigma
    cuda.syncthreads()




@cuda.jit
def GammaPePred_inGPU(mu_arr, sigma_arr, calcpe, rng_states):
    thread_id = cuda.grid(1)
    tx = cuda.threadIdx.x
    ty = cuda.blockIdx.x
    bw = cuda.blockDim.x
    pos = tx + ty * bw
    tmppe = 0
    if pos < calcpe.size:
        tmppe += cuda.random.xoroshiro128p_normal_float32(rng_states, thread_id) * sigma_arr[pos] + mu_arr[pos]
        calcpe[pos] = tmppe

    cuda.syncthreads()



def gammaSource(name, Etrue):
    start = time.time()
    path = "./data/"
    totpe = []
    filename = path + name + "_totpe.root"
    print("------> Reading " + filename)
    tmptotpe = loadSimTruth(filename)
    totpe.extend(tmptotpe)
    
    secBetaArr, secAntiBetaArr = loadPrmBeta("./data/"+name+"_J19.txt")

    maxN = 5000
    threadsperblock = 128
    blockspergrid = 40

    #st = time.time()
    #mu_arr, sigma_arr = np.zeros((1, maxN)), np.zeros((1, maxN))
    #index_arr = np.random.randint(0, 5000, size=maxN)
    #preSampling_inGPU[blockspergrid, threadsperblock](secBetaArr, secAntiBetaArr, index_arr, mu_arr, sigma_arr)
    #print(mu_arr)
    #print(sigma_arr)
    

    #mu_arr, sigma_arr = preSampling(secBetaArr, secAntiBetaArr)

    '''
    calcpe = np.array([0 for i in range(maxN)])
    ##blockspergrid = (tmp.size + (threadsperblock - 1)) // threadsperblock
    rng_states = create_xoroshiro128p_states(threadsperblock * blockspergrid, seed=1)
    GammaPePred_inGPU[blockspergrid, threadsperblock](mu_arr, sigma_arr, calcpe, rng_states)

    #calcpe, ctime = GammaPePred(secBetaArr, secAntiBetaArr)
    et = time.time()
    print("Total time consumed during prediction in GPU: %s s" %(et-st))


    totpe = np.array(totpe)
    calcpe = np.array(calcpe)
    low = totpe.min() - 50
    high = totpe.max() + 50


    f1 = TF1("f1", "gaus", low, high)
    f2 = TF1("f2", "gaus", low, high)
    h1 = TH1D(name+"h1", "", 100, low, high)
    h2 = TH1D(name+"h2", "", 100, low, high)
    f1.SetParameter(1, totpe.mean())
    f2.SetParameter(1, calcpe.mean())
    for i, j in zip(totpe, calcpe):
        h1.Fill(i)
        h2.Fill(j)
    h1.Fit(f1, "RE")
    h2.Fit(f2, "RE")


    sim_mean, sim_sigma   = f1.GetParameter(1), f1.GetParameter(2)
    calc_mean, calc_sigma = f2.GetParameter(1), f2.GetParameter(2)
    sim_mean_err, sim_sigma_err   = f1.GetParError(1), f1.GetParError(2)
    calc_mean_err, calc_sigma_err = f2.GetParError(1), f2.GetParError(2)

    end = time.time()
    print("Total time consumed in whole model: %s s" %(end-start))
    

    return sim_mean, sim_mean_err, sim_sigma, sim_sigma_err, calc_mean, calc_mean_err, calc_sigma, calc_sigma_err
    ''' 
    return 0
