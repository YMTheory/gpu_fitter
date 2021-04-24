import numpy as np
import gammaSource as gam

name = ["Cs137", "Mn54", "Ge68", "K40", "nH", "Co60", "AmBe", "nC12", "AmC"]
Etrue = [0.662, 0.835, 1.022, 1.461, 2.223, 2.506, 4.43, 4.94, 6.13]


def gamCol():
    simMeanArr      = [0 for i in range(9)]
    simMeanErrArr   = [0 for i in range(9)]
    simSigmaArr     = [0 for i in range(9)]
    simSigmaErrArr  = [0 for i in range(9)]
    calcMeanArr     = [0 for i in range(9)]
    calcMeanErrArr  = [0 for i in range(9)]
    calcSigmaArr    = [0 for i in range(9)]
    calcSigmaErrArr = [0 for i in range(9)]

    for i in range(9):
        simMeanArr[i], simMeanErrArr[i], simSigmaArr[i], simSigmaErrArr[i], \
        calcMeanArr[i], calcMeanErrArr[i], calcSigmaArr[i], calcSigmaErrArr[i] = gam.gammaSource(name[i], Etrue[i])

    scale = simMeanArr[5] / Etrue[5]
    simnonl, simnonl_err, simres, simres_err = [], [], [], []
    calcnonl, calcnonl_err, calcres, calcres_err = [], [], [], []
    for i in range(9):
        if name[i] == "Ge68" or name[i]=="Co60":
            simnonl.append(simMeanArr[i]/2/scale/Etrue[i])
            simnonl_err.append(simMeanErrArr[i]/2/scale/Etrue[i])
            calcnonl.append(calcMeanArr[i]/2/scale/Etrue[i])
            calcnonl_err.append(calcMeanErrArr[i]/2/scale/Etrue[i])
        else:
            simnonl.append(simMeanArr[i]/scale/Etrue[i])
            simnonl_err.append(simMeanErrArr[i]/scale/Etrue[i])
            calcnonl.append(calcMeanArr[i]/scale/Etrue[i])
            calcnonl_err.append(calcMeanErrArr[i]/scale/Etrue[i])

        simres.append(simSigmaArr[i]/simMeanArr[i])
        simres_err.append(np.sqrt(simSigmaErrArr[i]**2/simMeanArr[i]**2 + simMeanErrArr[i]**2*simSigmaArr[i]**2/simMeanArr[i]**4))
        calcres.append(calcSigmaArr[i]/calcMeanArr[i])
        calcres_err.append(np.sqrt(calcSigmaErrArr[i]**2/calcMeanArr[i]**2 + calcMeanErrArr[i]**2*calcSigmaArr[i]**2/calcMeanArr[i]**4))

    with open("test.txt", "w") as f:
        for i in range(9):
            f.write("%.5f %.6f %.5f %.6f %.5f %.6f %.5f %.6f" %(simnonl[i], simnonl_err[i], calcnonl[i], calcnonl_err[i], simres[i], simres_err[i], calcres[i], calcres_err[i]))
            f.write("\n")
    
    

