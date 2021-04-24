import numpy as np

def loadPEFile(filename):
    Earr, PEarr = [], []
    with open(filename) as f:
        for lines in f.readlines():
            line = lines.strip("\n")
            data = line.split(" ")
            Earr.append(float(data[0]))
            PEarr.append(float(data[1]))
    return Earr, PEarr

def loadResFile(filename):
    Earr, resol = [], []
    with open(filename) as f:
        for lines in f.readlines():
            line = lines.strip("\n")
            data = line.split(" ")
            Earr.append(float(data[0]))
            resol.append(float(data[3]))
    return Earr, resol

sctE, sctPE = loadPEFile("./data/sctPE.txt")
print(">>>>>>>>>> Load scintillation PE file <<<<<<<<<<")
cerE, cerPE = loadPEFile("./data/cerPE.txt") 
print(">>>>>>>>>> Load cerenkov PE file <<<<<<<<<<")
resolE, resol = loadResFile("./data/elecResol.txt")
print(">>>>>>>>>> Load resolution PE file <<<<<<<<<<")

from ROOT import TGraph
g1 = TGraph()
g2 = TGraph()
for i in range(len(sctE)):
    g1.SetPoint(i, sctE[i], sctPE[i]+cerPE[i])
for i in range(len(resolE)):
    g2.SetPoint(i, resolE[i], resol[i])


def getNPE(Etrue):
    return g1.Eval(Etrue, 0, "S")

def getSPE(Etrue):
    return g2.Eval(Etrue, 0, "S")


def getArray():
    sctE = np.array(sctE)
    sctPE = np.array(sctPE)
    cerPE = np.array(cerPE)
    resol = np.array(resol)
    return sctE, sctPE, cerPE, resol
