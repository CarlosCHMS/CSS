

import numpy
import matplotlib.pyplot as plt
from readTables import readTables
import sys
 
class solution():

    def __init__(self, meshFile, solFile):

        self.gamma = 1.4
        self.Rgas = 287.5


        rt = readTables(meshFile)
        
        x0 = rt.tabList[0]
        y0 = rt.tabList[1]
        
        self.x = numpy.zeros((x0.shape[0]-1, x0.shape[1]-1))
        self.y = numpy.zeros((y0.shape[0]-1, y0.shape[1]-1))
        
        for ii in range(0, x0.shape[0]-1):
            for jj in range(0, x0.shape[1]-1):
            
                self.x[ii, jj] = (x0[ii, jj] + x0[ii+1, jj] + x0[ii, jj+1] + x0[ii+1, jj+1])/4
                self.y[ii, jj] = (y0[ii, jj] + y0[ii+1, jj] + y0[ii, jj+1] + y0[ii+1, jj+1])/4
                
        rt = readTables(solFile)
        
        self.r = rt.tabList[0]
        self.ru = rt.tabList[1]
        self.rv = rt.tabList[2]
        self.rE = rt.tabList[3]

    def save(self, M, fileName):

        ff = open(fileName, "w")

        for ii in range(0, M.shape[0]):
            for jj in range(0, M.shape[1]):
                ff.write(" %.4f," % M[ii, jj])

            ff.write(" \n")

        ff.close()

        return None
            
    def calcPMT(self):
       
        u = self.ru/self.r
        v = self.rv/self.r
        E = self.rE/self.r
        
        RT = (E - (u**2 + v**2)/2)*(self.gamma - 1)
        self.T = RT/self.Rgas
        
        self.p = RT*self.r
        
        c = numpy.sqrt(self.gamma*RT)
        V = numpy.sqrt(u**2 + v**2)
        
        self.mach = V/c
        
        return None        
                
    
if __name__=="__main__":

    if len(sys.argv) < 1:
        print("Insert the input file")
        exit(0)
    
    path = sys.argv[1]
        
    s = solution(path+"mesh.csv", path+"solution.csv")        
        
    s.calcPMT()
    
    s.save(s.x, "./octave/x.csv")
    s.save(s.y, "./octave/y.csv")
    s.save(s.rv, "./octave/rv.csv")
    s.save(s.ru, "./octave/ru.csv")
    s.save(s.r, "./octave/r.csv")
    s.save(s.rE, "./octave/rE.csv")

