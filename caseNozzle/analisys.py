

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
        self.u = u        
        self.v = v

        self.entro = self.p/(self.r**self.gamma)

        self.H = E + self.p/self.r
        
        return None        

def levels(v, n):    

    max1 = v[0][0]
    min1 = v[0][0]
    for ii in range(0, v.shape[0]):
        for jj in range(0, v.shape[1]):
            max1 = max(v[ii][jj], max1)
            min1 = min(v[ii][jj], min1)
                            
    d = (max1-min1)/(n-1)
    levels = []
    for ii in range(0, n):
        levels.append(min1 + d*ii)
    
    return levels                
    
if __name__=="__main__":

    if len(sys.argv) < 1:
        print("Insert the input file")
        exit(0)
    
    path = sys.argv[1]
        
    s = solution(path+"mesh.csv", path+"./solution.csv")        
        
    s.calcPMT()
            
        
    plt.figure()
    plt.title("pressure")    
    plt.contourf(s.x, s.y, s.p)
    plt.axis("equal")
    plt.colorbar()
    plt.show()        
        
    plt.figure()
    plt.title("mach")
    plt.contourf(s.x, s.y, s.mach)
    plt.axis("equal")
    plt.colorbar()    
    plt.show()    

    plt.figure()
    plt.title("r")
    plt.contourf(s.x, s.y, s.r)
    plt.axis("equal")
    plt.colorbar()    
    plt.show()    
            
    plt.figure()
    plt.title("ru")
    plt.contourf(s.x, s.y, s.ru)
    plt.axis("equal")
    plt.colorbar()    
    plt.show()            
    
    plt.figure()
    plt.title("rv")
    plt.contourf(s.x, s.y, s.rv, levels=levels(s.rv, 10))
    plt.axis("equal")
    plt.colorbar()    
    plt.show()    

    plt.figure()
    plt.title("rE")
    plt.contourf(s.x, s.y, s.rE)
    plt.axis("equal")
    plt.colorbar()    
    plt.show() 

    plt.figure()
    plt.title("entropia")
    plt.contourf(s.x, s.y, s.entro)
    plt.axis("equal")
    plt.colorbar()    
    plt.show()

    plt.figure()
    plt.title("entalpia")
    plt.contourf(s.x, s.y, s.H)
    plt.axis("equal")
    plt.colorbar()    
    plt.show()

