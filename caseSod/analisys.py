

import numpy
import matplotlib.pyplot as plt
from readTables import readTables
import sys
import characteristics as ch
 
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
        
        return None        
                
    
if __name__=="__main__":

    if len(sys.argv) < 1:
        print("Insert the input file")
        exit(0)
    
    path = sys.argv[1]
        
    s = solution(path+"mesh.csv", path+"./solution.csv")        
        
    s.calcPMT()
        
    """
        
    plt.figure()
    plt.title("pressure")    
    plt.contourf(s.x, s.y, s.p)
    plt.axis("equal")
    plt.colorbar()
    plt.show()    
       
    """
        
    plt.figure()
    plt.title("mach")
    plt.contourf(s.x, s.y, s.mach)
    plt.axis("equal")
    plt.colorbar()    
    plt.show()    
    
        
    n = int(s.x.shape[1]/2)
    x = s.x[:, n]   
        
    char1 = ch.problem(xchange=25.0, p1=1e4, T1=240, p4=1e5, T4=300)
    charVar = char1.calcVar(0.02, x)

    v = 'rho'
    plt.figure()
    plt.plot(charVar['x'], charVar[v],'-b')
    plt.plot(x, s.r[:,n],'.r')
    plt.legend(['charac.', 'CFD'])
    plt.title(v)
    plt.grid(True)
    plt.show()
    
    v = 'u'
    plt.figure()
    plt.plot(charVar['x'], charVar[v],'-b')
    plt.plot(x, s.u[:,n],'.r')
    plt.legend(['charac.', 'CFD'])
    plt.title(v)
    plt.grid(True)
    plt.show()    
    
    v = 'p'
    plt.figure()
    plt.plot(charVar['x'], charVar[v],'-b')
    plt.plot(x, s.p[:,n],'.r')
    plt.legend(['charac.', 'CFD'])
    plt.title(v)
    plt.grid(True)
    plt.show()    
    
    v = 'mach'
    mach = charVar['u']/numpy.sqrt(1.4*charVar['p']/charVar['rho'])
    plt.figure()
    plt.plot(charVar['x'], mach,'-b')
    plt.plot(x, s.mach[:,n],'.r')
    plt.legend(['charac.', 'CFD'])
    plt.title(v)
    plt.grid(True)
    plt.show()    
