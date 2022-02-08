

import matplotlib.pyplot as plt
import numpy

class read():

    def __init__(self, fileName):
    
        ff = open(fileName, "r")
        
        ii = 0
        read = False
        
        self.x = []
        self.y = []
        
        for row in ff:
        
            if ii == 5:
                read = True
            
            if row[0] == ')':
                read = False
                
            if read:
                aux = row.split('\t')
                self.x.append(float(aux[0]))
                self.y.append(float(aux[1]))                
                            
            ii += 1
            
        ff.close()
        self.x = numpy.array(self.x)
        self.y = numpy.array(self.y)        
            
if __name__=="__main__":

    m = read("machxX")
    
    plt.figure()
    plt.plot(m.x, m.y)
    plt.plot(m.x, m.x*0 +  1.21021838)
    plt.show()
    
    p = read("pressurexX")
    
    plt.figure()
    plt.plot(p.x, p.y)
    plt.plot(p.x, p.x*0 +   2.84286270*1e5)
    plt.show()    
