


import numpy
import matplotlib.pyplot as plt


class basic():
     
    def plot(self):
    
        plt.figure()
    
        for ii in range(0, self.x.shape[0]):
            plt.plot(self.x[ii, :], self.y[ii, :], 'b')

        for ii in range(0, self.x.shape[1]):
            plt.plot(self.x[:, ii], self.y[:, ii], 'r')


        plt.axis('equal')        
        plt.show()
        
        return None
        
    def write(self):
    
        ff = open('mesh.csv', 'w')

        ff.write("2, \n")
        
        ff.write('0, %i, %i, \n' % self.x.shape)
        for ii in range(0, self.x.shape[0]):
            for jj in range(0, self.x.shape[1]):        
            
                ff.write('%.4f,' % self.x[ii, jj])
                
            ff.write('\n')
            
        ff.write('1, %i, %i, \n' % self.y.shape)
        for ii in range(0, self.y.shape[0]):
            for jj in range(0, self.y.shape[1]):        
            
                ff.write('%.4f,' % self.y[ii, jj])
                
            ff.write('\n')
            
        ff.close()    
        
class cylinder(basic):

    def __init__(self, r0, r1, n0, n1):
                
        self.x = numpy.zeros((n0, n1))
        self.y = numpy.zeros((n0, n1))
        
        t0 = 0
        
        dt = (numpy.pi/2 - 0)/(n0-1)
        q = (r1/r0)**(1/(n1-1))
        for ii in range(0, n0):
            t = ii*dt + t0            
            for jj in range(0, n1):
                r = r0*(q**jj)                    
                self.x[ii, jj] = -r*numpy.cos(t)
                self.y[ii, jj] = r*numpy.sin(t)
                            
            
                 
if __name__=='__main__':

    w = cylinder(0.5, 2.5, 80, 100) 
    w.plot()        
    w.write()
    
