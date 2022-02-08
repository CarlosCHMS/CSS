


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
            
                ff.write('%.10e,' % self.x[ii, jj])
                
            ff.write('\n')
            
        ff.write('1, %i, %i, \n' % self.y.shape)
        for ii in range(0, self.y.shape[0]):
            for jj in range(0, self.y.shape[1]):        
            
                ff.write('%.10e,' % self.y[ii, jj])
                
            ff.write('\n')
            
        ff.close()    


            
class nozzle(basic):

    def __init__(self, x0, x1, y0, n0, n1, n2, a):
    
                
        self.x = numpy.zeros((n0+n1, n2))
        self.y = numpy.zeros((n0+n1, n2))
        
        t0 = numpy.pi*a/180  
        
        b = 0.2
        
        dt = t0/n0
        r0 = 0.5
        for ii in range(0, n0):
            t = ii*dt
            aux = y0 + r0*(1 - numpy.cos(t))
            x = r0*numpy.sin(t)
            dy = aux/(n2-1)
            for jj in range(0, n2):                    
                self.x[ii, jj] = x*(1 + b*(aux**2 - (dy*jj)**2))
                self.y[ii, jj] = dy*jj
        
        x0 = r0*numpy.sin(t0)
        y01 = y0 + r0*(1 - numpy.cos(t0))
            
        dx = x1/(n1 - 1)
        for ii in range(0, n1):    
            aux = y01 + numpy.tan(t0)*dx*ii
            dy = aux/(n2-1)
            for jj in range(0, n2):                    
                self.x[ii + n0, jj] = (x0 + dx*ii)*(1 + b*(aux**2 - (dy*jj)**2))
                self.y[ii + n0, jj] = dy*jj
 
                 
if __name__=='__main__':

    w = nozzle(1, 1, 0.25, 10, 80, 20, 20)
    w.plot()        
    w.write()
    
