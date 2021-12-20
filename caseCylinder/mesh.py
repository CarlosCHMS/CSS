


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


class wedge(basic):

    def __init__(self, x0, x1, y0, n0, n1, n2, alpha):
    
        a = numpy.tan(alpha*numpy.pi/180)
    
        self.x0 = x0
        self.x1 = x1
        self.y0 = y0
                
        self.x = numpy.zeros((n0+n1, n2))
        self.y = numpy.zeros((n0+n1, n2))
        
        dx = self.x0/n0
        dy = y0/(n2-1)
        for ii in range(0, n0):
            for jj in range(0, n2):                    
                self.x[ii, jj] = dx*ii
                self.y[ii, jj] = dy*jj
            
        dx = self.x1/(n1 - 1)
        for ii in range(0, n1):
            xaux = dx*ii
            yaux = a*xaux
            dy = (y0 - yaux)/(n2-1)
            for jj in range(0, n2):                    
                self.x[ii + n0, jj] = x0 + xaux
                self.y[ii + n0, jj] = yaux + dy*jj
                


class tunel(basic):

    def __init__(self, x0, x1, y0, n0, n1, n2, a):
    
        self.a = a
    
        self.x0 = x0
        self.x1 = x1
        self.y0 = y0
                
        self.x = numpy.zeros((2*n0+n1, n2))
        self.y = numpy.zeros((2*n0+n1, n2))
        
        dx = self.x0/n0
        dy = y0/(n2-1)
        for ii in range(0, n0):
            for jj in range(0, n2):                    
                self.x[ii, jj] = dx*ii
                self.y[ii, jj] = dy*jj
            
        dx = self.x1/n1
        for ii in range(0, n1):
            xaux = dx*ii
            yaux = self.func(xaux)
            dy = (y0 - yaux)/(n2-1)
            for jj in range(0, n2):                    
                self.x[ii + n0, jj] = x0 + xaux
                self.y[ii + n0, jj] = yaux + dy*jj

        dx = self.x0/(n0-1)
        dy = y0/(n2-1)
        for ii in range(0, n0):
            for jj in range(0, n2):                    
                self.x[ii + n0 + n1, jj] = dx*ii + x0 + x1
                self.y[ii + n0 + n1, jj] = dy*jj
                
    def func(self, x):
    
        b = self.x1/2
        x = x - b
        y = numpy.sqrt(self.a**2 + b**2 - x**2) - self.a 
                
        return y        

class wedge2(basic):

    def __init__(self, x0, x1, x2, y0, n0, n1, n2, alpha):
    
        a = numpy.tan(alpha*numpy.pi/180)
    
        self.x0 = x0
        self.x1 = x1
        self.y0 = y0
                
        self.x = numpy.zeros((n0+n1+n0, n2))
        self.y = numpy.zeros((n0+n1+n0, n2))
        
        dx = self.x0/n0
        dy = y0/(n2-1)
        for ii in range(0, n0):
            for jj in range(0, n2):                    
                self.x[ii, jj] = dx*ii
                self.y[ii, jj] = dy*jj
            
        dx = self.x1/(n1)
        for ii in range(0, n1):
            xaux = dx*ii
            yaux = a*xaux
            dy = (y0 - yaux)/(n2-1)
            for jj in range(0, n2):                    
                self.x[ii + n0, jj] = x0 + xaux
                self.y[ii + n0, jj] = yaux + dy*jj
                
        dx = x2/(n0-1)
        dy = (y0 - yaux)/(n2-1)        
        for ii in range(0, n0):
            for jj in range(0, n2):                    
                self.x[ii + n0+n1, jj] = dx*ii + x0 + x1
                self.y[ii + n0+n1, jj] = yaux + dy*jj        
                
        
class nose(basic):

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

    #w = wedge2(0.25, 0.25, 0.25, 0.25, 40, 40, 40, 20)
    #w = wedge(0.25, 0.25, 0.5, 40, 40, 80, 20)
    #w = nozzle(1, 1, 0.25, 10, 80, 20, 20)
    w = nose(0.5, 2.5, 80, 100) 
    #w = tunel(1, 1, 1, 40, 40, 40, 2)       
    w.plot()        
    w.write()
    
