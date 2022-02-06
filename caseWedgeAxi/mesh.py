


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
                
        self.y += 0.001
                       
                 
if __name__=='__main__':

    w = wedge(0.1, 0.4, 0.5, 16, 64, 80, 10)
    w.plot()        
    w.write()
    
