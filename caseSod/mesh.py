


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



class tube(basic):

    def __init__(self, x0, y0, n0, n1):
    
                
        self.x = numpy.zeros((n0, n1))
        self.y = numpy.zeros((n0, n1))
        
        dx = x0/(n0-1)
        dy = y0/(n1-1)
        for ii in range(0, n0):
            for jj in range(0, n1):
                self.x[ii, jj] = dx*ii
                self.y[ii, jj] = dy*jj
                             
if __name__=='__main__':

    w = tube(50, 5, 100, 10)       
    w.plot()        
    w.write()
    
