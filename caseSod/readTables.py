 
 
 
import numpy
 
 
class readTables():

    def __init__(self, fileName):

        self.tabList = []

        ff = open(fileName, "r")
        
        ii = 0
        tabInit = 1
        newTab = True
        first = True
        tab = []
        for row in ff:
            aux = row.split(',')
            
            if first:
                first = False
            else:            
                if newTab:
                            
                    Nrow = int(aux[1])
                    Ncol = int(aux[2])
                    
                    tabInit += Nrow
                    tab = numpy.zeros((Nrow, Ncol))
                    kk = 0
                    newTab = False
                    
                else:            
                    for jj in range(0, Ncol):
                        tab[kk, jj] = float(aux[jj])

                    kk += 1
                                    
                    if kk == Nrow:
                        self.tabList.append(tab)                        
                        newTab = True
                    
            ii += 1

        return None
                
if __name__=="__main__":

    rt = readTables("solution.csv")
    print(rt.tabList[0])
    print()
    print(rt.tabList[1])
    print()
    print(rt.tabList[2])
    print()
    print(rt.tabList[3])
    
                
            
                
