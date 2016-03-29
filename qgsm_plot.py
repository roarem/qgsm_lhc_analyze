import numpy as np
import matplotlib.pyplot as plt

class QGSM_Plot:

  def __init__(self):
    #self.nrParts  = np.loadtxt("data/B_MULT",usecols=(2,),dtype=int)
    self.NPOMS, self.NPOMH = np.loadtxt("data/NPOM.dat",unpack=True,dtype=int)
    self.event    = []
    self.C        = []
    self.pT       = []
    self.one_pT   = []
    self.rapidity = []
    
    #with open("data/collected.out",'r') as infile:
    #  for i,line in enumerate(infile.readlines()):
    #    line = line.strip().split(' , ')
    #    self.C       .append(int  (line[0])) 
    #    self.pT      .append(float(line[1])) 
    #    self.one_pT  .append(float(line[2])) 
    #    self.rapidity.append(float(line[3])) 

    #    if i%100000==0:
    #      print(i)

  def dN_dnc(self):
    dN_dnc = np.bincount(self.nrParts)
    return dN_dnc

  def dN_dnpom(self):
    softpoms = np.bincount(self.NPOMS)
    hardpoms = np.bincount(self.NPOMH)
    return softpoms, hardpoms

if __name__=="__main__":
  test = QGSM_Plot()
  #dN_dnc = test.dN_dnc()
  softpoms, hardpoms = test.dN_dnpom() 
  fig,ax = plt.subplots()
  #ax.plot(dN_dnc)
  #ax.set_yscale('log')
  ax.plot(softpoms)
  ax.plot(hardpoms)


  plt.show()
