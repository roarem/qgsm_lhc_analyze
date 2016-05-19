import sys
import numpy as np
import matplotlib.pyplot as plt

class QGSM_Plot:

  def __init__(self):
    self.nrParts  = np.loadtxt("../data/B_MULT",usecols=(2,),dtype=int)
    self.NPOMS    = [] 
    self.NPOMH    = []
    self.event    = []
    self.C        = []
    self.p        = []
    self.pT       = []
    self.rapidity = []
    self.eta      = []
    
    with open("../data/collected.out",'r') as infile:
      infile = infile.readlines()
      self.total = len(infile)
      for i,line in enumerate(infile):
        line = line.strip().split(' ')
        self.event   .append(int  (line[0]))
        self.C       .append(int  (line[1])) 
        self.p       .append(float(line[2]))
        self.pT      .append(float(line[3])) 
        self.rapidity.append(float(line[4])) 
        self.eta     .append(float(line[5]))
        self.NPOMS   .append(int  (line[6]))
        self.NPOMH   .append(int  (line[7]))

        if i%(self.total//100)==0:
          self.msg("%i"%(i/self.total*100+1)+'%')
    print("Done reading, moving on...")

  def msg(self,text):
    text = "  "+text+chr(13)
    sys.stdout.write(text)
    sys.stdout.flush()


  def histogram(self,data,bins=100,range=None,weights=None):
    n, bin_edges = np.histogram(data,bins=bins,range=range,weights=weights)
    bin_centers = 0.5*(bin_edges[1:]+bin_edges[:-1])
    return n, bin_centers
  
  def dN_dnFnB(self):
    NPOM = np.zeros((len(self.nrParts),4))
    old_event = self.event[0]
    i = 0
    for event,eta in zip(enumerate(self.event),self.eta):
      if old_event != event[1]: 
        old_event = event[1]
        i += 1
        NPOM[i,2] = self.NPOMS[event[0]]
        NPOM[i,3] = self.NPOMH[event[0]]
      if eta>0:
        NPOM[i,0] +=1
      else:
        NPOM[i,1] +=1

    npoms = []
    npomh = []
    print("Finding indicies...")
    for i in range(2):
      npoms.append(np.where(NPOM[:,2]==i)[0])
      npomh.append(np.where(NPOM[:,3]==i)[0])
    print("Finding corresponding NPOMS NPOMH...")
    for i,indicies in enumerate(npoms):
      fig, ax = plt.subplots()
      for j,indexes in enumerate(npomh):
        print("NPOMS={}, NPOMH={}".format(i,j))
        npomshf = []
        npomshb = []
        print(len(indicies))        
        for count,k in enumerate(indicies):
          if count%10==0:
            self.msg(str(count))
          for l in indexes:
            if l == k:
              npomshf.append(NPOM[l,0])
              npomshb.append(NPOM[l,1])
        n,x = self.histogram(npomshf,bins=20,weights=npomsb)
        nf,xf = self.histogram(npomshf,bins=20)
        
        ax.plot(x,n/nf,label="NPOMH={}".format(j))
      plt.legend()
      plt.savefig("../graphs/npomsh/s{}.pdf".format(i))
    '''
    fig, ax = plt.subplots(ncols=2,nrows=1)
    n,x = self.histogram(NPOM[:,0],bins=20,weights=NPOM[:,1])
    nf,xf = self.histogram(NPOM[:,0],bins=20)
    label = 'All'
    ax[0].plot(x,n/nf,label=label,marker='*',linewidth=3)
    ax[1].plot(x,n/nf,label=label,marker='*',linewidth=3)
      
    for i in range(6):
      npomsindicies = np.where(NPOM[:,2]==i)[0]
      npomsf = []
      npomsb = []
      for j in npomsindicies:
        npomsf.append(NPOM[j,0])
        npomsb.append(NPOM[j,1])
      label = 'NPOMS={}'.format(i)    
      n,xi = self.histogram(npomsf,bins=20,weights=npomsb)
      nf,xif = self.histogram(npomsf,bins=20)
      ax[0].plot(x,n/nf,label=label,marker='o',linewidth=0.5)
    ax[0].legend()
    for i in range(3):
      npomhindicies = np.where(NPOM[:,3]==i)[0]
      npomhf = []
      npomhb = []
      for j in npomhindicies:
        npomhf.append(NPOM[j,0])
        npomhb.append(NPOM[j,1])
      label = 'NPOMH={}'.format(i)    
      n,xi = self.histogram(npomhf,bins=20,weights=npomhb)
      nf,xif = self.histogram(npomhf,bins=20)
      ax[1].plot(x,n/nf,label=label,marker='o',linewidth=0.5)
    ax[1].legend()         
    '''

  def dN_detaNPOMS(self):
    eta = [[[] for j in range(0,np.amax(self.NPOMH)+1)] for i in range(0,np.amax(self.NPOMS)+1)]
    print('Finding NPOMS and NPOMH matches...')
    for i in range(len(self.event)):
      if i%(self.total//100) == 0:
        self.msg('%i'%(i/self.total*100+1)+'%')
      eta[self.NPOMS[i]][self.NPOMH[i]].append(self.eta[i])
    print('Creating plots...')
    for npoms,i in enumerate(eta):
      fig, ax = plt.subplots(figsize=(20,13))
      for npomh,j in enumerate(i):
        label = 'NPOMH={}'.format(npomh)
        title = 'NPOMS = {}'.format(npoms)
        self.msg('NPOMS={},NPOMH={}'.format(npoms,npomh))
        n,x = self.histogram(np.abs(j),bins=200)
        ax.plot(x,n,linewidth=0.3,label=label)

      ax.set_yscale('log')
      ax.set_title(title)
      plt.legend()
      plt.savefig('../graphs/eta_NPOMS{}.pdf'.format(npoms))

  def dN_dnc(self,bins=500):
    dN_dnc,x = self.histogram(self.nrParts,bins=bins)
    fig,ax = plt.subplots()
    ax.plot(x,dN_dnc)
    ax.set_title('$\\frac{dN}{dn_c}$')

  def dN_dnpom(self):
    softpoms,x = self.histogram(self.NPOMS)
    hardpoms,x = self.histogram(self.NPOMH)
    fig, ax = plt.subplots(2)
    ax[0].plot(x,softpoms)
    ax[1].plot(x,hardpoms)
    ax[0].set_title('$\\frac{dN}{dNPOMS}$')
    ax[1].set_title('$\\frac{dN}{dNPOMH}$')

  def dN_deta(self,bins=500):
    n, x = self.histogram(self.eta,bins=bins)
    fig, ax = plt.subplots()
    ax.plot(x,n)
    ax.set_yscale('log')
    ax.set_title('$\\frac{dN}{d\eta}$')

  def dN_drap(self,bins=500):
    dN_drap,x = self.histogram(self.rapidity,bins=bins)
    fig, ax = plt.subplots()
    ax.plot(x,dN_drap)
    ax.set_yscale('log')
    ax.set_title('$\\frac{dN}{dy}$')

     
if __name__=="__main__":
  test = QGSM_Plot()
  #test.dN_deta()
  #test.dN_drap()
  test.dN_dnc()
  #test.dN_dnpom() 
  #test.dN_detaNPOMS()
  #test.dN_dnFnB()
  plt.show()
