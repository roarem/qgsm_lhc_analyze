import sys
import numpy as np

class QGSM_Distributions:

  def __init__(self,energy):

    print("Reading files...")
    path = "../data/rawData/"+str(energy)+"/data/"
    eventNr, nrParticles  = np.loadtxt(path+"B_MULT",dtype=int,usecols=(0,2),unpack=True)
    NPOM                  = np.loadtxt(path+"NPOM.dat",dtype=int)
    finalpr               = open(path+"finalpr.data",'r')
    
    self.path = path
    self.eventNr, self.nrParticles, self.finalpr, self.NPOM =\
        eventNr, nrParticles, finalpr, NPOM

  def closeFile(self):
    self.finalpr.close()

  def collectData(self):

    eventNr, nrParticles, finalpr, NPOM =\
        self.eventNr, self.nrParticles, self.finalpr, self.NPOM

    self.lineCount  = 0
    total           = len(eventNr)
    tol             = 10e-9
    inElasticEvents = 0
    elasticEvents   = 0
    startParticle   = 0
    endParticle     = 0

    ALL = []
    eventNr = eventNr[::]
    nrParticles = nrParticles[::]
    print("Collecting data...")
    for event,nrParts in zip(eventNr,nrParticles):

      if nrParts == 2:
        elasticEvents += 1
        finalpr.readline()
        finalpr.readline()
      else:
        inElasticEvents += 1

        for count in range(nrParts):
          parton = finalpr.readline()
          #[startParticle : startParticle+nrParts]:

          parton = list(map(float,parton.strip().split()))

          E   = parton[4]
          px  = parton[5]
          py  = parton[6]
          pz  = parton[7]
          C   = parton[13]

          if C != 0:
            p = np.sqrt(px**2 + py**2 + pz**2) 
            pTransverse = np.sqrt(px**2 + py**2)
            if pTransverse > 0.3 and pTransverse < 1.5:
              if abs(E-pz)<tol:
                rapidity = 20
              elif abs(E+pz)<tol:
                rapidity = 0
              else:
                rapidity = 0.5*np.log((E+pz)/(E-pz))

              if abs(p-pz)<tol:
                eta = 20
              elif abs(p+pz)<tol:
                eta = 0
              else:
                eta = 0.5*np.log((p+pz)/(p-pz))

              if np.abs(eta) < 1:
                NPOMS = NPOM[event-1,0]
                NPOMH = NPOM[event-1,1]


                self.lineCount += 1
                ALL.append([event,C,p,pTransverse,rapidity,eta,NPOMS,NPOMH])
            
      #startParticle += nrParts
      self.countedEvents = inElasticEvents

      if event%(total//100)==0:
        self.msg("%i"%(event/total*100)+'%')
     
    self.ALL = ALL
    print("")
   
  def msg(self,text):
    text = text+chr(13)
    sys.stdout.write(text)
    sys.stdout.flush()

  def writeAnalysis(self):
    print("Writing data...")
    fmts = ['%i','%i']
    output =[[self.lineCount,self.countedEvents]]
    np.savetxt('../data/meta.out',output,fmt=fmts,delimiter=' ')
    fmts = ['%i','%i','%.7e','%.7e','%.7e','%.7e','%i','%i']
    np.savetxt('../data/collected.out',self.ALL,fmt=fmts,delimiter=' ') 
    

if __name__=="__main__":
  energy = 7000
  test = QGSM_Distributions(energy)
  test.collectData()
  test.writeAnalysis()


