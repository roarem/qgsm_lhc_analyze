import numpy as np

class QGSM_Distributions:

  def __init__(self):

    print("Reading files...")
    eventNr, nrParticles  = np.loadtxt("data/B_MULT",dtype=int,usecols=(0,2),unpack=True)
    NPOMS, NPOMH          = np.loadtxt("data/NPOM.dat",unpack=True)
    finalpr               = open("data/finalpr.data",'r').readlines()
    print("Done reading files...")

    self.eventNr, self.nrParticles, self.NPOMS, self.NPOMH, self.finalpr =\
        eventNr, nrParticles, NPOMS, NPOMH, finalpr

  def collectData(self):

    eventNr, nrParticles, NPOMS, NPOMH, finalpr =\
        self.eventNr, self.nrParticles, self.NPOMS, self.NPOMH, self.finalpr

    inElasticEvents = 0
    elasticEvents   = 0
    startParticle   = 0
    endParticle     = 0

    ALL = []
    for event,nrParts in zip(eventNr,nrParticles):

      if nrParts == 2:
        inElasticEvents += 1
      else:
        elasticEvents += 1

        for parton in finalpr[startParticle : startParticle+nrParts]:

          parton = list(map(float,parton.strip().split()))

          E   = parton[4]
          px  = parton[5]
          py  = parton[6]
          pz  = parton[7]
          C   = parton[13]

          if C != 0:
            pTransverse = np.sqrt(px**2 + py**2)

            if (E-pz)==0:
              rapidity = 10
            elif (E+pz)==0:
              rapidity = 0
            else:
              rapidity = 0.5*np.log((E+pz)/(E-pz))
            
            ALL.append([C,pTransverse,1/pTransverse,rapidity])
      startParticle += nrParts

      if event%1000==0:
        print(event)
     
    self.ALL = ALL
   
    #self.charged          = charged
    #self.pT               = pT
    #self.pTWeight         = pTWeight
    #self.rapidity         = rapidity
    #self.inElasticEvents  = inElasticEvents
    #self.elasticEvents    = elasticEvents

  def writeAnalysis(self):
    fmts = ['%i','%i','%.7e','%.7e','%.7e']
    np.savetxt('data/collected.out',self.ALL,fmt=fmts,delimiter=' , ') 
    

if __name__=="__main__":

  test = QGSM_Distributions()
  test.collectData()
  test.writeAnalysis()


