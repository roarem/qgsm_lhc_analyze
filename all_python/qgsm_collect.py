from __future__ import print_function
import sys
import Histograms as hist
import ROOT
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
    self.bcorr              = []
    self.elasticEvents      = 0
    self.inElasticEvents    = 0
    self.countedEvents      = 0
    self.bcorr_linecount    = 0


    self.Nbins = 50
    start = -0.5
    stop  = self.Nbins + start

    self.nfnb_file = ROOT.TFile.Open('7000_4M_nfnb_all.root','recreate')
    
    H_all_input =[["all","ALL"],["nsd","Non-single diffraction"],\
                  ["etalim","|\eta|<1"],["ptcut","0.3<p_{T}<1.5"]]  
    H_all_treename = ["H","all"]
    H_div_input =[["all_div","DIV"],["nsd_div","Non-single diffraction"],\
                  ["etalim_div","|\eta|<1"],["ptcut_div","0.3<p_{T}<1.5"]]  
    H_div_treename = ["H","div"]
    H_nf_input =[["nf","NF"],["nsd_nf","Non-single diffraction"],\
                 ["etalim_nf","|\eta|<1"],["ptcut_nf","0.3<p_{T}<1.5"]]  
    H_nf_treename = ["H","nf"]

    self.trees = [hist.Histograms(H_all_treename,H_all_input,self.Nbins,start,stop),\
                  hist.Histograms(H_div_treename,H_div_input,self.Nbins,start,stop),\
                  hist.Histograms(H_nf_treename,H_nf_input,self.Nbins,start,stop)]    

  def closeFile(self):
    self.finalpr.close()

  def collectData(self,bcorr=False,nbnf=False):
    nbnf_all           = np.zeros(2)
    nbnf_nsd           = np.zeros(2)
    nbnf_etalim        = np.zeros(2)
    nbnf_ptcut         = np.zeros(2)
    nfnb_countFlag             = 0
    bcorr_counted              = 0
    nbnf_all_counted           = 0
    nbnf_nsd_counted           = 0
    nbnf_etalim_counted        = 0
    nbnf_ptcut_counted         = 0 
    self.nbnf_all_linecount    = 0
    self.nbnf_nsd_linecount    = 0
    self.nbnf_etalim_linecount = 0
    self.nbnf_ptcut_linecount  = 0 
    self.bcorr_Nevents         = 0      
    self.bcorr_Nevents_etalim  = 0 
    self.bcorr_Nevents_ptcut   = 0 
    self.nbnf_all_Nevents      = 0 
    self.nbnf_nsd_Nevents      = 0 
    self.nbnf_etalim_Nevents   = 0
    self.nbnf_ptcut_Nevents    = 0 
    nf_max                   = [0,0,0,0]

    eventNr, nrParticles, finalpr, NPOM =\
        self.eventNr, self.nrParticles, self.finalpr, self.NPOM

    self.lineCount      = 0
    total               = len(eventNr)
    tol                 = 10e-9
    idiag_vs_nrParts    = 0 
    countedFlag         = 0
    linjenummer         = 0
    nb = nf             = 0

    eventNr = eventNr[::]
    nrParticles = nrParticles[::]
    print("Collecting data...")

    for event,nrParts in zip(eventNr,nrParticles):
      for count in range(nrParts):
        store = False
        parton = finalpr.readline()
        linjenummer += 1
            
        parton = list(map(float,parton.strip().split()))

        E     = parton[4]
        px    = parton[5]
        py    = parton[6]
        pz    = parton[7]
        idiag = parton[10]
        C     = parton[13]
        
        if count == 0 and idiag == 4:
          self.elasticEvents += 1
          finalpr.readline()
          linjenummer += 1
          break

        else:
          if C != 0:
            p = np.sqrt(px**2 + py**2 + pz**2) 
            pTransverse = np.sqrt(px**2 + py**2)
            if abs(E-pz)<tol:
              rapidity = 10
            elif abs(E+pz)<tol:
              rapidity = 0
            else:
              rapidity = 0.5*np.log((E+pz)/(E-pz))

            if abs(p-pz)<tol:
              eta = 10
            elif abs(p+pz)<tol:
              eta = 0
            else:
              eta = 0.5*np.log((p+pz)/(p-pz))

            if bcorr:
              bcorr_countFlag = 1
              bcorr_counted   = 1
              if np.abs(eta)<1 and pTransverse > 0.05:
                bcorr_counted_etalim = 1
                if pTransverse > 0.3 and pTransverse < 1.5:
                  bcorr_counted_ptcut = 1
                  self.bcorr.append(np.asarray([event,eta]))
                  self.bcorr_linecount += 1
                  

            if nbnf:
              nbnf_index = self.eta_check(eta)
              nfnb_countFlag = 1
              nbnf_all[nbnf_index] += 1
              nbnf_all_counted = 1
              self.nbnf_all_linecount += 1

              if idiag != 1 and idiag != 6 and idiag != 10:
                nbnf_nsd_counted = 1
                nbnf_nsd[nbnf_index] += 1
                self.nbnf_nsd_linecount += 1

              if np.abs(eta) < 1:
                nbnf_etalim_counted = 1
                nbnf_etalim[nbnf_index] += 1
                self.nbnf_etalim_linecount += 1

                if np.abs(eta) > 0.2 and np.abs(eta) < 0.8:
                  if pTransverse > 0.3 and pTransverse < 1.5:
                    nbnf_ptcut[nbnf_index] += 1
                    nbnf_ptcut_counted = 1
                    self.nbnf_ptcut_linecount += 1

          
      if nfnb_countFlag:
        nfnb_countFlag = 0
        weighted = [True, True, False]
        nbnf_lists = np.asarray([nbnf_all,nbnf_nsd,nbnf_etalim,nbnf_ptcut])
        for i in range(len(self.trees)):
            self.trees[i].fillHistograms(nbnf_lists,weighted[i])
        nf_max[0]   = nbnf_all[0] if nf_max[0]<nbnf_all[0] else nf_max[0]
        nf_max[1]   = nbnf_nsd[0] if nf_max[1]<nbnf_nsd[0] else nf_max[1]
        nf_max[2]   = nbnf_etalim[0] if nf_max[2]<nbnf_etalim[0] else nf_max[2]
        nf_max[3]   = nbnf_ptcut[0] if nf_max[3]<nbnf_ptcut[0] else nf_max[3]

        nbnf_all           = np.zeros(2)
        nbnf_nsd           = np.zeros(2)
        nbnf_etalim        = np.zeros(2)
        nbnf_ptcut         = np.zeros(2)

      self.bcorr_Nevents        += bcorr_counted
      self.bcorr_Nevents_etalim += bcorr_counted_etalim
      self.bcorr_Nevents_ptcut  += bcorr_counted_ptcut
      self.nbnf_all_Nevents     += nbnf_all_counted
      self.nbnf_nsd_Nevents     += nbnf_nsd_counted
      self.nbnf_etalim_Nevents  += nbnf_etalim_counted
      self.nbnf_ptcut_Nevents   += nbnf_ptcut_counted

      bcorr_counted         = 0
      nbnf_all_counted      = 0
      nbnf_nsd_counted      = 0
      nbnf_etalim_counted   = 0
      nbnf_ptcut_counted    = 0 

      if event%(total//100)==0:
        self.msg("%i"%(float(event)/total*100)+'%')#+'  Event #%i'%(event))

    print("")
    if nbnf:
      for i in range(len(self.trees[1].getHistograms())):
        self.trees[1].getHistograms()[i].Divide(self.trees[2].getHistograms()[i])

      for i in range(len(self.trees)):
        self.trees[i].store_histograms(nf_max)
      
        #print("all: {}, nsd: {}, etalim: {}, ptcut: {}"\
        #.format(self.nbnf_all_Nevents,self.nbnf_nsd_Nevents,self.nbnf_etalim_Nevents,self.nbnf_ptcut_Nevents))
      self.nfnb_file.Write()
      self.nfnb_file.Close()
    self.bcorr = np.asarray(self.bcorr)
    print(self.bcorr_linecount)
    print(self.bcorr_Nevents,self.bcorr_Nevents_etalim,self.bcorr_Nevents_ptcut)
    np.savetxt("bcorr_temp.out",self.bcorr,delimiter=',')
     
  def create_nbnf_hists(self):
    self.nfnb_file.Write()
    self.nfnb_file.Close()

  def eta_check(self,eta):
    if eta > 0:
      return 0
    else:
      return 1

  def msg(self,text):
    text = "  "+text+chr(13)
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


