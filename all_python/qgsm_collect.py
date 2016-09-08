from __future__ import print_function
import sys
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
    self.ALL                = []
    self.bcorr              = []
    self.nbnf_all           = []
    self.nbnf_nsd           = []
    self.nbnf_etalim        = []
    self.nbnf_ptcut         = []
    self.elasticEvents      = 0
    self.inElasticEvents    = 0
    self.countedEvents      = 0
    self.bcorr_linecount    = 0

    self.nb_nf_all           = np.zeros(2)
    self.nb_nf_nsd           = np.zeros(2)
    self.nb_nf_etalim        = np.zeros(2)
    self.nb_nf_ptcut         = np.zeros(2)

    Nbins = 50
    start = 0
    stop  = 50

    self.H_all       = ROOT.TH1F("all","ALL",Nbins,start,stop) 
    self.H_nsd       = ROOT.TH1F("nsd","ALL",Nbins,start,stop) 
    self.H_etalim    = ROOT.TH1F("|\eta|<1","ALL",Nbins,start,stop) 
    self.H_ptcut     = ROOT.TH1F("0.3<p_{T}<1.5","ALL",Nbins,start,stop) 
    self.H_all_nf    = ROOT.TH1F("nf","ALL",Nbins,start,stop) 
    self.H_nsd_nf    = ROOT.TH1F("nf","ALL",Nbins,start,stop) 
    self.H_etalim_nf = ROOT.TH1F("nf","ALL",Nbins,start,stop) 
    self.H_ptcut_nf  = ROOT.TH1F("nf","ALL",Nbins,start,stop) 

  def closeFile(self):
    self.finalpr.close()

  def collectData(self,bcorr=False,nbnf=False):
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
    self.nbnf_all_Nevents      = 0 
    self.nbnf_nsd_Nevents      = 0 
    self.nbnf_etalim_Nevents   = 0
    self.nbnf_ptcut_Nevents    = 0 

    eventNr, nrParticles, finalpr, NPOM =\
        self.eventNr, self.nrParticles, self.finalpr, self.NPOM

    self.lineCount      = 0
    total               = len(eventNr)
    tol                 = 10e-9
    idiag_vs_nrParts    = 0 
    countedFlag         = 0
    linjenummer         = 0
    nb = nf             = 0

    eventNr = eventNr[::10]
    nrParticles = nrParticles[::10]
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
          #if nrParts == 2 and count == 0:
          #  idiag_vs_nrParts += 1 
          #  print(event, idiag, linjenummer)
          if C != 0:
            countFlag = 1
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

            nbnf_index = self.eta_check(eta)
            if bcorr:
              if np.abs(eta)<1 and pTransverse > 0.05:
                if pTransverse > 0.3 and pTransverse < 1.5:
                  bcorr_counted = 1
                  self.bcorr.append([event,eta])                  
                  self.bcorr_linecount += 1

            elif nbnf:
              #self.nbnf_all.append([event,eta])
              self.nb_nf_all[nbnf_index] += 1
              nbnf_all_counted = 1
              self.nbnf_all_linecount += 1

              if idiag != 1 or idiag != 6 or idiag != 10:
                nbnf_nsd_counted = 1
                #self.nbnf_nsd.append([event,eta])
                self.nb_nf_nsd[nbnf_index] += 1
                self.nbnf_nsd_linecount += 1

              if np.abs(eta) < 1:
                nbnf_etalim_counted = 1
                #self.nbnf_etalim.append([event,eta])
                self.nb_nf_etalim[nbnf_index] += 1
                self.nbnf_etalim_linecount += 1

                if np.abs(eta) > 0.2 and np.abs(eta) < 0.8:
                  if pTransverse > 0.3 and pTransverse < 1.5:
                    #self.nbnf_ptcut.append([event,eta]) 
                    self.nb_nf_ptcut[nbnf_index] += 1
                    nbnf_ptcut_counted = 1
                    self.nbnf_ptcut_linecount += 1

      
      if countFlag:
        self.H_all.Fill(self.nb_nf_all[1],self.nb_nf_all[0])
        self.H_nsd.Fill(self.nb_nf_nsd[1],self.nb_nf_nsd[0]) 
        self.H_etalim.Fill(self.nb_nf_etalim[1],self.nb_nf_etalim[0]) 
        self.H_ptcut.Fill(self.nb_nf_ptcut[1],self.nb_nf_ptcut[0]) 
        self.H_all_nf.Fill(self.nb_nf_all[0])
        self.H_nsd_nf.Fill(self.nb_nf_nsd[0]) 
        self.H_etalim_nf.Fill(self.nb_nf_etalim[0]) 
        self.H_ptcut_nf.Fill(self.nb_nf_ptcut[0]) 

        
       
      self.nb_nf_all           = np.zeros(2)
      self.nb_nf_nsd           = np.zeros(2)
      self.nb_nf_etalim        = np.zeros(2)
      self.nb_nf_ptcut         = np.zeros(2)

      self.bcorr_Nevents        += bcorr_counted
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
    #self.bcorr          = np.asarray(self.bcorr)
    #self.nbnf_all       = np.asarray(self.nbnf_all)
    #self.nbnf_nsd       = np.asarray(self.nbnf_nsd)
    #self.nbnf_etalim    = np.asarray(self.nbnf_etalim)
    #self.nbnf_ptcut     = np.asarray(self.nbnf_ptcut)
    #self.ALL = np.asarray(self.ALL)
    print("")

  def create_hists(self):
    c1 = ROOT.TCanvas( 'c1', 'Example with Formula', 200, 10, 700, 500 )

    self.H_all.Divide(self.H_all_nf)
    self.H_nsd.Divide(self.H_nsd_nf)
    self.H_etalim.Divide(self.H_etalim_nf)
    self.H_ptcut.Divide(self.H_ptcut_nf)

    self.H_all.SetMarkerStyle(20)
    self.H_all.SetMarkerColor(4)
    self.H_all.Draw("e1")
    self.H_nsd.SetMarkerStyle(21)
    self.H_nsd.SetMarkerColor(3)
    self.H_nsd.Draw("same e1")
    self.H_etalim.SetMarkerStyle(29)
    self.H_etalim.SetMarkerColor(5)
    self.H_etalim.Draw("same e1")
    self.H_ptcut.SetMarkerStyle(33)
    self.H_ptcut.SetMarkerColor(2)
    self.H_ptcut.Draw("same e1")

    leg = ROOT.TLegend(0.8,0.8,0.9,0.9)
    leg.SetFillColor(0)
    leg.AddEntry(self.H_all,"all", "p")
    leg.AddEntry(self.H_nsd,"nsd","p")
    leg.AddEntry(self.H_etalim,"$|\eta|<1$","p")
    leg.AddEntry(self.H_ptcut,"0.3<p_{T}<1.5","p")
    leg.Draw()
    c1.Update()
    raw_input()

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


