from __future__ import print_function
#import sys
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

    self.trees = [hist(H_all_treename,H_all_input,self.Nbins,start,stop),\
                  hist(H_div_treename,H_div_input,self.Nbins,start,stop),\
                  hist(H_nf_treename,H_nf_input,self.Nbins,start,stop)]    

    #self.H_all       = ROOT.TH1F("all","ALL",self.Nbins,start,stop) 
    #self.H_nsd       = ROOT.TH1F("nsd","ALL",self.Nbins,start,stop) 
    #self.H_etalim    = ROOT.TH1F("|\eta|<1","ALL",self.Nbins,start,stop) 
    #self.H_ptcut     = ROOT.TH1F("0.3<p_{T}<1.5","ALL",self.Nbins,start,stop) 
    #self.H_all_div   = ROOT.TH1F("all_div","ALL",self.Nbins,start,stop) 
    #self.H_nsd_div   = ROOT.TH1F("nsd_div","ALL",self.Nbins,start,stop) 
    #self.H_etalim_div= ROOT.TH1F("|\eta|<1 div","ALL",self.Nbins,start,stop) 
    #self.H_ptcut_div = ROOT.TH1F("0.3<p_{T}<1.5 div","ALL",self.Nbins,start,stop) 
    #self.H_all_nf    = ROOT.TH1F("nf_all","ALL",self.Nbins,start,stop) 
    #self.H_nsd_nf    = ROOT.TH1F("nf_nsd","ALL",self.Nbins,start,stop) 
    #self.H_etalim_nf = ROOT.TH1F("nf_etalim","ALL",self.Nbins,start,stop) 
    #self.H_ptcut_nf  = ROOT.TH1F("nf_ptcut","ALL",self.Nbins,start,stop) 


    #self.H = ROOT.TTree("H","all")
    #self.H.Branch("all",self.H_all)
    #self.H.Branch("nsd",self.H_nsd)
    #self.H.Branch("etalim",self.H_etalim)
    #self.H.Branch("ptcut",self.H_ptcut)

    #self.H_nf = ROOT.TTree("H_nf","nf")
    #self.H_nf.Branch("all_nf",self.H_all_nf)
    #self.H_nf.Branch("nsd_nf",self.H_nsd_nf)
    #self.H_nf.Branch("etalim_nf",self.H_etalim_nf)
    #self.H_nf.Branch("ptcut_nf",self.H_ptcut_nf)

    #self.H_div = ROOT.TTree("H_div","div")
    #self.H_div.Branch("all_div",self.H_all_div)
    #self.H_div.Branch("nsd_div",self.H_nsd_div)
    #self.H_div.Branch("etalim_div",self.H_etalim_div)
    #self.H_div.Branch("ptcut_div",self.H_ptcut_div)

  def closeFile(self):
    self.finalpr.close()

  def collectData(self,bcorr=False,nbnf=False):
    nbnf_all           = np.zeros(2)
    nbnf_nsd           = np.zeros(2)
    nbnf_etalim        = np.zeros(2)
    nbnf_ptcut         = np.zeros(2)
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
    self.nf_all                = 0
    self.nf_nsd                = 0
    self.nf_etalim             = 0
    self.nf_ptcut              = 0

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
              if np.abs(eta)<1 and pTransverse > 0.05:
                if pTransverse > 0.3 and pTransverse < 1.5:
                  bcorr_counted = 1
                  self.bcorr.append([event,eta])                  
                  self.bcorr_linecount += 1
                  

            elif nbnf:
              nbnf_index = self.eta_check(eta)
              nfnb_countFlag = 1
              self.nbnf_all[nbnf_index] += 1
              nbnf_all_counted = 1
              self.nbnf_all_linecount += 1

              if idiag != 1 and idiag != 6 and idiag != 10:
                nbnf_nsd_counted = 1
                self.nbnf_nsd[nbnf_index] += 1
                self.nbnf_nsd_linecount += 1

              if np.abs(eta) < 1:
                nbnf_etalim_counted = 1
                self.nbnf_etalim[nbnf_index] += 1
                self.nbnf_etalim_linecount += 1

                if np.abs(eta) > 0.2 and np.abs(eta) < 0.8:
                  if pTransverse > 0.3 and pTransverse < 1.5:
                    self.nbnf_ptcut[nbnf_index] += 1
                    nbnf_ptcut_counted = 1
                    self.nbnf_ptcut_linecount += 1

      
      if nfnb_countFlag:
        weighted = [True, True, False]
        nbnf = np.asarray([nbnf_all,nbnf_nsd,nbnf_etalim,nbnf_ptcut])
        for i in range(len(self.trees)):
            self.trees[i].fillHistograms(nbnf,weighted)

      self.nf_all            += self.nbnf_all[0]
      self.nf_nsd            += self.nbnf_nsd[0]
      self.nf_etalim         += self.nbnf_etalim[0]
      self.nf_ptcut          += self.nbnf_ptcut[0]
      nbnf_all           = np.zeros(2)
      nbnf_nsd           = np.zeros(2)
      nbnf_etalim        = np.zeros(2)
      nbnf_ptcut         = np.zeros(2)

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

    print("")
    if nbnf:
      self.H_all_div.Divide(self.H_all_nf)
      self.H_nsd_div.Divide(self.H_nsd_nf)
      self.H_etalim_div.Divide(self.H_etalim_nf)
      self.H_ptcut_div.Divide(self.H_ptcut_nf)
    
      print("all: {}, nsd: {}, etalim: {}, ptcut: {}"\
      .format(self.nbnf_all_Nevents,self.nbnf_nsd_Nevents,self.nbnf_etalim_Nevents,self.nbnf_ptcut_Nevents))

  def create_nbnf_hists(self):
    #self.nfnb_file = ROOT.TFile.Open('7000_4M_nfnb_all.root','recreate')
    #c1 = ROOT.TCanvas( 'c1', '<n_{B}(n_{F})>', 200, 10, 700, 500 )

    #self.H_all.Divide(self.H_all_nf)
    #self.H_nsd.Divide(self.H_nsd_nf)
    #self.H_etalim.Divide(self.H_etalim_nf)
    #self.H_ptcut.Divide(self.H_ptcut_nf)

    #self.H_all.SetMarkerStyle(20)
    #self.H_all.SetMarkerColor(4)
    #self.H_all.Draw("e1")

    #self.H_nsd.SetMarkerStyle(21)
    #self.H_nsd.SetMarkerColor(3)
    #self.H_nsd.Draw("same e1")

    #self.H_etalim.SetMarkerStyle(29)
    #self.H_etalim.SetMarkerColor(5)
    #self.H_etalim.Draw("same e1")

    #self.H_ptcut.SetMarkerStyle(33)
    #self.H_ptcut.SetMarkerColor(2)
    #self.H_ptcut.Draw("same e1")

    #leg = ROOT.TLegend(0.8,0.8,0.9,0.9)
    #leg.SetFillColor(0)
    #leg.AddEntry(self.H_all,"all","p")
    #leg.AddEntry(self.H_nsd,"nsd","p")
    #leg.AddEntry(self.H_etalim,"$|\eta|<1$","p")
    #leg.AddEntry(self.H_ptcut,"|\eta|<1 & 0.3<p_{T}<1.5","p")
    #leg.Draw()
    #c1.Update()
    
    self.nfnb_file.Write()
    self.nfnb_file.Close()
    #raw_input()

  def store_bin_data(self):
    histograms = np.zeros((self.Nbins,16))
    header = "nf_all,nf_nsd,nf_etalim,nf_ptcut,"+\
             "H_all,H_nsd,H_etalim,H_ptcut,H_all_div,H_nsd_div,"+\
             "H_etalim_div,H_ptcut_div,H_all_nf,H_nsd_nf,H_etalim_nf,H_ptcut_nf"
    nf_all    = np.linspace(0,self.nf_all,self.Nbins)
    nf_nsd    = np.linspace(0,self.nf_nsd,self.Nbins)
    nf_etalim = np.linspace(0,self.nf_etalim,self.Nbins)
    nf_ptcut  = np.linspace(0,self.nf_ptcut,self.Nbins)
    for i in range(self.Nbins):
      histograms[i,0] = nf_all[i]
      histograms[i,1] = nf_nsd[i]
      histograms[i,2] = nf_etalim[i]
      histograms[i,3] = nf_ptcut[i]
      histograms[i,4] = self.H_all.GetBinContent(i) 
      histograms[i,5] = self.H_nsd.GetBinContent(i) 
      histograms[i,6] = self.H_etalim.GetBinContent(i) 
      histograms[i,7] = self.H_ptcut.GetBinContent(i) 

      histograms[i,8] = self.H_all_div.GetBinContent(i) 
      histograms[i,9] = self.H_nsd_div.GetBinContent(i) 
      histograms[i,10] = self.H_etalim_div.GetBinContent(i) 
      histograms[i,11] = self.H_ptcut_div.GetBinContent(i) 

      histograms[i,12]  = self.H_all_nf.GetBinContent(i) 
      histograms[i,13]  = self.H_nsd_nf.GetBinContent(i) 
      histograms[i,14] = self.H_etalim_nf.GetBinContent(i) 
      histograms[i,15] = self.H_ptcut_nf.GetBinContent(i) 

    np.savetxt("all_histograms.out",histograms,delimiter=',',header=header)
    
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


