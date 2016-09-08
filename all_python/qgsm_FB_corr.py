from __future__ import print_function
import ROOT
import numpy as np
import matplotlib.pyplot as plt

class FB:
    def __init__(self,energy):
        self.energy = energy

    def b_corr_collect_file(self):
        print("Reading from meta.out and collected.out ...")
        self.linecount, self.Nevents= np.loadtxt('../data/meta.out',dtype=np.int32)
        self.event, self.pT, self.eta = \
                np.loadtxt("../data/collected.out",unpack=True,usecols=[0,3,5])
        print("Arrays put in memory")

    def b_corr_collect_mem(self):
        self.event      = bcorr[:,0]
        self.eta        = bcorr[:,2]
        self.linecount  = bcorr_linecount
        self.Nevents    = bcorr_countedEvents
    
    def nB_nF_collect_mem(self,all_linecount,nsd_linecount,etalim_linecount,ptcut_linecount,\
                               all_Nevents,nsd_Nevents,etalim_Nevents,ptcut_Nevents,\
                               all_list,nsd_list,etalim_list,ptcut_list):
        self.nbnf_all_linecount     = all_linecount
        self.nbnf_nsd_linecount     = nsd_linecount
        self.nbnf_etalim_linecount  = etalim_linecount
        self.nbnf_ptcut_linecount   = ptcut_linecount
        self.nbnf_all_Nevents       = all_Nevents
        self.nbnf_nsd_Nevents       = nsd_Nevents
        self.nbnf_etalim_Nevents    = etalim_Nevents
        self.nbnf_ptcut_Nevents     = ptcut_Nevents 
        self.nbnf_all_event         = all_list[:,0] 
        self.nbnf_nsd_event         = nsd_list[:,0] 
        self.nbnf_etalim_event      = etalim_list[:,0]    
        self.nbnf_ptcut_event       = ptcut_list[:,0]    
        self.nbnf_all_eta           = all_list[:,1]     
        self.nbnf_nsd_eta           = nsd_list[:,1]     
        self.nbnf_etalim_eta        = etalim_list[:,1]     
        self.nbnf_ptcut_eta         = ptcut_list[:,1]   

        Nbins = 50
        start = 0
        stop  = 50
        c1 = ROOT.TCanvas( 'c1', 'Example with Formula', 200, 10, 700, 500 )

        print("Starting histogram for all")
        H_all_in = ROOT.TH1F("all","ALL",Nbins,start,stop) 
        H_all = self.nB_nF(self.nbnf_all_linecount,self.nbnf_all_Nevents,\
                           self.nbnf_all_event,self.nbnf_all_eta,H_all_in)
        H_all.SetMarkerStyle(20)
        H_all.SetMarkerColor(4)
        H_all.Draw("e1")

        print("Starting histogram for nsd")
        H_nsd_in = ROOT.TH1F("nsd","ALL",Nbins,start,stop)
        H_nsd = self.nB_nF(self.nbnf_nsd_linecount,self.nbnf_nsd_Nevents,\
                           self.nbnf_nsd_event,self.nbnf_nsd_eta,H_nsd_in)
        H_nsd.SetMarkerStyle(21)
        H_nsd.SetMarkerColor(3)
        H_nsd.Draw("same e1")

        print("Starting histogram for etalim")
        H_etalim_in = ROOT.TH1F("etalim","ALL",Nbins,start,stop)
        H_etalim = self.nB_nF(self.nbnf_etalim_linecount,self.nbnf_etalim_Nevents,\
                                self.nbnf_etalim_event,self.nbnf_etalim_eta,H_etalim_in)
        H_etalim.SetMarkerStyle(29)
        H_etalim.SetMarkerColor(5)
        H_etalim.Draw("same e1")

        print("Starting histogram for cuts")
        H_cuts_in = ROOT.TH1F("cuts","ALL",Nbins,start,stop) 
        H_cuts = self.nB_nF(self.nbnf_ptcut_linecount,self.nbnf_ptcut_Nevents,\
                            self.nbnf_ptcut_event,self.nbnf_ptcut_eta,H_cuts_in)
        H_cuts.SetMarkerStyle(33)
        H_cuts.SetMarkerColor(2)
        H_cuts.Draw("same e1")

        leg = ROOT.TLegend(0.8,0.8,0.9,0.9)
        leg.SetFillColor(0)
        leg.AddEntry(H_all,"all", "p")
        leg.AddEntry(H_nsd,"nsd","p")
        leg.AddEntry(H_etalim,"$|\eta|<1$","p")
        leg.AddEntry(H_cuts,"0.3<p_{T}<1.5","p")
        leg.Draw()
        c1.Update()
        raw_input()

    def nB_nF(self,linecount,Nevents,event,eta,H):
        old_event = event[0]
        Nbins = 50
        start = 0
        stop  = 50
        nb = nf = 0 
        #c1 = ROOT.TCanvas( 'c1', 'Example with Formula', 200, 10, 700, 500 )
        #H = ROOT.TH1F("nbdnf", "ALL", Nbins,start,stop)
        Hnf = ROOT.TH1F("dnf", "dnf dist",Nbins,start,stop)
        
        for line in range(linecount):
            if(old_event != event[line]):
                old_event = event[line]
                H.Fill(nb,nf)
                Hnf.Fill(nf)
                nb = nf = 0 
            
            #if np.abs(self.eta[line]) > 0.2 and np.abs(self.eta[line]) < 0.8: 
            if eta[line] < 0:
                nb += 1
            else:
                nf += 1

        H.Divide(Hnf)
        return H
        
    def b_corr_count(self,Nevents,bcorr,linecount):
        print("Starting counting ...")
        event = bcorr[:,0]
        eta   = bcorr[:,1]
        Nevents = float(Nevents)
        self.b_corr = [] 
        Ngaps  = [7,3,2,1]
        delta_eta_gap = [0.1,0.2,0.2,0.0]
        i = 0 
        print("number of events {}".format(Nevents))
        for delta_eta,measure in zip(delta_eta_gap,[0.2,0.4,0.6,0.8]):
            eta_lower_lim = 0
            eta_upper_lim = measure

            for j in range(Ngaps[i]):
                nF = nB = nBnF = nFnF = nb = nf = 0
                old_event = 1

                for line in range(linecount):
                    if(old_event != event[line]):
                        old_event = event[line]
                        nBnF += nb*nf
                        nFnF += nf*nf
                        nB   += nb
                        nF   += nf
                        nb = nf = 0
                    if (eta[line] >= eta_lower_lim and eta[line] <= eta_upper_lim):
                        nf += 1
                    elif (eta[line] >= -eta_upper_lim and eta[line] <= -eta_lower_lim):
                        nb += 1
                nBnF += nb*nf
                nFnF += nf*nf
                nB   += nb
                nF   += nf

                b_corr_numb = (nBnF - nB*nF/Nevents)/(nFnF - nF*nF/Nevents)
                self.b_corr.append(b_corr_numb)
                print("----------")
                print(eta_lower_lim," - ", eta_upper_lim)
                print("nB: {}, nF: {}, nFnF: {}, nBnF: {},".format(nB,nF,nFnF,nBnF))
                #print(eta_lower_lim," - ", eta_upper_lim)
                print(b_corr_numb)
                print("----------")

                eta_lower_lim += delta_eta
                eta_upper_lim += delta_eta
            print("\n")
            i+=1

        print("Saving results to file b_corr.py.out ...")
        np.savetxt('b_corr.py.out',self.b_corr)

    def b_corr_plot(self,from_file=False):

        if from_file:
          self.b_corr = np.loadtxt("b_corr.py.out")
    
        fig, ax = plt.subplots()
        if self.energy==900:
            experimental = np.asarray([0.212,0.203,0.193,0.182,0.172,0.163,\
                    0.159,0.335,0.300,0.274,0.406,0.368,0.452])
            ax.set_xlim(-0.1,1.3)
            ax.set_ylim(0,0.7)
        if self.energy==7000:
            experimental = np.asarray([0.366,0.358,0.345,0.334,0.327,0.316,\
                    0.311,0.521,0.487,0.463,0.598,0.564,0.643])
            ax.set_xlim(-0.1,1.3)
            ax.set_ylim(0,0.75)

        ax.plot([0.0,0.2,0.4,0.6,0.8,1.0,1.2],self.b_corr[:7],'bo-',markersize=15,label='$\delta\eta = 0.2$')
        ax.plot([0.0,0.2,0.4,0.6,0.8,1.0,1.2],experimental[:7],'ro-',markersize=15)
        
        ax.plot([0.0,0.4,0.8],self.b_corr[7:10],'b*-',markersize=15,label='$\delta\eta = 0.4$')
        ax.plot([0.0,0.4,0.8],experimental[7:10],'r*-',markersize=15)
        
        ax.plot([0.0,0.4],self.b_corr[10:12],'b^-',markersize=15,label='$\delta\eta = 0.6$')
        ax.plot([0.0,0.4],experimental[10:12],'r^-',markersize=15)
        
        ax.plot([0],self.b_corr[12],'bs-',markersize=15,label='$\delta\eta = 0.8$')
        ax.plot([0],experimental[12],'rs-',markersize=15)
        
        ax.set_title('Blue simulated, red experimental')
        plt.legend()
        
        plt.show()

    def b_corr_compare(self):
        b_corr_py = np.loadtxt('b_corr.py.out')
        b_corr_c  = np.loadtxt('../c++/build/b_corr.out')
        print(b_corr_py-b_corr_c)

if __name__=="__main__":
    energy = 7000
    FBs = FB(energy)
    FBs.b_corr_collect()
    FBs.b_corr_count()
    FBs.b_corr_plot()
    FBs.b_corr_compare()
