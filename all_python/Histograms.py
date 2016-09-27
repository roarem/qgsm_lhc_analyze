import ROOT
import numpy as np

'''
hists_input = [["all","ALL"],["nsd","Non-single diffraction"],["etalim","|\eta|<1"],["ptcut","0.3<p_{T}<1.5"]]
Nbins = 50
start = -0.5
stop  = Nbins + start
tree_name = ["Tree_name","Tree_title"]
H_base = Histograms(tree_name,hists_input,Nbins,start,stop)
'''

class Histograms:
    
    def __init__(self,tree_name=[None,None],hists_input=[None],Nbins=0,start=0,stop=0):
        histograms = [] 

        for hist in hists_input:
            histograms.append(ROOT.TH1F(hist[0],hist[1],Nbins,start,stop))

        self.Tree = ROOT.TTree(tree_name[0],tree_name[1])

        for hist,names in zip(histograms,hists_input):
            self.Tree.Branch(names[0],hist)
   
        self.NHistograms    = len(histograms)
        self.histograms     = histograms
        self.Nbins          = Nbins
        self.start          = start
        self.stop           = stop
        self.hists_input    = hists_input

    def getHistograms(self):
        return self.histograms

    def fillHistograms(self,nbnf=None,weighted=None):
        if weighted: 
            for i in range(self.NHistograms):
                self.histograms[i].Fill(nbnf[i,1],nbnf[i,0])
        else:
            for i in range(self.NHistograms):
                self.histograms[i].Fill(nbnf[i,0])

        self.Tree.Fill()

    def store_histograms(self,nf_total=None):
        NColumns = np.zeros((self.Nbins,len(nf_total)+self.NHistograms))
        header = ""
        nf_xaxis = []
        print(nf_total)
        for i in range(len(nf_total)):
            header += 'nf_'+self.hists_input[i][0]+','+self.hists_input[i][0]+','
            nf_xaxis.append(np.linspace(0,nf_total[i],self.Nbins))
        NXaxis = len(nf_xaxis)
        for i in range(self.Nbins):
            for j in range(self.NHistograms):
                NColumns[i,j*2] = nf_xaxis[j][i]
                NColumns[i,j*2+1] = self.histograms[j].GetBinContent(i)
            #for j in range(NXaxis-1,NXaxis+self.NHistograms-1):
        np.savetxt(self.hists_input[0][1]+".out",NColumns,delimiter=',',header=header)
