import ROOT

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
        self.histograms = []
        for hist in hists_input:
            self.histograms.append(ROOT.TH1F(hist[0],hist[1],Nbins,start,stop)

        self.NHistograms = len(self.histograms)
        self.Tree = ROOT.TTree(tree_name[0],tree_name[1])

        for hist,names in zip(self.histograms,hists_input):
            self.Tree.Branch(names[0],hist)
   

    def getHistograms(self):
        return self.histograms

    def fillHistograms(self,nbnf=None,weighted=True):
        if weighted: 
            for i in range(self.NHistograms):
                self.histograms[i].Fill(nbnf[i,1],nbnf[i,0])
        else:
            for i in range(self.NHistograms):
                self.histograms[i].Fill(nbnf[i,0])

        self.Tree.Fill()
