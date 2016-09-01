import numpy as np
import matplotlib.pyplot as plt

class FB:
    def __init__(self,energy):
        self.energy = energy

    def b_corr_collect_file(self):
        print("Reading from meta.out and collected.out ...")
        self.linecount, self.Nevents= np.loadtxt('../data/meta.out',dtype=np.int32)
        self.event, self.pT, self.eta = \
                np.loadtxt("../data/collected.out",unpack=True,usecols=[0,3,5])#,dtype=np.int32)
        print("Arrays put in memory")

    def b_corr_collect_mem(self,ALL,linecount,Nevents):
        self.event      = ALL[:,0]
        self.pT         = ALL[:,1]
        self.eta        = ALL[:,2]
        self.linecount  = linecount
        self.Nevents    = Nevents
        print(len(self.event), len(self.pT), len(self.eta))

    def b_corr_count(self):
        print("Starting counting ...")
        self.b_corr = [] 
        Ngaps  = [7,3,2,1]
        delta_eta_gap = [0.1,0.2,0.2,0.0]
        i = 0 
        for delta_eta,measure in zip(delta_eta_gap,[0.2,0.4,0.6,0.8]):
            eta_lower_lim = 0
            eta_upper_lim = measure

            for j in range(Ngaps[i]):
                nF = nB = nBnF = nFnF = nb = nf = 0
                old_event = 1

                for line in range(self.linecount):
                    #if (self.pT[line] > 0.3 and self.pT[line] < 1.5):
                    if(old_event != self.event[line]):
                        old_event = self.event[line]
                        nBnF += nb*nf
                        nFnF += nf*nf
                        nB   += nb
                        nF   += nf
                        nb = nf = 0

                    if (self.eta[line] >= eta_lower_lim and self.eta[line] <= eta_upper_lim):
                        nf += 1
                    elif (self.eta[line] >= -eta_upper_lim and self.eta[line] <= -eta_lower_lim):
                        nb += 1

                nBnF += nb*nf
                nFnF += nf*nf
                nB   += nb
                nF   += nf

                b_corr_numb = (nBnF - nB*nF/self.Nevents)/(nFnF - nF*nF/self.Nevents)
                self.b_corr.append(b_corr_numb)
                print("----------")
                print(eta_lower_lim," - ", eta_upper_lim)
                print(b_corr_numb)
                print("----------")

                eta_lower_lim += delta_eta
                eta_upper_lim += delta_eta
            print("\n")
            i+=1

        print("Saving results to file b_corr.py.out ...")
        np.savetxt('b_corr.py.out',self.b_corr)

    def b_corr_plot(self):
    
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
