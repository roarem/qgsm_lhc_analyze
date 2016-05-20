#include <fstream>
#include <iostream>
#include <vector>

void rootTest()
{
  int Nbins = 100;
  int start = -20;
  int stop = 20;
  //TFile *test = new TFile("test.root","RECREATE");
  TH1F *HA   = new TH1F("dnf0"  ,"tot dist",Nbins,start,stop);
  TH1F *Hnf  = new TH1F("dnf"   ,"dnf dist",Nbins,start,stop);
  TH1F *H0nf = new TH1F("dnf1"  ,"H0 dist" ,Nbins,start,stop);
  TH1F *H2nf = new TH1F("dnf2"  ,"H2 dist" ,Nbins,start,stop);
  TH1F *H4nf = new TH1F("dnf3"  ,"H4 dist" ,Nbins,start,stop);
  TH1F *H6nf = new TH1F("dnf4"  ,"H6 dist" ,Nbins,start,stop);
  TH1F *H8nf = new TH1F("dnf5"  ,"H8 dist" ,Nbins,start,stop);
  TH1F *H0nb = new TH1F("nbdnf1","H0 dist" ,Nbins,start,stop);
  TH1F *H2nb = new TH1F("nbdnf2","H2 dist" ,Nbins,start,stop);
  TH1F *H4nb = new TH1F("nbdnf3","H4 dist" ,Nbins,start,stop);
  TH1F *H6nb = new TH1F("nbdnf4","H6 dist" ,Nbins,start,stop);
  TH1F *H8nb = new TH1F("nbdnf5","H8 dist" ,Nbins,start,stop);
  std::vector<int> NPOMSs; 
  std::vector<int> NPOMHs; 
  std::vector<int> nf; 
  std::vector<int> nb; 
  std::string line;
  std::string file = "NPOM.out";
  std::ifstream f(file.c_str());
  int npoms, npomh, nfs, nbs;
  int N = 0;
  while (std::getline(f, line)){
    std::istringstream ss(line);
    ss >> npoms >> npomh >> nfs >> nbs;
    NPOMSs.push_back(npoms);
    NPOMHs.push_back(npomh);
    nf.push_back(nfs);
    nb.push_back(nbs);
    N++;
  }
  int S = 4; 
  for (int i=0 ; i<N ; i++){
    std::cout<<nf[i]<<", ";
    //if (NPOMSs[i] == S){
      HA->Fill(nf[i],nb[i]);
      Hnf->Fill(nf[i]);
      if (NPOMHs[i] == 0){
        H0nf->Fill(nf[i]); 
        H0nb->Fill(nf[i],nb[i]); 
      }
      else if (NPOMHs[i]==2){
        H2nf->Fill(nf[i]); 
        H2nb->Fill(nf[i],nb[i]); 
      }
      else if (NPOMHs[i]==4){
        H4nf->Fill(nf[i]); 
        H4nb->Fill(nf[i],nb[i]); 
      }
      else if (NPOMHs[i]==6){
        H6nf->Fill(nf[i]); 
        H6nb->Fill(nf[i],nb[i]); 
      }
      else if (NPOMHs[i]==8){
        H8nf->Fill(nf[i]); 
        H8nb->Fill(nf[i],nb[i]); 
      }
    //}
  }
  std::cout<<std::endl;
  //test->Write();
  HA->Draw("B");
}

