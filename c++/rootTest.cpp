#include <fstream>
#include <iostream>
#include <vector>

void rootTest();
void larissa_root(TLegend& leg);

void rootTest()
{
  const Int_t Nbins = 50;
  const Int_t start = 0;
  const Int_t stop  = 50;

  Int_t logStart;
  if(start>0){logStart = log(start);}
  else{logStart = 0;}
  const Double_t logStop  = log(stop);
  const Double_t logStep  = (logStop-logStart)/Nbins;
  Double_t edges[Nbins+1];
  Double_t tempStep = logStart;
  for (int i=0 ; i<Nbins+1 ; i++){
    edges[i] = exp(tempStep);
    tempStep+= logStep;
  }
  
  //TFile *test = new TFile("test.root","RECREATE");
  
  TH1F *HA   = new TH1F("dnf0"  ,"ALL",Nbins,start,stop);
  TH1F *Hnf  = new TH1F("dnf"   ,"dnf dist",Nbins,start,stop);
//  TH1F *H0nf = new TH1F("dnf1"  ,"H0 dist" ,Nbins,edges);
//  TH1F *H2nf = new TH1F("dnf2"  ,"H2 dist" ,Nbins,edges);
//  TH1F *H4nf = new TH1F("dnf3"  ,"H4 dist" ,Nbins,edges);
//  TH1F *H6nf = new TH1F("dnf4"  ,"H6 dist" ,Nbins,edges);
//  TH1F *H8nf = new TH1F("dnf5"  ,"H8 dist" ,Nbins,edges);
//  TH1F *H0nb = new TH1F("nbdnf1","H0 dist" ,Nbins,edges);
//  TH1F *H2nb = new TH1F("nbdnf2","H2 dist" ,Nbins,edges);
//  TH1F *H4nb = new TH1F("nbdnf3","H4 dist" ,Nbins,edges);
//  TH1F *H6nb = new TH1F("nbdnf4","H6 dist" ,Nbins,edges);
//  TH1F *H8nb = new TH1F("nbdnf5","H8 dist" ,Nbins,edges);
 
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
    ss >> nfs >> nbs >> npoms >> npomh;
    NPOMSs.push_back(npoms);
    NPOMHs.push_back(npomh);
    nf.push_back(nfs);
    nb.push_back(nbs);
    N++;
  }
  
  int S = 4; 
  for (int i=0 ; i<N ; i++){
    //if (NPOMSs[i] == S){
      HA->Fill(nf[i],nb[i]);
      Hnf->Fill(nf[i]);
      //if (NPOMHs[i] == 0){
      //  H0nf->Fill(nf[i]); 
      //  H0nb->Fill(nf[i],nb[i]); 
      //}
      //else if (NPOMHs[i]==2){
      //  H2nf->Fill(nf[i]); 
      //  H2nb->Fill(nf[i],nb[i]); 
      //}
      //else if (NPOMHs[i]==4){
      //  H4nf->Fill(nf[i]); 
      //  H4nb->Fill(nf[i],nb[i]); 
      //}
      //else if (NPOMHs[i]==6){
      //  H6nf->Fill(nf[i]); 
      //  H6nb->Fill(nf[i],nb[i]); 
      //}
      //else if (NPOMHs[i]==8){
      //  H8nf->Fill(nf[i]); 
      //  H8nb->Fill(nf[i],nb[i]); 
      //}
    //}
  }

  //test->Write();
  HA->Divide(Hnf);
  HA->SetMarkerStyle(4);
  HA->SetMarkerSize(2);
  //HA->SetLineWidth(1);
  HA->Draw("HIST P");
  leg = new TLegend(0.8,0.8,0.9,0.9);
  leg->AddEntry(HA,"<nb(nf)>", "p");
  larissa_root(*leg);
  leg->Draw();
  /*
  H0nb->Divide(H0nf);
  H2nb->Divide(H2nf);
  H4nb->Divide(H4nf);
  H6nb->Divide(H6nf);
  H8nb->Divide(H8nf);

  H0nb->SetLineColor(2);
  H2nb->SetLineColor(3);
  H4nb->SetLineColor(4);
  H6nb->SetLineColor(5);
  H8nb->SetLineColor(6);
  H0nb->SetLineWidth(3);
  H2nb->SetLineWidth(3);
  H4nb->SetLineWidth(3);
  H6nb->SetLineWidth(3);
  H8nb->SetLineWidth(3);

  H0nb->Draw("hist,c,same");
  H2nb->Draw("hist,c,same");
  H4nb->Draw("hist,c,same");
  H6nb->Draw("hist,c,same");
  H8nb->Draw("hist,c,same");

  legend = new TLegend(0.8,0.2,0.9,0.4);
  legend->AddEntry(H0nb,"H0");
  legend->AddEntry(H2nb,"H2");
  legend->AddEntry(H4nb,"H4");
  legend->AddEntry(H6nb,"H6");
  legend->AddEntry(H8nb,"H8");
  legend->AddEntry(HA,"Total");
  legend->Draw();
  */
}

void larissa_root(TLegend &leg)
{
  const Int_t Nbins = 50;
  const Int_t start = 0;
  const Int_t stop  = 50;

  Int_t logStart;
  if(start>0){logStart = log(start);}
  else{logStart = 0;}
  const Double_t logStop  = log(stop);
  const Double_t logStep  = (logStop-logStart)/Nbins;
  Double_t edges[Nbins+1];
  Double_t tempStep = logStart;
  for (int i=0 ; i<Nbins+1 ; i++){
    edges[i] = exp(tempStep);
    tempStep+= logStep;
  }

  TH1F *HA  = new TH1F("dnf", "ALL", Nbins,edges);
  TH1F *HB  = new TH1F("dnf", "ALL", Nbins,edges);

  std::vector<Double_t> ratio;
  std::vector<Double_t> eratio;
  std::vector<Double_t> dnfnb;
  std::vector<Double_t> dnf;
  std::vector<Double_t> X;

  std::string line;
  std::string file = "fort.16_13t";
  std::ifstream f(file.c_str());
  Double_t a,b,c,ratios, eratios;
  int N = 0;
  while (std::getline(f, line)){
    std::istringstream ss(line);
    ss >> a >> b >> c >> ratios >> eratios;

    X.push_back(N); 
    dnfnb.push_back(a);
    dnf.push_back(b);
    ratio.push_back(ratios);
    eratio.push_back(eratios);
    N++;
  }

  TGraph *gr1 = new TGraph(N,&(X[0]),&(ratio[0]));
  TGraph *gr2 = new TGraph(N,&(X[0]),&(eratio[0]));

  for(int i=0 ; i<N ; i++){
    HA->Fill(dnfnb[i]);
    //HB->Fill(dnf[i]);
  }


  gr1->SetMarkerStyle(5);
  gr1->SetMarkerSize(2);
  gr1->Draw("P same");
  //gr2->SetMarkerStyle(3);
  //gr2->SetMarkerSize(2);
  //gr2->Draw("P same");
  
  //HA->Divide(HB);
  //HA->Draw("C");
  //HA->SetMarkerStyle(20);
  //HB->Draw("H1 same");

  //leg = new TLegend(0.8,0.8,0.9,0.9);
  //leg->SetFillColor(0);
  leg.AddEntry(gr1,"ratio", "p");
  //leg.AddEntry(gr2,"eratio","p");
  //leg->Draw();
}
