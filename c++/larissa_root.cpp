#include <fstream>
#include <vector>

void larissa_root()
{
  const Int_t Nbins = 500;
  const Int_t start = 0;
  const Int_t stop  = 500;

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
    HB->Fill(dnf[i]);
  }

  gr2->SetMarkerStyle(20);
  gr1->Draw();
  gr2->Draw("P same");
  
  //HA->Divide(HB);
  //HA->Draw("C");
  //HA->SetMarkerStyle(20);
  //HB->Draw("H1 same");

  leg = new TLegend(0.8,0.8,0.9,0.9);
  leg->SetFillColor(0);
  leg->AddEntry(gr1,"ratio","l");
  leg->AddEntry(gr2,"eratio","p");
  leg->Draw();
}
