#include "Histogram.h"

Histogram::Histogram(std::string metafilepath, std::string filepath){
  my_metafilepath = metafilepath;
  my_filepath     = filepath;
  ReadIn();
}

void Histogram::ReadIn()
{
  int     i = 0;
  int     lineCount,Nevents,event,C,NPOMS,NPOMH;
  double  p,pT,rap,eta;

  std::string line;
  std::ifstream meta(my_metafilepath.c_str());
  std::getline(meta,line);
  std::istringstream ss(line);
  ss >> lineCount >> Nevents;
  my_lineCount = lineCount;
  my_Nevents   = Nevents;

  Eigen::VectorXi events(lineCount);
  Eigen::VectorXi Cs    (lineCount);
  Eigen::VectorXi NPOMSs(lineCount);
  Eigen::VectorXi NPOMHs(lineCount);
  Eigen::VectorXd ps    (lineCount);
  Eigen::VectorXd pTs   (lineCount);
  Eigen::VectorXd raps  (lineCount);
  Eigen::VectorXd etas  (lineCount);
  std::ifstream f(my_filepath.c_str());

  while (std::getline(f, line)){
    std::istringstream ss(line);
    ss >> event >> C >> p >> pT >> rap >> eta >> NPOMS >> NPOMH;
    events(i) = event;
    Cs    (i) = C;
    ps    (i) = p;
    pTs   (i) = pT;
    raps  (i) = rap;
    etas  (i) = eta;
    NPOMSs(i) = NPOMS;
    NPOMHs(i) = NPOMH;
    i++;
  }
  my_events = events;
  my_Cs     = Cs;
  my_NPOMSs = NPOMSs;
  my_NPOMHs = NPOMHs;
  my_ps     = ps;
  my_pTs    = pTs;
  my_raps   = raps;
  my_etas   = etas;
}


void Histogram::EtaMax ()
{
  const int Neta = my_etas.size();
  for (int i = 0 ; i < Neta ; i++)
    if (my_etas[i] > my_maxEta){
      my_maxEta = my_etas(i);
    }
}

void Histogram::nfMax ()
{
  const int Nf = my_NPOM.rows();
  for (int i=0 ; i<Nf ; i++){
    if (my_NPOM(i,0) > my_maxNF){
      my_maxNF = my_NPOM(i,0);
    }
  }
}

void Histogram::Rapidity(int Nbins)
{
  //std::vector<int> bin (Nbins);
  Eigen::VectorXi bin (Nbins);
  const double rap_max  =  10;
  const double rap_min  = -10;
  double bin_i    = rap_min; 
  const double bin_size = (rap_max - rap_min)/Nbins;
  int i = 0; 
  while (bin_i < rap_max){
    for (int j = 0 ; j<my_lineCount ;j++){
      if (my_raps(j) > bin_i and my_raps(j) < bin_i+bin_size)
        bin(i) += 1;
    }
    i++;
    bin_i += bin_size;
  }
  std::string name = "raps.out";
//  Write(bin,name);
}

void Histogram::NPOM_count()
{
  int old_event = -1;
  int i = -1;
  Eigen::MatrixXd NPOM (my_Nevents,4);
  for(int j = 0; j < my_lineCount ; j++){
    if (old_event != my_events(j)){
      old_event = my_events(j);
      i++;
      NPOM(i,2) = my_NPOMSs(j);
      NPOM(i,3) = my_NPOMHs(j);
    }
    if (my_etas(j)>0)
      NPOM(i,0) += 1;
    else
      NPOM(i,1) += 1;
  }
  my_NPOM = NPOM;
  std::string name = "NPOM.out";
  WriteMat(name, NPOM);
}

void Histogram::NPOMS_NPOMH (int NPS_max, int NPH_max, int Nbins)
{
  my_maxEta = 15;
  double eta_step = my_maxEta/Nbins;
  for (int l=0 ; l<NPS_max ; l++){
    cout << l << endl;
    int old_event = -1;
    int j = -1;
    double eta = 0;
    Eigen::MatrixXd bins (NPH_max,Nbins);
    bins.setZero();
    for (int i = 0 ; i < my_lineCount ; i++){
      if (old_event != my_events(i)){
        old_event = my_events(i);
        j++;
      }
      eta = 0;
      for (int k=0 ; k < Nbins ; k++){
        if (my_etas(i) > eta and my_etas(i) < eta+eta_step){
          if (my_NPOM(j,2)==l and my_NPOM(j,3)<NPH_max)
           bins(my_NPOM(j,3),k) += 1;
        }
        eta += eta_step;
      }
    }
    std::string filename = std::to_string(l)+".txt";
    WriteMat(filename, bins);
  }
}

void Histogram::NFNB (int Nbins)
{
  Eigen::VectorXd nf_bin(Nbins);
  Eigen::VectorXd nfH_bin(Nbins);
  Eigen::VectorXd nbdnfH_bin(Nbins);
  Eigen::VectorXd nbdnf_bin(Nbins);
  Eigen::VectorXd bins(Nbins);
  Eigen::VectorXd bin_edges(Nbins);
  const int Nevents = my_NPOM.rows();
  const double bin_step = (double)my_maxNF/(double)Nbins;
  cout << my_maxNF << "/"<<Nbins << " = " <<bin_step << endl;
  double bin_i = 0;
  std::string path = "../../counted/";
  for (int S=0 ; S<9 ; S++){
    nf_bin.setZero();
    nbdnf_bin.setZero();
    for (int H=0 ; H<9 ; H+=2){
      bin_i = 0;
      nbdnfH_bin.setZero();
      nfH_bin.setZero();
      bins.setZero();
      for (int i=0 ; i<Nbins ; i++){ 
        for (int j=0 ; j<Nevents ;j++){
          if (my_NPOM(j,0) > bin_i and my_NPOM(j,0) < bin_i+bin_step){
            if (my_NPOM(j,2)==S and my_NPOM(j,3)==H){
              nbdnfH_bin(i) += my_NPOM(j,0)*my_NPOM(j,1);
              nfH_bin(i)    += my_NPOM(j,0);
              nf_bin(i)     += my_NPOM(j,0);
              nbdnf_bin(i)  += my_NPOM(j,0)*my_NPOM(j,1);
            }
          }
        }
        bin_edges(i) = bin_i+bin_step/2;
        bin_i += bin_step;
        if(nfH_bin(i)== 0)
          bins(i) = 0;
        else
          bins(i) = nbdnfH_bin(i)/nfH_bin(i);
      }
      //std::string name1 = path+"S"+std::to_string(S)
      //                        +"/H"+std::to_string(H)+"dnf.out";
      //std::string name2 = path+"S"+std::to_string(S)
      //                        +"/H"+std::to_string(H)+"nbdnf.out";
      std::string name3 = path+"S"+std::to_string(S)
                              +"/H"+std::to_string(H)+"nbdnf_dnf.out";
      //WriteVec(name1,nfH_bin);
      //WriteVec(name2,nfnb_bin);
      WriteVec(name3,bins);
    }
    bins.setZero();
    for(int i=0 ; i<Nbins ; i++)
      bins(i) = nbdnf_bin(i)/nf_bin(i);
    std::string name4 = path+"S"+std::to_string(S)+"/tot_sum.out";
    WriteVec(name4,bins); 
    std::string name5 = path+"S"+std::to_string(S)+"/bin_edges.out";
    WriteVec(name5,bin_edges);
  }
}

void Histogram::WriteVec(std::string filename, Eigen::VectorXd& bins)
{
  std::ofstream file(filename);
  if (file.is_open())
  {
    file << bins << '\n';
  }
}

void Histogram::WriteMat(std::string filename, Eigen::MatrixXd& bins)
{
  std::ofstream file(filename);
  if (file.is_open())
  {
    file << bins << '\n';
  }
}

