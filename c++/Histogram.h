#pragma once
#include "Eigen/Core"
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <iostream>
using std::cout;
using std::endl;
class Histogram {
  public:
    Histogram       (std::string metafilepath,std::string filepath);
    void ReadIn     ();
    void EtaMax     ();
    void nfMax      ();
    void Rapidity   (int Nbins);
    void NPOM_count ();
    void NFNB       (int Nbins);
    void NPOMS_NPOMH(int NPS_max, int NPH_max,int Nbins);
    void Write      (std::vector<int>& bin, std::string name);
    void WriteVec   (std::string filename, Eigen::VectorXd& bins);
    void WriteMat   (std::string filename, Eigen::MatrixXd& bins);

    double get_maxEta (){return my_maxEta;}

  protected:
    int                 my_Nevents;
    int                 my_lineCount; 
    int                 my_maxNF = 0;
    double              my_maxEta = 0;

    std::string         my_filepath;
    std::string         my_metafilepath;
    Eigen::VectorXi     my_events;
    Eigen::VectorXi     my_Cs;
    Eigen::VectorXi     my_NPOMSs;
    Eigen::VectorXi     my_NPOMHs;
    Eigen::VectorXd     my_ps;
    Eigen::VectorXd     my_pTs;
    Eigen::VectorXd     my_raps;
    Eigen::VectorXd     my_etas;
    Eigen::MatrixXd     my_NPOM;
    //std::vector<int>    my_events;
    //std::vector<int>    my_Cs;
    //std::vector<int>    my_NPOMSs;
    //std::vector<int>    my_NPOMHs;
    //std::vector<double> my_ps;
    //std::vector<double> my_pTs;
    //std::vector<double> my_raps;
    //std::vector<double> my_etas;
    //std::vector<std::vector<int>> my_NPOM;
};
