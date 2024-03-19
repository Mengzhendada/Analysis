#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <istream>
#include <vector>
#include <cmath>
#include <ios>
#include <iosfwd>
#include <iomanip>
#include <streambuf>

#include "ROOT/RDataFrame.hxx"
#include "TCanvas.h"
#include "TChain.h"
#include "TFile.h"
#include "Math/Vector3D.h"
#include "Math/Vector4D.h"
#include "Math/VectorUtil.h"
#include "ROOT/RVec.hxx"
#include "TVector3.h"

bool root_file_exists(std::string rootfile) {
  bool found_good_file = false;
  if (!gSystem->AccessPathName(rootfile.c_str())) {
    TFile file(rootfile.c_str());
    if (file.IsZombie()) {
      std::cout << rootfile << " is a zombie.\n";
      std::cout
      << " Check that the file before running this script.\n";
      return false;
      // return;
    } else {
      std::cout << " using : " << rootfile << "\n";
      return true;
    }
  }
  return false;
}
//This script reads 2d histos in Acceptance/*root and replot them with numbers
void Check_acceptance(std::string file_name="file_name"){

  //std::cout<<" Enter file name "<<std::endl;
  //std::cin>>file_name;
  //file_name = "acceptance_solid_SIDIS_He3_electron_201701_1e7_output_final";
  //file_name = "acceptance_solid_SIDIS_He3_pim_201701_1e7_output_final";
  file_name = "acceptance_solid_SIDIS_He3_pip_201701_1e7_output_final";
  //file_name = "acceptance_solid_SIDIS_He3_km_201701_1e7_output_final";
  //file_name = "acceptance_solid_SIDIS_He3_kp_201701_1e7_output_final";
  
  std::string full_file_name="Acceptance/"+file_name+".root";
  if(root_file_exists(full_file_name.c_str())){
    TFile *file = new TFile(full_file_name.c_str()); 
    TH2F* h_forward = (TH2F *) file->Get("acceptance_ThetaP_forwardangle");
    TH2F* h_large = (TH2F *) file->Get("acceptance_ThetaP_largeangle");



    //TH1F* h_forward_1D = (TH1F*) h_forward->ProjectionX();
    //TH1F* h_large_1D = (TH1F*) h_large->ProjectionX();

    TCanvas *c_forward = new TCanvas("","",1200,1200);
    gStyle->SetOptTitle(0);
    c_forward->SetLeftMargin(0.15);
    c_forward->SetBottomMargin(0.2);
    h_forward->Draw("colz");
    c_forward->SaveAs("results/Acceptance_forward.pdf");
    TCanvas *c_large = new TCanvas("","",1200,1200);
    gStyle->SetOptTitle(0);
    c_large->SetLeftMargin(0.15);
    c_large->SetBottomMargin(0.2);
    h_large->Draw("colz");
    c_large->SaveAs("results/Acceptance_large.pdf");
  }
}

