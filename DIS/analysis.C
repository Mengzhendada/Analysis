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

//This script checks the output from the DIS generator and plot the rates vs theta and momentum 2024Jan Shuo


using Pvec3D = ROOT::Math::XYZVector;
using Pvec4D = ROOT::Math::PxPyPzMVector;
//auto p_hadron = [](double px, double py, double pz) {
//  return Pvec4D{px , py , pz , M_hadron};
//};
auto hadron_momentum = [](double px,double py,double pz){
  TVector3 v(px,py,pz);
  return v;
};

bool root_file_exists(std::string rootfile) {
  bool found_good_file = false;
  if (!gSystem->AccessPathName(rootfile.c_str())) {
    TFile file(rootfile.c_str());
    if (file.IsZombie()) {
      std::cout << rootfile << " is a zombie.\n";
      std::cout
      << " Did your replay finish?  Check that the it is done before running this script.\n";
      return false;
      // return;
    } else {
      std::cout << " using : " << rootfile << "\n";
      return true;
    }
  }
  return false;
}


void analysis(std::string file_name="file_name"){
  std::cout<<" Enter file name "<<std::endl;
  std::cin>>file_name;

  if(root_file_exists(file_name.c_str())){

    ROOT::RDataFrame d_raw("T",file_name);
    std::cout<<"raw counts"<<*d_raw.Count()<<std::endl;
    auto d_nocut = d_raw
    //.Define("theta_deg",,{"theta"})
    .Define("hadron_P",hadron_momentum,{"px","py","pz"})
    .Define("momentum","sqrt(hadron_P.Dot(hadron_P))")
    ;
    auto d_cut = d_nocut
    .Filter("W>2")
    .Filter("Q2>1")
    ;
    std::cout<<"counts after cuts"<<*d_cut.Count()<<std::endl;
    //TH1* h_theta = new TH1("","",100,0,30);
    auto h_theta = d_cut.Histo1D({"","",100,0,30},"theta","rate");
    std::cout<<" counts from histo" <<h_theta->Integral()<<std::endl;
    gStyle->SetOptStat(0);
    TCanvas *c_theta = new TCanvas();
    //h_theta->GetXaxis()->SetRangeUser(15,25);
    h_theta->GetXaxis()->SetTitle("#theta (deg)");
    h_theta->GetYaxis()->SetTitle("rate(Hz)");
    h_theta->DrawCopy("hist");
    c_theta->SaveAs("theta.pdf");
    auto h_momentum = d_cut.Histo1D({"","",100,0,10},"momentum","rate");
    std::cout<<" counts from histo" <<h_momentum->Integral()<<std::endl;
    gStyle->SetOptStat(0);
    TCanvas *c_momentum = new TCanvas();
    //h_momentum->GetXaxis()->SetRange(35,75);
    h_momentum->GetXaxis()->SetTitle("momentum (GeV/c)");
    h_momentum->GetYaxis()->SetTitle("rate(Hz)");
    h_momentum->DrawCopy("hist");
    c_momentum->SaveAs("momentum.pdf");
    auto h_theta_momentum_raw = d_nocut.Histo2D({"","",100,0,7,100,0,30},"momentum","theta","rate");
    gStyle->SetOptStat(0);
    TCanvas *c_theta_momentum_raw = new TCanvas();
    c_theta_momentum_raw->SetRightMargin(0.1);
    c_theta_momentum_raw->SetLeftMargin(0.15);
    c_theta_momentum_raw->SetBottomMargin(0.15);
    //h_theta_momentum->GetXaxis()->SetRangeUser(5,10);
    h_theta_momentum_raw->GetYaxis()->SetRangeUser(5,25);
    h_theta_momentum_raw->GetXaxis()->SetTitleSize(0.05);
    h_theta_momentum_raw->GetYaxis()->SetTitleSize(0.05);
    h_theta_momentum_raw->GetXaxis()->SetTitle("momentum (GeV/c)");
    h_theta_momentum_raw->GetYaxis()->SetTitle("#theta (deg)");
    h_theta_momentum_raw->DrawCopy("colz");
    c_theta_momentum_raw->SaveAs("theta_momentum_raw.pdf");

    auto h_theta_momentum = d_cut.Histo2D({"","",100,0,7,100,0,30},"momentum","theta","rate");
    gStyle->SetOptStat(0);
    TCanvas *c_theta_momentum = new TCanvas();
    c_theta_momentum->SetRightMargin(0.1);
    c_theta_momentum->SetLeftMargin(0.15);
    c_theta_momentum->SetBottomMargin(0.15);
    //h_theta_momentum->GetXaxis()->SetRangeUser(5,10);
    //h_theta_momentum->GetYaxis()->SetRangeUser(14,25);
    h_theta_momentum->GetXaxis()->SetTitleSize(0.05);
    h_theta_momentum->GetYaxis()->SetTitleSize(0.05);
    h_theta_momentum->GetXaxis()->SetTitle("momentum (GeV/c)");
    h_theta_momentum->GetYaxis()->SetTitle("#theta (deg)");
    h_theta_momentum->DrawCopy("colz");
    c_theta_momentum->SaveAs("theta_momentum.pdf");
    
    auto h_theta_momentum_langle = d_cut.Histo2D({"","",100,0,7,100,0,30},"momentum","theta","rate");
    gStyle->SetOptStat(0);
    TCanvas *c_theta_momentum_langle = new TCanvas();
    c_theta_momentum_langle->SetRightMargin(0.1);
    c_theta_momentum_langle->SetLeftMargin(0.15);
    c_theta_momentum_langle->SetBottomMargin(0.15);
    //h_theta_momentum_langle->GetXaxis()->SetRangeUser(5,10);
    h_theta_momentum_langle->GetYaxis()->SetRangeUser(14,25);
    h_theta_momentum_langle->GetXaxis()->SetTitleSize(0.05);
    h_theta_momentum_langle->GetYaxis()->SetTitleSize(0.05);
    h_theta_momentum_langle->GetXaxis()->SetTitle("momentum (GeV/c)");
    h_theta_momentum_langle->GetYaxis()->SetTitle("#theta (deg)");
    h_theta_momentum_langle->DrawCopy("colz");
    c_theta_momentum_langle->SaveAs("theta_momentum.pdf");
    
    TCanvas *c_Q2_x = new TCanvas();
    auto h_Q2_x=d_cut.Histo2D({"","",100,0,1,100,0,10},"x","Q2","rate");
    //h_Q2_x->GetXaxis()->SetRangeUser(5,10);
    h_Q2_x->GetXaxis()->SetTitle("x");
    h_Q2_x->GetYaxis()->SetTitle("Q^2 (GeV^2)");
    h_Q2_x->DrawCopy("colz");
    c_Q2_x->SaveAs("Q2_x.pdf");

  }
}
