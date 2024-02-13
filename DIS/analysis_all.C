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
using Pvec3D = ROOT::Math::XYZVector;
using Pvec4D = ROOT::Math::PxPyPzMVector;
//auto p_hadron = [](double px, double py, double pz) {
//  return Pvec4D{px , py , pz , M_hadron};
//};
auto hadron_momentum = [](double px,double py,double pz){
  TVector3 v(px,py,pz);
  return v;
};


void analysis_all(){

    
    ROOT::RDataFrame d_CT14NLO("T","DIS/gen_H2_11GeV_1.root");
    ROOT::RDataFrame d_CT14LO("T","DIS/gen_H2_11GeV_2.root");
    ROOT::RDataFrame d_CT18LO("T","DIS/gen_H2_11GeV_3.root");
    ROOT::RDataFrame d_CT18NLO("T","DIS/gen_H2_11GeV_4.root");
    ROOT::RDataFrame d_CJ15NLO("T","DIS/gen_H2_11GeV_5.root");
    //ROOT::RDataFrame d_CJ15LO("T","DIS/gen_H2_11GeV_6.root");
    std::cout<<"raw CT14NLO counts"<<*d_CT14NLO.Count()<<std::endl;
    std::cout<<"raw CT14LO counts"<<*d_CT14LO.Count()<<std::endl;
    std::cout<<"raw CT18NLO counts"<<*d_CT18NLO.Count()<<std::endl;
    std::cout<<"raw CT18LO counts"<<*d_CT18LO.Count()<<std::endl;
    std::cout<<"raw CJ15NLO counts"<<*d_CJ15NLO.Count()<<std::endl;
    //std::cout<<"raw CJ15LO counts"<<*d_CJ15LO.Count()<<std::endl;
    auto d_cut_CT14NLO = d_CT14NLO
    //.Define("theta_deg",,{"theta"})
    .Define("hadron_P",hadron_momentum,{"px","py","pz"})
    .Define("momentum","sqrt(hadron_P.Dot(hadron_P))")
    .Filter("W>2")
    .Filter("Q2>1")
    ;
    auto d_cut_CT14LO = d_CT14LO
    //.Define("theta_deg",,{"theta"})
    .Define("hadron_P",hadron_momentum,{"px","py","pz"})
    .Define("momentum","sqrt(hadron_P.Dot(hadron_P))")
    .Filter("W>2")
    .Filter("Q2>1")
    ;
    auto d_cut_CT18NLO = d_CT18NLO
    //.Define("theta_deg",,{"theta"})
    .Define("hadron_P",hadron_momentum,{"px","py","pz"})
    .Define("momentum","sqrt(hadron_P.Dot(hadron_P))")
    .Filter("W>2")
    .Filter("Q2>1")
    ;
    auto d_cut_CT18LO = d_CT18LO
    //.Define("theta_deg",,{"theta"})
    .Define("hadron_P",hadron_momentum,{"px","py","pz"})
    .Define("momentum","sqrt(hadron_P.Dot(hadron_P))")
    .Filter("W>2")
    .Filter("Q2>1")
    ;
    auto d_cut_CJ15NLO = d_CJ15NLO
    //.Define("theta_deg",,{"theta"})
    .Define("hadron_P",hadron_momentum,{"px","py","pz"})
    .Define("momentum","sqrt(hadron_P.Dot(hadron_P))")
    .Filter("W>2")
    .Filter("Q2>1")
    ;
    //std::cout<<"counts after cuts"<<*d_cut.Count()<<std::endl;
    //TH1* h_theta = new TH1("","",100,0,30);
    auto h_theta_CT14NLO = d_cut_CT14NLO.Histo1D({"CT14NLO","",100,0,30},"theta","rate");
    auto h_theta_CT14LO = d_cut_CT14LO.Histo1D({"CT14LO","",100,0,30},"theta","rate");
    auto h_theta_CT18NLO = d_cut_CT18NLO.Histo1D({"CT18NLO","",100,0,30},"theta","rate");
    auto h_theta_CT18LO = d_cut_CT18LO.Histo1D({"CT18LO","",100,0,30},"theta","rate");
    auto h_theta_CJ15NLO = d_cut_CJ15NLO.Histo1D({"CJ15NLO","",100,0,30},"theta","rate");
    //std::cout<<" counts from histo" <<h_theta->Integral()<<std::endl;
    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);
    TCanvas *c_theta = new TCanvas();
    h_theta_CT14NLO->GetXaxis()->SetRange(35,75);
    h_theta_CT14NLO->GetXaxis()->SetTitle("#theta (deg)");
    h_theta_CT14NLO->GetYaxis()->SetTitle("rate(Hz)");
    h_theta_CT14NLO->SetLineColor(1);
    h_theta_CT14LO->SetLineColor(4);
    h_theta_CT18LO->SetLineColor(2);
    h_theta_CT18NLO->SetLineColor(6);
    h_theta_CJ15NLO->SetLineColor(9);
    h_theta_CT14NLO->SetTitle("CT14NLO");
    h_theta_CT14LO->SetTitle("CT14LO");
    h_theta_CT18NLO->SetTitle("CT18NLO");
    h_theta_CT18LO->SetTitle("CT18LO");
    h_theta_CJ15NLO->SetTitle("CJ15NLO");
    h_theta_CT14NLO->DrawCopy("hist");
    h_theta_CT14LO->DrawCopy("hist same");
    h_theta_CT18NLO->DrawCopy("hist same");
    h_theta_CT18LO->DrawCopy("hist same");
    h_theta_CJ15NLO->DrawCopy("hist same");
    c_theta->BuildLegend(0.65,0.65,0.9,0.9);
    c_theta->SaveAs("theta_compare.pdf");
    
    auto h_momentum_CT14NLO = d_cut_CT14NLO.Histo1D({"CT14NLO","",100,0,10},"momentum","rate");
    auto h_momentum_CT14LO = d_cut_CT14LO.Histo1D({"CT14LO","",100,0,10},"momentum","rate");
    auto h_momentum_CT18NLO = d_cut_CT18NLO.Histo1D({"CT18NLO","",100,0,10},"momentum","rate");
    auto h_momentum_CT18LO = d_cut_CT18LO.Histo1D({"CT18LO","",100,0,10},"momentum","rate");
    auto h_momentum_CJ15NLO = d_cut_CJ15NLO.Histo1D({"CJ15NLO","",100,0,10},"momentum","rate");
    //std::cout<<" counts from histo" <<h_momentum->Integral()<<std::endl;
    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);
    TCanvas *c_momentum = new TCanvas();
    //h_momentum_CT14NLO->GetXaxis()->SetRange(35,75);
    h_momentum_CT14NLO->GetXaxis()->SetTitle("#momentum (deg)");
    h_momentum_CT14NLO->GetYaxis()->SetTitle("rate(Hz)");
    h_momentum_CT14NLO->SetLineColor(1);
    h_momentum_CT14LO->SetLineColor(4);
    h_momentum_CT18LO->SetLineColor(2);
    h_momentum_CT18NLO->SetLineColor(6);
    h_momentum_CJ15NLO->SetLineColor(9);
    h_momentum_CT14NLO->SetTitle("CT14NLO");
    h_momentum_CT14LO->SetTitle("CT14LO");
    h_momentum_CT18NLO->SetTitle("CT18NLO");
    h_momentum_CT18LO->SetTitle("CT18LO");
    h_momentum_CJ15NLO->SetTitle("CJ15NLO");
    h_momentum_CT14NLO->DrawCopy("hist");
    h_momentum_CT14LO->DrawCopy("hist same");
    h_momentum_CT18NLO->DrawCopy("hist same");
    h_momentum_CT18LO->DrawCopy("hist same");
    h_momentum_CJ15NLO->DrawCopy("hist same");
    c_momentum->BuildLegend(0.15,0.6,0.35,0.9);
    c_momentum->SaveAs("momentum_compare.pdf");
}
