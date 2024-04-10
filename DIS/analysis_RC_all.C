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
const double M_particle = 0.00051;
auto p_particle = [](double px, double py, double pz) {
  //return Pvec4D{px , py , pz , M_particle};
  double E = px*px+py*py+pz*pz+M_particle*M_particle;
  return TLorentzVector{px,py,pz,E};
};
auto particle_momentum = [](double px,double py,double pz){
  TVector3 v(px,py,pz);
  return v;
};
TFile * file_e = new TFile("Acceptance/acceptance_solid_SIDIS_He3_electron_201701_1e7_output_final.root", "r");
TH2F * acc_FA_e = (TH2F *) file_e->Get("acceptance_ThetaP_forwardangle");
TH2F * acc_LA_e = (TH2F *) file_e->Get("acceptance_ThetaP_largeangle");
double thetamin = 8.0;
double GetAcceptance_e(const TLorentzVector p){//Get electron acceptance
  double theta = p.Theta() / M_PI * 180.0;
  if (theta < thetamin || theta > 30.0) return 0; //theta > 180
  const char * detector = "all";
  double mom = p.P();
  double acc = 0;
  if (strcmp(detector, "FA") == 0 || strcmp(detector, "all") == 0)
    acc += acc_FA_e->GetBinContent(acc_FA_e->GetXaxis()->FindBin(theta), acc_FA_e->GetYaxis()->FindBin(mom));
  if ( (strcmp(detector, "LA") == 0 || strcmp(detector, "all") == 0))
  //if (mom > 3.5 && (strcmp(detector, "LA") == 0 || strcmp(detector, "all") == 0))
    acc += acc_LA_e->GetBinContent(acc_LA_e->GetXaxis()->FindBin(theta), acc_LA_e->GetYaxis()->FindBin(mom));
  //if (theta > thetamin && theta < 8.0 && mom > 2.0) return 0.5; //this line should be deleted
  return acc;
}


void analysis_RC_all(){

    
    ROOT::RDataFrame d_CT14NLO("T","DIS/gen_H2_11GeV_1.root");
    ROOT::RDataFrame d_CT14LO("T","DIS/gen_H2_11GeV_2.root");
    ROOT::RDataFrame d_CT18LO("T","DIS/gen_H2_11GeV_3.root");
    ROOT::RDataFrame d_CT18NLO("T","DIS/gen_H2_11GeV_4.root");
    ROOT::RDataFrame d_CJ15NLO("T","DIS/gen_H2_11GeV_5.root");
    auto d_cut_CT14NLO = d_CT14NLO
    //.Define("theta_deg",,{"theta"})
    .Define("particle_4Vector",p_particle,{"px","py","pz"})
    .Define("acc",GetAcceptance_e,{"particle_4Vector"})
    .Define("acc_rate","acc*rate")
    .Define("particle_P",particle_momentum,{"px","py","pz"})
    .Define("momentum","sqrt(particle_P.Dot(particle_P))")
    .Filter("W>2")
    .Filter("Q2>1")
    ;
    auto d_cut_CT14LO = d_CT14LO
    //.Define("theta_deg",,{"theta"})
    .Define("particle_4Vector",p_particle,{"px","py","pz"})
    .Define("acc",GetAcceptance_e,{"particle_4Vector"})
    .Define("acc_rate","acc*rate")
    .Define("particle_P",particle_momentum,{"px","py","pz"})
    .Define("momentum","sqrt(particle_P.Dot(particle_P))")
    .Filter("W>2")
    .Filter("Q2>1")
    ;
    auto d_cut_CT18NLO = d_CT18NLO
    //.Define("theta_deg",,{"theta"})
    .Define("particle_4Vector",p_particle,{"px","py","pz"})
    .Define("acc",GetAcceptance_e,{"particle_4Vector"})
    .Define("acc_rate","acc*rate")
    .Define("particle_P",particle_momentum,{"px","py","pz"})
    .Define("momentum","sqrt(particle_P.Dot(particle_P))")
    .Filter("W>2")
    .Filter("Q2>1")
    ;
    auto d_cut_CT18LO = d_CT18LO
    //.Define("theta_deg",,{"theta"})
    .Define("particle_4Vector",p_particle,{"px","py","pz"})
    .Define("acc",GetAcceptance_e,{"particle_4Vector"})
    .Define("acc_rate","acc*rate")
    .Define("particle_P",particle_momentum,{"px","py","pz"})
    .Define("momentum","sqrt(particle_P.Dot(particle_P))")
    .Filter("W>2")
    .Filter("Q2>1")
    ;
    //For RC case
    ROOT::RDataFrame d_RC_CT14NLO("T","DIS/gen_H2_11GeV_1RC.root");
    ROOT::RDataFrame d_RC_CT14LO("T","DIS/gen_H2_11GeV_2RC.root");
    //ROOT::RDataFrame d_RC_CT18LO("T","DIS/gen_H2_11GeV_3RC.root");
    //ROOT::RDataFrame d_RC_CT18NLO("T","DIS/gen_H2_11GeV_4RC.root");
    //ROOT::RDataFrame d_RC_CJ15NLO("T","DIS/gen_H2_11GeV_5RC.root");
    auto d_RC_cut_CT14NLO = d_RC_CT14NLO
    //.Define("theta_deg",,{"theta"})
    .Define("particle_4Vector",p_particle,{"px","py","pz"})
    .Define("acc",GetAcceptance_e,{"particle_4Vector"})
    .Define("acc_rate","acc*rate")
    .Define("particle_P",particle_momentum,{"px","py","pz"})
    .Define("momentum","sqrt(particle_P.Dot(particle_P))")
    .Filter("W>2")
    .Filter("Q2>1")
    ;
    auto d_RC_cut_CT14LO = d_RC_CT14LO
    //.Define("theta_deg",,{"theta"})
    .Define("particle_4Vector",p_particle,{"px","py","pz"})
    .Define("acc",GetAcceptance_e,{"particle_4Vector"})
    .Define("acc_rate","acc*rate")
    .Define("particle_P",particle_momentum,{"px","py","pz"})
    .Define("momentum","sqrt(particle_P.Dot(particle_P))")
    .Filter("W>2")
    .Filter("Q2>1")
    ;
    //auto d_RC_cut_CT18NLO = d_RC_CT18NLO
    ////.Define("theta_deg",,{"theta"})
    //.Define("particle_P",particle_momentum,{"px","py","pz"})
    //.Define("momentum","sqrt(particle_P.Dot(particle_P))")
    //.Filter("W>2")
    //.Filter("Q2>1")
    //;
    //auto d_RC_cut_CT18LO = d_RC_CT18LO
    ////.Define("theta_deg",,{"theta"})
    //.Define("particle_P",particle_momentum,{"px","py","pz"})
    //.Define("momentum","sqrt(particle_P.Dot(particle_P))")
    //.Filter("W>2")
    //.Filter("Q2>1")
    //;
    //std::cout<<"counts after cuts"<<*d_cut.Count()<<std::endl;
    //TH1* h_theta = new TH1("","",20,0,30);
    auto h_theta_CT14NLO = d_cut_CT14NLO.Histo1D({"","CT14NLO",20,0,30},"theta","acc_rate");
    auto h_theta_CT14LO = d_cut_CT14LO.Histo1D({"","CT14LO",20,0,30},"theta","acc_rate");
    auto h_theta_CT18NLO = d_cut_CT18NLO.Histo1D({"","CT18NLO",20,0,30},"theta","acc_rate");
    auto h_theta_CT18LO = d_cut_CT18LO.Histo1D({"","CT18LO",20,0,30},"theta","acc_rate");
    //std::cout<<" counts from histo" <<h_theta->Integral()<<std::endl;
    
    //for RC case
    auto h_RC_theta_CT14NLO = d_RC_cut_CT14NLO.Histo1D({"","CT14NLO",20,0,30},"theta","acc_rate");
    auto h_RC_theta_CT14LO = d_RC_cut_CT14LO.Histo1D({"","CT14LO",20,0,30},"theta","acc_rate");
    //auto h_RC_theta_CT18NLO = d_RC_cut_CT18NLO.Histo1D({"","CT18NLO",20,0,30},"theta","acc_rate");
    //auto h_RC_theta_CT18LO = d_RC_cut_CT18LO.Histo1D({"","CT18LO",20,0,30},"theta","acc_rate");
    
    //take ratio
    TH1 *h_theta_CT14NLO_ratio = (TH1*)h_theta_CT14NLO->Clone("h_theta_CT14NLO_ratio");
    h_theta_CT14NLO_ratio->Divide(h_theta_CT14NLO.GetPtr(),h_RC_theta_CT14NLO.GetPtr());
    TH1 *h_theta_CT14LO_ratio = (TH1*)h_theta_CT14LO->Clone("h_theta_CT14LO_ratio");
    h_theta_CT14LO_ratio->Divide(h_theta_CT14LO.GetPtr(),h_RC_theta_CT14LO.GetPtr());

    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gStyle->SetTitleFontSize(.1);
    gStyle->SetLabelSize(.07,"XY");
    TCanvas *c_theta = new TCanvas();
    c_theta->SetRightMargin(0.1);
    c_theta->SetLeftMargin(0.15);
    c_theta->SetBottomMargin(0.15);
    h_theta_CT14NLO_ratio->GetYaxis()->SetRangeUser(0.8,1.2);
    //h_theta_CT14NLO->GetXaxis()->SetRange(35,75);
    h_theta_CT14NLO_ratio->GetXaxis()->SetTitle("#theta (deg)");
    h_theta_CT14NLO_ratio->GetYaxis()->SetTitle("RC");
    h_theta_CT14NLO_ratio->GetXaxis()->SetTitleSize(0.05);
    h_theta_CT14NLO_ratio->GetYaxis()->SetTitleSize(0.05);
    h_theta_CT14NLO_ratio->DrawCopy("L");
    h_theta_CT14LO_ratio->SetLineColor(2);
    h_theta_CT14LO_ratio->DrawCopy("L same");
    //h_theta_CT18NLO->DrawCopy("L same");
    //h_theta_CT18LO->DrawCopy("L same");
    c_theta->BuildLegend(0.75,0.75,0.9,0.9);
    c_theta->SaveAs("theta_RC_compare.pdf");
    
    auto h_momentum_CT14NLO = d_cut_CT14NLO.Histo1D({"","CT14NLO",20,0,10},"momentum","rate");
    auto h_momentum_CT14LO = d_cut_CT14LO.Histo1D({"","CT14LO",20,0,10},"momentum","rate");
    auto h_momentum_CT18NLO = d_cut_CT18NLO.Histo1D({"","CT18NLO",20,0,10},"momentum","rate");
    auto h_momentum_CT18LO = d_cut_CT18LO.Histo1D({"","CT18LO",20,0,10},"momentum","rate");
    //std::cout<<" counts from histo" <<h_momentum->Integral()<<std::endl;
    //for RC case
    auto h_RC_momentum_CT14NLO = d_RC_cut_CT14NLO.Histo1D({"","CT14NLO",20,0,10},"momentum","rate");
    auto h_RC_momentum_CT14LO = d_RC_cut_CT14LO.Histo1D({"","CT14LO",20,0,10},"momentum","rate");
    //auto h_RC_momentum_CT18NLO = d_RC_cut_CT18NLO.Histo1D({"","CT18NLO",20,0,10},"momentum","rate");
    //auto h_RC_momentum_CT18LO = d_RC_cut_CT18LO.Histo1D({"","CT18LO",20,0,10},"momentum","rate");
    
    //take ratio
    TH1 *h_momentum_CT14NLO_ratio = (TH1*)h_momentum_CT14NLO->Clone("h_momentum_CT14NLO_ratio");
    h_momentum_CT14NLO_ratio->Divide(h_momentum_CT14NLO.GetPtr(),h_RC_momentum_CT14NLO.GetPtr());
    TH1 *h_momentum_CT14LO_ratio = (TH1*)h_momentum_CT14LO->Clone("h_momentum_CT14LO_ratio");
    h_momentum_CT14LO_ratio->Divide(h_momentum_CT14LO.GetPtr(),h_RC_momentum_CT14LO.GetPtr());
    
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    TCanvas *c_momentum = new TCanvas();
    //h_momentum_CT14NLO->GetXaxis()->SetRange(35,75);
    h_momentum_CT14NLO_ratio->GetXaxis()->SetTitle("momentum (GeV/c)");
    h_momentum_CT14NLO_ratio->GetYaxis()->SetTitle("RC");
    h_momentum_CT14NLO_ratio->GetYaxis()->SetRangeUser(0.8,1.2);
    h_momentum_CT14NLO_ratio->DrawCopy("L");
    h_momentum_CT14LO_ratio->SetLineColor(2);
    h_momentum_CT14LO_ratio->DrawCopy("L same");
    //h_momentum_CT18NLO->DrawCopy("L same");
    //h_momentum_CT18LO->DrawCopy("L same");
    c_momentum->BuildLegend(0.75,0.75,0.9,0.9);
    c_momentum->SaveAs("momentum_RCcompare.pdf");
}
