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
    auto d_cut_CJ15NLO = d_CJ15NLO
    //.Define("theta_deg",,{"theta"})
    .Define("particle_4Vector",p_particle,{"px","py","pz"})
    .Define("acc",GetAcceptance_e,{"particle_4Vector"})
    .Define("acc_rate","acc*rate")
    .Define("particle_P",particle_momentum,{"px","py","pz"})
    .Define("momentum","sqrt(particle_P.Dot(particle_P))")
    .Filter("W>2")
    .Filter("Q2>1")
    ;
    //std::cout<<"counts after cuts"<<*d_cut.Count()<<std::endl;
    TH1D* h_theta_ave = new TH1D("","",100,0,30);
    auto h_theta_CT14NLO = d_cut_CT14NLO.Histo1D({"CT14NLO","",100,0,30},"theta","acc_rate");
    auto h_theta_CT14LO = d_cut_CT14LO.Histo1D({"CT14LO","",100,0,30},"theta","acc_rate");
    auto h_theta_CT18NLO = d_cut_CT18NLO.Histo1D({"CT18NLO","",100,0,30},"theta","acc_rate");
    auto h_theta_CT18LO = d_cut_CT18LO.Histo1D({"CT18LO","",100,0,30},"theta","acc_rate");
    auto h_theta_CJ15NLO = d_cut_CJ15NLO.Histo1D({"CJ15NLO","",100,0,30},"theta","acc_rate");
    
    //std::cout<<" counts from histo" <<h_theta->Integral()<<std::endl;
    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);
    gStyle->SetTitleFontSize(.1);
    gStyle->SetLabelSize(.07,"XY");
    
    TCanvas *c_theta = new TCanvas();
    c_theta->SetRightMargin(0.1);
    c_theta->SetLeftMargin(0.15);
    c_theta->SetBottomMargin(0.15);
    h_theta_CT14NLO->GetXaxis()->SetRange(35,75);
    h_theta_CT14NLO->GetXaxis()->SetTitle("#theta (deg)");
    h_theta_CT14NLO->GetYaxis()->SetTitle("rate(Hz)");
    h_theta_CT14NLO->GetXaxis()->SetTitleSize(0.05);
    h_theta_CT14NLO->GetYaxis()->SetTitleSize(0.05);
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
    //h_theta_ave->SetTitle("average");
    //h_theta_ave->DrawCopy("hist same");
    c_theta->BuildLegend(0.65,0.65,0.9,0.9);
    c_theta->SaveAs("theta_compare.pdf");
    
    h_theta_ave->Add(h_theta_CT14NLO.GetPtr());
    h_theta_ave->Add(h_theta_CT14LO.GetPtr());
    h_theta_ave->Add(h_theta_CT18LO.GetPtr());
    h_theta_ave->Add(h_theta_CT18NLO.GetPtr());
    h_theta_ave->Add(h_theta_CJ15NLO.GetPtr());
    h_theta_ave->Scale(1/5.0);
    TCanvas *c_ratio = new TCanvas();
    TH1D *h_theta_CT14NLO_ratio=(TH1D*)h_theta_CT14NLO->Clone("h_theta_CT14NLO_ratio");
    h_theta_CT14NLO_ratio->Divide(h_theta_CT14NLO.GetPtr(),h_theta_ave);
    TH1D *h_theta_CT14LO_ratio=(TH1D*)h_theta_CT14LO->Clone("h_theta_CT14LO_ratio");
    h_theta_CT14LO_ratio->Divide(h_theta_CT14LO.GetPtr(),h_theta_ave);
    TH1D *h_theta_CT18NLO_ratio=(TH1D*)h_theta_CT18NLO->Clone("h_theta_CT18NLO_ratio");
    h_theta_CT18NLO_ratio->Divide(h_theta_CT18NLO.GetPtr(),h_theta_ave);
    TH1D *h_theta_CT18LO_ratio=(TH1D*)h_theta_CT18LO->Clone("h_theta_CT18LO_ratio");
    h_theta_CT18LO_ratio->Divide(h_theta_CT18LO.GetPtr(),h_theta_ave);
    TH1D *h_theta_CJ15NLO_ratio=(TH1D*)h_theta_CJ15NLO->Clone("h_theta_CJ15NLO_ratio");
    h_theta_CJ15NLO_ratio->Divide(h_theta_CJ15NLO.GetPtr(),h_theta_ave);
    
    h_theta_CT14NLO_ratio->SetTitle("CT14NLO");
    h_theta_CT14LO_ratio->SetTitle("CT14LO");
    h_theta_CT18NLO_ratio->SetTitle("CT18NLO");
    h_theta_CT18LO_ratio->SetTitle("CT18LO");
    h_theta_CJ15NLO_ratio->SetTitle("CJ15NLO");
    h_theta_CT14NLO_ratio->DrawCopy("L");
    h_theta_CT14LO_ratio->DrawCopy("L same");
    h_theta_CT18NLO_ratio->DrawCopy("L same");
    h_theta_CT18LO_ratio->DrawCopy("L same");
    h_theta_CJ15NLO_ratio->DrawCopy("L same");
    c_ratio->BuildLegend();
    c_ratio->SaveAs("theta_ratio.pdf");

    auto h_momentum_CT14NLO = d_cut_CT14NLO.Histo1D({"CT14NLO","",100,0,10},"momentum","acc_rate");
    auto h_momentum_CT14LO = d_cut_CT14LO.Histo1D({"CT14LO","",100,0,10},"momentum","acc_rate");
    auto h_momentum_CT18NLO = d_cut_CT18NLO.Histo1D({"CT18NLO","",100,0,10},"momentum","acc_rate");
    auto h_momentum_CT18LO = d_cut_CT18LO.Histo1D({"CT18LO","",100,0,10},"momentum","acc_rate");
    auto h_momentum_CJ15NLO = d_cut_CJ15NLO.Histo1D({"CJ15NLO","",100,0,10},"momentum","acc_rate");
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
