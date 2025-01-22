#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <ctime>
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TFile.h"
#include "TTree.h"
#include "TDirectoryFile.h"
#include "TCanvas.h"
#include "TImage.h"
#include "TStyle.h"
#include "ILDStyle_Standalone.cc"
#include <TMultiGraph.h>
#include <string>
#include <TProfile.h>


// ----- COMMAND ---> root -l preselection.cc

void GetcategorisationTree(TFile* file, const std::string& objectname, TTree*& ttree_file, int& TrueCat, int& RecoCatBasic, int& RecoCatAdv);
void GetobservablesTree(TFile* file, const std::string& objectname, TTree*& ttree_file, double& Evis, double& pT, double& Evis_nonIso, double& pT_nonIso, double& missE, double& misspT, double& mInv);
TProfile* PlotObservable(TTree*& tree_file, const std::string& variable, const std::string& histName, int bins, double xMin, double xMax, const std::string& xTitle = "", const std::string& yTitle = "Events", int fillColor = kRed);
void PlotStackedHistograms(const std::string& variable, TProfile* h_signal, TProfile* h_3f, TProfile* h_6f, TProfile* h_5f, TProfile* h_2f, TProfile* h_4fbkg);

string fileFormat = ".pdf";
string outspec = "";

// ---- input root files ----
string background_2f = "../../ROOTFiles/PreSelection2f.root";
string background_3f = "../../ROOTFiles/PreSelection3f.root";
string background_4f = "../../ROOTFiles/PreSelection4f_bkg.root";
string background_6f = "../../ROOTFiles/PreSelection6f.root";
string ww_4f = "../../ROOTFiles/PreSelection4f_signal.root";
string background_5f = "../../ROOTFiles/PreSelection5f.root";
string outfolder = ".";

void plotting(){
  TStyle* MyStyle = new  TStyle("MyStyle", "My Style");
  ILDStyle(MyStyle);
  MyStyle->cd();
  gROOT->ForceStyle();
  gStyle->ls();
  gStyle->SetPalette(kBird);
  
  // ----------- ACCESS ROOT FILE ------------
  TFile* input_2fbackground = TFile::Open(background_2f.c_str());
  TFile* input_3fbackground = TFile::Open(background_3f.c_str());
  TFile* input_4fbackground = TFile::Open(background_4f.c_str());
  TFile* input_6fbackground = TFile::Open(background_6f.c_str());
  TFile* input_5fbackground = TFile::Open(background_5f.c_str());
  TFile* input_ww = TFile::Open(ww_4f.c_str());

  if(!input_ww || input_ww->IsZombie()){
    std::cerr << "Failed to open file" << std::endl;
    return;
  }

  // Map to store histograms
  std::map<std::string, TProfile*> histograms;

	if(1){
    // --------  TREES  ------------
    // observablesTree
    TTree *tree, *tree_4fbkg, *tree_2f, *tree_3f, *tree_6f, *tree_5f;
    // categorisationTree
    TTree *categorisationTree4f, *categorisationTree4f_bkg, *categorisationTree2f, *categorisationTree3f, *categorisationTree6f, *categorisationTree5f;

    // ------ VARIABLES ------------
    double Evis, pT, mInv;
    double Evis_nonIso, pT_nonIso;
    double missE, misspT;

    int TrueCat, RecoCatBasic, RecoCatAdv;

    // ------- Get Categorisation Integers -----------
    GetcategorisationTree(input_ww, "observablesTree", categorisationTree4f, TrueCat, RecoCatBasic, RecoCatAdv);
    GetcategorisationTree(input_2fbackground, "observablesTree", categorisationTree2f, TrueCat, RecoCatBasic, RecoCatAdv);
    GetcategorisationTree(input_3fbackground, "observablesTree", categorisationTree3f, TrueCat, RecoCatBasic, RecoCatAdv);
    GetcategorisationTree(input_4fbackground, "observablesTree", categorisationTree4f_bkg, TrueCat, RecoCatBasic, RecoCatAdv);
    GetcategorisationTree(input_6fbackground, "observablesTree", categorisationTree6f, TrueCat, RecoCatBasic, RecoCatAdv);
    GetcategorisationTree(input_5fbackground, "observablesTree", categorisationTree5f, TrueCat, RecoCatBasic, RecoCatAdv);

    // ------- Get Observables -----------
    GetobservablesTree(input_ww, "observablesTree", tree, Evis, pT, Evis_nonIso, pT_nonIso, missE, misspT, mInv);
    GetobservablesTree(input_2fbackground, "observablesTree", tree_2f, Evis, pT, Evis_nonIso, pT_nonIso, missE, misspT, mInv);
    GetobservablesTree(input_3fbackground, "observablesTree", tree_3f, Evis, pT, Evis_nonIso, pT_nonIso, missE, misspT, mInv);
    GetobservablesTree(input_4fbackground, "observablesTree", tree_4fbkg, Evis, pT, Evis_nonIso, pT_nonIso, missE, misspT, mInv);
    GetobservablesTree(input_6fbackground, "observablesTree", tree_6f, Evis, pT, Evis_nonIso, pT_nonIso, missE, misspT, mInv);
    GetobservablesTree(input_5fbackground, "observablesTree", tree_5f, Evis, pT, Evis_nonIso, pT_nonIso, missE, misspT, mInv);

    // Print values of categorisationTree
    for (int i = 0; i < categorisationTree4f->GetEntries(); i++) {
      categorisationTree4f->GetEntry(i); 
    }
    for (int i = 0; i < categorisationTree3f->GetEntries(); i++) {
      categorisationTree3f->GetEntry(i);
      //std::cout << "3f - TrueCat: " << TrueCat << ", RecoCatBasic: " << RecoCatBasic << ", RecoCatAdv: " << RecoCatAdv << std::endl;
    }
    for (int i = 0; i < categorisationTree6f->GetEntries(); i++) {
      categorisationTree6f->GetEntry(i);
      //std::cout << "6f - TrueCat: " << TrueCat << ", RecoCatBasic: " << RecoCatBasic << ", RecoCatAdv: " << RecoCatAdv << std::endl;
    }
    // test if the values are correctly retrieved
    if(TrueCat == 0){
      std::cout << "Successfully retrieved Categorisation value." << std::endl;
      // Process the tree, e.g., loop over entries
    } else {
      std::cerr << "Failed to retrieve Categorisation value." << std::endl;
    }

    if (tree && tree_3f && tree_6f && tree_2f && tree_5f && tree_4fbkg) {
      std::cout << "Successfully retrieved TTrees." << std::endl;
      // Process the tree, e.g., loop over entries
    } else {
      std::cerr << "Failed to retrieve TTree." << std::endl;
    }

    // ------- PLOTTING ------------
    // --------- 2 FERMION ------------
    histograms["h_evis_2f"] = PlotObservable(tree_2f, "Evis", "h_evis_2f", 401, -0.5, 400.5, "Visible Energy [GeV]", "Events", kYellow);
    histograms["h_pt_2f"] = PlotObservable(tree_2f, "pT", "h_pt_2f", 101, -0.5, 100.5, "pT [GeV]", "Events", kYellow);
    histograms["h_minv_2f"] = PlotObservable(tree_2f, "mInv", "h_minv_2f", 401, -0.5, 400.5, "M_{inv} [GeV]", "Events", kYellow);
    histograms["h_pt_nonIso_2f"] = PlotObservable(tree_2f, "pT_nonIso", "h_pt_nonIso_2f", 101, -0.5, 100.5, "pT [GeV]", "Events", kYellow);
    histograms["h_evis_nonIso_2f"] = PlotObservable(tree_2f, "Evis_nonIso", "h_evis_nonIso_2f", 401, -0.5, 400.5, "Visible Energy [GeV]", "Events", kYellow);
    histograms["h_misspT_2f"] = PlotObservable(tree_2f, "misspT", "h_misspT_2f", 101, -0.5, 100.5, "Missing pT [GeV]", "Events", kYellow);
    histograms["h_missE_2f"] = PlotObservable(tree_2f, "missE", "h_missE_2f", 401, -0.5, 400.5, "Missing Energy [GeV]", "Events", kYellow);

    /* TH1F* h_pt_2f = PlotObservable(tree_2f, "pT", "h_pt_2f", 101, -0.5, 100.5, "pT [GeV]", "Events", kYellow);
    TH1F* h_minv_2f = PlotObservable(tree_2f, "mInv", "h_minv_2f", 401, -0.5, 00.5, "M_{inv} [GeV]", "Events", kYellow);
    TH1F* h_pt_nonIso_2f = PlotObservable(tree_2f, "pT_nonIso", "h_pt_nonIso_2f", 101, -0.5, 100.5, "pT [GeV]", "Events", kYellow);
    TH1F* h_evis_nonIso_2f = PlotObservable(tree_2f, "Evis_nonIso", "h_evis_nonIso_2f", 401, -0.5, 00.5, "Visible Energy [GeV]", "Events", kYellow);
    TH1F* h_misspT_2f = PlotObservable(tree_2f, "misspT", "h_misspT_2f", 101, -0.5, 100.5, "Missing pT [GeV]", "Events", kYellow);
    TH1F* h_missE_2f = PlotObservable(tree_2f, "missE", "h_missE_2f", 401, -0.5, 00.5, "Missing Energy [GeV]", "Events", kYellow); */

    TProfile* h_evis_2f = histograms["h_evis_2f"];
    TProfile* h_pt_2f = histograms["h_pt_2f"];
    TProfile* h_minv_2f = histograms["h_minv_2f"];
    TProfile* h_pt_nonIso_2f = histograms["h_pt_nonIso_2f"];
    TProfile* h_evis_nonIso_2f = histograms["h_evis_nonIso_2f"];
    TProfile* h_misspT_2f = histograms["h_misspT_2f"];
    TProfile* h_missE_2f = histograms["h_missE_2f"];
    
    // --------- 3 FERMION ------------
    histograms["h_evis_bkg"] = PlotObservable(tree_3f, "Evis", "h_evis_bkg", 401, -0.5, 400.5, "Visible Energy [GeV]", "Events", kBlue);
    histograms["h_pt_bkg"] = PlotObservable(tree_3f, "pT", "h_pt_bkg", 101, -0.5, 100.5, "pT [GeV]", "Events", kBlue);
    histograms["h_minv_bkg"] = PlotObservable(tree_3f, "mInv", "h_minv_bkg", 401, -0.5, 400.5, "M_{inv} [GeV]", "Events", kBlue);
    histograms["h_pt_nonIso_bkg"] = PlotObservable(tree_3f, "pT_nonIso", "h_pt_nonIso_bkg", 101, -0.5, 100.5, "pT [GeV]", "Events", kBlue);
    histograms["h_evis_nonIso_3f"] = PlotObservable(tree_3f, "Evis_nonIso", "h_evis_nonIso_3f", 401, -0.5, 400.5, "Visible Energy [GeV]", "Events", kBlue);
    histograms["h_misspT_bkg"] = PlotObservable(tree_3f, "misspT", "h_misspT_bkg", 101, -0.5, 100.5, "Missing pT [GeV]", "Events", kBlue);
    histograms["h_missE_bkg"] = PlotObservable(tree_3f, "missE", "h_missE_bkg", 401, -0.5, 400.5, "Missing Energy [GeV]", "Events", kBlue);

    TProfile* h_evis_bkg = histograms["h_evis_bkg"];
    TProfile* h_pt_bkg = histograms["h_pt_bkg"];
    TProfile* h_minv_bkg = histograms["h_minv_bkg"];
    TProfile* h_pt_nonIso_bkg = histograms["h_pt_nonIso_bkg"];
    TProfile* h_evis_nonIso_3f = histograms["h_evis_nonIso_3f"];
    TProfile* h_misspT_bkg = histograms["h_misspT_bkg"];
    TProfile* h_missE_bkg = histograms["h_missE_bkg"];
    
    // --------- 6 FERMION ------------
    histograms["h_evis_6fbkg"] = PlotObservable(tree_6f, "Evis", "h_evis_6fbkg", 401, -0.5, 400.5, "Visible Energy [GeV]", "Events", kRed);
    histograms["h_pt_6fbkg"] = PlotObservable(tree_6f, "pT", "h_pt_6fbkg", 101, -0.5, 100.5, "pT [GeV]", "Events", kRed);
    histograms["h_minv_6fbkg"] = PlotObservable(tree_6f, "mInv", "h_minv_6fbkg", 401, -0.5, 400.5, "M_{inv} [GeV]", "Events", kRed);
    histograms["h_pt_nonIso_6fbkg"] = PlotObservable(tree_6f, "pT_nonIso", "h_pt_nonIso_6fbkg", 101, -0.5, 100.5, "pT [GeV]", "Events", kRed);
    histograms["h_evis_nonIso_6fbkg"] = PlotObservable(tree_6f, "Evis_nonIso", "h_evis_nonIso_6fbkg", 401, -0.5, 400.5, "Visible Energy [GeV]", "Events", kRed);
    histograms["h_misspT_6fbkg"] = PlotObservable(tree_6f, "misspT", "h_misspT_6fbkg", 101, -0.5, 100.5, "Missing pT [GeV]", "Events", kRed);
    histograms["h_missE_6fbkg"] = PlotObservable(tree_6f, "missE", "h_missE_6fbkg", 401, -0.5, 400.5, "Missing Energy [GeV]", "Events", kRed);

    TProfile* h_evis_6fbkg = histograms["h_evis_6fbkg"];
    TProfile* h_pt_6fbkg = histograms["h_pt_6fbkg"];
    TProfile* h_minv_6fbkg = histograms["h_minv_6fbkg"];
    TProfile* h_pt_nonIso_6fbkg = histograms["h_pt_nonIso_6fbkg"];
    TProfile* h_evis_nonIso_6fbkg = histograms["h_evis_nonIso_6fbkg"];
    TProfile* h_misspT_6fbkg = histograms["h_misspT_6fbkg"];
    TProfile* h_missE_6fbkg = histograms["h_missE_6fbkg"];

    // --------- 4 FERMION BACKGROUND ------------
    histograms["h_evis_4fbkg"] = PlotObservable(tree_4fbkg, "Evis", "h_evis_4fbkg", 401, -0.5, 400.5, "Visible Energy [GeV]", "Events", kPink+10);
    histograms["h_pt_4fbkg"] = PlotObservable(tree_4fbkg, "pT", "h_pt_4fbkg", 101, -0.5, 100.5, "pT [GeV]", "Events", kPink+10);
    histograms["h_minv_4fbkg"] = PlotObservable(tree_4fbkg, "mInv", "h_minv_4fbkg", 401, -0.5, 400.5, "M_{inv} [GeV]", "Events", kPink+10);
    histograms["h_pt_nonIso_4fbkg"] = PlotObservable(tree_4fbkg, "pT_nonIso", "h_pt_nonIso_4fbkg", 101, -0.5, 100.5, "pT [GeV]", "Events", kPink+10);
    histograms["h_evis_nonIso_4fbkg"] = PlotObservable(tree_4fbkg, "Evis_nonIso", "h_evis_nonIso_4fbkg", 401, -0.5, 400.5, "Visible Energy [GeV]", "Events", kPink+10);
    histograms["h_misspT_4fbkg"] = PlotObservable(tree_4fbkg, "misspT", "h_misspT_4fbkg", 101, -0.5, 100.5, "Missing pT [GeV]", "Events", kPink+10);
    histograms["h_missE_4fbkg"] = PlotObservable(tree_4fbkg, "missE", "h_missE_4fbkg", 401, -0.5, 400.5, "Missing Energy [GeV]", "Events", kPink+10);

    TProfile* h_evis_4fbkg = histograms["h_evis_4fbkg"];
    TProfile* h_pt_4fbkg = histograms["h_pt_4fbkg"];
    TProfile* h_minv_4fbkg = histograms["h_minv_4fbkg"];
    TProfile* h_pt_nonIso_4fbkg = histograms["h_pt_nonIso_4fbkg"];
    TProfile* h_evis_nonIso_4fbkg = histograms["h_evis_nonIso_4fbkg"];
    TProfile* h_misspT_4fbkg = histograms["h_misspT_4fbkg"];
    TProfile* h_missE_4fbkg = histograms["h_missE_4fbkg"];
    // --------- 5 FERMION ------------
    histograms["h_evis_5f"] = PlotObservable(tree_5f, "Evis", "h_evis_5f", 401, -0.5, 400.5, "Visible Energy [GeV]", "Events", kGray);
    histograms["h_pt_5f"] = PlotObservable(tree_5f, "pT", "h_pt_5f", 101, -0.5, 100.5, "pT [GeV]", "Events", kGray);
    histograms["h_minv_5f"] = PlotObservable(tree_5f, "mInv", "h_minv_5f", 401, -0.5, 400.5, "M_{inv} [GeV]", "Events", kGray);
    histograms["h_pt_nonIso_5f"] = PlotObservable(tree_5f, "pT_nonIso", "h_pt_nonIso_5f", 101, -0.5, 100.5, "pT [GeV]", "Events", kGray);
    histograms["h_evis_nonIso_5f"] = PlotObservable(tree_5f, "Evis_nonIso", "h_evis_nonIso_5f", 401, -0.5, 400.5, "Visible Energy [GeV]", "Events", kGray);
    histograms["h_misspT_5f"] = PlotObservable(tree_5f, "misspT", "h_misspT_5f", 101, -0.5, 100.5, "Missing pT [GeV]", "Events", kGray);
    histograms["h_missE_5f"] = PlotObservable(tree_5f, "missE", "h_missE_5f", 401, -0.5, 400.5, "Missing Energy [GeV]", "Events", kGray);

    TProfile* h_evis_5f = histograms["h_evis_5f"];
    TProfile* h_pt_5f = histograms["h_pt_5f"];
    TProfile* h_minv_5f = histograms["h_minv_5f"];
    TProfile* h_pt_nonIso_5f = histograms["h_pt_nonIso_5f"];
    TProfile* h_evis_nonIso_5f = histograms["h_evis_nonIso_5f"];
    TProfile* h_misspT_5f = histograms["h_misspT_5f"];
    TProfile* h_missE_5f = histograms["h_missE_5f"];
    // --------- 4 FERMION SIGNAL ------------
    histograms["h_evis"] = PlotObservable(tree, "Evis", "h_evis", 401, -0.5, 400.5, "Visible Energy [GeV]", "Events", kGreen);
    histograms["h_pt"] = PlotObservable(tree, "pT", "h_pt", 101, -0.5, 100.5, "pT [GeV]", "Events", kGreen);
    histograms["h_minv"] = PlotObservable(tree, "mInv", "h_minv", 401, -0.5, 400.5, "M_{inv} [GeV]", "Events", kGreen);
    histograms["h_pt_nonIso"] = PlotObservable(tree, "pT_nonIso", "h_pt_nonIso", 101, -0.5, 100.5, "pT [GeV]", "Events", kGreen);
    histograms["h_evis_nonIso"] = PlotObservable(tree, "Evis_nonIso", "h_evis_nonIso", 401, -0.5, 400.5, "Visible Energy [GeV]", "Events", kGreen);
    histograms["h_misspT"] = PlotObservable(tree, "misspT", "h_misspT", 101, -0.5, 100.5, "Missing pT [GeV]", "Events", kGreen);
    histograms["h_missE"] = PlotObservable(tree, "missE", "h_missE", 401, -0.5, 400.5, "Missing Energy [GeV]", "Events", kGreen);

    TProfile* h_evis = histograms["h_evis"];
    TProfile* h_pt = histograms["h_pt"];
    TProfile* h_minv = histograms["h_minv"];
    TProfile* h_pt_nonIso = histograms["h_pt_nonIso"];
    TProfile* h_evis_nonIso = histograms["h_evis_nonIso"];
    TProfile* h_misspT = histograms["h_misspT"];
    TProfile* h_missE = histograms["h_missE"];

    // ---- STACKED PLOTS ----
    PlotStackedHistograms("Evis", histograms["h_evis"], histograms["h_evis_bkg"], histograms["h_evis_6fbkg"], histograms["h_evis_5f"], histograms["h_evis_2f"], histograms["h_evis_4fbkg"]);
    PlotStackedHistograms("pT", histograms["h_pt"], histograms["h_pt_bkg"], histograms["h_pt_6fbkg"], histograms["h_pt_5f"], histograms["h_pt_2f"], histograms["h_pt_4fbkg"]);
    PlotStackedHistograms("mInv", histograms["h_minv"], histograms["h_minv_bkg"], histograms["h_minv_6fbkg"], histograms["h_minv_5f"], histograms["h_minv_2f"], histograms["h_minv_4fbkg"]);
    PlotStackedHistograms("pT_nonIso", histograms["h_pt_nonIso"], histograms["h_pt_nonIso_bkg"], histograms["h_pt_nonIso_6fbkg"], histograms["h_pt_nonIso_5f"], histograms["h_pt_nonIso_2f"], histograms["h_pt_nonIso_4fbkg"]);
    PlotStackedHistograms("misspT", histograms["h_misspT"], histograms["h_misspT_bkg"], histograms["h_misspT_6fbkg"], histograms["h_misspT_5f"], histograms["h_misspT_2f"], histograms["h_misspT_4fbkg"]);
    PlotStackedHistograms("missE", histograms["h_missE"], histograms["h_missE_bkg"], histograms["h_missE_6fbkg"], histograms["h_missE_5f"], histograms["h_missE_2f"], histograms["h_missE_4fbkg"]);
    PlotStackedHistograms("Evis_nonIso", histograms["h_evis_nonIso"], histograms["h_evis_nonIso_3f"], histograms["h_evis_nonIso_6fbkg"], histograms["h_evis_nonIso_5f"], histograms["h_evis_nonIso_2f"], histograms["h_evis_nonIso_4fbkg"]);

    /* PlotStackedHistograms("pT", h_pt, h_pt_bkg, h_pt_6fbkg, h_pt_5f, h_pt_2f, h_pt_4fbkg);
    PlotStackedHistograms("mInv", h_minv, h_minv_bkg, h_minv_6fbkg, h_minv_5f, h_minv_2f, h_minv_4fbkg);
    PlotStackedHistograms("pT_nonIso", h_pt_nonIso, h_pt_nonIso_bkg, h_pt_nonIso_6fbkg, h_pt_nonIso_5f, h_pt_nonIso_2f, h_pt_nonIso_4fbkg);
    PlotStackedHistograms("misspT", h_misspT, h_misspT_bkg, h_misspT_6fbkg, h_misspT_5f, h_misspT_2f, h_misspT_4fbkg);
    PlotStackedHistograms("missE", h_missE, h_missE_bkg, h_missE_6fbkg, h_missE_5f, h_missE_2f, h_missE_4fbkg);
    PlotStackedHistograms("Evis_nonIso", h_evis_nonIso, h_evis_nonIso_3f, h_evis_nonIso_6fbkg, h_evis_nonIso_5f, h_evis_nonIso_2f, h_evis_nonIso_4fbkg); */

    // --------- confusion matrix summed ------------
    TH2F* sum_cm = nullptr;
    if (categorisationTree4f && categorisationTree3f && categorisationTree6f && categorisationTree2f && categorisationTree5f && categorisationTree4f_bkg) {
      sum_cm = new TH2F("sum_cm", "Summed Confusion Matrix", 7, -0.5, 6.5, 7, -0.5, 6.5);
      sum_cm->SetXTitle("True Category");
      sum_cm->SetYTitle("Reconstructed Category");
      sum_cm->SetZTitle("Events");

      // -------- CONFUSION MATRIX REBUILT ------------

      // Loop over entries in the tree to fill the confusion matrix
      for (int i = 0; i < categorisationTree4f->GetEntries(); i++) {
        categorisationTree4f->GetEntry(i);
        sum_cm->Fill(TrueCat, RecoCatAdv);
      }
      for (int i = 0; i < categorisationTree3f->GetEntries(); i++) {
        categorisationTree3f->GetEntry(i);
        sum_cm->Fill(TrueCat, RecoCatAdv);
      }
      for (int i = 0; i < categorisationTree6f->GetEntries(); i++) {
        categorisationTree6f->GetEntry(i);
        sum_cm->Fill(TrueCat, RecoCatAdv);
      }
      for (int i = 0; i < categorisationTree2f->GetEntries(); i++) {
        categorisationTree2f->GetEntry(i);
        sum_cm->Fill(TrueCat, RecoCatAdv);
      }
      for (int i = 0; i < categorisationTree5f->GetEntries(); i++) {
        categorisationTree5f->GetEntry(i);
        sum_cm->Fill(TrueCat, RecoCatAdv);
      }
      for (int i = 0; i < categorisationTree4f_bkg->GetEntries(); i++) {
        categorisationTree4f_bkg->GetEntry(i);
        sum_cm->Fill(TrueCat, RecoCatAdv);
      }
     
      int sumEventsRecoCatAdv3 = 0;
      std::vector<TTree*> categorisationTrees = {categorisationTree3f, categorisationTree4f_bkg, categorisationTree2f, categorisationTree6f, categorisationTree4f, categorisationTree5f};
      
      for (auto& tree : categorisationTrees) {
        for (int i = 0; i < tree->GetEntries(); i++) {
          tree->GetEntry(i);
          if (RecoCatAdv == 3) {
          sumEventsRecoCatAdv3++;
          }
        }
      }
      std::cout << "Total number of events in RecoCatAdv == 3: " << sumEventsRecoCatAdv3 << std::endl;

      int nBinsX = sum_cm->GetNbinsX();
      int nBinsY = sum_cm->GetNbinsY();

      /* std::cout << "Bins with zero content:" << std::endl;
      for (int x = 1; x <= nBinsX; ++x) { // Bins start from 1 in ROOT
        for (int y = 1; y <= nBinsY; ++y) {
          if (sum_cm->GetBinContent(x, y) == 0) {
            std::cout << "Bin (" << x << ", " << y << ") is zero." << std::endl;
          }
        }
      } */
      // ------------------------------------------------------

      // Create a canvas to draw the confusion matrix
      TCanvas *canvas_cm = new TCanvas("canvas_cm", "Confusion Matrix", 800, 600);
      sum_cm->Draw("COLZ");
      canvas_cm->Update();
      int entries= 0;
      for(int i = 0; i < 7; i++){
        sum_cm->GetXaxis()->SetBinLabel(i+1, std::to_string(i).c_str());
        sum_cm->GetYaxis()->SetBinLabel(i+1, std::to_string(i).c_str());
        entries += sum_cm->GetBinContent(i + 1, 4); // Accumulate values in the 4th row (reconstruction of sl muon)
      }
      std::cout << "Total Entries (2f,3f,4f,5f,6f) in RecoCatAdv == 3 (muon) : " << entries << std::endl;
    }
	}
  
  // Function to check and write histograms if not already present
  auto write_histogram = [](TFile* file, const std::string& hist_name, TH1* hist) {
    if (file->Get(hist_name.c_str()) == nullptr) { // Check if histogram exists
      file->cd();
      hist->Write();
      std::cout << "Writing histogram: " << hist_name << std::endl;
    }else{
      std::cout << "Histogram already exists: " << hist_name << std::endl;
    }
  };
 
  // Store TProfiles in their respective files
  TFile* output_2f = TFile::Open(background_2f.c_str(), "UPDATE");
  TFile* output_3f = TFile::Open(background_3f.c_str(), "UPDATE");
  TFile* output_4f = TFile::Open(background_4f.c_str(), "UPDATE");
  TFile* output_6f = TFile::Open(background_6f.c_str(), "UPDATE");
  TFile* output_5f = TFile::Open(background_5f.c_str(), "UPDATE");
  TFile* output_ww = TFile::Open(ww_4f.c_str(), "UPDATE");

  // Write histograms if they don't already exist
  write_histogram(output_ww, "h_evis", histograms["h_evis"]);
  write_histogram(output_ww, "h_pt", histograms["h_pt"]);
  write_histogram(output_ww, "h_minv", histograms["h_minv"]);
  write_histogram(output_ww, "h_pt_nonIso", histograms["h_pt_nonIso"]);
  write_histogram(output_ww, "h_evis_nonIso", histograms["h_evis_nonIso"]);
  write_histogram(output_ww, "h_misspT", histograms["h_misspT"]);
  write_histogram(output_ww, "h_missE", histograms["h_missE"]);

  write_histogram(output_2f, "h_evis_2f", histograms["h_evis_2f"]);
  write_histogram(output_2f, "h_pt_2f", histograms["h_pt_2f"]);
  write_histogram(output_2f, "h_minv_2f", histograms["h_minv_2f"]);
  write_histogram(output_2f, "h_pt_nonIso_2f", histograms["h_pt_nonIso_2f"]);
  write_histogram(output_2f, "h_evis_nonIso_2f", histograms["h_evis_nonIso_2f"]);
  write_histogram(output_2f, "h_misspT_2f", histograms["h_misspT_2f"]);
  write_histogram(output_2f, "h_missE_2f", histograms["h_missE_2f"]);

  write_histogram(output_3f, "h_evis_bkg", histograms["h_evis_bkg"]);
  write_histogram(output_3f, "h_pt_bkg", histograms["h_pt_bkg"]);
  write_histogram(output_3f, "h_minv_bkg", histograms["h_minv_bkg"]);
  write_histogram(output_3f, "h_pt_nonIso_bkg", histograms["h_pt_nonIso_bkg"]);
  write_histogram(output_3f, "h_evis_nonIso_3f", histograms["h_evis_nonIso_3f"]);
  write_histogram(output_3f, "h_misspT_bkg", histograms["h_misspT_bkg"]);
  write_histogram(output_3f, "h_missE_bkg", histograms["h_missE_bkg"]);

  write_histogram(output_6f, "h_evis_6fbkg", histograms["h_evis_6fbkg"]);
  write_histogram(output_6f, "h_pt_6fbkg", histograms["h_pt_6fbkg"]);
  write_histogram(output_6f, "h_minv_6fbkg", histograms["h_minv_6fbkg"]);
  write_histogram(output_6f, "h_pt_nonIso_6fbkg", histograms["h_pt_nonIso_6fbkg"]);
  write_histogram(output_6f, "h_evis_nonIso_6fbkg", histograms["h_evis_nonIso_6fbkg"]);
  write_histogram(output_6f, "h_misspT_6fbkg", histograms["h_misspT_6fbkg"]);
  write_histogram(output_6f, "h_missE_6fbkg", histograms["h_missE_6fbkg"]);

  write_histogram(output_5f, "h_evis_5f", histograms["h_evis_5f"]);
  write_histogram(output_5f, "h_pt_5f", histograms["h_pt_5f"]);
  write_histogram(output_5f, "h_minv_5f", histograms["h_minv_5f"]);
  write_histogram(output_5f, "h_pt_nonIso_5f", histograms["h_pt_nonIso_5f"]);
  write_histogram(output_5f, "h_evis_nonIso_5f", histograms["h_evis_nonIso_5f"]);
  write_histogram(output_5f, "h_misspT_5f", histograms["h_misspT_5f"]);
  write_histogram(output_5f, "h_missE_5f", histograms["h_missE_5f"]);

  write_histogram(output_4f, "h_evis_4fbkg", histograms["h_evis_4fbkg"]);
  write_histogram(output_4f, "h_pt_4fbkg", histograms["h_pt_4fbkg"]);
  write_histogram(output_4f, "h_minv_4fbkg", histograms["h_minv_4fbkg"]);
  write_histogram(output_4f, "h_pt_nonIso_4fbkg", histograms["h_pt_nonIso_4fbkg"]);
  write_histogram(output_4f, "h_evis_nonIso_4fbkg", histograms["h_evis_nonIso_4fbkg"]);
  write_histogram(output_4f, "h_misspT_4fbkg", histograms["h_misspT_4fbkg"]);
  write_histogram(output_4f, "h_missE_4fbkg", histograms["h_missE_4fbkg"]);

  // Close the output files
  output_2f->Close();
  output_3f->Close();
  output_4f->Close();
  output_6f->Close();
  output_5f->Close();
  output_ww->Close();
}
void GetcategorisationTree(TFile* file, const std::string& objectname, TTree*& ttree_file, int& TrueCat, int& RecoCatBasic, int& RecoCatAdv){
  // Retrieve the TTree from the file
  file->GetObject(objectname.c_str(), ttree_file);
  if (!ttree_file) {
    std::cerr << "Error: TTree not found in the file." << std::endl;
    return;
  }

  // Set up branches to read data into variables by reference
  ttree_file->SetBranchAddress("TrueCat", &TrueCat);
  ttree_file->SetBranchAddress("RecoCatBasic", &RecoCatBasic);
  ttree_file->SetBranchAddress("RecoCatAdv", &RecoCatAdv);
 
}
void GetobservablesTree(TFile* file, const std::string& objectname, TTree*& ttree_file, double& Evis, double& pT, double& Evis_nonIso, double& pT_nonIso, double& missE, double& misspT, double& mInv){
  // Retrieve the TTree from the file
  file->GetObject(objectname.c_str(), ttree_file);
  if (!ttree_file) {
    std::cerr << "Error: TTree not found in the file." << std::endl;
    return;
  }

  // Set up branches to read data into variables by reference
  ttree_file->SetBranchAddress("Evis", &Evis);
  ttree_file->SetBranchAddress("pT", &pT);
  ttree_file->SetBranchAddress("Evis_nonIso", &Evis_nonIso);
  ttree_file->SetBranchAddress("pT_nonIso", &pT_nonIso);
  ttree_file->SetBranchAddress("missE", &missE);
  ttree_file->SetBranchAddress("misspT", &misspT);
  ttree_file->SetBranchAddress("mInv", &mInv);
}
TProfile* PlotObservable(TTree*& tree_file, const std::string& variable, const std::string& histName, int bins, double xMin, double xMax, const std::string& xTitle = "", const std::string& yTitle = "Events", int fillColor = kRed){
  gStyle->SetOptStat(11);
  // Define the draw command with the histogram name and binning
  std::string drawCmd = variable + ">>" + histName + "(" + std::to_string(bins) + "," + std::to_string(xMin) + "," + std::to_string(xMax) + ")";
  
  std::stringstream category; category << "RecoCatAdv == 3";

  // Draw the histogram from the tree
  tree_file->Draw(drawCmd.c_str() , category.str().c_str());

  // Retrieve the histogram from the canvas
  TProfile* hist = (TProfile*)gPad->GetPrimitive(histName.c_str());
  if (!hist) {
      std::cerr << "Error: Histogram could not be created." << std::endl;
      return nullptr;
  }

  // Customize histogram properties
  hist->SetXTitle(xTitle.c_str());
  hist->SetYTitle(yTitle.c_str());
  hist->SetFillColor(fillColor);

  // Return the histogram in case further manipulation is needed
  return hist;
}
void PlotStackedHistograms(const std::string& variable, TProfile* h_signal, TProfile* h_3f, TProfile* h_6f, TProfile* h_5f, TProfile* h_2f, TProfile* h_4fbkg){
  //Create a string for the observable and axis titles
  std::string Observable = "Stacked " + variable;
  std::string xAxis = variable + " [GeV]";

  // Create a THStack to hold and stack the histograms
  auto hs = new THStack("hs", Observable.c_str());

  // Add histograms to the stack
  hs->Add(h_signal);
  hs->Add(h_2f);
  hs->Add(h_3f);
  hs->Add(h_6f);
  hs->Add(h_5f);
  hs->Add(h_4fbkg);
  
  // Create a canvas to draw the stack
  TCanvas *canvas = new TCanvas(Observable.c_str(), Observable.c_str(), 800, 600);
  
  // Draw the stacked histograms
  hs->Draw("HIST");

  // Customize axis titles
  hs->GetXaxis()->SetTitle(xAxis.c_str());
  hs->GetYaxis()->SetTitle("Events");

  // Create and customize the legend
  auto legend = new TLegend(0.1, 0.75, 0.35, 0.85);
  legend->AddEntry(h_signal, "4f Signal", "f");
  legend->AddEntry(h_4fbkg, "4f Background", "f");
  legend->AddEntry(h_2f, "2f Background", "f");
  legend->AddEntry(h_3f, "3f Background", "f");
  legend->AddEntry(h_6f, "6f Background", "f");
  legend->AddEntry(h_5f, "5f Background", "f");
  legend->Draw();

  // Set the overall title for the stack
  hs->SetTitle(Observable.c_str());

  // Add "ILD work in progress" text to the canvas
  TLatex text;
  text.SetTextSize(0.035);          // Set text size
  text.SetTextAlign(34);           // Align at top-left corner
  text.DrawLatexNDC(0.95, 0.93935, "ILD #font[52]{work in progress}");
  
  // Update the canvas to display everything
  canvas->Update();
}