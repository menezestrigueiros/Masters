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
#include "nlohmann/json.hpp"

using jsonf = nlohmann::json;

// ----- COMMAND ---> root -l readJSON.cc
void GetcategorisationTree(TFile* file, const std::string& objectname, TTree*& ttree_file, int& TrueCat, int& RecoCatBasic, int& RecoCatAdv);
void GetobservablesTree(TFile* file, const std::string& objectname, TTree*& ttree_file, double& Evis, double& pT, double& Evis_nonIso, double& pT_nonIso, double& missE, double& misspT, double& mInv, double& weights);

string fileFormat = ".pdf";
string outspec = "";
// ---- input root files ----

string ww = "../../../workarea/ROOTFiles/PreSelection4f_bkg.root";
string outfolder = ".";

void readJSON(){
	// Path to the JSON file
	const std::string jsonFileName = "cuts.json";

	// Open the JSON file
	std::ifstream jsonFile(jsonFileName);
	if (!jsonFile.is_open()) {
		std::cerr << "Error: Could not open file " << jsonFileName << std::endl;
		return;
	}

	// Parse the JSON file
    nlohmann::json jsonData;
    try {
        jsonFile >> jsonData;
    } catch (const nlohmann::json::parse_error& error) {
        std::cerr << "JSON parsing error: " << error.what() << std::endl;
        return;
    }

	// Access JSON data
	// Visible Energy
    double EvisCut_min = jsonData["PreSelection"]["Evis"][0];
	double EvisCut_max = jsonData["PreSelection"]["Evis"][1];
	// pT
	double pTCut = jsonData["PreSelection"]["pT"];
	// Missing pT
	double misspTCut = jsonData["PreSelection"]["misspT"];
	// Missing Energy
	double missECut_min = jsonData["PreSelection"]["missE"][0];
	double missECut_max = jsonData["PreSelection"]["missE"][1];
	// Invariant Mass
	double mInvCut_min = jsonData["PreSelection"]["mInv"][0];
	double mInvCut_max = jsonData["PreSelection"]["mInv"][1];
	// pT (non-isolated)
	double pT_nonIsoCut = jsonData["PreSelection"]["pT_nonIso"];
	// Description
    std::string description = jsonData["description"];

	// Print values
    std::cout << "Description: " << description << std::endl;
	std::cout << "EvisCut_min: " << EvisCut_min << std::endl;
	std::cout << "EvisCut_max: " << EvisCut_max << std::endl;
	std::cout << "pTCut: " << pTCut << std::endl;

	// ----------- ACCESS ROOT FILE ------------
	TFile* input_ww = TFile::Open(ww.c_str());

	if(!input_ww || input_ww->IsZombie()){
		std::cerr << "Failed to open file" << std::endl;
		return;
	}

	// Retrieve the histogram (e.g., "VisibleEnergy")
	TProfile* h_scaled_visible_energy = (TProfile*)input_ww->Get("scaled_h_evis_4fbkg");
	if (!h_scaled_visible_energy) {
		std::cerr << "Failed to retrieve histogram: VisibleEnergy" << std::endl;
		return;
	}
	TProfile* h_scaled_pT = (TProfile*)input_ww->Get("scaled_h_minv_4fbkg");
	if (!h_scaled_pT) {
		std::cerr << "Failed to retrieve histogram: pT" << std::endl;
		return;
	}
	TProfile* h_scaled_pT_nonIso = (TProfile*)input_ww->Get("scaled_h_pt_nonIso_4fbkg");
	if (!h_scaled_pT_nonIso) {
		std::cerr << "Failed to retrieve histogram: pT_nonIso" << std::endl;
		return;
	}
	TProfile* h_scaled_misspT = (TProfile*)input_ww->Get("scaled_h_misspT_4fbkg");
	if (!h_scaled_misspT) {
		std::cerr << "Failed to retrieve histogram: misspT" << std::endl;
		return;
	}
	TProfile* h_scaled_missE = (TProfile*)input_ww->Get("scaled_h_missE_4fbkg");
	if (!h_scaled_missE) {
		std::cerr << "Failed to retrieve histogram: missE" << std::endl;
		return;
	}
	TProfile* h_scaled_mInv = (TProfile*)input_ww->Get("scaled_h_minv_4fbkg");
	if (!h_scaled_mInv) {
		std::cerr << "Failed to retrieve histogram: mInv" << std::endl;
		return;
	}


	if(1){
		// --------  TREES  ------------
    	TTree *tree, *tree_2f, *tree_3f, *tree_6f, *tree_5f;
		TTree *categorisationTree4f, *categorisationTree4f_bkg, *categorisationTree2f, *categorisationTree3f, *categorisationTree6f, *categorisationTree5f;
		// ------ VARIABLES ------------
		double Evis, pT, mInv;
		double Evis_nonIso, pT_nonIso;
		double missE, misspT;

		int TrueCat, RecoCatBasic, RecoCatAdv;
		double event_weights;

		// ------- Get Observables -----------
		GetobservablesTree(input_ww, "observablesTree", tree, Evis, pT, Evis_nonIso, pT_nonIso, missE, misspT, mInv, event_weights);
		GetcategorisationTree(input_ww, "observablesTree", categorisationTree4f, TrueCat, RecoCatBasic, RecoCatAdv);
		

		// Output file and new TTree
    	TFile *outputFile = TFile::Open("filtered.root", "RECREATE");
    	TTree *filteredTree = tree->CloneTree(0); // Clone structure, no entries

		//debug
		if (tree ) {
      		std::cout << "Successfully retrieved TTree." << std::endl;
      		// Process the tree, e.g., loop over entries
    	} else {
      		std::cerr << "Failed to retrieve TTree." << std::endl;
    	}

		// Initialize counters for total entries and cuts
    	double totalEntries = std::round(h_scaled_visible_energy->GetSumOfWeights());
		//int totalEntries = tree->GetEntries("RecoCatAdv == 3"); if we want a specific bin or a set of
		int _nPassed1 = 0, _nPassed2 = 0, _nPassed3 = 0, _nPassed4 = 0, _nPassed5 = 0, _nPassed6 = 0;
		int _nEvents = 0, _nEvents_Weighted = 0;

		double _nPassed1_Weighted = 0, _nPassed2_Weighted = 0, _nPassed3_Weighted = 0, _nPassed4_Weighted = 0, _nPassed5_Weighted = 0, _nPassed6_Weighted = 0;
		
		// Loop through all entries in the TTree
		int nEntries = tree->GetEntries();
		for (int i = 0; i < nEntries; i++) {
			tree->GetEntry(i);
			
			//std::cout << "Event weight: " << event_weights << std::endl;

			// Filter entries based on RecoCatAdv == 3
			if (RecoCatAdv != 3) continue;

			// Increment total events
			_nEvents++;

			// Apply the first cut on Evis (visible energy)
			if (Evis < EvisCut_min || Evis > EvisCut_max) continue;  // Skip event if it fails this cut
			_nPassed1++;  // Event survived cut
			_nPassed1_Weighted += event_weights;  // Add the weighted count

			// Apply the second cut on mInv (invariant mass)
			if (mInv < mInvCut_min || mInv > mInvCut_max) continue;  // Skip event if it fails this cut
			_nPassed2++;  // Event survived both cuts (Evis and mInv)
			_nPassed2_Weighted += event_weights;  // Add the weighted count

			// Apply the third cut on pT_nonIso (non-isolated transverse momentum)
			if (pT_nonIso < pT_nonIsoCut) continue;  // Skip event if it fails this cut
			_nPassed3++;  // Event survived all previous cuts (Evis, mInv, pT_nonIso)
			_nPassed3_Weighted += event_weights;  // Add the weighted count

			// Apply the fourth cut on pT (transverse momentum)
			if (pT < pTCut) continue;  // Skip event if it fails this cut
			_nPassed4++;  // Event survived all previous cuts (Evis, mInv, pT_nonIso, pT)
			_nPassed4_Weighted += event_weights;  // Add the weighted count

			// Apply the fifth cut on misspT (missing transverse momentum)
			if (misspT < misspTCut) continue;  // Skip event if it fails this cut
			_nPassed5++;  // Event survived all previous cuts (Evis, mInv, pT_nonIso, pT, misspT)
			_nPassed5_Weighted += event_weights;  // Add the weighted count

			// Apply the sixth cut on missE (missing energy)
			if (missE < missECut_min || missE > missECut_max) continue;  // Skip event if it fails this cut
			_nPassed6++;  // Event survived all previous cuts (Evis, mInv, pT_nonIso, pT, misspT, missE)
			_nPassed6_Weighted += event_weights;  // Add the weighted count
		}

		double total_weight_sum = 0;

		for (int i = 0; i < nEntries; i++) {
			tree->GetEntry(i);
			// Filter entries based on RecoCatAdv == 3
			if (RecoCatAdv != 3) continue;
    		
   			total_weight_sum += event_weights;

		}
		std::cout << "Total weight sum: " << total_weight_sum << std::endl;

		// Print the cut flow table
		std::cout << "Cut Flow Summary:" << std::endl;
		std::cout << "--------------------------------" << std::endl;
		std::cout << "Total Weighted entries: " << std::round(h_scaled_visible_energy->GetSumOfWeights()) << std::endl;
		std::cout << "Total entries: " << nEntries << std::endl;
		std::cout << "Entries in RecoCatAdv = 3: " << _nEvents << " (weighted: " << total_weight_sum << ")" << std::endl;
		std::cout << "After Evis Cut  (120.0 < Evis < 240.0 [GeV])  : " << _nPassed1 << " (weighted: " << _nPassed1_Weighted << ")" << std::endl;
		std::cout << "After mInv Cut (60.0 < M < 110.0 [GeV])       : " << _nPassed2 << " (weighted: " << _nPassed2_Weighted << ")" << std::endl;
		std::cout << "After pT_nonIso Cut (pT_nonIso > 20.0 [GeV])  : " << _nPassed3 << " (weighted: " << std::round(_nPassed3_Weighted) << ")" << std::endl;
		std::cout << "After pT Cut (pT > 10.0 [GeV])                : " << _nPassed4 << " (weighted: " << std::round(_nPassed4_Weighted) << ")" << std::endl;
		std::cout << "After misspT Cut (misspT > 10.0 [GeV])        : " << _nPassed5 << " (weighted: " << std::round(_nPassed5_Weighted) << ")" << std::endl;
		std::cout << "After missE Cut (20.0 < missE < 120.0 [GeV])  : " << _nPassed6 << " (weighted: " << std::round(_nPassed6_Weighted) << ")" << std::endl;
		std::cout << "--------------------------------" << std::endl;
		std::cout << "-> Efficiency: " << (double)_nPassed6_Weighted/total_weight_sum*100 << " %" << std::endl;

		filteredTree->Write();
		outputFile->Close();
		input_ww->Close();
	}
	

}
void GetobservablesTree(TFile* file, const std::string& objectname, TTree*& ttree_file, double& Evis, double& pT, double& Evis_nonIso, double& pT_nonIso, double& missE, double& misspT, double& mInv, double& weights){
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
  ttree_file->SetBranchAddress("event_weights", &weights);
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