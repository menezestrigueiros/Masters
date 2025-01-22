#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <string>
#include "TFile.h"
#include "TTree.h"
#include "TStyle.h"
#include "TProfile.h"

// ----- COMMAND ---> root -l preselection.cc

void GetEventInfoTree(TFile* file, const std::string& objectname, TTree*& ttree_file, int& TrueCat, int& RecoCatBasic, int& RecoCatAdv, float& xSec, float& BPol0, float& BPol1);
void GetobservablesTree(TFile* file, const std::string& objectname, TTree*& ttree_file, double& Evis, double& pT, double& Evis_nonIso, double& pT_nonIso, double& missE, double& misspT, double& mInv);
bool query_histogram(const std::map<std::string, TH1*>& histograms, const std::string& histName, TH1*& histVariable);

void PolWeights(){
  // Input ROOT file

  std::string ww_4f = "../../../workarea/ROOTFiles/PreSelection4f_bkg.root";

  TFile* input_ww = TFile::Open(ww_4f.c_str(), "UPDATE");
  if (!input_ww || input_ww->IsZombie()) {
    std::cerr << "Failed to open file" << std::endl;
    return;
  }
  
  // Map to store histograms by their names
  std::map<std::string, TH1*> histograms;

  // Iterate over all keys in the ROOT file
  TIter next(input_ww->GetListOfKeys());
  TKey* key;
  while ((key = (TKey*)next())) {
    // Get the class name of the object
    TClass* cl = TClass::GetClass(key->GetClassName());
    if (!cl) continue;

    // Check if the object is a histogram (or derived from TH1)
    if (cl->InheritsFrom("TH1")) {
      // Retrieve the histogram
      TH1* hist = (TH1*)key->ReadObj();
      if (hist) {
          histograms[hist->GetName()] = hist;
          std::cout << "Loaded histogram: " << hist->GetName() << std::endl;
      }
    }
  }

  double event_weights = 0.0;
  // Create the new branch in the TTree


  // Vector to store histogram pointers
  std::vector<TH1*> histogramPointers = {
    nullptr, // Placeholder for h_evis_6fbkg
    nullptr, // Placeholder for h_pt_6fbkg
    nullptr, // Placeholder for h_minv_6fbkg
    nullptr, // Placeholder for h_pt_nonIso_6fbkg
    nullptr,  // Placeholder for h_evis_nonIso_6fbkg
    nullptr, // Placeholder for h_misspT_6fbkg
    nullptr,  // Placeholder for h_missE_6fbkg
  };

  // Corresponding histogram names
  std::vector<std::string> histogramNames_6f = {
    "h_evis_6fbkg",
    "h_pt_6fbkg",
    "h_minv_6fbkg",
    "h_pt_nonIso_6fbkg",
    "h_evis_nonIso_6fbkg",
    "h_misspT_6fbkg",
    "h_missE_6fbkg"
  };
  std::vector<std::string> histogramNames_5f = {
    "h_evis_5f",
    "h_pt_5f",
    "h_minv_5f",
    "h_pt_nonIso_5f",
    "h_evis_nonIso_5f",
    "h_misspT_5f",
    "h_missE_5f"
  };
  std::vector<std::string> histogramNames_2f = {
    "h_evis_2f",
    "h_pt_2f",
    "h_minv_2f",
    "h_pt_nonIso_2f",
    "h_evis_nonIso_2f",
    "h_misspT_2f",
    "h_missE_2f"
  };
  std::vector<std::string> histogramNames_3f = {
    "h_evis_bkg",
    "h_pt_bkg",
    "h_minv_bkg",
    "h_pt_nonIso_bkg",
    "h_evis_nonIso_3f",
    "h_misspT_bkg",
    "h_missE_bkg"
  };
  std::vector<std::string> histogramNames_4f = {
    "h_evis",
    "h_pt",
    "h_minv",
    "h_pt_nonIso",
    "h_evis_nonIso",
    "h_misspT",
    "h_missE"
  };
  std::vector<std::string> histogramNames_4fbkg = {
    "h_evis_4fbkg",
    "h_pt_4fbkg",
    "h_minv_4fbkg",
    "h_pt_nonIso_4fbkg",
    "h_evis_nonIso_4fbkg",
    "h_misspT_4fbkg",
    "h_missE_4fbkg"
  };
  // Vector to store all histogram name vectors
  std::vector<std::vector<std::string>> fullSM = {
    histogramNames_6f,
    histogramNames_5f,
    histogramNames_2f,
    histogramNames_3f,
    histogramNames_4f,
    histogramNames_4fbkg
  };


  for (size_t i = 0; i < histogramNames_4fbkg.size(); ++i) {
    query_histogram(histograms, histogramNames_4fbkg[i], histogramPointers[i]);
  }

  // TTree and variables
  TTree* categorisationTree4f;
  int TrueCat, RecoCatBasic, RecoCatAdv;
  float BPol0, BPol1, xSec;
  double Evis, pT, Evis_nonIso, pT_nonIso, missE, misspT, mInv;

  std::string process;
  std::string* process_ptr = &process;

  const double Pol0[] = {-.8, .8, -.8, .8};
  const double Pol1[] = { .3, -.3, -.3, .3};
  const double Lumi[] = {900, 900, 100, 100};

  GetEventInfoTree(input_ww, "observablesTree", categorisationTree4f, TrueCat, RecoCatBasic, RecoCatAdv, xSec, BPol0, BPol1);
  GetobservablesTree(input_ww, "observablesTree", categorisationTree4f, Evis, pT, Evis_nonIso, pT_nonIso, missE, misspT, mInv);
  categorisationTree4f->SetBranchAddress("proc", &process_ptr);
  
  // Check if the TTree exists
  if (!categorisationTree4f) {
    std::cerr << "Error: Failed to retrieve TTree." << std::endl;
    return;
  }

  

  // Clone the original tree
  TTree *clonedTree = categorisationTree4f->CloneTree(0);
  clonedTree->Branch("event_weights", &event_weights, "event_weights/D"); // 'D' for double type


  // Map to count events per process and polarization
  std::map<std::string, std::map<std::string, int>> event_counts;
  std::map<std::string, std::map<std::string, double>> weights; // {process: {polarization: weight}}
  std::map<std::string, std::map<std::string, float>> xSec_map; // {process: {polarization: xSec}}
  std::set<std::string> unique_processes;

  // Create a map to store BPol0 and BPol1 for each polarization
  std::map<std::string, std::pair<int, int>> polarization_map = {
    {"LR", {-1, 1}},
    {"RL", {1, -1}},
    {"LL", {-1, -1}},
    {"RR", {1, 1}},
    {"gamma-gamma", {0, 0}},
    {"right-gamma", {1, 0}},
    {"left-gamma", {-1, 0}},
    {"gamma-right", {0, 1}},
    {"gamma-left", {0, -1}},
    //{"Unknown", {0, 0}}
  };

  // Loop through the tree to count events
  int nEntries = categorisationTree4f->GetEntries();
  for (int j = 0; j < nEntries; ++j) {
    categorisationTree4f->GetEntry(j);

    // Determine polarization configuration
    std::string polarization;
    if (BPol0 == -1 && BPol1 == 1) polarization = "LR";
    else if (BPol0 == 1 && BPol1 == -1) polarization = "RL";
    else if (BPol0 == -1 && BPol1 == -1) polarization = "LL";
    else if (BPol0 == 1 && BPol1 == 1) polarization = "RR";
    // including photon interations
    else if(BPol0 == 0 && BPol1 == 0) polarization = "gamma-gamma";
    else if(BPol0 == 1 && BPol1 == 0) polarization = "right-gamma";
    else if(BPol0 == -1 && BPol1 == 0) polarization = "left-gamma";
    else if(BPol0 == 0 && BPol1 == 1) polarization = "gamma-right";
    else if(BPol0 == 0 && BPol1 == -1) polarization = "gamma-left";
    else polarization = "Unknown";

    // Keep track of unique processes
    unique_processes.insert(process);

    // Increment count for the process and polarization
    event_counts[process][polarization]++;

    // Store the xSec for the process-polarization combination
    xSec_map[process][polarization] = xSec;
  }
  std::cout << "Entries before scaling: " << histogramPointers[0]->GetEntries() << std::endl;

  // Second loop to calculate weights
  for (const auto& process_entry : event_counts) {
    const std::string& process_name = process_entry.first;
    const auto& polarizations = process_entry.second;

    // Loop through each polarization for the current process
    for (const auto& pol_entry : polarizations) {
      const std::string& polarization = pol_entry.first;
      int total_events_for_process_polarization = pol_entry.second;  // Events for this specific polarization
      cout << "------------------------------------------------------------------------------"  << endl;
      std::cout << "Total Events for Process: " << process_name << " Polarization: " << polarization
                << " = " << total_events_for_process_polarization << std::endl;

      // Retrieve the cross section and polarization factors
      double xSecProc = xSec_map[process_name][polarization]; // Cross section for the process-polarization pair

      // Retrieve BPol0 and BPol1 values from the map based on polarization
      int BeamPolarization0 = polarization_map[polarization].first;
      int BeamPolarization1 = polarization_map[polarization].second;

      // Calculate polarization factor and expected events
      double PolFac = (1 + Pol0[0] * BeamPolarization0) / 2. * (1 + Pol1[0] * BeamPolarization1) / 2.;
      //cout << "  (DEBUG) ---BeamPolarization0: " << BeamPolarization0 << " and BeamPolarization1: " << BeamPolarization1 << std::endl;
      double nEvt_expec = xSecProc * Lumi[0] * PolFac; // Expected number of events

      cout << "  (DEBUG) ---xSecProc: " << xSecProc << std::endl;
      cout << "  (DEBUG) ---PolFac: " << PolFac << std::endl;
      cout << "  (DEBUG) ---nEvt_exp: " << nEvt_expec << std::endl;
      

      // Calculate the weight based on the number of events for the current process-polarization pair
      double weight = nEvt_expec / total_events_for_process_polarization;
      std::cout << "Process: " << process_name << " Polarization: " << polarization
                << " Weight: " << weight << std::endl;
      cout << "------------------------------------------------------------------------------"  << endl;

      // Store the calculated weight for this specific process-polarization pair
      weights[process_name][polarization] = weight;     
    }
  }
  // Reset the histogram to clear existing entries
  histogramPointers[0]->Reset();
  histogramPointers[1]->Reset();
  histogramPointers[2]->Reset();
  histogramPointers[3]->Reset();
  histogramPointers[4]->Reset();
  histogramPointers[5]->Reset();
  histogramPointers[6]->Reset();

  // Third loop: Fill the histogram with weights
  for (int j = 0; j < nEntries; j++) {
    categorisationTree4f->GetEntry(j); // Load entry

    // Determine polarization configuration
    std::string polarization;
    if (BPol0 == -1 && BPol1 == 1) polarization = "LR";
    else if (BPol0 == 1 && BPol1 == -1) polarization = "RL";
    else if (BPol0 == -1 && BPol1 == -1) polarization = "LL";
    else if (BPol0 == 1 && BPol1 == 1) polarization = "RR";
    // Including photon interactions
    else if(BPol0 == 0 && BPol1 == 0) polarization = "gamma-gamma";
    else if(BPol0 == 1 && BPol1 == 0) polarization = "right-gamma";
    else if(BPol0 == -1 && BPol1 == 0) polarization = "left-gamma";
    else if(BPol0 == 0 && BPol1 == 1) polarization = "gamma-right";
    else if(BPol0 == 0 && BPol1 == -1) polarization = "gamma-left";
    else polarization = "Unknown";

    // Skip unknown configurations
    if (polarization == "Unknown") continue;

    // Retrieve the weight for this process-polarization pair
    double weight = weights[process][polarization];

    // Fill the histogram with "Visible Energy" and apply the weight
    histogramPointers[0]->Fill(Evis, weight);
    histogramPointers[1]->Fill(pT, weight);
    histogramPointers[2]->Fill(mInv, weight);
    histogramPointers[3]->Fill(pT_nonIso, weight);
    histogramPointers[4]->Fill(Evis_nonIso, weight);
    histogramPointers[5]->Fill(misspT, weight);
    histogramPointers[6]->Fill(missE, weight);

    // Set the weight for this event
    event_weights = weight;

    // Fill the TTree with the new weight
    clonedTree->Fill();

  }
  input_ww->cd();
  clonedTree->Write("observablesTree", TObject::kOverwrite);  // Overwrite the TTree
  //clonedTree->Write("weights", TObject::kOverwrite);  // Overwrite the TTree



  // -------------------------------------------------------------------------------------------------
  // Print sum of weights to verify if total number of expected events closely match for every histogram.
  std::cout << "(DEBUG) --- expected number of events rounded" << std::endl;
  double sum_of_weights = std::round(histogramPointers[0]->GetSumOfWeights());
  Double_t underflow_weight1 = histogramPointers[0]->GetBinContent(0);  // Underflow bin (bin 0)
  Double_t overflow_weight1 = histogramPointers[0]->GetBinContent(histogramPointers[0]->GetNbinsX() + 1);  // Overflow bin (last bin + 1)
  
  sum_of_weights += underflow_weight1 + overflow_weight1;
  std::cout << "Sum of weights in 'h_evis': " << sum_of_weights << std::endl;

  double sum_of_weightss = std::round(histogramPointers[1]->GetSumOfWeights());
  Double_t underflow_weight2 = histogramPointers[1]->GetBinContent(0);  // Underflow bin (bin 0)
  Double_t overflow_weight2 = histogramPointers[1]->GetBinContent(histogramPointers[1]->GetNbinsX() + 1);  // Overflow bin (last bin + 1)

  sum_of_weightss += underflow_weight2 + overflow_weight2;
  std::cout << "Sum of weights in 'h_pt': " << std::round(sum_of_weightss) << std::endl;

  double sum_of_weightsss = std::round(histogramPointers[2]->GetSumOfWeights());
  std::cout << "Sum of weights in 'h_minv: " << sum_of_weightsss << std::endl;

  double sum_of_weightssw = std::round(histogramPointers[3]->GetSumOfWeights());
  Double_t underflow_weight3 = histogramPointers[3]->GetBinContent(0);  // Underflow bin (bin 0)
  Double_t overflow_weight3 = histogramPointers[3]->GetBinContent(histogramPointers[3]->GetNbinsX() + 1);  // Overflow bin (last bin + 1)

  sum_of_weightssw += underflow_weight3 + overflow_weight3;
  std::cout << "Sum of weights in 'h_pt_nonIso': " << std::round(sum_of_weightssw) << std::endl;

  double sum_of_weightssww = std::round(histogramPointers[4]->GetSumOfWeights());
  Double_t underflow_weight4 = histogramPointers[4]->GetBinContent(0);  // Underflow bin (bin 0)
  Double_t overflow_weight4 = histogramPointers[4]->GetBinContent(histogramPointers[4]->GetNbinsX() + 1);  // Overflow bin (last bin + 1)

  sum_of_weightssww += + underflow_weight4 + overflow_weight4;
  std::cout << "Sum of weights in 'h_evis_nonIso': " << std::round(sum_of_weightssww) << std::endl;
  std::cout << "----------------------------------------------------" << std::endl;

  // -------------------------------------------------------------------------------------------------
  
  // Create a new histogram for the scaled version
  std::string newName;
  // Loop through the histograms and clone them
  for (size_t i = 0; i < histogramPointers.size(); ++i) {
    if (histogramPointers[i]) { // Ensure the pointer is valid
      // Create a new name by prepending "scaled_" to the original name
      newName = "scaled_" + histogramNames_4fbkg[i];

      // Clone the histogram and assign the new name
      TH1F* scaledHist = (TH1F*)histogramPointers[i]->Clone(newName.c_str());
      
      if (scaledHist) {
            std::cout << "Cloned and renamed histogram: " << newName << std::endl;
            if (input_ww->Get(newName.c_str()) != nullptr) {
              std::cout << "Histogram '" << newName << "' already exists in the file. Overwrite histogram." << std::endl;
              scaledHist->Write(newName.c_str(), TObject::kOverwrite);  // Overwrite the histogram
            } else {
              // Write the new scaled histogram
              scaledHist->Write(newName.c_str());
            }
        } else {
            std::cerr << "Failed to clone histogram: " << histogramNames_4fbkg[i] << std::endl;
        }
    } else {
      std::cerr << "Invalid histogram pointer for: " << histogramNames_4fbkg[i] << std::endl;
    }
  }
  
  input_ww->Close();

  // Print results
  std::cout << "==========================================" << std::endl;
  std::cout << "           Event Counts by Process        " << std::endl;
  std::cout << "==========================================" << std::endl;
  std::cout << "Total Entries: " <<  nEntries << std::endl;

  for (const auto& process_entry : event_counts) {
    const std::string& process_name = process_entry.first;
    
    const auto& polarizations = process_entry.second;
    int total_events_for_process = 0;
    for (const auto& pol_entry : process_entry.second) {
        total_events_for_process += pol_entry.second; // Sum all events for the process
    }

    std::cout << "Total events for process " << process_name << ": " << total_events_for_process << std::endl;

    std::cout << "Process: " << process_name << std::endl;
    for (const auto& pol_entry : polarizations) {
      std::cout << "  - " << pol_entry.first << ": " << pol_entry.second << " events" << std::endl;
    }
  }

  std::cout << "==========================================" << std::endl;
  std::cout << "        Process-Polarisation Weights      " << std::endl;
  std::cout << "==========================================" << std::endl;

  for (const auto& proc : unique_processes) {
    std::cout << "Process: " << proc << std::endl;
    for (const auto& pol_weight : weights[proc]) {
      std::cout << "  " << pol_weight.first << ": Weight = " << pol_weight.second << std::endl;
    }
    std::cout << "------------------------------------------" << std::endl;
  }  
}

void GetEventInfoTree(TFile* file, const std::string& objectname, TTree*& ttree_file, int& TrueCat, int& RecoCatBasic, int& RecoCatAdv, float& xSec, float& BPol0, float& BPol1){
  file->GetObject(objectname.c_str(), ttree_file);
  if (!ttree_file) {
    std::cerr << "Error: TTree not found in the file." << std::endl;
    return;
  }
  // Set up branches to read data into variables by reference
  ttree_file->SetBranchAddress("TrueCat", &TrueCat);
  ttree_file->SetBranchAddress("RecoCatBasic", &RecoCatBasic);
  ttree_file->SetBranchAddress("RecoCatAdv", &RecoCatAdv);
  ttree_file->SetBranchAddress("xSection", &xSec);
  ttree_file->SetBranchAddress("beamPol0", &BPol0);
  ttree_file->SetBranchAddress("beamPol1", &BPol1);
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
// Function to retrieve a histogram from the map and assign it to a variable
bool query_histogram(const std::map<std::string, TH1*>& histograms, const std::string& histName, TH1*& histVariable) {
    auto it = histograms.find(histName);
    if (it != histograms.end()) {
        histVariable = it->second; // Assign the histogram to the user-provided pointer
        if (histVariable) {
          std::cout << "Histogram: " << histVariable->GetName() << std::endl;
          std::cout << "Entries: " << histVariable->GetEntries() << std::endl;
          return true; // Success
        }
    } else {
        std::cerr << "Histogram '" << histName << "' not found in the map!" << std::endl;
    }
    histVariable = nullptr; // Ensure the variable is null if not found
    return false; // Failure
}