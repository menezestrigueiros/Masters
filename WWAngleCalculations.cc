#include "../include/WWAngleCalculationProcessor.h"
#include "../include/WWTools.h"
#include <iostream>
#include <cmath>
#include <vector>
#include <array> 
#include <algorithm>
#include <sstream>
#include "TStyle.h"


#ifdef MARLIN_USE_AIDA
#include <marlin/AIDAProcessor.h>
#include <AIDA/IHistogramFactory.h>
#include <AIDA/ICloud1D.h>
//#include <AIDA/IHistogram1D.h>
#endif
#include "IMPL/LCEventImpl.h" 
#include "TrueJet.h"
#include "IMPL/LCCollectionVec.h"
#include "IMPL/SimCalorimeterHitImpl.h"
#include "IMPL/CalorimeterHitImpl.h"
#include "IMPL/MCParticleImpl.h" 
#include "IMPL/TrackerHitImpl.h" 
#include "IMPL/TrackImpl.h" 
#include "IMPL/ClusterImpl.h" 
#include "IMPL/ReconstructedParticleImpl.h" 
#include <CLHEP/Vector/LorentzVector.h>
#include "JetFitObject.h"
#include "IMPL/ParticleIDImpl.h" 
#include "IMPL/LCFlagImpl.h" 
#include "MCTree.h"
#include "IMPL/LCRelationImpl.h"
#include "UTIL/LCRelationNavigator.h"
#include "UTIL/PIDHandler.h"
#include <Math/GenVector/AxisAngle.h>


#include <EVENT/LCCollection.h>
#include <EVENT/ReconstructedParticle.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TText.h>

#include <Math/Vector4D.h>
#include <Math/Vector3D.h>
#include "Math/GenVector/Rotation3D.h"
#include <Math/GenVector/VectorUtil.h>
#include "Math/GenVector/PxPyPzE4D.h"

#include "Math/GenVector/DisplacementVector3D.h"
#include "Math/GenVector/GenVectorIO.h"

using namespace lcio ;
using namespace marlin ;
using namespace std;
using namespace WWTools;

vector <double> legendre_recursive(const double& x , const int& n,const vector<int> moments );

WWAngleCalculationProcessor  aWWAngleCalculationProcessor ;


WWAngleCalculationProcessor::WWAngleCalculationProcessor() : Processor("WWAngleCalculationProcessor") {
  
  // modify processor description
  _description = "WW Angle Calculation " ;
  
   
  // register steering parameters: name, description, class-variable, default value

  registerInputCollection( LCIO::MCPARTICLE, 
			   "MCParticles" , 
			   "Name of the MCParticle collection"  ,
			   _MCParColName ,
			   std::string("MCParticlesSkimmed") ) ;
  registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE, 
			   "PandoraPFOs" , 
			   "Name of the ReconstructedParticle collection"  ,
			   _PandoraPFOsCol ,
			   std::string("PandoraPFOs") ) ;
}


void WWAngleCalculationProcessor::init() { 
  //gStyle->SetPaintTextFormat("4.3f");
  //gStyle->SetPalette(1);


  _TTreeFile = new TFile("WWAngle.root", "RECREATE");
  gStyle->SetNumberContours(256);

  _star_phi  = new TH1F(" #phi * of lepton from W+", " ; #phi_{l} * ; Events", 100, -3.3, 3.3);
  _star_phi_minus  = new TH1F(" #phi * of lepton from W-", " ; #phi_{l} * ; Events", 100, -3.3, 3.3);


  _star_theta  = new TH1F(" cos(#theta_{lepton}) *", " ;  cos(#theta_{lepton}) * ; Events", 40, -1.1, 1.1);
  _cos_prodw = new TH1F(" Full simulation cos(#theta_{W}) ", " ; cos(#theta_{w}) ; Events", 100, -1.1, 1.1);

  _hTest = new TH1F(" WW explicit cos(#theta_{W})", " ; cos(#theta_{w}) ; Events", 100, -1.1, 1.1);

  // usually a good idea to
  printParameters() ;
  
  _nRun = 0 ;
  _nEvt = 0 ;
  
}

void WWAngleCalculationProcessor::processRunHeader( LCRunHeader*  /*run*/) {
  _nRun++ ;
  _nEvt = 0;
}

double WWAngleCalculationProcessor::CosineTheta(std::vector<double>& v1){
  double mag = std::sqrt( v1[0]*v1[0] + v1[1]*v1[1] + v1[2]*v1[2]);
  double cosine = v1[2]/mag;
  return cosine;
}


void WWAngleCalculationProcessor::processEvent( LCEvent * evt ){ 

  _nEvt ++;


  streamlog_out(DEBUG) 
	<< " processing event " << evt->getEventNumber() 
	<< "  in run "          << evt->getRunNumber() 
	<< std::endl ;

  LCCollection *col_muon{}, *col_electron{}, *col_tau{}, *col_mcparticles{}, *col_pfo{}, *col_jets{};
  
  // fill histogram from LCIO data :
  try{//col_muon    = evt->getCollection(_MuonColName);
  //col_electron    = evt->getCollection(_ElectronColName);
  //col_tau         = evt->getCollection(_TauColName);
  col_mcparticles = evt->getCollection(_MCParColName);
  col_pfo         = evt->getCollection(_PandoraPFOsCol);
  //col_jets = evt->getCollection(_PFOsMinusIsolatedObjetcs);
  }
  catch(DataNotAvailableException &e){
    streamlog_out(DEBUG) << "Input collections not found - skipping event " << _nEvt << std::endl;
    return;
  }

  MCParticle* mc = dynamic_cast <MCParticle*> (col_mcparticles->getElementAt(4)); //electron after ISR
  //MCParticle* mc1 = dynamic_cast <MCParticle*> (col_mcparticles->getElementAt(0)); before ISR

  std::vector<MCParticle*> daughters = mc->getDaughters();

  int incoming_pdg = abs(mc->getPDG()); // get PDG of incoming e- or e+ after ISR
  //int daughters_pdg = 0;

  int lWcharge = 0; //leptonic W charge (Reconstructed from daughters particles)
  int hWcharge = 0; //hadronic W charge (Reconstructed from daughters particles)
  int wcharge = 0; //explicit W charge (PDG == abs(24))


  ////////// W BOSONS ///////////////
  PxPyPzEVector lv_W1(0, 0, 0, 0); //Explicit W
  PxPyPzEVector lv_lW(0, 0, 0, 0); // Non Explicit Leptonic W
  PxPyPzEVector lv_hW(0, 0, 0, 0); // Non Explicit hadronic W

  ////////// LEPTONS/QUARKS //////////
  PxPyPzEVector lv_leptons_w1(0, 0, 0, 0); //leptons from explicit W
  PxPyPzEVector lv_leptons(0, 0, 0, 0); //leptons from reconstructed W
  PxPyPzEVector lv_quarks(0, 0, 0, 0); //quarks from reconstructed W

  PxPyPzEVector transf_lv(0, 0, 0, 0); //general transformed lorentz vector

  _nW.push_back(0);
  _nNonW.push_back(0);
 

  std::vector<int> general_daughters_pdgs = {}; //array with all daughters from e+e-
  std::vector<int> boson_daughters_pdg_array{}; //daughters that come from a W or Z boson

  std::vector<double> leptonic_w(4); //0,0,0,0
  std::vector<double> hadronic_w(4); //0,0,0,0
  std::vector<double> vector_w(4);
  std::vector<double> boost(4);

  std::vector<MCParticle*> boson_daughters = {};

  ROOT::Math::XYZVector v1(0, 0, 1);

  MCParticle* quarks; //MC information of each particle in hadronic decay
  MCParticle* leptons; //MC information of each particle in a leptonic decay

  // TEST ATTEMPTS************** For every explicit W (PDG == 24)
  {
    for(int i = 0; i < col_mcparticles->getNumberOfElements(); i++){
      MCParticle* mc1 = static_cast<MCParticle*> (col_mcparticles->getElementAt(i));
      if ( mc1->getPDG() != -24 ) continue;
      _hTest->Fill(mc1->getMomentum()[2] / std::hypot(mc1->getMomentum()[0], mc1->getMomentum()[1], mc1->getMomentum()[2]));
    }
  }

  if(incoming_pdg == 11){
    if(daughters.size() <= 3){
      streamlog_out(MESSAGE) << " ////////// [2/3 DAUGHTERS] in _nEvt " << _nEvt << " //////////" << endl;
      for ( long unsigned int i = 0; i < daughters.size(); i++){ 
        int daughters_pdg = abs(daughters[i]->getPDG());

        if (daughters_pdg == 24){ //Explicit MC W boson 4v and charge
          MCParticle* wboson = daughters[i];
          lv_W1 = PxPyPzEVector(wboson->getMomentum()[0], wboson->getMomentum()[1], wboson->getMomentum()[2], wboson->getEnergy());
          wcharge = wboson->getCharge();

          lv_W1.GetCoordinates(vector_w.begin(), vector_w.end());      

          if(vector_w[0] != 0){
            _W++;
            { //for debugging
              streamlog_out(MESSAGE) << "[ON SHELL] W 4v in " << _nEvt << ": [ "; 
              for (double k = 0; k < 4; k++){
                streamlog_out(MESSAGE) << vector_w[k] << ", ";
              }
              streamlog_out(MESSAGE) << "]" << endl;
              streamlog_out(MESSAGE) << "[ON SHELL] W CHARGE in _nEvt " << _nEvt << ": " << wcharge << endl;
            }

            //////////////////////// W- SELECTION //////////////////
            double prod_cosine = CosineTheta(vector_w);
            if(wcharge == -1) _cos_prodw->Fill(prod_cosine);

            //////////// LEPTONIC DECAY ///////////
            boson_daughters = wboson->getDaughters();
            for ( long unsigned int l = 0; l < boson_daughters.size(); l++){
              MCParticle* lepton_from_w = boson_daughters[l];
              int lepton_pdg = abs(boson_daughters.at(l)->getPDG());

              if(lepton_pdg == 11 || lepton_pdg == 13 || lepton_pdg == 15){ //if it`s a electron, muon or tau
                lv_leptons_w1 = PxPyPzEVector(lepton_from_w->getMomentum()[0], lepton_from_w->getMomentum()[1], lepton_from_w->getMomentum()[2], lepton_from_w->getEnergy());
                
                /////////// DECAY ANGLES //////////
                transf_lv = starVector(lv_W1, lv_leptons_w1, v1);
                transf_lv.GetCoordinates(boost.begin(), boost.end());
                double cos_star_theta = CosineTheta(boost);

                /////// W- or W+ //////
                if(wcharge == -1){
                  _Wleptonic++; //explicit W- that decay leptonically
                  _star_theta->Fill(cos_star_theta);
                  _star_phi_minus->Fill(transf_lv.Phi());
                }
                if(wcharge == 1){
                  _Wleptonic1++; //explicit W+ that decay leptonically
                  _star_phi->Fill(transf_lv.Phi());
                }
              }

              if((lepton_pdg == 1 || lepton_pdg == 3 || lepton_pdg == 5)){
                if(wcharge == -1) _Whadronic++; //explicit W- that decay hadronically
                if(wcharge == 1) _Whadronic1++;//explicit W+ that decay hadronically
              } 
            }
          }
        }

        if (daughters_pdg != 24){ // Non-explicit W boson

          if(daughters_pdg < 10){ 
            quarks = daughters[i];
            hWcharge = -wcharge; //opposite charge to the explicit W charge 
            lv_hW += PxPyPzEVector(quarks->getMomentum()[0], quarks->getMomentum()[1], quarks->getMomentum()[2], quarks->getEnergy()); //w 4v
            lv_quarks = PxPyPzEVector(quarks->getMomentum()[0], quarks->getMomentum()[1], quarks->getMomentum()[2], quarks->getEnergy()); //leptons 4v
          }else if(daughters_pdg > 10){
            leptons = daughters[i];
            lWcharge = -wcharge;
            lv_lW += PxPyPzEVector(leptons->getMomentum()[0], leptons->getMomentum()[1], leptons->getMomentum()[2], leptons->getEnergy());
            if(daughters_pdg == 11 || daughters_pdg == 13 || daughters_pdg == 15) lv_leptons = PxPyPzEVector(leptons->getMomentum()[0], leptons->getMomentum()[1], leptons->getMomentum()[2], leptons->getEnergy());
          }
        }
        lv_lW.GetCoordinates(leptonic_w.begin(), leptonic_w.end()); //LEPTONIC W
        lv_hW.GetCoordinates(hadronic_w.begin(), hadronic_w.end()); //HADRONIC W

      }

      if(hadronic_w[0] != 0){
        _W++;
        { //for debbuging
          streamlog_out(MESSAGE) << "[OFFSHELL] W 4v in " << _nEvt << ": [ "; 
          for (double k = 0; k < 4; k++){
            streamlog_out(MESSAGE) << hadronic_w[k] << ", ";
          }
          streamlog_out(MESSAGE) << "]" << endl;
          streamlog_out(MESSAGE) << "[OFFSHELL] hadronic W CHARGE in _nEvt " << _nEvt << ": " << hWcharge << endl;
        }

        /////// W- or W+ //////
        if(hWcharge == -1){
          _Whadronic++; // non explicit W- that decay hadronically
          double prod_cosine = CosineTheta(hadronic_w);
          _cos_prodw->Fill(prod_cosine);     
        }
        if(hWcharge == 1) _Whadronic1++ ; // non explicit W- that decay hadronically
      }
      if(leptonic_w[0] != 0){
        _W++;
        { //for debugging
          streamlog_out(MESSAGE) << "[OFFSHELL] W 4v in " <<_nEvt << ": [ "; 
          for (double k = 0; k < 4; k++){
            streamlog_out(MESSAGE) << leptonic_w[k] << ", ";
          }
          streamlog_out(MESSAGE) << "]" << endl;
          streamlog_out(MESSAGE) << "[OFFSHELL] W LEPTON CHARGE in _nEvt " << _nEvt << ": " << lWcharge << endl;
        }

        //Invariant Mass
        double inv_mass = lv_lW.M();
        streamlog_out(MESSAGE) << "[OFFSHELL] W INV MASS " << inv_mass << endl;
        
        /////////// W REST FRAME //////////
        transf_lv = starVector(lv_lW, lv_leptons, v1);
        transf_lv.GetCoordinates(boost.begin(), boost.end());

        { //for debbuging
          streamlog_out(MESSAGE) << "[LEP] transformed 4v in " << _nEvt << ": [ "; 
          for (double k = 0; k < 4; k++){
            streamlog_out(MESSAGE) << boost[k] << ", ";
          }
          streamlog_out(MESSAGE) << "]" << endl;
        }
        
        /////// W- or W+ //////
        if(lWcharge == -1){
          _Wleptonic++; 

          /////////// W PRODUCTION ANGLE //////////
          double prod_cosine = CosineTheta(leptonic_w);

          /////////// DECAY ANGLES //////////
          transf_lv.GetCoordinates(boost.begin(), boost.end());
          double cos_star_theta = CosineTheta(boost);

          /////////// PLOTS //////////
          _cos_prodw->Fill(prod_cosine);
          _star_theta->Fill(cos_star_theta);
          _star_phi_minus->Fill(transf_lv.Phi());
        }
        if(lWcharge == 1){
          _Wleptonic1++;
          _star_phi->Fill(transf_lv.Phi());
        }
      }
    
    }else if (daughters.size() == 4){ //non-explicit W procedure
      streamlog_out(MESSAGE) << "/////////////// [4 DAUGHTERS] in _nEvt " << _nEvt << " ///////////////" << endl;
      for(long unsigned int x= 0; x < daughters.size(); x++){
        int daughters_pdg = abs(daughters[x]->getPDG()); //get PDG of each daughter
        if(daughters_pdg < 10){ //hadronic decay
          quarks = daughters[x];
          //hWcharge += daughters[x]->getCharge(); if uncommented, the charge sum is 0 for hW due to 'quarks' and 'leptons' same definition (Ibelieve)
          lv_hW += PxPyPzEVector(quarks->getMomentum()[0], quarks->getMomentum()[1], quarks->getMomentum()[2], quarks->getEnergy());
          lv_quarks = PxPyPzEVector(quarks->getMomentum()[0], quarks->getMomentum()[1], quarks->getMomentum()[2], quarks->getEnergy()); //leptons 4v
        }
        if(daughters_pdg > 10){ //leptonic decay
          leptons = daughters[x];
          lWcharge += daughters[x]->getCharge();
          hWcharge = -lWcharge;
          lv_lW += PxPyPzEVector(leptons->getMomentum()[0], leptons->getMomentum()[1], leptons->getMomentum()[2], leptons->getEnergy());
          if(daughters_pdg == 11 || daughters_pdg == 13 || daughters_pdg == 15) lv_leptons = PxPyPzEVector(leptons->getMomentum()[0], leptons->getMomentum()[1], leptons->getMomentum()[2], leptons->getEnergy());
        }
        lv_lW.GetCoordinates(leptonic_w.begin(), leptonic_w.end());//for debbuging
        lv_hW.GetCoordinates(hadronic_w.begin(), hadronic_w.end());//for debbuging
      }

      if (leptonic_w[0] != 0){
        _W++;
        { //for debbugging
          streamlog_out(MESSAGE) << "[4 DAUGHTERS] W 4v in " <<_nEvt << ": [ "; 
          for (double k = 0; k < 4; k++){
            streamlog_out(MESSAGE) << leptonic_w[k] << ", ";
          }
          streamlog_out(MESSAGE) << "]" << endl;
          streamlog_out(DEBUG) << "[CHARGE] LEPTONIC W in _nEvt " << _nEvt << ": " << lWcharge << endl;
        }
        
        double inv_mass = lv_lW.M();
        streamlog_out(MESSAGE) << "[4 DAUGHTERS] W INV MASS: " << inv_mass << endl;

        /////////// W REST FRAME //////////
        transf_lv = starVector(lv_lW, lv_leptons, v1);
        transf_lv.GetCoordinates(boost.begin(), boost.end()); //W REST FRAME 

        { //for debugging
          streamlog_out(MESSAGE) << "[LEPTON] TRANSF 4v in " << _nEvt << ": [ "; 
          for (double k = 0; k < 4; k++){
            streamlog_out(MESSAGE) << boost[k] << ", ";
          }
          streamlog_out(MESSAGE) << "]" << endl;
        }

        /////// W- or W+ //////
        if(lWcharge == -1){
          _Wleptonic++;

          /////////// W PRODUCTION ANGLE //////////
          double prod_cosine = CosineTheta(leptonic_w);

          /////////// DECAY ANGLES //////////
          transf_lv.GetCoordinates(boost.begin(), boost.end()); 
          double cos_star_theta = CosineTheta(boost);

          ///////// PLOTS //////////
          _cos_prodw->Fill(prod_cosine);
          _star_theta->Fill(cos_star_theta);
          _star_phi_minus->Fill(transf_lv.Phi());
        }
        if(lWcharge == 1){
          _Wleptonic1++;
          _star_phi->Fill(transf_lv.Phi());
        }
      }
      if (hadronic_w[0] != 0){
        _W++;
        { //debugging
          streamlog_out(DEBUG) << "[4 DAUGHTERS] W 4v in " << _nEvt << ": [ "; 
          for (double k = 0; k < 4; k++){
            streamlog_out(DEBUG) << hadronic_w[k] << ", ";
          }
          streamlog_out(DEBUG) << "]" << endl;
        }
        
        /////// W- or W+ //////
        if(hWcharge == -1){ 
          _Whadronic++;
          double prod_cosine = CosineTheta(hadronic_w);
          _cos_prodw->Fill(prod_cosine);     
        }
        if(hWcharge == 1) _Whadronic1++ ;
        
      } 
      if (lWcharge != 0)streamlog_out(MESSAGE) << "[CHARGE] LEPTONIC W in _nEvt " << _nEvt << ": " << lWcharge << endl;
      if (hWcharge != 0)streamlog_out(MESSAGE) << "[CHARGE] HADRONIC W in _nEvt " << _nEvt << ": " << hWcharge << endl; 
    }
  }
}

void WWAngleCalculationProcessor::check( LCEvent *  /*evt*/ ) { 
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}


void WWAngleCalculationProcessor::end(){ 

  _TTreeFile->cd();
  //Histograms
  _star_phi->Write();
  _star_phi_minus->Write();
  _star_theta->Write();
  _cos_prodw->Write();
  _hTest->Write();

  _nW[_W]++;
  _nNonW[_nonW]++;
  

  streamlog_out(MESSAGE) << "W- that decay hadronically: " << _Whadronic << std::endl; 
  streamlog_out(MESSAGE) << "W- that decay leptonically: " << _Wleptonic << std::endl; 
  streamlog_out(MESSAGE) << "W+ that decay hadronically: " << _Whadronic1 << std::endl; 
  streamlog_out(MESSAGE) << "W+ that decay leptonically: " << _Wleptonic1 << std::endl; 

  streamlog_out(MESSAGE) << "TOTAL W-: " << _Wleptonic + _Whadronic << std::endl; 
  streamlog_out(MESSAGE) << "TOTAL W+: " << _Wleptonic1 + _Whadronic1 << std::endl; 

  streamlog_out(MESSAGE) << "How many total W bosons? --> " << _W << std::endl; 

  

  _TTreeFile->Close();
  // 	    << " processed " << _nEvt << " events in " << _nRun << " runs "
  // 	    << std::endl ;
}



