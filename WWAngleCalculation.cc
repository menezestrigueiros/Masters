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
#include "UTIL/LCRelationNavigator.h"
#include "UTIL/PIDHandler.h"


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

#include "../../../ILDStyle_Standalone.cc"


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
  registerInputCollection(LCIO::RECONSTRUCTEDPARTICLE,
          "IsolatedElectrons",
          "Name of the ReconstructedParticle collection",
          _ElectronColName,
          std::string("IsolatedElectrons"));
  registerInputCollection(LCIO::RECONSTRUCTEDPARTICLE,
          "IsolatedMuons",
          "Name of the ReconstructedParticle collection",
          _MuonColName,
          std::string("IsolatedMuons"));
  registerInputCollection(LCIO::RECONSTRUCTEDPARTICLE,
          "IsolatedTaus",
          "Name of the ReconstructedParticle collection",
          _TauColName,
          std::string("IsolatedTaus"));
  registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
          "IsolatedPhotons",
          "Name of the ReconstructedParticle collection",
          _PhotonColName,
          std::string("IsolatedPhotons"));
  registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
          "OutputCollectionPFOsOverlayCut",
          "Name of the PFOsMinusOverlayCut collection",
          _PFOsWithoutOverlay,
          std::string("PFOsWithoutOverlay"));
  registerInputCollection( LCIO::LCRELATION,
          "RecoMCTruthLink",
          "Name of the RecoMCTruthLink input collection"  ,
          _recoMCTruthLink,
          std::string("RecoMCTruthLink")
          );
  registerInputCollection( LCIO::LCRELATION,
          "MCTruthRecoLink",
          "Name of the MCTruthRecoLink input collection"  ,
          _mcTruthRecoLink,
          std::string("MCTruthRecoLink")
          );
    registerInputCollection(LCIO::RECONSTRUCTEDPARTICLE,
          "PFOsminusphoton",
          "Name of the ReconstructedParticle collection",
          _PFOsMinusIsolatedObjetcs,
          std::string("PFOsminusphoton")
          );
  registerOutputCollection(LCIO::RECONSTRUCTEDPARTICLE,
          "OutputPfoCollection",
          "Name of output muon PFOs collection with ISR particles removed",
          _OutputPFOsCollection,
          std::string("RemainingPFOsWithoutISR")
          );
  registerInputCollection(LCIO::RECONSTRUCTEDPARTICLE,
          "JetOut",
          "Name of the ReconstructedParticle collection",
          _JetOut,
          std::string("JetOut")
          );
  registerInputCollection(LCIO::RECONSTRUCTEDPARTICLE,
          "JetPFOs",
          "Name of the ReconstructedParticle collection",
          _JetPFOs,
          std::string("JetPFOs")
          );
  registerInputCollection(LCIO::RECONSTRUCTEDPARTICLE,
          "JetOut1",
          "Name of the ReconstructedParticle collection",
          _JetOuts,
          std::string("JetOut1")
          );
  registerInputCollection(LCIO::RECONSTRUCTEDPARTICLE,
          "JetPFOs1",
          "Name of the ReconstructedParticle collection",
          _JetPFOss,
          std::string("JetPFOs1")
          );
  registerProcessorParameter("TTreeFileName",
          "Name of the root file in which the TTree with the 4 event observables is stored; if left empty no root file is created; default: WWCategorisationTTreeFile.root",
          _TTreeFileName,
          std::string("WWAngle.root"));       
}


void WWAngleCalculationProcessor::init() { 
  //gStyle->SetPaintTextFormat("4.3f");
  //gStyle->SetPalette(1);

  TStyle* MyStyle = new  TStyle("MyStyle", "My Style");
  ILDStyle(MyStyle);
  MyStyle->cd();
  gROOT->ForceStyle();
  gStyle->ls();
  gStyle->SetOptStat(1111);

  //logarithmic binning for energy vs cosine(theta) for IsoPhotons
  double _minBinX = 0;
  double _maxBinX = 1; 

  double _minBinY = 0;  // for 1 GeV = 10^0 = 1 GeV
  double _maxBinY = 2.18;  // for 150 GeV = 10^2.176 GeV

  int _nBinsY = 300;
  int _nBinsX = 300; // choose a adequate number here

  double* histbinsX  = new double[_nBinsX+1];
  double* histbinsY  = new double[_nBinsY+1];

  for (int i=0; i<_nBinsY+1; i++) histbinsY [i] = pow( 10, _minBinY + (_maxBinY-_minBinY)*i/(float)(_nBinsY));
  for (int i=0; i<_nBinsX+1; i++) histbinsX [i] = _minBinX + (_maxBinX-_minBinX)*i/(float)(_nBinsX);


  if(_TTreeFileName != ""){

    _doTT = true;
    _TTreeFile = new TFile(_TTreeFileName.c_str(), "RECREATE");

    _TTreeFile->mkdir("MCTruthID");
    _TTreeFile->mkdir("RecoID");

    _TTreeFile->cd();
    _WPropertiesTTree = new TTree("WWAngle", "Tree of all WW properties");
    /* gStyle->SetNumberContours(256); */

    _WPropertiesTTree->Branch("nW", &_nW, "nW/F");
    _WPropertiesTTree->Branch("mInv",&_mInv,"mInv/F");
    _WPropertiesTTree->Branch("mInvJets", &_mInvJets, "mInvJets/F");
    _WPropertiesTTree->Branch("charge", &_charge, "charge/F");

    //MC: W+ and W- angles
    _MC_StarPhi_Wplus  = new TH1F(" #phi * of lepton from W^{+}", " ; #phi_{l} * ;nEvents", 100, -3.3, 3.3);
    _MC_StarPhi_Wminus  = new TH1F(" #phi * of lepton from W^{-}", " ; #phi_{l} * ;nEvents", 100, -3.3, 3.3);
    _MC_decayAngle_lep  = new TH1F(" cos(#theta_{#mu^{-}}) *", " ;  cos(#theta_{#mu^{-}}) * ;nEvents", 40, -1.1, 1.1);
    _MC_decayAngle_muplus  = new TH1F(" cos(#theta_{#mu^{+}}) *", " ;  cos(#theta_{#mu^{+}}) * ;nEvents", 40, -1.1, 1.1);
    
    _MC_prodAngle_Wminus = new TH1F("cos(#theta_{W^{-}}) ", " ; cos(#theta_{w^{-}}) ;nEvents", 100, -1.1, 1.1);
    _MC_prodAngle_Wplus = new TH1F("cos(#theta_{W^{+}}) ", " ; cos(#theta_{w^{+}}) ;nEvents", 100, -1.1, 1.1);

    _hTest = new TH1F("explicit W cos(#theta_{W})", " ; cos(#theta_{w}) ; Events", 100, -1.1, 1.1);
    
    // RECO
    _isophotons = new TH2D("isophotons_distribution", " ; cos(#theta) ; Energy [GeV]", _nBinsX, histbinsX, _nBinsY, histbinsY);
    
    //invariant mass 
    _Wminus_mass = new TH1F("reco W- inv mass ", " ; M_{W^{-}}_{inv} [GeV] ; Events", 200, 0., 200.);
    _Wplus_mass = new TH1F("reco W+ inv mass ", " ; M_{W^{+}}_{inv} [GeV] ; Events", 200, 0., 200.);
    _Wmassjets = new TH1F("reco W inv mass from jets no overlay", " ; M_{inv} [GeV] ; Events", 200, 0., 200.);
    _Wmassjetss = new TH1F("reco W inv mass from jets", " ; M_{inv} [GeV] ; Events", 200, 0., 200.);

    _charge_distribution = new TH1F("charge distribution", " ; Charge ; Events", 100, -10., 10);
    //prod angles
    _reco_cos_prodw = new TH1F("Reconstruced cos(#theta_{W-}) ", " ; cos(#theta_{w^{-}}) ;nEvents", 100, -1.1, 1.1);
    _reco_cos_prodwplus = new TH1F("Reconstruced cos(#theta_{W^{+}}) ", " ; cos(#theta_{w^{+}}) ;nEvents", 100, -1.1, 1.1);
    //decay angles
    //w- decay angle
    _reco_startheta = new TH1F("Reconstructed cos(#theta_{#mu^{-}}) *", " ;  cos(#theta_{#mu^{-}}) * ;nEvents", 40, -1.1, 1.1);
    _reco_starphi = new TH1F("Reconstructed #phi * of #mu^{-} from W^{-}", " ; #phi_{#mu^{-}} * ;nEvents", 100, -3.3, 3.3);
    //w+ decay angle
    _reco_startheta_plus = new TH1F("Reconstructed cos(#theta_{#mu^{+}}) *", " ;  cos(#theta_{#mu^{+}}) * ;nEvents", 40, -1.1, 1.1);
    _reco_starphi_plus = new TH1F("Reconstructed #phi * of #mu^{+} from W^{+}", " ; #phi_{#mu^{+}} * ;nEvents", 100, -3.3, 3.3);
  }


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

double WWAngleCalculationProcessor::Azimuth(std::vector<double>& v1){
  double phi = atan2(v1[0], v1[1]);
  return phi;
}

void WWAngleCalculationProcessor::processEvent( LCEvent * evt ){ 



  _nEvt ++;

  //extract categorisation parameters
  LCParameters& pars = evt->parameters();
  auto parTrueCat = pars.getIntVal("WWCategorisation.TrueCat");
  auto parRecoCatAdvanced = pars.getIntVal("WWCategorisation.RecoCatAdvanced");  

  streamlog_out(DEBUG) 
	<< " processing event " << evt->getEventNumber() 
	<< "  in run "          << evt->getRunNumber() 
	<< std::endl ;

  LCCollection *col_muon{}, *col_electron{}, *col_tau{}, *col_mcparticles{}, *col_pfo{}, *col_photons{}, *col_recoMCTruth{}, *col_pfominusphoton{}, *col_pfooverlaycut{}, *col_jets{}, *col_jetpfos{}, *col_jets1{}, *col_jetpfos1{};

  IMPL::LCCollectionVec* OutputPFOsCollection(NULL);
  OutputPFOsCollection = new IMPL::LCCollectionVec(LCIO::RECONSTRUCTEDPARTICLE);
  OutputPFOsCollection->setSubset( true ); // This is a subset of a collection
  
  // fill histogram from LCIO data :
  try{col_muon    = evt->getCollection(_MuonColName);
  col_electron    = evt->getCollection(_ElectronColName);
  col_tau         = evt->getCollection(_TauColName);
  col_mcparticles = evt->getCollection(_MCParColName);
  col_pfo         = evt->getCollection(_PandoraPFOsCol);
  col_photons     = evt->getCollection(_PhotonColName);
  col_pfominusphoton = evt->getCollection(_PFOsMinusIsolatedObjetcs);
  col_recoMCTruth = evt->getCollection(_recoMCTruthLink);
  col_pfooverlaycut = evt->getCollection(_PFOsWithoutOverlay);
  col_jets = evt->getCollection(_JetOut);
  col_jetpfos = evt->getCollection(_JetPFOs);

  /* col_jets1 = evt->getCollection(_JetOuts);
  col_jetpfos1 = evt->getCollection(_JetPFOss); */

  //LCRelationNavigator RecoMCParticleNav( col_recoMCTruth );
  LCRelationNavigator MCParticleRecoNav( evt->getCollection( _mcTruthRecoLink ) );

  }
  catch(DataNotAvailableException &e){
    streamlog_out(MESSAGE) << "Input collections not found - skipping event " << _nEvt << std::endl;
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

  streamlog_out(DEBUG) << "OverlayCut jetPFOs (subset): " << col_pfooverlaycut->getNumberOfElements() << endl;
  streamlog_out(DEBUG) << "Photon jetPFOs: " << col_pfominusphoton->getNumberOfElements() << endl;


  ////////// W BOSONS ///////////////
  PxPyPzEVector lv_W1(0, 0, 0, 0); //Explicit W
  PxPyPzEVector lv_lW(0, 0, 0, 0); // Non Explicit Leptonic W
  PxPyPzEVector lv_hW(0, 0, 0, 0); // Non Explicit hadronic W

  ////////// LEPTONS/QUARKS //////////
  PxPyPzEVector lv_leptons_w1(0, 0, 0, 0); //leptons from explicit W
  PxPyPzEVector lv_leptons(0, 0, 0, 0); //leptons from reconstructed W
  PxPyPzEVector lv_quarks(0, 0, 0, 0); //quarks from reconstructed W

  PxPyPzEVector transf_lv(0, 0, 0, 0); //general transformed lorentz vector

  PxPyPzEVector photon_lv(0, 0, 0, 0); //general transformed lorentz vector of photons



  _nW.push_back(0);
  _nNonW.push_back(0);
 

  std::vector<int> general_daughters_pdgs = {}; //array with all daughters from e+e-
  std::vector<int> boson_daughters_pdg_array{}; //daughters that come from a W or Z boson

  std::vector<double> leptonic_w(4); //0,0,0,0
  std::vector<double> hadronic_w(4); //0,0,0,0
  std::vector<double> vector_w(4);
  std::vector<double> boost(4);
/*   std::vector<double> isophotons(4);*/
  std::vector<MCParticle*> boson_daughters = {};

  ROOT::Math::XYZVector v1(0, 0, 1);

  MCParticle* quarks; //MC information of each particle in hadronic decay
  MCParticle* leptons; //MC information of each particle in a leptonic decay

  // TEST ATTEMPTS************** For every explicit W (PDG == 24)
  if(parTrueCat == 3 && parRecoCatAdvanced == 3){
    {
      for(int i = 0; i < col_mcparticles->getNumberOfElements(); i++){
        MCParticle* mc1 = static_cast<MCParticle*> (col_mcparticles->getElementAt(i));
        if ( mc1->getPDG() != -24 ) continue;
        _hTest->Fill(mc1->getMomentum()[2] / std::hypot(mc1->getMomentum()[0], mc1->getMomentum()[1], mc1->getMomentum()[2]));
      }
    }

    howmany_muons++;
  
    if(incoming_pdg == 11){
      if(daughters.size() <= 3){
        streamlog_out(DEBUG) << " ////////// [2/3 DAUGHTERS] in _nEvt " << _nEvt << " //////////" << endl;
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
                streamlog_out(DEBUG) << "[ON SHELL] W 4v in " << _nEvt << ": [ "; 
                for (double k = 0; k < 4; k++){
                  streamlog_out(DEBUG) << vector_w[k] << ", ";
                }
                streamlog_out(DEBUG) << "]" << endl;
                streamlog_out(DEBUG) << "[ON SHELL] W CHARGE in _nEvt " << _nEvt << ": " << wcharge << endl;
              }

              //////////////////////// W- SELECTION //////////////////
              double prod_cosine = CosineTheta(vector_w);
              /* if(wcharge == -1) _MC_prodAngle_Wminus->Fill(prod_cosine); */

              //////////// LEPTONIC DECAY ///////////
              boson_daughters = wboson->getDaughters();
              for ( long unsigned int l = 0; l < boson_daughters.size(); l++){
                MCParticle* lepton_from_w = boson_daughters[l];
                int lepton_pdg = abs(boson_daughters.at(l)->getPDG());

                if(lepton_pdg == 11 || lepton_pdg == 13 || lepton_pdg == 15){ //if it's a electron, muon or tau
                  lv_leptons_w1 = PxPyPzEVector(lepton_from_w->getMomentum()[0], lepton_from_w->getMomentum()[1], lepton_from_w->getMomentum()[2], lepton_from_w->getEnergy());
                  

                  /////////// DECAY ANGLES //////////
                 
                  transf_lv = starVector(lv_W1, lv_leptons_w1, v1);
                  transf_lv.GetCoordinates(boost.begin(), boost.end());
                  double cos_star_theta = CosineTheta(boost);

                  /////// W- or W+ //////
                  if(wcharge == -1){
                    _Wleptonic++; //explicit W- that decay leptonically
                    _MC_prodAngle_Wminus->Fill(prod_cosine); //prod angle of W-
                    _MC_decayAngle_lep->Fill(cos_star_theta);
                    _MC_StarPhi_Wminus->Fill(transf_lv.Phi());
                  }
                  if(wcharge == 1){
                    _Wleptonic1++; //explicit W+ that decay leptonically
                    _MC_prodAngle_Wplus->Fill(-prod_cosine); //prod angle of W+
                    _MC_decayAngle_muplus->Fill(cos_star_theta);
                    _MC_StarPhi_Wplus->Fill(transf_lv.Phi());
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
            streamlog_out(DEBUG) << "[OFFSHELL] W 4v in " << _nEvt << ": [ "; 
            for (double k = 0; k < 4; k++){
              streamlog_out(DEBUG) << hadronic_w[k] << ", ";
            }
            streamlog_out(DEBUG) << "]" << endl;
            streamlog_out(DEBUG) << "[OFFSHELL] hadronic W CHARGE in _nEvt " << _nEvt << ": " << hWcharge << endl;
          }

          /////// W- or W+ //////
          if(hWcharge == -1){
            _Whadronic++; // non explicit W- that decay hadronically
            double prod_cosine = CosineTheta(hadronic_w);
           /*  _MC_prodAngle_Wminus->Fill(prod_cosine) */;     
          }
          if(hWcharge == 1) _Whadronic1++ ; // non explicit W- that decay hadronically
        }
        if(leptonic_w[0] != 0){
          _W++;
          { //for debugging
            streamlog_out(DEBUG) << "[OFFSHELL] W 4v in " <<_nEvt << ": [ "; 
            for (double k = 0; k < 4; k++){
              streamlog_out(DEBUG) << leptonic_w[k] << ", ";
            }
            streamlog_out(DEBUG) << "]" << endl;
            streamlog_out(DEBUG) << "[OFFSHELL] W LEPTON CHARGE in _nEvt " << _nEvt << ": " << lWcharge << endl;
          }

          //Invariant Mass
          double inv_mass = lv_lW.M();
          streamlog_out(DEBUG) << "[OFFSHELL] W INV MASS " << inv_mass << endl;
          
          /////////// W REST FRAME //////////
          transf_lv = starVector(lv_lW, lv_leptons, v1);
          transf_lv.GetCoordinates(boost.begin(), boost.end());

          { //for debbuging
            streamlog_out(DEBUG) << "[LEP] transformed 4v in " << _nEvt << ": [ "; 
            for (double k = 0; k < 4; k++){
              streamlog_out(DEBUG) << boost[k] << ", ";
            }
            streamlog_out(DEBUG) << "]" << endl;
          }
          
          /////////// W PRODUCTION ANGLE //////////
          double prod_cosine = CosineTheta(leptonic_w);

          /////// W- or W+ //////
          if(lWcharge == -1){
            _Wleptonic++; 
            /////////// DECAY ANGLES //////////
            transf_lv.GetCoordinates(boost.begin(), boost.end());
            double cos_star_theta = CosineTheta(boost);

            /////////// PLOTS //////////
            _MC_prodAngle_Wminus->Fill(prod_cosine); //W- production angle
            _MC_decayAngle_lep->Fill(cos_star_theta); // lepton decay angle
            _MC_StarPhi_Wminus->Fill(transf_lv.Phi()); //lepton azimuthal angle
          }
          if(lWcharge == 1){
            _Wleptonic1++;

            /////////// DECAY ANGLES //////////
            transf_lv.GetCoordinates(boost.begin(), boost.end());
            double cos_star_theta = CosineTheta(boost);

            _MC_prodAngle_Wplus->Fill(-prod_cosine); //W+ production angle
            _MC_decayAngle_muplus->Fill(cos_star_theta);
            _MC_StarPhi_Wplus->Fill(transf_lv.Phi());
          }
        }
      
      }else if (daughters.size() == 4){ //non-explicit W procedure
        streamlog_out(DEBUG) << "/////////////// [4 DAUGHTERS] in _nEvt " << _nEvt << " ///////////////" << endl;
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
            streamlog_out(DEBUG) << "[4 DAUGHTERS] W 4v in " <<_nEvt << ": [ "; 
            for (double k = 0; k < 4; k++){
              streamlog_out(DEBUG) << leptonic_w[k] << ", ";
            }
            streamlog_out(DEBUG) << "]" << endl;
            streamlog_out(DEBUG) << "[CHARGE] LEPTONIC W in _nEvt " << _nEvt << ": " << lWcharge << endl;
          }
          
          double inv_mass = lv_lW.M();
          streamlog_out(DEBUG) << "[4 DAUGHTERS] W INV MASS: " << inv_mass << endl;

          /////////// W REST FRAME //////////
          transf_lv = starVector(lv_lW, lv_leptons, v1);
          transf_lv.GetCoordinates(boost.begin(), boost.end()); //W REST FRAME 

          { //for debugging
            streamlog_out(DEBUG) << "[LEPTON] TRANSF 4v in " << _nEvt << ": [ "; 
            for (double k = 0; k < 4; k++){
              streamlog_out(DEBUG) << boost[k] << ", ";
            }
            streamlog_out(DEBUG) << "]" << endl;
          }

          /////////// W PRODUCTION ANGLE //////////
          double prod_cosine = CosineTheta(leptonic_w);

          /////// W- or W+ //////
          if(lWcharge == -1){
            _Wleptonic++;

            /////////// DECAY ANGLES //////////
            transf_lv.GetCoordinates(boost.begin(), boost.end()); 
            double cos_star_theta = CosineTheta(boost);

            ///////// PLOTS //////////
            _MC_prodAngle_Wminus->Fill(prod_cosine); //W- production angle
            _MC_decayAngle_lep->Fill(cos_star_theta);
            _MC_StarPhi_Wminus->Fill(transf_lv.Phi());
          }
          if(lWcharge == 1){
            _Wleptonic1++;

            /////////// DECAY ANGLES //////////
            transf_lv.GetCoordinates(boost.begin(), boost.end()); 
            double cos_star_theta = CosineTheta(boost);

            _MC_prodAngle_Wplus->Fill(-prod_cosine); //prod angle of W+
            _MC_decayAngle_muplus->Fill(cos_star_theta);
            _MC_StarPhi_Wplus->Fill(transf_lv.Phi());
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
            /* _MC_prodAngle_Wminus->Fill(prod_cosine);  */    
          }
          if(hWcharge == 1) _Whadronic1++ ;
          
        } 
        if (lWcharge != 0)streamlog_out(DEBUG) << "[CHARGE] LEPTONIC W in _nEvt " << _nEvt << ": " << lWcharge << endl;
        if (hWcharge != 0)streamlog_out(DEBUG) << "[CHARGE] HADRONIC W in _nEvt " << _nEvt << ": " << hWcharge << endl; 
      }
    }
  }


  //////////////// RECO //////////////

  std::vector<double> recoboost(4);

  // get number of isolated leptons and photons -> used for basic categorisation
  int n_muons = col_muon->getNumberOfElements();
  int n_photons = col_photons->getNumberOfElements();
  int n_electrons = col_electron->getNumberOfElements();
  int n_taus = col_tau->getNumberOfElements();
  int n_pfos = col_pfo->getNumberOfElements();

  LCCollection* cols[] = {col_muon, col_photons};
  int n_parts[] = { n_muons, n_photons};

  LCCollection* cols1[] = {col_electron, col_muon, col_tau, col_photons};
  int n_parts1[] = {n_electrons, n_muons, n_taus, n_photons};
  PxPyPzEVector lv_IsoCumul(0, 0, 0, 0);  // cumulative 4-momentum for all isolated particles


  if(parTrueCat == 3 && parRecoCatAdvanced == 3){

    PxPyPzEVector lv_all_iso_pfos(0, 0, 0, 0); //4-momentum for all isolated leptons
    PxPyPzEVector IsolatedPhotons(0, 0, 0, 0);

    //////////////////////////////// for reco ////////////////////////////////
    
    double recow_charge = 0;
    double cosine_photons = 0;
    double chargephotons = 0;
    double event_charge = 0;

    // loop through isolated muons and photons
    double highest_momentum = 0;
    double highest_momentum_charge = 0;

    std::vector<double> isomuon(4);
    PxPyPzEVector lv_iso_muon(0, 0, 0, 0);
    //charge summation and getting lorentz vector of muon with highest momentum.
    for(int i = 0; i < 2; ++i){ 
      for (int j = 0; j < n_parts[i]; ++j){
        EVENT::ReconstructedParticle* isopar = dynamic_cast <EVENT::ReconstructedParticle*>(cols[i]->getElementAt(j));
        if (i == 0){
          double momtot = sqrt(isopar->getMomentum()[0]*isopar->getMomentum()[0] + isopar->getMomentum()[1]*isopar->getMomentum()[1] + isopar->getMomentum()[2]*isopar->getMomentum()[2]);
          if(momtot > highest_momentum){
            highest_momentum = momtot;
            highest_momentum_charge = isopar->getCharge();
            lv_iso_muon = PxPyPzEVector(isopar->getMomentum()[0], isopar->getMomentum()[1], isopar->getMomentum()[2], isopar->getEnergy());
          }
          recow_charge = highest_momentum_charge;
          
        } else{
          chargephotons += isopar->getCharge();
        }
      }
    }
    //get lorentz vector of muon with highest momentum
    lv_iso_muon.GetCoordinates(isomuon.begin(), isomuon.end());

    for (int i = 0; i < 4; ++i){  // loop through isolated species
    for (int j = 0; j < n_parts1[i]; ++j){ // loop through all isolated particles
    EVENT::ReconstructedParticle* par = dynamic_cast <EVENT::ReconstructedParticle*>(cols1[i]->getElementAt(j));
    lv_IsoCumul += PxPyPzEVector(par->getMomentum()[0], par->getMomentum()[1], par->getMomentum()[2], par->getEnergy()); // sum 4-mom for all particles
    }
  }  

    event_charge = recow_charge + chargephotons;
    //debugging purposes -> two muons on a event would give total charge sum diff from -1 and 1. 
    if(n_muons == 2) _muon_count++;
    _charge = event_charge;

  
    //check with RecoMCTruthLink wether this picks actual ISR most of the time
    // Create an LCRelationNavigator for the LCRelations between ReconstructedParticle and MCParticle
    LCRelationNavigator pfo2mc(col_recoMCTruth);

    std::vector<double> isophotons(4);
    PxPyPzEVector lv_isopho(0, 0, 0, 0);

    // Loop through all isolated photons
    for (int l = 0; l < n_parts[1]; ++l) {
      EVENT::ReconstructedParticle* isopho = dynamic_cast<EVENT::ReconstructedParticle*>(cols[1]->getElementAt(l));
      lv_isopho += PxPyPzEVector(isopho->getMomentum()[0], isopho->getMomentum()[1], isopho->getMomentum()[2], isopho->getEnergy());
      lv_isopho.GetCoordinates(isophotons.begin(), isophotons.end());
      cosine_photons = CosineTheta(isophotons);


      // Get the related MCParticles for the current ReconstructedParticle
      const EVENT::LCObjectVec& mcParticles = pfo2mc.getRelatedToObjects(isopho);

      // Loop through the related MCParticles
      for (const auto& mcParticle : mcParticles) {
        // Do something with the MCParticle
        // For example, you can access its properties using dynamic_cast
        EVENT::MCParticle* mc1 = dynamic_cast<EVENT::MCParticle*>(mcParticle);
        if (mc1) {
          totalphotons++;
          // Access the properties of the MCParticle
          int pdg = mc1->getPDG();
          if(pdg == 22){
            actualphotons++;
            //continue; //take this out of the collection
          }
        }
      }
    }
    

    //remove muon and - if found, the ISR candidate  - from the PFO list
    //not sure if this is working as intended... when I try to cluster w / FastJet the PFOs, it gives me errors... skipping this atm
    for(int m = 0; m < n_pfos; ++m){
      EVENT::ReconstructedParticle* pfo = dynamic_cast<EVENT::ReconstructedParticle*>(col_pfo->getElementAt(m));
      LCRelationNavigator pfo2mc(col_recoMCTruth);
      const EVENT::LCObjectVec& mcParticle = pfo2mc.getRelatedToObjects(pfo);
      if(mcParticle.empty()) continue;
      EVENT::MCParticle* mc2 = dynamic_cast<EVENT::MCParticle*>(mcParticle[0]);
      int pfo_type = abs(mc2->getPDG());
      if(pfo_type == 13 || (pfo_type == 22 && (cosine_photons > 0.9 && isophotons[3] > 10))){
        continue; // Skip muons and ISR photons
      }
      
      // Process the remaining particles
      // ...
      OutputPFOsCollection->addElement(pfo);
    }


    //this is a new collection with only jet pfos
    evt->addCollection(OutputPFOsCollection, _OutputPFOsCollection.c_str());

    // new collection
    const EVENT::LCCollection* PFOs = evt->getCollection(_OutputPFOsCollection);

    //After clustering, we have the jet pfos
    PxPyPzEVector jet4v(0, 0, 0, 0);
    std::vector<double> jetv(4);
    std::vector<double> jetfv(4);
    int jets = 0;
    double jet_charge = 0;
    for(int j = 0; j < col_jets->getNumberOfElements() ; ++j){
      jets++;
      EVENT::ReconstructedParticle* jet = dynamic_cast<EVENT::ReconstructedParticle*>(col_jets->getElementAt(j));
      jet4v += PxPyPzEVector(jet->getMomentum()[0], jet->getMomentum()[1], jet->getMomentum()[2], jet->getEnergy());
    }

    for (int h = 0; h < col_jetpfos->getNumberOfElements(); h++){
      EVENT::ReconstructedParticle* jetpfo = dynamic_cast<EVENT::ReconstructedParticle*>(col_jetpfos->getElementAt(h));
      jet_charge += jetpfo->getCharge();
    }
    jet4v.GetCoordinates(jetv.begin(), jetv.end());
    
    // inv mass jets 
    double inv_massjet = jet4v.M(); // inv mass of jets pfosminusphotos minus overlay removal
    streamlog_out(DEBUG) << "NJets: " << jets << std::endl;
    _charge_distribution->Fill(jet_charge);

   
    
    //invariant mass
    _mInvJets = jet4v.M();
    _Wmassjets->Fill(abs(inv_massjet));

    //4-momentum for all PFOsminusphotons (except those deemed overlay from the categorisation)
    PxPyPzEVector lv_all_pfosmp(0, 0, 0, 0);
    for(int r = 0; r < col_pfooverlaycut->getNumberOfElements(); ++r){
      EVENT::ReconstructedParticle* par = dynamic_cast <EVENT::ReconstructedParticle*>(col_pfooverlaycut->getElementAt(r));
      lv_all_pfosmp += PxPyPzEVector(par->getMomentum()[0], par->getMomentum()[1], par->getMomentum()[2], par->getEnergy());
    }
    std::vector<double> reco_hadronic_w(4);
    lv_all_pfosmp.GetCoordinates(reco_hadronic_w.begin(), reco_hadronic_w.end());


    // take into account beam crossing angle, get missing momentum and energy
    PxPyPzEVector initial_conditions = PxPyPzEVector(1.75, 0, 0, 250); 
    
    //initial - all pfos (IsoParticles + JetPFOs --> excluding overlay)
    PxPyPzEVector missing_neutrino4v_general = initial_conditions - (lv_all_pfosmp + lv_IsoCumul);
    //now we have W boson 4vector
    PxPyPzEVector reco_wboson = missing_neutrino4v_general + lv_iso_muon;

    
    //4v leptonically decayed W
    std::vector<double> reco_w(4);
    reco_wboson.GetCoordinates(reco_w.begin(), reco_w.end());

    //invariant mass
    double inv_mass = reco_wboson.M();
    _mInv = reco_wboson.M();

    PxPyPzEVector reco_transf_lv(0, 0, 0, 0); //general transformed lorentz vector

    /////////// W REST FRAME //////////
    
    reco_transf_lv = starVector(reco_wboson, lv_iso_muon, v1);
    reco_transf_lv.GetCoordinates(recoboost.begin(), recoboost.end()); 
    double beta = lv_iso_muon.Beta();

    streamlog_out(MESSAGE) << "BOOSTED VECTOR BETA: " << beta << std::endl;

    //leptonically decayed W
    if(reco_w[0] != 0){

      double cos_star_theta = CosineTheta(reco_w);
      //double muonstarphi = reco_transf_lv.Phi();
      //double muonstarphi = Azimuth(reco_w);

      if(event_charge == -1){
        _wminus++;
        _Wminus_mass->Fill(abs(inv_mass));

        /////////// W PRODUCTION ANGLE //////////
        _reco_cos_prodw->Fill(cos_star_theta);

        /////////// DECAY ANGLES //////////
        reco_transf_lv.GetCoordinates(recoboost.begin(), recoboost.end()); 

        if(recoboost[1] == 0){ //for debbuging of zero-ed four-momentas


          streamlog_out(MESSAGE) << "[RECO W] 4v in " << _nEvt << ": [ "; 
          for (double i = 0; i < 4; i++){
            streamlog_out(MESSAGE) << reco_w[i] << ", ";
          }
          streamlog_out(MESSAGE) << "]" << endl;

          streamlog_out(MESSAGE) << "[MUON] 4v in " << _nEvt << ": [ "; 
          for (double i = 0; i < 4; i++){
            streamlog_out(MESSAGE) << isomuon[i] << ", ";
          }
          streamlog_out(MESSAGE) << "]" << endl;

          streamlog_out(MESSAGE) << "[TRANSFORMED MUON] transformed 4v in " << _nEvt << ": [ "; 
          for (double k = 0; k < 4; k++){
            streamlog_out(MESSAGE) << recoboost[k] << ", ";
          }
          streamlog_out(MESSAGE) << "]" << endl;
        }
        
        // polar angle
        double cosprod_mu = CosineTheta(recoboost);
        _reco_startheta->Fill(cosprod_mu);

        // azimuthal angle and folding
        double muonstarphi = reco_transf_lv.Phi();     
        _reco_starphi->Fill(muonstarphi);
      }
      if(event_charge == 1){
        _not_wminus++;

        _Wplus_mass->Fill(abs(inv_mass));

        /////////// W PRODUCTION ANGLE //////////
        _reco_cos_prodwplus->Fill(-cos_star_theta);

        /////////// DECAY ANGLES //////////
        reco_transf_lv.GetCoordinates(recoboost.begin(), recoboost.end()); 

        // polar angle
        double cosprod_mu = CosineTheta(recoboost);
        _reco_startheta_plus->Fill(cosprod_mu);

        // azimuthal angle and folding
        double muonstarphi = reco_transf_lv.Phi();
        _reco_starphi_plus->Fill(muonstarphi);

      }
    }
    //hadronically decayed W
    if(jetv[0] != 0){
      if(jet_charge == -1) _wplus++;   
    }
  }

  // store in TTree and/or confusion matrix
  if (_doTT) _WPropertiesTTree->Fill();


  _nW[_W]++;
  _nNonW[_nonW]++;
}


void WWAngleCalculationProcessor::check( LCEvent *  /*evt*/ ) { 
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}


void WWAngleCalculationProcessor::end(){ 

  //for debbugging
  TCanvas* canvas = new TCanvas("Energy vs CosTheta semileptonic #mu", "Energy vs SinTheta semileptonic #mu");
  Angular(canvas, _isophotons);
  TCanvas* canvasone = new TCanvas;
  InvMass(canvasone, _Wmassjets, _Wmassjetss);


 if (_doTT){
    _TTreeFile->cd();
    _WPropertiesTTree->Write();
  }
  
  //MC
  _TTreeFile->cd("MCTruthID");
  _MC_prodAngle_Wminus->Write();
  _MC_prodAngle_Wplus->Write();
  _MC_decayAngle_lep->Write();
  _MC_decayAngle_muplus->Write();
  _MC_StarPhi_Wplus->Write();
  _MC_StarPhi_Wminus->Write();

  _hTest->Write();

  //RECO
  _TTreeFile->cd("RecoID");

  _Wminus_mass->Write();
  _Wplus_mass->Write();

  _reco_cos_prodw->Write();
  _reco_cos_prodwplus->Write();
  _isophotons->Write();

  _reco_starphi->Write();
  _reco_starphi_plus->Write();

  canvas->Write();
  canvasone->Write();

  _reco_startheta->Write();
  _reco_startheta_plus->Write();
  _charge_distribution->Write();
  
  if (_doTT) _TTreeFile->Close();
  
  //double ratio = actualphotons/totalphotons;
  streamlog_out(DEBUG) << "Ratio of actual photons: " << actualphotons/totalphotons * 100.0 << std::endl;
  streamlog_out(DEBUG) << "actual photons: " << actualphotons << std::endl;
  streamlog_out(DEBUG) << "total photons: " << totalphotons << std::endl;
  //MC DEBBUGGING
  streamlog_out(MESSAGE) << "MC W- that decay hadronically: " << _Whadronic << std::endl; 
  streamlog_out(MESSAGE) << "MC W- that decay leptonically: " << _Wleptonic << std::endl; 
  streamlog_out(MESSAGE) << "MC W+ that decay hadronically: " << _Whadronic1 << std::endl; 
  streamlog_out(MESSAGE) << "MC W+ that decay leptonically: " << _Wleptonic1 << std::endl; 
  streamlog_out(MESSAGE) << "MC TOTAL W-: " << _Wleptonic + _Whadronic << std::endl; 
  streamlog_out(MESSAGE) << "MC TOTAL W+: " << _Wleptonic1 + _Whadronic1 << std::endl; 
  streamlog_out(MESSAGE) << "MC: TOTAL W BOSONS:" << _W << std::endl; 
  streamlog_out(MESSAGE) << "MC: How many muons channel? --> " << howmany_muons << std::endl;

  //RECO DEBBUGGING
  streamlog_out(MESSAGE) << "RECO: W- leptonically: " << _wminus << std::endl;
  streamlog_out(MESSAGE) << "RECO: W+ leptonically: " << _not_wminus << std::endl;

  streamlog_out(MESSAGE) << "RECO: W- hadronically: " << _wplus << std::endl;

  streamlog_out(MESSAGE) << "Two Muons: " << _muon_count << std::endl;
  

  streamlog_out(MESSAGE) << "In event " << _nEvt << ": Kept MuPFOs = " << nKeptMuons << " VS removed number of PFOs = " << actualphotons << " VS total number of MuPFOs = " << howmany_muons << std::endl;


  

  _TTreeFile->Close();
  // 	    << " processed " << _nEvt << " events in " << _nRun << " runs "
  // 	    << std::endl ;
}

void WWAngleCalculationProcessor::Angular( TCanvas* canvas, TH2* histo){
  TH2D* histo1 = (TH2D*)histo->Clone();
  canvas->SetLogy();
  canvas->SetLogz();
  
  histo1->SetMinimum(0);
  histo1->Draw("COLZ");
  

  canvas->Update();
  std::stringstream s; s << "./" << histo1->GetName() << ".pdf";
  canvas->Print(s.str().c_str());
}
void WWAngleCalculationProcessor::InvMass( TCanvas* canvas, TH1* histo1, TH1* histo2){
  TH1F* histo = (TH1F*)histo1->Clone();
  TH1F* histos = (TH1F*)histo2->Clone();

  histo->SetLineColor(kRed);
  histos->SetLineColor(kBlue);

  histo->Draw();
  histos->Draw("SAME");

  canvas->Update();
}
