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

    _WPropertiesTTree->Branch("mInv",&_mInv,"mInv/F");
    _WPropertiesTTree->Branch("mInvJets", &_mInvJets, "mInvJets/F");
    _WPropertiesTTree->Branch("charge", &_charge, "charge/F");

    //MC: W+ and W- angles
    _MC_StarPhi_Wplus  = new TH1F(" #phi * of lepton from W^{+}", " ; #phi_{l} * ; Events", 100, -3.3, 3.3);
    _MC_StarPhi_Wminus  = new TH1F(" #phi * of lepton from W^{-}", " ; #phi_{l} * ; Events", 100, -3.3, 3.3);
    _MC_decayAngle_lep  = new TH1F(" cos(#theta_{lepton}) *", " ;  cos(#theta_{lepton}) * ; Events", 40, -1.1, 1.1);

    _MC_prodAngle_Wminus = new TH1F("cos(#theta_{W^{-}}) ", " ; cos(#theta_{w^{-}}) ;nEvents", 100, -1.1, 1.1);
    _MC_prodAngle_Wplus = new TH1F("cos(#theta_{W^{+}}) ", " ; cos(#theta_{w^{+}}) ;nEvents", 100, -1.1, 1.1);

    _hTest = new TH1F("explicit W cos(#theta_{W})", " ; cos(#theta_{w}) ; Events", 100, -1.1, 1.1);
    
    // RECO
    _isophotons = new TH2D("isophotons_distribution", " ; cos(#theta) ; Energy [GeV]", _nBinsX, histbinsX, _nBinsY, histbinsY);
    
    //invariant mass 
    _Wminus_mass = new TH1F("reco W- inv mass ", " ; M_{W^{-}}_{inv} [GeV] ; Events", 200, 0., 200.);
    _Wplus_mass = new TH1F("reco W+ inv mass ", " ; M_{W^{+}}_{inv} [GeV] ; Events", 200, 0., 200.);
    _Wmassjets = new TH1F("reco W inv mass from jets", " ; M_{inv} [GeV] ; Events", 200, 0., 200.);

    _charge_distribution = new TH1F("charge distribution", " ; Charge ; Events", 100, -10., 10);
    //prod angles
    _reco_cos_prodw = new TH1F("Reconstruced cos(#theta_{W-}) ", " ; cos(#theta_{w^{-}}) ;nEvents", 100, -1.1, 1.1);
    _reco_cos_prodwplus = new TH1F("Reconstruced cos(#theta_{W^{+}}) ", " ; cos(#theta_{w^{+}}) ;nEvents", 100, -1.1, 1.1);
    //decay angles
    //w- decay angle
    _reco_startheta = new TH1F("Reconstructed cos(#theta_{#mu^{-}}) *", " ;  cos(#theta_{#mu^{-}}) * ;nEvents", 100, -1.1, 1.1);
    _reco_starphi = new TH1F("Reconstructed #phi * of #mu^{-} from W^{-}", " ; #phi_{#mu^{-}} * ;nEvents", 100, -3.3, 3.3);
    //w+ decay angle
    _reco_startheta_plus = new TH1F("Reconstructed cos(#theta_{#mu^{+}}) *", " ;  cos(#theta_{#mu^{+}}) * ;nEvents", 100, -1.1, 1.1);
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

  LCCollection *col_muon{}, *col_electron{}, *col_tau{}, *col_mcparticles{}, *col_pfo{}, *col_photons{}, *col_recoMCTruth{}, *col_pfominusphoton{}, *col_pfooverlaycut{}, *col_jets{}, *col_jetpfos{};

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

  //LCRelationNavigator RecoMCParticleNav( col_recoMCTruth );
  LCRelationNavigator MCParticleRecoNav( evt->getCollection( _mcTruthRecoLink ) );

  }
  catch(DataNotAvailableException &e){
    streamlog_out(MESSAGE) << "Input collections not found - skipping event " << _nEvt << std::endl;
    return;
  }

  


  ROOT::Math::XYZVector v1(0, 0, 1);

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
    double inv_massjet = jet4v.M();
    streamlog_out(DEBUG) << "NJets: " << jets << std::endl;
    _charge_distribution->Fill(jet_charge);
    
    //invariant mass
    _mInvJets = jet4v.M();
    _Wmassjets->Fill(abs(inv_massjet));


    //4-momentum for all PFOsminusphotons (except those deemed overlay from the categorisation)
    PxPyPzEVector lv_all_pfosmp(0, 0, 0, 0);
    for(int r = 0; r < col_pfominusphoton->getNumberOfElements(); ++r){
      EVENT::ReconstructedParticle* par = dynamic_cast <EVENT::ReconstructedParticle*>(col_pfominusphoton->getElementAt(r));
      lv_all_pfosmp += PxPyPzEVector(par->getMomentum()[0], par->getMomentum()[1], par->getMomentum()[2], par->getEnergy());
    }
    std::vector<double> reco_hadronic_w(4);
    lv_all_pfosmp.GetCoordinates(reco_hadronic_w.begin(), reco_hadronic_w.end());


    // take into account beam crossing angle, get missing momentum and energy
    PxPyPzEVector initial_conditions = PxPyPzEVector(1.75, 0, 0, 250); 
    
    //initial - all pfos (IsoParticles + JetPFOs --> excluding overlay)
    PxPyPzEVector missing_neutrino4v_general = initial_conditions - (lv_all_pfosmp + lv_IsoCumul) ;
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

         //for debbuging


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

          streamlog_out(MESSAGE) << "[TRANSFORMED MUON] 4v in " << _nEvt << ": [ "; 
          for (double k = 0; k < 4; k++){
            streamlog_out(MESSAGE) << recoboost[k] << ", ";
          }
          streamlog_out(MESSAGE) << "]" << endl;
        
        
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
        _reco_cos_prodwplus->Fill(cos_star_theta);

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


}


void WWAngleCalculationProcessor::check( LCEvent *  /*evt*/ ) { 
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}


void WWAngleCalculationProcessor::end(){ 

  //for debbugging
  TCanvas* canvas = new TCanvas("Energy vs CosTheta semileptonic #mu", "Energy vs SinTheta semileptonic #mu");
  Angular(canvas, _isophotons);


 if (_doTT){
    _TTreeFile->cd();
    _WPropertiesTTree->Write();
  }
  
  //MC
  _TTreeFile->cd("MCTruthID");
  _MC_prodAngle_Wminus->Write();
  _MC_prodAngle_Wplus->Write();
  _MC_decayAngle_lep->Write();
  _MC_StarPhi_Wplus->Write();
  _MC_StarPhi_Wminus->Write();

  _hTest->Write();

  //RECO
  _TTreeFile->cd("RecoID");
  
  _Wminus_mass->Write();
  _Wplus_mass->Write();

  _Wmassjets->Write();
  _reco_cos_prodw->Write();
  _isophotons->Write();

  _reco_starphi->Write();
  _reco_starphi_plus->Write();

  canvas->Write();
  _reco_startheta->Write();
  _reco_startheta_plus->Write();
  _charge_distribution->Write();
  
  if (_doTT) _TTreeFile->Close();
  
  //double ratio = actualphotons/totalphotons;
 /*  streamlog_out(DEBUG) << "Ratio of actual photons: " << actualphotons/totalphotons * 100.0 << std::endl;
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
  streamlog_out(MESSAGE) << "MC: How many muons channel? --> " << howmany_muons << std::endl; */

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
