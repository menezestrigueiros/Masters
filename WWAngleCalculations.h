#ifndef WWAngleCalculationProcessor_h
#define WWAngleCalculationProcessor_h 1

#include "marlin/Processor.h"
#include "../include/WWTools.h"
#include "lcio.h"
#include <string>
#include <vector>

#include "IMPL/MCParticleImpl.h" 
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TApplication.h"
#include "TFile.h"
#include "TLorentzVector.h"

using namespace lcio ;
using namespace marlin ;
using namespace WWTools;



/** Processor calculates MC decay angles from W boson.
 *  PROTOCOL
 * 1) Goes through MC Particles
 * 
 * 2) [ON SHELL] W has PDG abs(24) ---> Get four-momentum && Charge
 *  -  Production angle of W- (CosineTheta)
 *  -  LEPTONIC and HADRONIC decay discrimination (Directly from W daughters).
 *  -  LEPTONIC DECAY FROM ---> W- or W+ (cross check by counting) ---> For decay angles
 *  -  HADRONIC DECAY FROM ---> W- or W+ (cross check by counting)
 *  -  DECAY ANGLES ---> Use starVector ---> *theta and *phi
 * 
 * 3) [OFF SHELL] W reconstructed from it's daughters ---> Get four-momentum && Charge
 *  -  Production angle of W- (CosineTheta)
 *  -  LEPTONIC and HADRONIC decay discrimination (pdg > or < than 10) 
 *  -  LEPTONIC FROM ---> W- or W+ (cross check by counting) ---> Use for decay angles
 *  -  HADRONIC FROM ---> W- or W+ (cross check by counting)
 *  -  DECAY ANGLES ---> Use starVector -> *theta and *phi
 * 
 *  @author A.Silva, DESY
 */
class WWAngleCalculationProcessor: public Processor {
  
 public:
  
  virtual Processor*  newProcessor() { return new WWAngleCalculationProcessor; }
  
  
  WWAngleCalculationProcessor() ;

  WWAngleCalculationProcessor(const WWAngleCalculationProcessor&) = delete;
  WWAngleCalculationProcessor& operator=(const WWAngleCalculationProcessor&) = delete;
  
  /** Called at the begin of the job before anything is read.
   * Use to initialize the processor, e.g. book histograms.
   */
  virtual void init() ;
  
  /** Called for every run.
   */
  virtual void processRunHeader( LCRunHeader* run ) ;

  double magnitude(std::vector<double>& v1);

  
  /** Called for every event - the working horse.
   */
  virtual void processEvent( LCEvent * evt ) ; 
 
  
  virtual void check( LCEvent * evt ) ; 
  
  
  /** Called after data processing for clean up.
   */
  virtual void end() ;
 
  
 protected:

 TH1F* _star_phi{};
 TH1F* _star_phi_minus{};
 TH1F* _star_phi_total{};
 TH1F* _star_theta{};
 TH1F* _cos_prodw{};

 TH1F* _hBohdan{};

 //collection names
  std::string _MCParColName{};
  std::string _ElectronColName{};
  std::string _TauColName{};
  std::string _MuonColName{};
  std::string _PhotonColName{};
  std::string _jetcolName{};
  std::string _PandoraPFOsCol{};
  std::string _PFOsMinusIsolatedObjetcs{};
  std::string _LikelihoodPIDMethod = "dEdxPIDv2";
  std::string _TTreeFileName{};

  TFile* _TTreeFile{};

  std::vector<long long int> _nW{};
  std::vector<long long int> _nNonW{};
  int _W = 0;
  int _Whadronic = 0;
  int _Whadronic1 = 0;
  int _Wleptonic = 0;
  int _Wleptonic1 = 0;
  int _nonW = 0;


  TApplication _application = TApplication("app", 0, nullptr);
  int _nRun{};
  int _nEvt{};
  
} ;

#endif


