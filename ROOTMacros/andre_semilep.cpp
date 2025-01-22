#include <vector>
#include <fstream>
#include <iostream>


using namespace std;


void andre_semilep(){
    // Open the slcio file
    fstream File;

    IO::LCReader* lcReader = IOIMPL::LCFactory::getInstance()->createLCReader();
    lcReader->setReadCollectionNames({"PandoraPFOs", "MCParticlesSkimmed", "IsolatedMuons", "IsolatedElectrons", "IsolatedPhotons","IsolatedTaus"});

    //std::vector<std::string> files = {"./rv02-02-01.sv02-02-01.mILD_l5_o2_v02.E250-SetA.I500084.P4f_ww_sl.eR.pL.n000.d_mini_dstm_15669_0.slcio"}

    lcReader->open("./rv02-02-01.sv02-02-01.mILD_l5_o2_v02.E250-SetA.I500084.P4f_ww_sl.eR.pL.n000.d_mini_dstm_15669_0.slcio");

    int i, j;
    int sumCount = 0;

    TH1F* histo_number_muons = new TH1F("histo_number_muons", " ;  Number of Isolated Muons ; Events", 5, -0.5, 4.5);
    TH1F* histo_number_photons = new TH1F("  ", " ;  Number of Isolated Photons ; Events", 5, -0.5, 4.5);
    TH1F* histo_number_electrons = new TH1F("   ", " ;  Number of Isolated Electrons ; Events", 5, -0.5, 4.5);
    TH1F* histo_number_taus = new TH1F("   ", " ;  Number of Isolated Taus ; Events", 5, -0.5, 4.5);

    TH1F* histo_number_leptons = new TH1F("    ", " ;  Number of Isolated leptons ; Events", 5, -0.5, 4.5);



    EVENT::LCEvent* event;
    int nEvents = 0;
    int maxEvent = 9999999;
    
    // isolated muons
    int iso_muon_zero = 0;
    int iso_muon_one = 0;
    //isolated photons
    int iso_photon_zero = 0;
    int iso_photon_one = 0;
    //isolated electrons
    int iso_electrons_zero = 0;
    int iso_electrons_one = 0;

    std::vector<int> iso_muon_number = {iso_muon_zero, iso_muon_one};
    std::vector<int> iso_photon_number = {iso_photon_zero, iso_photon_one};
    std::vector<int> iso_electrons_number = {iso_electrons_zero, iso_electrons_one};

    while( ( event = lcReader->readNextEvent() ) != nullptr && (nEvents++) < maxEvent ){
	    //printing status every 500 events
      if (nEvents % 500 == 0){
        std::cout<<"Analysing event #"<<nEvents<<"/"<<min(lcReader->getNumberOfEvents(), maxEvent)<<std::endl;
      }

      //reconstructed particles
      EVENT::LCCollection* MuonsCollection = event->getCollection("IsolatedMuons");
      EVENT::LCCollection* PhotonsCollection = event->getCollection("IsolatedPhotons");
      EVENT::LCCollection* ElectronsCollection = event->getCollection("IsolatedElectrons");
      EVENT::LCCollection* TausCollection = event->getCollection("IsolatedTaus");

      //mcparticles
      EVENT::LCCollection* MCParticles = event->getCollection("MCParticlesSkimmed");
      

      int n_muons = MuonsCollection->getNumberOfElements(); //gets size of collection
      int n_photons = PhotonsCollection->getNumberOfElements();
      int n_electrons = ElectronsCollection->getNumberOfElements();
      int n_taus = TausCollection->getNumberOfElements();

      int n_leptons = n_electrons + n_muons + n_taus;
      

      //leptons individually
      histo_number_muons->Fill(n_muons);
      histo_number_electrons->Fill(n_electrons);
      histo_number_taus->Fill(n_taus);
      //Photons
      histo_number_photons->Fill(n_photons);
      //Leptons
      histo_number_leptons->Fill(n_leptons);
    }  
    

    lcReader->close(); 

    TCanvas* canvas1 = new TCanvas("canvas1","Muon Number");
    histo_number_muons->Draw("LEGO");
    histo_number_muons->SetLineColor(1);
    int a = histo_number_muons->GetBinContent(1, maxEvent);
    int b = histo_number_muons->GetBinContent(2, maxEvent);
    int c = histo_number_muons->GetBinContent(3, maxEvent);
    int d = histo_number_muons->GetBinContent(4, maxEvent);
    std::cout << "Events w/ 0 Isolated Muons: " << a << std::endl;
    std::cout << "Events w/ 1 Isolated Muons: " << b << std::endl;
    std::cout << "Events w/ 2 Isolated Muons: " << c << std::endl;
    std::cout << "Events w/ 3 Isolated Muons: " << d << std::endl;

    TCanvas* canvas2 = new TCanvas("canvas2", "Tau Number");
    histo_number_taus->Draw("LEGO");
    histo_number_taus->SetLineColor(1);
    int aa = histo_number_taus->GetBinContent(1, maxEvent);
    int bb = histo_number_taus->GetBinContent(2, maxEvent);
    int cc = histo_number_taus->GetBinContent(3, maxEvent);
    int dd = histo_number_taus->GetBinContent(4, maxEvent);
    int ee = histo_number_taus->GetBinContent(5, maxEvent);
    std::cout << "Events w/ 0 Isolated Taus: " << aa << std::endl;
    std::cout << "Events w/ 1 Isolated Taus: " << bb << std::endl;
    std::cout << "Events w/ 2 Isolated Taus: " << cc << std::endl;
    std::cout << "Events w/ 3 Isolated Taus: " << dd << std::endl;
    std::cout << "Events w/ 4 Isolated Taus: " << ee << std::endl;


    TCanvas* canvas3 = new TCanvas("canvas3","Electron Number");
    histo_number_electrons->Draw("LEGO");
    histo_number_electrons->SetLineColor(1);
    int aaa = histo_number_electrons->GetBinContent(1, maxEvent);
    int bbb = histo_number_electrons->GetBinContent(2, maxEvent);
    int ccc = histo_number_electrons->GetBinContent(3, maxEvent);
    int ddd = histo_number_electrons->GetBinContent(4, maxEvent);
    std::cout << "Events w/ 0 Isolated Electrons: " << aaa << std::endl;
    std::cout << "Events w/ 1 Isolated Electrons: " << bbb << std::endl;
    std::cout << "Events w/ 2 Isolated Electrons: " << ccc << std::endl;
    std::cout << "Events w/ 3 Isolated Electrons: " << ddd << std::endl;

    //Classifying Events by Isolated Lepton Numbers
    TCanvas* canvas4 = new TCanvas("canvas4","Isolated Lepton Number");
    histo_number_leptons->Draw("LEGO");
    histo_number_leptons->SetLineColor(1);
    int zero_lepton = histo_number_leptons->GetBinContent(1, maxEvent);
    int one_lepton = histo_number_leptons->GetBinContent(2, maxEvent);
    int two_lepton = histo_number_leptons->GetBinContent(3, maxEvent);
    int three_lepton = histo_number_leptons->GetBinContent(4, maxEvent);
    int four_lepton = histo_number_leptons->GetBinContent(5, maxEvent);
    std::cout << "Events w/ 0 Isolated Leptons: " << zero_lepton << std::endl;
    std::cout << "Events w/ 1 Isolated Leptons: " << one_lepton << std::endl;
    std::cout << "Events w/ 2 Isolated Leptons: " << two_lepton<< std::endl;
    std::cout << "Events w/ 3 Isolated Leptons: " << three_lepton << std::endl;
    std::cout << "Events w/ 4 Isolated Leptons: " << four_lepton << std::endl;








    //histo_invariant_mass_WZ->SetLineColor(6);
  
    /*TLegend *leg = new TLegend(0.5, 0.6, 0.8, 0.8);
    leg->SetBorderSize(0);
    leg->AddEntry(histo_invariant_mass, "Invariant W Di-jet", "l");
    leg->AddEntry(histo_reference_W, "Reference Invariant W Di-jet", "l");
    leg->AddEntry(histo_invariant_mass_Z, "Invariant Z Di-jet", "l");
    leg->AddEntry(histo_reference_Z, "Reference Invariant Z Di-jet", "l");
    leg->Draw();*/

    //TCanvas* canvas6 = new TCanvas("canvas6","Mass distribuition");
    //histo_masses_WZ->Draw();





    // you can try this. It doesn't look pretty, but my goal is to show you what you can modify in style, there are more of course..
   /* TCanvas* canvas2 = new TCanvas("canvas2","W Events Invariant Mass");
    histo_total->SetLineColor(9);
    histo_total->SetLineWidth(2);
    histo_total->SetMarkerStyle(21);
    histo_total->SetMarkerColor(9);
    histo_total->SetMarkerSize(1.);
    //histo_total->Draw("PL");
    canvas2->SetLogy();
    canvas2->SetGridx();
    //canvas2->SetGridy();
    canvas2->BuildLegend();
    canvas2->Update(); */

   /* TCanvas* canvas4 = new TCanvas("canvas4","Z Events Invariant Mass");
    histo_total_Z->SetLineColor(4);
    histo_total_Z->SetLineWidth(2);
    histo_total_Z->SetMarkerStyle(21);
    histo_total_Z->SetMarkerColor(6);
    histo_total_Z->SetMarkerSize(1.);
    //histo_total_Z->Draw("PL");
    canvas4->SetLogy();
    canvas4->SetGridx();
    //canvas2->SetGridy();
    canvas4->BuildLegend();
    canvas4->Update();*/


}

