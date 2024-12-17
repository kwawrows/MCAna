///////////////////////////////////////////////////////////////////////
// Class:       EventAnaTree
// Plugin Type: analyzer (art v3_05_01)
// File:        EventAnaTree_module.cc
// Generated Dec 17 2024 by Klaudia Wawrowska
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "art_root_io/RootIOPolicy.h"
#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"

#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Persistency/Common/FindManyP.h" // for associations
#include "canvas/Persistency/Common/PtrVector.h"
#include "canvas/Persistency/Common/FindOneP.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/WireGeo.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_vectors.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCFlux.h"
#include "lardataobj/MCBase/MCShower.h"
#include "lardataobj/MCBase/MCTrack.h"
#include "lardataobj/MCBase/MCStep.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "lardataobj/Simulation/AuxDetSimChannel.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RawData/raw.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/PhotonBackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "larreco/SpacePointSolver/Solver.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "lardataobj/RawData/RDTimeStamp.h"
#include "larreco/RecoAlg/TrackMomentumCalculator.h"
#include "lardataobj/AnalysisBase/CosmicTag.h"
#include "lardataobj/AnalysisBase/FlashMatch.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "lardataobj/AnalysisBase/BackTrackerMatchingData.h"
#include "larevt/SpaceChargeServices/SpaceChargeService.h"


//Additional framework
#include "TTree.h"
#include <string>
#include "TFile.h"
#include "TDirectory.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TProfile.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TMath.h"
#include "TGraphErrors.h"
#include "TMinuit.h"
#include "TString.h"
#include "TTimeStamp.h"
#include "TVectorD.h"
#include "TCanvas.h"
#include "TFrame.h"
#include "TLine.h"
#include "TAxis.h"
#include "TTimeStamp.h"
#include "TVector3.h"
#include <vector>
#include <fstream>
#include "TPaveStats.h"
#include <iostream>
#include <string>
#include "math.h"
#include "stdio.h"
#include <iterator>
#include "TGeoMatrix.h"

namespace dune {
  class EventAnaTree;
}


class dune::EventAnaTree : public art::EDAnalyzer {
public:
  explicit EventAnaTree(fhicl::ParameterSet const& p);

  EventAnaTree(EventAnaTree const&) = delete;
  EventAnaTree(EventAnaTree&&) = delete;
  EventAnaTree& operator=(EventAnaTree const&) = delete;
  EventAnaTree& operator=(EventAnaTree&&) = delete;

  // function initializations
  void beginJob() override;
  void endJob() override;
  void analyze(const art::Event& e) override;

private:

  //--- Function definitions ---
  void ResetVariables();
  long unsigned int WhichGeneratorType(int TrID); 
  void FillMaps(std::map<int, simb::MCParticle> &MyMap, art::FindManyP<simb::MCParticle> Assn, art::ValidHandle<std::vector<simb::MCTruth>> Hand); 
  bool InMap(int TrID, std::map< int, simb::MCParticle> ParMap);
  bool InMap(int TrID, std::map<int, float> TrackIDMap); 

  //--- Root Tree ---
  TTree *fTree; 


  // --- Fcl Configurables ---
  bool fSaveGenieInfo;
  std::string fTruthLabel; //which module produced simulation
  std::string fHitLabel;
  std::string fGenieLabel;
  std::string fGeantLabel; //G4 label 
  std::vector<std::string> fLabels; //generator label vector 
  
  // ===========================================================
  //                          EVENT VARIABLES 
  // ============================================================

  unsigned int fEvent; 
  unsigned int fRun;
  unsigned int fSubRun;


  //Map for associating MC particles to their truth label 
  std::vector<std::map<int, simb::MCParticle>> GeneratorParticles = {}; 
  std::vector<int> fNGenParts; // number of particles per generator type 

  //Genie / generator info
  unsigned int fgenie_no_primaries; 
  std::vector<int>    fgenie_primaries_pdg; 
  std::vector<float>  fgenie_Eng; 
  std::vector<float>  fgenie_Px;
  std::vector<float>  fgenie_Py;
  std::vector<float>  fgenie_Pz;
  std::vector<float>  fgenie_P;
  std::vector<int>    fgenie_status_code; 

  std::vector<int>    fnuPDG_truth; 
  std::vector<int>    fccnc_truth; // 0 = CC, 1 = NC
  std::vector<int>    fmode_truth; // 0 = QE/El, 1=RES, 2=DIS, 3= Coherent production
  std::vector<float>  fhitnuc_truth;//hit nucleon 
  std::vector<float>  fnu_vx_truth; // neutrino vertex x 
  std::vector<float>  fnu_vy_truth; // neutrino vertex y 
  std::vector<float>  fnu_vz_truth; // neutrino vertex z 
  std::vector<float>  fnu_dcosx_truth; 
  std::vector<float>  fnu_dcosy_truth; 
  std::vector<float>  fnu_dcosz_truth; 
  std::vector<float>  flep_mom_truth;
  std::vector<float>  flep_dcosx_truth;
  std::vector<float>  flep_dcosy_truth;
  std::vector<float>  flep_dcosz_truth;


  //Geant/truth info
  unsigned int fnPrimaries; // no. of primary particles 
  unsigned int fnGeantParticles; //total number of geant particles 
  std::vector<int>     fTrackId;
  std::vector<float>   fMother; 
  std::vector<float>   fEng;
  std::vector<float>   fEkin;
  std::vector<float>   fMass;
  std::vector<int>     fPdg;
  std::vector<float>   fP;
  std::vector<float>   fPx;
  std::vector<float>   fPy;
  std::vector<float>   fPz;
  std::vector<int>     fND; //# daughters 
  std::vector<double>  fstartX;
  std::vector<double>  fstartY;
  std::vector<double>  fstartZ;
  std::vector<double>  fendX;
  std::vector<double>  fendY;
  std::vector<double>  fendZ;

  //TPC info
  unsigned int fnHits; //number of hits 
  unsigned int fnColHits; // # of hits on collection plane 
  std::vector<int>     fhit_tpc;
  std::vector<int>     fhit_channel;
  std::vector<float>   fhit_time;
  std::vector<float>   fhit_SADC;
  std::vector<float>   fhit_wire;
  std::vector<float>   fhit_charge;
  std::vector<int>     fhit_plane;
  std::vector<float>   fhit_width; //width RMS in time
  std::vector<int>     fhit_TOT; // hit width in ticks
  std::vector<float>   fhit_trueX;
  std::vector<float>   fhit_trueY;
  std::vector<float>   fhit_trueZ;
  std::vector<int>     fhit_clusterId; 

  //collection info
  std::vector<int>     fcolhit_tpc;
  std::vector<int>     fcolhit_channel;
  std::vector<float>   fcolhit_time;
  std::vector<float>   fcolhit_SADC;
  std::vector<float>   fcolhit_wire;
  std::vector<float>   fcolhit_charge;
  std::vector<float>   fcolhit_width; 
  std::vector<float>   fcolhit_trueX;
  std::vector<float>   fcolhit_trueY;
  std::vector<float>   fcolhit_trueZ;
  std::vector<int>     fcolhit_label;

  //Declare services
  art::ServiceHandle<geo::Geometry> geo;
  art::ServiceHandle<cheat::PhotonBackTrackerService> pbt_serv;
  art::ServiceHandle<cheat::BackTrackerService> bt_serv; // MC cheat for mapping back to truth 
  art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
};

// Constructor with config. parameters 
dune::EventAnaTree::EventAnaTree(fhicl::ParameterSet const& p) :
  EDAnalyzer{p}//,
{
  //Module label configs
  fTruthLabel          = p.get<std::string>("TruthLabel");
  fHitLabel            = p.get<std::string>("HitLabel");
  fGenieLabel          = p.get<std::string>("GenieLabel");
  fGeantLabel          = p.get<std::string>("GeantLabel");
  fSaveGenieInfo       = p.get<bool>("SaveGenieInfo", 0);
  fLabels              = p.get<std::vector<std::string>>("GeneratorLabelVector");
}


void dune::EventAnaTree::ResetVariables()
{

  fEvent = fRun = fSubRun = -1;
  fgenie_no_primaries = 0;
  fnPrimaries = fnGeantParticles = 0;
  fnHits = fnColHits = 0;

  fNGenParts = {}; 
  GeneratorParticles = {}; 
    
  fgenie_primaries_pdg.clear();
  fgenie_Eng.clear();
  fgenie_Px.clear();
  fgenie_Py.clear();
  fgenie_Pz.clear();
  fgenie_P.clear();
  fgenie_status_code.clear();

  fnuPDG_truth.clear();
  fccnc_truth.clear();
  fmode_truth.clear();
  fhitnuc_truth.clear();
  fnu_vx_truth.clear();
  fnu_vy_truth.clear();
  fnu_vz_truth.clear();
  fnu_dcosx_truth.clear();
  fnu_dcosy_truth.clear();
  fnu_dcosz_truth.clear();
  flep_mom_truth.clear();
  flep_dcosx_truth.clear();
  flep_dcosy_truth.clear();
  flep_dcosz_truth.clear();


  fTrackId.clear();
  fMother.clear();
  fEng.clear();
  fEkin.clear();
  fMass.clear();
  fPdg.clear();
  fP.clear();
  fPx.clear();
  fPy.clear();
  fPz.clear();
  fND.clear();
  fstartX.clear();
  fstartY.clear();
  fstartZ.clear();
  fendX.clear();
  fendY.clear();
  fendZ.clear();


  fhit_tpc.clear();
  fhit_channel.clear();
  fhit_time.clear();
  fhit_SADC.clear();
  fhit_wire.clear();
  fhit_charge.clear();
  fhit_plane.clear();
  fhit_width.clear();
  fhit_TOT.clear();
  fhit_trueX.clear();
  fhit_trueY.clear();
  fhit_trueZ.clear();
  fcolhit_tpc.clear();
  fcolhit_channel.clear();
  fcolhit_time.clear();
  fcolhit_SADC.clear();
  fcolhit_wire.clear();
  fcolhit_charge.clear();
  fcolhit_width.clear();
  fcolhit_trueX.clear();
  fcolhit_trueY.clear();
  fcolhit_trueZ.clear();
  fcolhit_label.clear();
 
  
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///
///                                       ANALYZER
///
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void dune::EventAnaTree::analyze(art::Event const& e) 
{
  ResetVariables();
 
  auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(e);
  
  fEvent   = e.id().event();
  fRun     = e.run();
  fSubRun  = e.subRun();

  bool isMC = !e.isRealData();
  
  //Prepare for particle tagging 
  std::vector<std::set<int>> trackids = {};
  std::map<int, simb::MCParticle> ThisGeneratorParts;

  if (isMC){

    //Create maps for ID tracking 
    const sim::ParticleList &PartList = pi_serv->ParticleList();
    fnGeantParticles = PartList.size(); // total number of particles in event 

    //Loop over all signal and background labels and collect track IDs 
    for (size_t l = 0; l < fLabels.size(); l++){

      GeneratorParticles.push_back(ThisGeneratorParts); // For each label insert empty list

      art::Handle<std::vector<simb::MCTruth>> ThisHandle;
      e.getByLabel(fLabels[l], ThisHandle);

      if (ThisHandle)
        {
          auto LabelHandle = e.getValidHandle<std::vector<simb::MCTruth>>(fLabels[l]); //valid generator handle
          art::FindManyP<simb::MCParticle> Assn(LabelHandle, e, fGeantLabel); // Assign labels to MC particles 
          FillMaps(GeneratorParticles[l], Assn, LabelHandle);                 // Fill empty list with previously assigned particles

          fNGenParts.push_back(GeneratorParticles[l].size()); // save number of particles for this generator 

          //Empty set for track IDs for current generator 
          std::set<int> ThisGeneratorIDs = {}; 

          //Loop over particles for generator [l] & populate the set with track IDs  
          if (!GeneratorParticles[l].empty()){

            for (const auto& part : GeneratorParticles[l]){
              trackids.push_back(ThisGeneratorIDs);
              trackids[l].insert(part.first);
            }
          } else{
            //Handle case of particle list is empty 
            fNGenParts.push_back(0); 
            trackids.push_back(ThisGeneratorIDs); 
          }
        }
    }


    if ( fSaveGenieInfo ) {

      art::Ptr<simb::MCTruth> mctruth;  

      // ----------------------- Genie:: neutrino truth info ---------------------------------------------
      art::ValidHandle<std::vector <simb::MCTruth> > mctruthListHandle = e.getValidHandle< std::vector<simb::MCTruth> >(fGenieLabel);
      std::vector< art::Ptr<simb::MCTruth> > mclist;

      if (mctruthListHandle.isValid()){
        art::fill_ptr_vector(mclist, mctruthListHandle);
      }

      if (!mclist.empty()){
        mctruth = mclist[0];
      }
      if (mctruth->NeutrinoSet() ){
        fgenie_no_primaries = mctruth->NParticles();
    
        for (size_t iPart = 0; iPart < fgenie_no_primaries; ++iPart){

          const simb::MCParticle& part(mctruth->GetParticle(iPart));
          fgenie_primaries_pdg.push_back( part.PdgCode() );
          fgenie_Eng.push_back( part.E() );
          fgenie_Px.push_back( part.Px() );
          fgenie_Py.push_back( part.Py() );
          fgenie_Pz.push_back( part.Pz() );
          fgenie_P.push_back( part.P() );
          fgenie_status_code.push_back( part.StatusCode() );
            
        } // loop over genie particles 
      } // if (mctruth->NeutrinoSet()) 

      //Save neutrino interaction information 31.01.2022
      if (mclist.size() > 0) {
        int neutrino_i = 0;
        for (unsigned int iList = 0; iList < mclist.size(); ++iList){
          if ( mclist[iList]->NeutrinoSet() ){
            fnuPDG_truth.push_back( mclist[iList]->GetNeutrino().Nu().PdgCode() );
            fccnc_truth.push_back( mclist[iList]->GetNeutrino().CCNC() );
            fmode_truth.push_back( mclist[iList]->GetNeutrino().Mode() );
            fhitnuc_truth.push_back( mclist[iList]->GetNeutrino().HitNuc() );
            fnu_vx_truth.push_back( mclist[iList]->GetNeutrino().Nu().Vx() );
            fnu_vy_truth.push_back( mclist[iList]->GetNeutrino().Nu().Vy() );
            fnu_vz_truth.push_back( mclist[iList]->GetNeutrino().Nu().Vz() );

            if ( mclist[iList]->GetNeutrino().Nu().P() ){
              fnu_dcosx_truth.push_back( mclist[iList]->GetNeutrino().Nu().Px()/(mclist[iList]->GetNeutrino().Nu().P()) );
              fnu_dcosy_truth.push_back( mclist[iList]->GetNeutrino().Nu().Py()/(mclist[iList]->GetNeutrino().Nu().P()) );
              fnu_dcosz_truth.push_back( mclist[iList]->GetNeutrino().Nu().Pz()/(mclist[iList]->GetNeutrino().Nu().P()) );
            } //angular info for neutrino

            flep_mom_truth.push_back( mclist[iList]->GetNeutrino().Lepton().P() );
            if ( mclist[iList]->GetNeutrino().Lepton().P()){
              flep_dcosx_truth.push_back(  mclist[iList]->GetNeutrino().Lepton().Px()/(mclist[iList]->GetNeutrino().Lepton().P()) );
              flep_dcosy_truth.push_back(  mclist[iList]->GetNeutrino().Lepton().Py()/(mclist[iList]->GetNeutrino().Lepton().P()) );
              flep_dcosz_truth.push_back(  mclist[iList]->GetNeutrino().Lepton().Pz()/(mclist[iList]->GetNeutrino().Lepton().P()) );
            } //angular info for lepton 
             
            neutrino_i++;
          }// mclist is NeutrinoSet()
        }  // loop over MC records 
      } //at least one MC record 

    } // if SaveGenieInfo

  
    // ----------------------- GEANT :: truth list of MC particles -----------------------------------

    art::ValidHandle<std::vector <simb::MCParticle> > mcParticles = e.getValidHandle<std::vector <simb::MCParticle> >(fTruthLabel);
    if (mcParticles.isValid()){
      fnPrimaries = mcParticles->size();
      
      //Get MC particle list from G4 
      for (unsigned int t = 0; t < mcParticles->size(); ++t){ 
        const simb::MCParticle trueParticle = mcParticles->at(t);

        fTrackId.push_back( trueParticle.TrackId());
        fMother.push_back( trueParticle.Mother());
        fEng.push_back( trueParticle.E() );
        fPdg.push_back( trueParticle.PdgCode());
        fEkin.push_back( trueParticle.E() - trueParticle.Mass() );
        fMass.push_back( trueParticle.Mass() );
        fP.push_back( trueParticle.P()) ;
        fPx.push_back( trueParticle.Px());
        fPy.push_back( trueParticle.Py());
        fPz.push_back( trueParticle.Pz());
        fND.push_back( trueParticle.NumberDaughters());
        fstartX.push_back( trueParticle.Vx());
        fstartY.push_back( trueParticle.Vy());
        fstartZ.push_back( trueParticle.Vz());
        fendX.push_back( trueParticle.EndPosition()[0]);
        fendY.push_back( trueParticle.EndPosition()[1]);
        fendZ.push_back( trueParticle.EndPosition()[2]);

      }
    } //if MC particles are valid 

    //----------------------- LarSoft ::  hit list -------------------------- 
    //art::ValidHandle<std::vector <recob::Hit> > recoHits = e.getValidHandle<std::vector <recob::Hit> >(fHitLabel);
    std::vector< art::Ptr<recob::Hit> >recoHits; 
    auto RecoHitHandle = e.getValidHandle<std::vector<recob::Hit>>(fHitLabel); // art handle for the reco hits 
    if(RecoHitHandle){ art::fill_ptr_vector(recoHits, RecoHitHandle); }

    //Run over the reco hits 
    if ( RecoHitHandle ){
      fnHits = recoHits.size();
      fnColHits = 0;
      
      for (unsigned int hit = 0; hit < recoHits.size(); ++hit){
        //const recob::Hit thisHit = recoHits->at(hit);
        const art::Ptr<recob::Hit> thisHit = recoHits.at(hit);

        fhit_tpc.push_back(thisHit->WireID().TPC);
        fhit_channel.push_back(thisHit->Channel());
        fhit_time.push_back(thisHit->PeakTime());
        fhit_SADC.push_back(thisHit->SummedADC());
        fhit_wire.push_back(thisHit->WireID().Wire);
        fhit_charge.push_back(thisHit->Integral());
        fhit_plane.push_back(thisHit->View());
        fhit_width.push_back( std::abs( thisHit->PeakTimePlusRMS() - thisHit->PeakTimeMinusRMS()) ); 
        fhit_TOT.push_back(thisHit->EndTick() - thisHit->StartTick()); // time over threshold in ticks 
        // 1 tick = 500 ns (2 MHz sampling frequency)

        // extracting physical position of the hit in detector geometry  {
        std::vector<const sim::IDE*> ides;
        try{ ides = bt_serv->HitToSimIDEs_Ps(clockData,thisHit); }
        catch(...){}
        if (ides.size() > 0){
          //getting physical coordinates of hit from ides 
          std::vector<double> xyz = bt_serv->SimIDEsToXYZ(ides);
          fhit_trueX.push_back(xyz[0]);
          fhit_trueY.push_back(xyz[1]);
          fhit_trueZ.push_back(xyz[2]);
        }


        if (thisHit->View() == 2){  //collection plane hits only
          fnColHits +=1;   
          fcolhit_tpc.push_back(thisHit->WireID().TPC);
          fcolhit_channel.push_back(thisHit->Channel());
          fcolhit_time.push_back(thisHit->PeakTime());
          fcolhit_SADC.push_back(thisHit->SummedADC());
          fcolhit_wire.push_back(thisHit->WireID().Wire);
          fcolhit_charge.push_back(thisHit->Integral());
          fcolhit_width.push_back( std::abs( thisHit->PeakTimePlusRMS() - thisHit->PeakTimeMinusRMS()) );
          std::vector<const sim::IDE*> ColIDES;
          try{
            ColIDES = bt_serv->HitToSimIDEs_Ps(clockData,thisHit);
          }
          catch(...){}
          if (ColIDES.size() > 0){
            std::vector<double> colxyz = bt_serv->SimIDEsToXYZ(ColIDES);
            fcolhit_trueX.push_back(colxyz[0]);
            fcolhit_trueY.push_back(colxyz[1]);
            fcolhit_trueZ.push_back(colxyz[2]);
          } // hit XYZ : ides > 0 
          
          //tag the collection hits 
          std::vector<sim::TrackIDE> hitTrkIDE = bt_serv-> HitToTrackIDEs(clockData, thisHit); 
          double MainTrkID = -1; 
          double TopEFrac  = 0; 
          
          //Loop over IDEs 
          for (size_t ideL = 0; ideL < hitTrkIDE.size(); ++ideL)
            {
              if (hitTrkIDE[ideL].energyFrac > TopEFrac)
                {
                  TopEFrac = hitTrkIDE[ideL].energyFrac;
                  MainTrkID = abs(hitTrkIDE[ideL].trackID); // HERE
                }
            }
          
          // Default value for the generator label index
          int labelIndex = -1;
          
          // Check if MainTrkID exists in the map and get its generator label index
          if (MainTrkID != -1)
            {
              labelIndex = WhichGeneratorType( MainTrkID);
            }
          fcolhit_label.push_back(labelIndex); // which generator produced this hit: 0 = marley, 1+.. = radiologicals, -1 = noise (or something went wrong)
        } // col hits 
      } // run over all hits
    } //if reco hits valid     
  }//if is MC
  

  fTree->Fill();
  
} // analyze function end 

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///
///                                       BEGIN JOB
///
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void dune::EventAnaTree::beginJob()
{

  art::ServiceHandle< art::TFileService> tfs;
  fTree = tfs->make<TTree>("tree", "Analyser Output Tree"); 

  fTree->Branch("event",&fEvent,"event/i");
  fTree->Branch("run",&fRun,"run/i");
  fTree->Branch("subrun",&fSubRun,"subrun/i");


  if (fSaveGenieInfo){

    fTree->Branch("genie_no_primaries",&fgenie_no_primaries,"genie_no_primaries/i");
    fTree->Branch("genie_primaries_pdg",&fgenie_primaries_pdg);
    fTree->Branch("genie_Eng",&fgenie_Eng);
    fTree->Branch("genie_Px",&fgenie_Px);
    fTree->Branch("genie_Py",&fgenie_Py);
    fTree->Branch("genie_Pz",&fgenie_Pz);
    fTree->Branch("genie_P",&fgenie_P);
    fTree->Branch("genie_status_code",&fgenie_status_code); 

    fTree->Branch("nuPDG_truth",&fnuPDG_truth);
    fTree->Branch("ccnc_truth",&fccnc_truth);
    fTree->Branch("mode_truth",&fmode_truth);
    fTree->Branch("hitnuc_truth",&fhitnuc_truth);
    fTree->Branch("nu_vx_truth",&fnu_vx_truth);
    fTree->Branch("nu_vy_truth",&fnu_vy_truth);
    fTree->Branch("nu_vz_truth",&fnu_vz_truth);
    fTree->Branch("nu_dcosx_truth",&fnu_dcosx_truth);
    fTree->Branch("nu_dcosy_truth",&fnu_dcosy_truth);
    fTree->Branch("nu_dcosz_truth",&fnu_dcosz_truth);
    fTree->Branch("lep_mom_truth",&flep_mom_truth);
    fTree->Branch("lep_dcosx_truth",&flep_dcosx_truth);
    fTree->Branch("lep_dcosy_truth",&flep_dcosy_truth);
    fTree->Branch("lep_dcosz_truth",&flep_dcosz_truth);


  }

  fTree->Branch("NGenParts", &fNGenParts); // number of g4 particles per generator 
  fTree->Branch("nGeantParticles",&fnGeantParticles,"nGeantParticles/i"); 
  fTree->Branch("nPrimaries",&fnPrimaries,"nPrimaries/i");
  fTree->Branch("TrackId",&fTrackId);
  fTree->Branch("Mother",&fMother);
  fTree->Branch("Pdg",&fPdg);
  fTree->Branch("Eng",&fEng);
  fTree->Branch("Ekin",&fEkin);
  fTree->Branch("Mass",&fMass);
  fTree->Branch("P",&fP);
  fTree->Branch("Px",&fPx);
  fTree->Branch("Py",&fPy);
  fTree->Branch("Pz",&fPz);
  fTree->Branch("ND",&fND);
  fTree->Branch("startX",&fstartX);
  fTree->Branch("startY",&fstartY);
  fTree->Branch("startZ",&fstartZ);
  fTree->Branch("endX",&fendX);
  fTree->Branch("endY",&fendY);
  fTree->Branch("endZ",&fendZ);

  fTree->Branch("nHits",&fnHits,"nHits/i");
  fTree->Branch("nColHits",&fnColHits,"nColHits/i");
  fTree->Branch("hit_tpc",&fhit_tpc);
  fTree->Branch("hit_channel",&fhit_channel);
  fTree->Branch("hit_time",&fhit_time);
  fTree->Branch("hit_SADC",&fhit_SADC);
  fTree->Branch("hit_wire",&fhit_wire);
  fTree->Branch("hit_charge",&fhit_charge);
  fTree->Branch("hit_plane",&fhit_plane);
  fTree->Branch("hit_width",&fhit_width);
  fTree->Branch("hit_TOT",&fhit_TOT);
  fTree->Branch("hit_trueX",&fhit_trueX);
  fTree->Branch("hit_trueY",&fhit_trueY);
  fTree->Branch("hit_trueZ",&fhit_trueZ);
  fTree->Branch("colhit_tpc",&fcolhit_tpc);
  fTree->Branch("colhit_channel",&fcolhit_channel);
  fTree->Branch("colhit_time",&fcolhit_time);
  fTree->Branch("colhit_SADC",&fcolhit_SADC);
  fTree->Branch("colhit_wire",&fcolhit_wire);
  fTree->Branch("colhit_charge",&fcolhit_charge);
  fTree->Branch("colhit_width",&fcolhit_width);
  fTree->Branch("colhit_trueX",&fcolhit_trueX);
  fTree->Branch("colhit_trueY",&fcolhit_trueY);
  fTree->Branch("colhit_trueZ",&fcolhit_trueZ);
  fTree->Branch("colhit_label",&fcolhit_label);
  
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///
///                                       END JOB
///
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void dune::EventAnaTree::endJob()
{
    
}


//Function which returns the type of generator particle which produced a given MC truth 
long unsigned int dune::EventAnaTree::WhichGeneratorType(int TrID)
{
  for (long unsigned int i = 0; i < fLabels.size(); i++)
    {
      if (InMap(TrID, GeneratorParticles[i]))
        {
          return i + 1;
        }
    }
  return 0;
}

 
//Function to fill a map with MCparticles from a given MCTruth 
void dune::EventAnaTree::FillMaps(std::map<int, simb::MCParticle> &MyMap, art::FindManyP<simb::MCParticle> Assn, art::ValidHandle<std::vector<simb::MCTruth>> Hand)
{
  for (size_t L1 = 0; L1 < Hand->size(); ++L1)
    {
      for (size_t L2 = 0; L2 < Assn.at(L1).size(); ++L2)
        {
          const simb::MCParticle ThisPar = (*Assn.at(L1).at(L2));
          MyMap[abs(ThisPar.TrackId())] = ThisPar;
        }
    }
  return;
}

// Function to check if a given trackID is in a given map 
bool dune::EventAnaTree::InMap( int TrID, std::map<int, simb::MCParticle> ParMap)
{
  std::map<int, simb::MCParticle>::iterator Particle;
  Particle = ParMap.find(TrID);
  if (Particle != ParMap.end())
    {
      return true;
    }
  else
    return false;
}

bool dune::EventAnaTree::InMap( int TrID, std::map<int, float> TrackIDMap)
{
  std::map<int, float>::iterator Particle;
  Particle = TrackIDMap.find(TrID);
  if (Particle != TrackIDMap.end())
    {
      return true;
    }
  else
    return false;
}


DEFINE_ART_MODULE(dune::EventAnaTree)
