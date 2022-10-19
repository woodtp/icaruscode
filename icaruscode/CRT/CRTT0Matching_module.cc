/////////////////////////////////////////////////////////////////////////////
/// Class:       CRTT0Matching
/// Module Type: producer
/// File:        CRTT0Matching_module.cc
///
/// Author:         Thomas Brooks
/// E-mail address: tbrooks@fnal.gov
///
/// Modified from CRTT0Matching by Thomas Warburton.
/////////////////////////////////////////////////////////////////////////////

// sbndcode includes
#include "sbnobj/Common/CRT/CRTHit.hh"
#include "icaruscode/CRT/CRTUtils/CRTT0MatchAlg.h"

// Framework includes
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Event.h" 
#include "canvas/Persistency/Common/Ptr.h" 
#include "canvas/Persistency/Common/PtrVector.h" 
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"

#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include <memory>
#include <iostream>
#include <map>
#include <iterator>
#include <algorithm>

// LArSoft
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcorealg/Geometry/PlaneGeo.h"
#include "larcorealg/Geometry/WireGeo.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Track.h"
#include "larcore/CoreUtils/ServiceUtil.h" // lar::providerFrom()
#include "lardata/Utilities/AssociationUtil.h"
#include "lardata/DetectorInfoServices/LArPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardataobj/RawData/ExternalTrigger.h"
#include "larcoreobj/SimpleTypesAndConstants/PhysicalConstants.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "sbnobj/Common/Trigger/ExtraTriggerInfo.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "lardataobj/RecoBase/PFParticleMetadata.h"
// ROOT
#include "TVector3.h"
#include "TH1.h"
#include "TH2.h"
#include "TVector3.h"
#include "TTree.h"

namespace icarus {
  
  class CRTT0Matching : public art::EDProducer {
  public:

    explicit CRTT0Matching(fhicl::ParameterSet const & p);

    // The destructor generated by the compiler is fine for classes
    // without bare pointers or other resource use.

    // Plugins should not be copied or assigned.
    CRTT0Matching(CRTT0Matching const &) = delete;
    CRTT0Matching(CRTT0Matching &&) = delete;
    CRTT0Matching & operator = (CRTT0Matching const &) = delete; 
    CRTT0Matching & operator = (CRTT0Matching &&) = delete;

    // Required functions.
    void produce(art::Event & e) override;

    // Selected optional functions.
    void beginJob() override;

    void endJob() override;

    void reconfigure(fhicl::ParameterSet const & p);

  private:

    // Params got from fcl file.......
    //    art::InputTag fTpcTrackModuleLabel; ///< name of track producer
    std::vector<art::InputTag> fTpcTrackModuleLabel; ///< name of track producer
    std::vector<art::InputTag> fPFParticleLabel; ///< labels for source of PFParticle
    art::InputTag              fCrtHitModuleLabel;   ///< name of crt producer
    art::InputTag              fTriggerLabel;        ///< labels for trigger
    CRTT0MatchAlg              t0Alg;


    geo::GeometryCore const* fGeometryService;   ///< pointer to Geometry provider
    icarus::crt::CRTCommonUtils* fCrtutils;
    //  CRTCommonUtils* fCrtutils;

    TTree* fTree;
    int fEvent;        ///< number of the event being processed
    int fRun;          ///< number of the run being processed
    int fSubRun;       ///< number of the sub-run being processed
    vector<int>            fCrtRegion;    //CRT hit region code
    vector<double>            fDCA;    //CRT hit region code
    vector<double>            fDOL;    //CRT hit region code
    vector<double>            fT0;    //CRT hit region code
    vector<double>            ftpcx;    // x-position from tpc
    vector<double>            ftpcy;    // y-position from tpc
    vector<double>            ftpcz;    // z-position from tpc
    vector<double>            fcrtx;    // x-position from crt
    vector<double>            fcrty;    // y-position from crt
    vector<double>            fcrtz;    // z-position from crt
    vector<double>         fpandorat0;    ///< Track T0 based on Pandora (Cathode Crossing Track)
    vector<int>            fcryo;         ///< cryo number
    vector<int>            fntracks;      ///< total number of tracks
 
    //add trigger data product vars
    unsigned int m_gate_type;
    std::string  m_gate_name;
    uint64_t     m_trigger_timestamp;
    uint64_t     m_gate_start_timestamp;
    uint64_t     m_trigger_gate_diff;
    uint64_t     m_gate_crt_diff;


    // Histograms
    std::map<std::string, TH1F*> hDCA;
    std::map<std::string, TH1F*> hMatchDCA;
    std::map<std::string, TH1F*> hNoMatchDCA;

    std::map<std::string, TH1F*> hDoL;
    std::map<std::string, TH1F*> hMatchDoL;
    std::map<std::string, TH1F*> hNoMatchDoL;

    std::map<std::string, TH1F*> hT0;
    std::map<std::string, TH1F*> hMatchT0;
    std::map<std::string, TH1F*> hNoMatchT0;
  }; // class CRTT0Matching


  CRTT0Matching::CRTT0Matching(fhicl::ParameterSet const & p)
    : EDProducer(p), t0Alg(p.get<fhicl::ParameterSet>("T0Alg"))
    , fCrtutils(new icarus::crt::CRTCommonUtils())
      // Initialize member data here, if know don't want to reconfigure on the fly
  {

    // Call appropriate produces<>() functions here.
    produces< std::vector<anab::T0>                  >();
    produces< art::Assns<recob::Track , anab::T0>    >();
    produces< art::Assns<sbn::crt::CRTHit, anab::T0> >();

    fGeometryService = lar::providerFrom<geo::Geometry>();
    reconfigure(p);

  } // CRTT0Matching()


  void CRTT0Matching::reconfigure(fhicl::ParameterSet const & p)
  {
    fTpcTrackModuleLabel = p.get< std::vector<art::InputTag>>("TpcTrackModuleLabel", {"pandoraTrackGausCryoE"});
    fCrtHitModuleLabel   = p.get<art::InputTag> ("CrtHitModuleLabel", "crthit"); 
    fPFParticleLabel    =  p.get< std::vector<art::InputTag> >("PFParticleLabel",             {""});
    fTriggerLabel        = p.get<art::InputTag>("TriggerLabel","daqTrigger");
  } // CRTT0Matching::reconfigure()


  void CRTT0Matching::beginJob()
  {
    // Access tfileservice to handle creating and writing histograms
    art::ServiceHandle<art::TFileService> tfs;

    fTree = tfs->make<TTree>("matchTree","CRTHit - TPC track matching analysis");

    fTree->Branch("Event",           &fEvent,           "Event/I");
    fTree->Branch("SubRun",          &fSubRun,          "SubRun/I");
    fTree->Branch("Run",             &fRun,             "Run/I");
    fTree->Branch("crtRegion",    "std::vector<int>",   &fCrtRegion);
    fTree->Branch("DCA",    "std::vector<double>",   &fDCA);
    fTree->Branch("DOL",    "std::vector<double>",   &fDOL);
    fTree->Branch("t0",     "std::vector<double>",   &fT0);
    fTree->Branch("tpcx",     "std::vector<double>",   &ftpcx);
    fTree->Branch("tpcy",     "std::vector<double>",   &ftpcy);
    fTree->Branch("tpcz",     "std::vector<double>",   &ftpcz);
    fTree->Branch("crtx",     "std::vector<double>",   &fcrtx);
    fTree->Branch("crty",     "std::vector<double>",   &fcrty);
    fTree->Branch("crtz",     "std::vector<double>",   &fcrtz);
    fTree->Branch("cryo",            "std::vector<int>",      &fcryo);
    fTree->Branch("ntracks",            "std::vector<int>",      &fntracks);
    fTree->Branch("pandorat0",       "std::vector<double>", &fpandorat0);
    fTree->Branch("gate_type", &m_gate_type, "gate_type/b");
    fTree->Branch("gate_name", &m_gate_name);
    fTree->Branch("trigger_timestamp", &m_trigger_timestamp, "trigger_timestamp/l");
    fTree->Branch("gate_start_timestamp", &m_gate_start_timestamp, "gate_start_timestamp/l");
    fTree->Branch("trigger_gate_diff", &m_trigger_gate_diff, "trigger_gate_diff/l");
    fTree->Branch("gate_crt_diff",&m_gate_crt_diff, "gate_crt_diff/l");

    for(int i = 30; i < 50 + 1; i++){
      std::string tagger = "All";
      if (i >= 35 && i < 40) continue;
      if (i==48 || i==49) continue;
      // if(i < ){
      tagger = fCrtutils->GetRegionNameFromNum(i);//fCrtGeo.GetTagger(i).name;
      //    std::cout << "tagger: " << tagger.c_str() << std::endl;
      hDCA[tagger]     = tfs->make<TH1F>(Form("DCA_%s", tagger.c_str()),        "", 50, 0, 100);
      hDoL[tagger]     = tfs->make<TH1F>(Form("DoL_%s", tagger.c_str()),        "", 100, 0, 0.25);
      hT0[tagger]      = tfs->make<TH1F>(Form("T0_%s", tagger.c_str()),        "", 600, -3000, 3000);
    }

  } // CRTT0Matching::beginJob()

  void CRTT0Matching::produce(art::Event & event)
  {
    fDCA.clear();
    fDOL.clear();
    fT0.clear();
    fcrtx.clear();
    fcrty.clear();
    fcrtz.clear();
    ftpcx.clear();
    ftpcy.clear();
    ftpcz.clear();
    fCrtRegion.clear();
    fpandorat0.clear();
    fcryo.clear();
    fntracks.clear();
    fEvent  = event.id().event();
    fRun    = event.run();
    fSubRun = event.subRun();

    // Create anab::T0 objects and make association with recob::Track
    std::unique_ptr< std::vector<anab::T0> > T0col( new std::vector<anab::T0>);
    std::unique_ptr< art::Assns<recob::Track, anab::T0> > Trackassn( new art::Assns<recob::Track, anab::T0>);
    std::unique_ptr< art::Assns <sbn::crt::CRTHit, anab::T0> > t0_crthit_assn( new art::Assns<sbn::crt::CRTHit, anab::T0> );

    //add trigger info
    if( !fTriggerLabel.empty() ) {

      art::Handle<sbn::ExtraTriggerInfo> trigger_handle;
      event.getByLabel( fTriggerLabel, trigger_handle );
      if( trigger_handle.isValid() ) {
	sbn::triggerSource bit = trigger_handle->sourceType;
	m_gate_type            = (unsigned int)bit;
	m_gate_name            = bitName(bit);
	m_trigger_timestamp    = trigger_handle->triggerTimestamp;
	m_gate_start_timestamp = trigger_handle->beamGateTimestamp;
	m_trigger_gate_diff    = trigger_handle->triggerTimestamp - trigger_handle->beamGateTimestamp;

      }
      else{
	mf::LogError("CRTT0Matching:") << "No raw::Trigger associated to label: " << fTriggerLabel.label() << "\n" ;
      }
    }
    else {
      mf::LogError("CRTT0Matching:") << "Trigger Data product " << fTriggerLabel.label() << " not found!\n" ;
    }

    // Retrieve CRT hit list
    art::Handle<std::vector<sbn::crt::CRTHit>> crtListHandle;
    std::vector<art::Ptr<sbn::crt::CRTHit>> crtList;
    if(event.getByLabel(fCrtHitModuleLabel, crtListHandle))
      art::fill_ptr_vector(crtList, crtListHandle);

    std::vector<sbn::crt::CRTHit> crtHits;
    for (auto const& crtHit : crtList){
      crtHits.push_back(*crtHit);
    }

    // Retrieve track list
    for(const auto& trackLabel : fTpcTrackModuleLabel){

      auto it = &trackLabel - fTpcTrackModuleLabel.data();

      art::Handle< std::vector<recob::Track> > trackListHandle;
      std::vector<art::Ptr<recob::Track> > trackList;
      if (event.getByLabel(trackLabel,trackListHandle))
      	art::fill_ptr_vector(trackList, trackListHandle);   

      //  if (event.getByLabel(fTpcTrackModuleLabel,trackListHandle))
      //    std::cout << "crtlabel: " << fCrtHitModuleLabel << " , Tpctrklabel: " << fTpcTrackModuleLabel << std::endl;

      mf::LogInfo("CRTT0Matching")
	<<"Number of reconstructed tracks = "<<trackList.size()<<"\n"
	<<"Number of CRT hits = "<<crtList.size();

      fntracks.push_back(trackList.size());
      //Get PFParticles
      auto pfpListHandle = event.getValidHandle<std::vector<recob::PFParticle> >(fPFParticleLabel[it]);
      if (!pfpListHandle.isValid()) continue;

      //Get PFParticle-Track association
      art::FindManyP<recob::PFParticle> fmpfp(trackListHandle, event, trackLabel);

      //Get T0-PFParticle association
      art::FindManyP<anab::T0> fmt0pandora(pfpListHandle, event, fPFParticleLabel[it]);

      
      if (trackListHandle.isValid() && crtListHandle.isValid() ){
	
	auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(event);
	art::FindManyP<recob::Hit> findManyHits(trackListHandle, event, trackLabel);

	// Loop over all the reconstructed tracks 
	for(size_t track_i = 0; track_i < trackList.size(); track_i++) {

	  double t0 = -99999999;
	  //Find PFParticle for track i
	  //art::Ptr::key() gives the index in the vector
	  auto pfps = fmpfp.at(trackList[track_i]->ID());

	  if (!pfps.empty()){
	    //Find T0 for PFParticle
	    auto t0s = fmt0pandora.at(pfps[0].key());
	    if (!t0s.empty()){

	      t0 = t0s[0]->Time();   //Get T0
	    }
	  }

	  std::vector<art::Ptr<recob::Hit>> hits = findManyHits.at(trackList[track_i]->ID());
	  if (hits.size() == 0) continue;
	  int const cryoNumber = hits[0]->WireID().Cryostat;
	  // std::pair<double, double> matchedTime = t0Alg.T0AndDCAFromCRTHits(detProp, *trackList[track_i], crtHits, event);
	  matchCand closest = t0Alg.GetClosestCRTHit(detProp, *trackList[track_i], hits, crtHits,  m_gate_start_timestamp, false);
	  // std::vector <matchCand> closestvec = t0Alg.GetClosestCRTHit(detProp, *trackList[track_i], crtHits, event);
	  // matchCand closest = closestvec.back();	  

	  if(closest.dca >=0 ){
	    mf::LogInfo("CRTT0Matching")
	      <<"Matched time = "<<closest.t0<<" [us] to track "<<trackList[track_i]->ID()<<" with DCA = "<<closest.dca;
	    T0col->push_back(anab::T0(closest.t0*1e3, trackList[track_i]->ID(),  closest.thishit.plane, (int)closest.extrapLen, closest.dca));
	    util::CreateAssn(*this, event, *T0col, trackList[track_i], *Trackassn);
	    
	    //std::cout << "---------------------- line #156 "  << std::endl;
	    double sin_angle = -99999;
	    if(closest.dca != -99999){
	      auto start = trackList[track_i]->Vertex<TVector3>();
	      //auto end   = trackList[track_i]->End<TVector3>();
	      hDCA[closest.thishit.tagger]->Fill(closest.dca);
	      //hDCA["All"]->Fill(closest.dca);
	      sin_angle = closest.dca/closest.extrapLen;
	      hDoL[closest.thishit.tagger]->Fill(sin_angle);
	      // hDoL["All"]->Fill(sin_angle);
	      hT0[closest.thishit.tagger]->Fill(closest.t0);
	      fDCA.push_back(closest.dca);
	      fDOL.push_back(sin_angle);
	      fT0.push_back(closest.t0);
	      fcryo.push_back(cryoNumber);
	      fpandorat0.push_back(t0);
	      ftpcx.push_back(start.X());
	      ftpcy.push_back(start.Y());
	      ftpcz.push_back(start.Z());
	      fcrtx.push_back(closest.thishit.x_pos);
	      fcrty.push_back(closest.thishit.y_pos);
	      fcrtz.push_back(closest.thishit.z_pos);
	      // fCrtRegion.push_back(closest.thishit.tagger);
	      fCrtRegion.push_back(fCrtutils->AuxDetRegionNameToNum(closest.thishit.tagger));	    
	    }
	      //find this CRThit in the collection
	    // note this does not work (both the loop and the assoc !!)
	    unsigned CRThitIndex = std::numeric_limits<unsigned>::max();
	    for (int ic=0; ic<(int)crtList.size(); ++ic){
	      if (crtList[ic]->ts0_ns==closest.thishit.ts0_ns && crtList[ic]->z_pos==closest.thishit.z_pos && crtList[ic]->peshit==closest.thishit.peshit)
		CRThitIndex=ic;
	      // std::cout << "CRThitIndex: \t" << CRThitIndex << std::endl;
	    }
	    if (CRThitIndex != std::numeric_limits<unsigned>::max()){
	      //  std::cout <<"CRThitIndex: " << CRThitIndex << "  passed......: \t"  << std::endl;
	      util::CreateAssn(*this, event, *T0col, crtList[CRThitIndex], *t0_crthit_assn);
	    }
	  } // DCA check
	  
	} // Loop over tracks  
	
      } // Validity check

    } // all track labels in a vector 
    event.put(std::move(T0col));
    event.put(std::move(Trackassn));
    event.put(std::move(t0_crthit_assn));
    fTree->Fill();
  } // CRTT0Matching::produce()


  void CRTT0Matching::endJob()
  {

  } // CRTT0Matching::endJob()


  DEFINE_ART_MODULE(CRTT0Matching)

} // sbnd namespace

namespace {

}
