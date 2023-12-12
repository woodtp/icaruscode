////////////////////////////////////////////////////////////////////////
// Class:       SliceHitPurityAna
// Plugin Type: analyzer (Unknown Unknown)
// File:        SliceHitPurityAna_module.cc
//
// Generated at Tue Oct 31 12:56:33 2023 by Anthony Wood using cetskelgen
// from  version .
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art_root_io/TFileService.h"

// #include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// #include "lardataobj/AnalysisBase/BackTrackerMatchingData.h"
// #include "lardataobj/MCBase/MCParticleLite.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/RecoBase/Vertex.h"

#include "TTree.h"

#include <cstring>

namespace SliceHitPurity {
  class SliceHitPurityAna : public art::EDAnalyzer {
  public:
    explicit SliceHitPurityAna(fhicl::ParameterSet const& p);
    // The compiler-generated destructor is fine for non-base
    // classes without bare pointers or other resource use.

    // Plugins should not be copied or assigned.
    SliceHitPurityAna(SliceHitPurityAna const&) = delete;
    SliceHitPurityAna(SliceHitPurityAna&&) = delete;
    SliceHitPurityAna& operator=(SliceHitPurityAna const&) = delete;
    SliceHitPurityAna& operator=(SliceHitPurityAna&&) = delete;

    // Required functions.
    void analyze(art::Event const& e) override;

  private:
    enum WirePlane { IND1 = 0, IND2 = 1, COLL = 2 };
    static constexpr const char* fModuleName = "SliceHitPurityAna";
    static constexpr std::array<const char*, 2> fCryostats{"E", "W"};
    inline static const std::string fSliceLabel = "pandoraGausCryo";
    inline static const std::string fHitLabel = "cluster3DCryo";
    // inline static const std::string fMCLabel = "mcassociationsGausCryo";
    constexpr static double fZGap = 200.; // for avoiding the z-Gap, e.g., std::abs(z) < fZGap

    art::Handle<std::vector<recob::Slice>> fSliceHandle;

    std::unique_ptr<art::FindManyP<recob::Hit>> fFindManyHits;
    // std::unique_ptr<art::FindMany<simb::MCParticle, anab::BackTrackerHitMatchingData>>
    //   fFindManyMCParticles;
    std::unique_ptr<art::FindManyP<recob::PFParticle>> fFindManyPFParticles;
    std::unique_ptr<art::FindManyP<recob::Vertex>> fFindManyVertices;

    unsigned int fRun;
    unsigned int fSubRun;
    unsigned int fEvent;
    unsigned int fCryo;
    unsigned int fSliceNumber;
    unsigned int fNoHitsInSlice;
    unsigned int fMixCandidatesInSlice;
    unsigned int fTotalNoSlices;
    unsigned int fTotalNoSlicesWithMixing;

    double fSliceVtxX, fSliceVtxY, fSliceVtxZ;

    TTree* fTree;

    void InitTTree();
    void SetupFindManyPointers(art::Event const& evt, const std::string_view cryo);
  };

  SliceHitPurityAna::SliceHitPurityAna(fhicl::ParameterSet const& p)
    : EDAnalyzer{p}
    , fFindManyHits(nullptr)
    , fFindManyPFParticles(nullptr)
    , fFindManyVertices(nullptr)
    , fTotalNoSlices(0)
    , fTotalNoSlicesWithMixing(0)
  {
    InitTTree();
  }

  void SliceHitPurityAna::InitTTree()
  {
    art::ServiceHandle<art::TFileService> tfs;
    fTree = tfs->make<TTree>("sliceHitPurityAna", "Tree by SliceHitPurityAna_module");
    fTree->Branch("run", &fRun, "run/I");
    fTree->Branch("subrun", &fSubRun, "subrun/I");
    fTree->Branch("event", &fEvent, "event/I");
    fTree->Branch("cryostat", &fCryo, "cryostat/I");
    fTree->Branch("slice", &fSliceNumber, "slice/I");
    fTree->Branch("slice_vtx_x", &fSliceVtxX, "slice_vtx_x/D");
    fTree->Branch("slice_vtx_y", &fSliceVtxY, "slice_vtx_y/D");
    fTree->Branch("slice_vtx_z", &fSliceVtxZ, "slice_vtx_z/D");
    fTree->Branch("no_hits", &fNoHitsInSlice, "no_hits/I");
    fTree->Branch("no_mixing_candidates", &fMixCandidatesInSlice, "no_mixing_candidates/I");
  }

  void SliceHitPurityAna::SetupFindManyPointers(art::Event const& evt, const std::string_view cryo)
  {
    std::string label = fSliceLabel + cryo.data();
    // std::string mcLabel = fMCLabel + cryo.data();

    evt.getByLabel(fSliceLabel + cryo.data(), fSliceHandle);
    if (!fSliceHandle.isValid()) {
      mf::LogWarning(fModuleName) << "Unabled to locate recob::Slice in "
                                  << fSliceLabel + cryo.data() << '\n';
    }

    art::Handle<std::vector<recob::Hit>> hitHandle;
    evt.getByLabel(fHitLabel + cryo.data(), hitHandle);
    if (!hitHandle.isValid()) {
      mf::LogWarning(fModuleName) << "Unabled to locate recob::Hit in " << fHitLabel + cryo.data()
                                  << '\n';
    }

    art::Handle<std::vector<recob::PFParticle>> pfpHandle;
    evt.getByLabel(fSliceLabel + cryo.data(), pfpHandle);
    if (!pfpHandle.isValid()) {
      mf::LogWarning(fModuleName) << "Unabled to locate recob::Hit in " << fHitLabel + cryo.data()
                                  << '\n';
    }

    // fFindManyMCParticles =
    //   std::make_unique<art::FindMany<simb::MCParticle, anab::BackTrackerHitMatchingData>>(
    //     hitHandle, evt, mcLabel);

    fFindManyHits = std::make_unique<art::FindManyP<recob::Hit>>(fSliceHandle, evt, label);
    fFindManyVertices = std::make_unique<art::FindManyP<recob::Vertex>>(pfpHandle, evt, label);
    fFindManyPFParticles =
      std::make_unique<art::FindManyP<recob::PFParticle>>(fSliceHandle, evt, label);
  }

  void SliceHitPurityAna::analyze(art::Event const& evt)
  {
    fRun = evt.id().run();
    fSubRun = evt.id().subRun();
    fEvent = evt.id().event();

    int noSlices = 0;
    int noSlicesWithMixCand = 0;
    for (auto const& cryo : fCryostats) {
      SetupFindManyPointers(evt, cryo);

      auto const& slices = *fSliceHandle;

      for (auto const& slc : slices) {
        const int sliceId = slc.ID();

        auto const& pfps = fFindManyPFParticles->at(sliceId);
        if (pfps.empty()) { continue; }
        auto const& prim =
          std::find_if(pfps.begin(), pfps.end(), [](const art::Ptr<recob::PFParticle>& pfp) {
            return pfp->IsPrimary();
          });
        auto const& vertices = fFindManyVertices->at(prim->key());

        if (vertices.size() != 1) { continue; }

        auto const vtx = TVector3(
          vertices[0]->position().X(), vertices[0]->position().Y(), vertices[0]->position().Z());

        if (std::abs(vtx.Z()) < fZGap) { continue; }

        auto const& hits = fFindManyHits->at(sliceId);

        if (hits.empty()) {
          mf::LogInfo(fModuleName) << "No hits in slice " << sliceId << '\n';
          continue;
        }

        auto const noHits = hits.size();

        unsigned int mixCandidates = 0;
        for (auto prevIt = hits.begin(); std::distance(prevIt, hits.end()) >= 3; ++prevIt) {

          auto const it = prevIt + 1;
          auto const nextIt = prevIt + 2;

          auto const& hit = *it;
          auto const& prevHit = *prevIt;
          auto const& nextHit = *nextIt;

          auto const& wireID = hit->WireID();
          const int planeId = wireID.asPlaneID().deepestIndex();

          if (planeId != WirePlane::IND1) { continue; }

          const int tpcId = wireID.asTPCID().deepestIndex();
          const int prevTPCId = prevHit->WireID().asTPCID().deepestIndex();
          const int lastTPCId = nextHit->WireID().asTPCID().deepestIndex();

          if ((tpcId != prevTPCId) && (prevTPCId == lastTPCId) &&
              (prevHit->Channel() + 1 == nextHit->Channel())) {
            mixCandidates++;
          }

        } // Hits

        noSlices++;
        if (mixCandidates > 0) { noSlicesWithMixCand++; }

        fCryo = strcmp(cryo, "E") == 0 ? 0 : 1;
        fSliceNumber = sliceId;
        fSliceVtxX = vtx.X();
        fSliceVtxY = vtx.Y();
        fSliceVtxZ = vtx.Z();
        fNoHitsInSlice = noHits;
        fMixCandidatesInSlice = mixCandidates;

        fTree->Fill();

        const double candsPerHits = mixCandidates / (double)noHits;

        mf::LogInfo(fModuleName) << "Slice Number " << sliceId << '\n'
                                 << "mix candidates " << mixCandidates << " / " << noHits << " = "
                                 << candsPerHits << '\n'
                                 << "Slice Vtx (" << vtx.X() << ", " << vtx.Y() << ", " << vtx.Z()
                                 << ")\n";
      } // Slices
    }   // East/West Cryostat

    fTotalNoSlices += noSlices;
    fTotalNoSlicesWithMixing += noSlicesWithMixCand;

    if (fTotalNoSlices > 0) {
      mf::LogInfo(fModuleName) << " Slices with mix / Total Slices = " << fTotalNoSlicesWithMixing
                               << " / " << fTotalNoSlices << " = "
                               << fTotalNoSlicesWithMixing / (double)fTotalNoSlices << '\n';
    }
    else {
      mf::LogInfo(fModuleName) << "No valid slices so far.\n";
    }
  } // end Analyze()

} // namespace SliceHitPurity

DEFINE_ART_MODULE(SliceHitPurity::SliceHitPurityAna)
