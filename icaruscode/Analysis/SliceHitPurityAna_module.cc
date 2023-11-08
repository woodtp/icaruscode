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

#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "lardataobj/MCBase/MCParticleLite.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Slice.h"

#include <iostream>

#include "TTree.h"

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
    // Declare member data here.
    const std::string fModuleName = "SliceHitPurityAna";
    const std::array<const char*, 2> fCryostats{"E", "W"};
    const std::string fSliceLabel = "pandoraGausCryo";
    const std::string fHitLabel = "cluster3DCryo";
    const std::string fMCLabel = "mcassociationsGausCryo";

    std::vector<int> fRuns;
    std::vector<int> fSubRuns;
    std::vector<int> fEvents;
    std::vector<int> fSliceNumbers;
    std::vector<int> fNHitsInSlice;
    std::vector<int> fUnmatchedHits;
    std::vector<int> fMatchedHits;
    std::vector<int> fMaxMCParticleMatches;
    std::vector<int> fMultiMatchedHits;

    std::vector<double> fUnmatchedHitFraction;
    std::vector<double> fMatchedHitFraction;
    std::vector<double> fMultiMatchedHitFraction;

    int fRun, fSubRun, fEvent;

    TTree* fTree;

    void InitTTree();
  };

  SliceHitPurityAna::SliceHitPurityAna(fhicl::ParameterSet const& p) : EDAnalyzer{p}
  {
    // Call appropriate consumes<>() for any products to be retrieved by this
    // TODO: what?
    InitTTree();
  }

  void SliceHitPurityAna::InitTTree()
  {
    art::ServiceHandle<art::TFileService> tfs;
    fTree = tfs->make<TTree>("sliceHitPurityAna", "Tree by SliceHitPurityAna_module");
    // fTree->Branch("run", "std::vector<int>", &fRuns);
    // fTree->Branch("subrun", "std::vector<int>", &fSubRuns);
    // fTree->Branch("event", "std::vector<int>", &fEvents);
    fTree->Branch("run", &fRun, "run/I");
    fTree->Branch("subrun", &fSubRun, "subrun/I");
    fTree->Branch("event", &fEvent, "event/I");

    fTree->Branch("slice", "std::vector<int>", &fSliceNumbers);
    fTree->Branch("nhits", "std::vector<int>", &fNHitsInSlice);
    fTree->Branch("unmatched_hits", "std::vector<int>", &fUnmatchedHits);
    fTree->Branch("matched_hits", "std::vector<int>", &fMatchedHits);
    fTree->Branch("max_mcparticle_matches", "std::vector<int>", &fMaxMCParticleMatches);
    fTree->Branch("multiple_matched_hits", "std::vector<int>", &fMultiMatchedHits);
    fTree->Branch("unmatched_fraction", "std::vector<double>", &fUnmatchedHitFraction);
    fTree->Branch("matched_fraction", "std::vector<double>", &fMatchedHitFraction);
    fTree->Branch("multimatched_fraction", "std::vector<double>", &fMultiMatchedHitFraction);

    // fTree->Branch("pdg", "std::vector<int>", &fMatchedPDG);
    // fTree->Branch("vertex_x", "std::vector<double>", &fMatchedVx);
    // fTree->Branch("vertex_y", "std::vector<double>", &fMatchedVy);
    // fTree->Branch("vertex_z", "std::vector<double>", &fMatchedVz);
    //
    // fTree->Branch("end_x", "std::vector<double>", &fMatchedEndX);
    // fTree->Branch("end_y", "std::vector<double>", &fMatchedEndY);
    // fTree->Branch("end_z", "std::vector<double>", &fMatchedEndZ);
  }

  void SliceHitPurityAna::analyze(art::Event const& evt)
  {

    for (auto const& cryo : fCryostats) {
      std::string label = fSliceLabel + cryo;
      std::string mcLabel = fMCLabel + cryo;

      auto sliceHandles = evt.getMany<std::vector<recob::Slice>>();

      art::Handle<std::vector<recob::Hit>> hitHandle;
      evt.getByLabel(fHitLabel + cryo, hitHandle);
      if (!hitHandle.isValid()) {
        mf::LogWarning(fModuleName)
          << "Unabled to locate recob::Hit in 'cluster3DCryo" << cryo << "'\n";
      }
      art::FindMany<simb::MCParticle> fmp(hitHandle, evt, mcLabel);

      for (auto const& slcHandle : sliceHandles) {
        if (!slcHandle.isValid()) continue;

        art::FindManyP<recob::Hit> fmh(slcHandle, evt, label);

        for (auto const& slc : *slcHandle) {
          const int sliceID = slc.ID();
          auto const hits = fmh.at(sliceID);

          const double nHits = hits.size();

          unsigned int noMatches = 0;
          unsigned int hitsMatchedToOneMCP = 0;
          unsigned int hitsMatchedToMultipleMCP = 0;
          unsigned int maxMCParticleMatches = 0;
          for (auto const& hit : hits) {
            auto const parts = fmp.at(hit.key());

            const unsigned int nMCP = parts.size();

            if (nMCP > maxMCParticleMatches) maxMCParticleMatches = nMCP;

            if (nMCP == 0) { noMatches++; }
            else if (nMCP == 1) {
              hitsMatchedToOneMCP++;
            }
            else {
              if (nMCP > 2) {
                std::cout << "Cryo " << cryo << '\n';
                std::cout << "Slice " << sliceID << '\n';
                std::cout << "nMCP " << nMCP << '\n';
                for (auto const& p : parts) {
                  std::cout << "TrackID " << p->TrackId() << " PDG " << p->PdgCode() << " (" << p->Vx() << ", " << p->Vy() << ", "
                            << p->Vz() << ") ";
                  std::cout << "(" << p->EndX() << ", " << p->EndY() << ", " << p->EndZ()
                            << ")\n\n";
                }
              }
              hitsMatchedToMultipleMCP++;
            }
          }
          const double unmatchedFraction = nHits > 0 ? noMatches / nHits : -1.;
          const double matchedFraction = nHits > 0 ? hitsMatchedToOneMCP / nHits : -1.;
          const double multimatchedFraction = nHits > 0 ? hitsMatchedToMultipleMCP / nHits : -1.;

          fSliceNumbers.push_back(sliceID);
          fNHitsInSlice.push_back(nHits);

          fUnmatchedHits.push_back(noMatches);
          fMatchedHits.push_back(hitsMatchedToOneMCP);
          fMultiMatchedHits.push_back(hitsMatchedToMultipleMCP);

          fUnmatchedHitFraction.push_back(unmatchedFraction);
          fMatchedHitFraction.push_back(matchedFraction);
          fMultiMatchedHitFraction.push_back(multimatchedFraction);

          fMaxMCParticleMatches.push_back(maxMCParticleMatches);
        }
      }
    }

    fRun = evt.id().run();
    fSubRun = evt.id().subRun();
    fEvent = evt.id().event();
    fTree->Fill();

    mf::LogInfo(fModuleName) << "Done." << std::endl;
  }
} // namespace SliceHitPurity

DEFINE_ART_MODULE(SliceHitPurity::SliceHitPurityAna)
