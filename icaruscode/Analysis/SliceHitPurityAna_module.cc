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

#include "lardataobj/AnalysisBase/BackTrackerMatchingData.h"
#include "lardataobj/MCBase/MCParticleLite.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Slice.h"

#include <iostream>
#include <larcoreobj/SimpleTypesAndConstants/geo_types.h>
#include <lardataobj/RecoBase/PFParticle.h>
#include <lardataobj/RecoBase/Vertex.h>
#include <map>
#include <set>
#include <string_view>
#include <unordered_map>
#include <unordered_set>

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
    enum WirePlane { IND1 = 0, IND2 = 1, COLL = 2 };
    static constexpr std::string_view fModuleName = "SliceHitPurityAna";
    const std::array<std::string_view, 2> fCryostats{"E", "W"};
    const std::string fSliceLabel = "pandoraGausCryo";
    const std::string fHitLabel = "cluster3DCryo";
    const std::string fMCLabel = "mcassociationsGausCryo";

    std::vector<int> fRuns;
    std::vector<int> fSubRuns;
    std::vector<int> fEvents;
    std::vector<int> fSliceNumbers;
    std::vector<bool> fIsCosmicSlice;

    std::vector<int> fNHitsInSlice;
    std::vector<int> fUnmatchedHits;
    std::vector<int> fDuplicates;
    std::vector<int> fMatchedHits;
    std::vector<int> fMaxMCParticleMatches;
    std::vector<int> fMultiMatchedHits;

    std::vector<double> fUnmatchedHitFraction;
    std::vector<double> fMatchedHitFraction;
    std::vector<double> fMultiMatchedHitFractionContainingDupes;
    std::vector<double> fMultiMatchedHitFraction;
    std::unordered_map<int, int> fPdgSetMap;

    int fRun, fSubRun, fEvent;

    TTree* fTree;

    void InitTTree();
    void DuplicateStudy(
      const std::vector<art::Ptr<recob::Hit>>& hits,
      const art::FindMany<simb::MCParticle, anab::BackTrackerHitMatchingData>& fmp,
      const int sliceID);
    static bool IsProbablyACosmic(const simb::MCParticle& particle);
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
    fTree->Branch("run", &fRun, "run/I");
    fTree->Branch("subrun", &fSubRun, "subrun/I");
    fTree->Branch("event", &fEvent, "event/I");

    fTree->Branch("slice", "std::vector<int>", &fSliceNumbers);
    fTree->Branch("is_cosmic_slice", "std::vector<bool>", &fIsCosmicSlice);
    fTree->Branch("nhits", "std::vector<int>", &fNHitsInSlice);
    fTree->Branch("unmatched_hits", "std::vector<int>", &fUnmatchedHits);
    fTree->Branch("matched_hits", "std::vector<int>", &fMatchedHits);
    fTree->Branch("n_duplicates", "std::vector<int>", &fDuplicates);
    fTree->Branch("max_mcparticle_matches", "std::vector<int>", &fMaxMCParticleMatches);
    fTree->Branch("multiple_matched_hits", "std::vector<int>", &fMultiMatchedHits);
    fTree->Branch("unmatched_fraction", "std::vector<double>", &fUnmatchedHitFraction);
    fTree->Branch("matched_fraction", "std::vector<double>", &fMatchedHitFraction);
    fTree->Branch("multimatched_fraction", "std::vector<double>", &fMultiMatchedHitFraction);
    fTree->Branch("multimatched_fraction_containing_dupes",
                  "std::vector<double>",
                  &fMultiMatchedHitFractionContainingDupes);
  }

  bool SliceHitPurityAna::IsProbablyACosmic(const simb::MCParticle& particle)
  {
    const bool isMuon = std::abs(particle.PdgCode());
    const bool generatedOutsideVolume = std::fabs(particle.Vx()) > 360. ||
                                        std::fabs(particle.Vz()) > 900. || particle.Vy() > 135. ||
                                        particle.Vy() < -185.;

    return isMuon && generatedOutsideVolume;
  }

  void SliceHitPurityAna::DuplicateStudy(
    const std::vector<art::Ptr<recob::Hit>>& hits,
    const art::FindMany<simb::MCParticle, anab::BackTrackerHitMatchingData>& fmp,
    const int sliceID)
  {
    auto const nHits = hits.size();
    int hit_count = 0;

    unsigned int noMatches = 0;
    unsigned int hitsMatchedToOneMCP = 0;
    unsigned int hitsMatchedToMultipleMCP = 0;
    unsigned int maxMCParticleMatches = 0;
    unsigned int duplicates = 0;

    std::unordered_set<int> mcpIds;
    std::vector<bool> cosmicTags;

    for (auto const& hit : hits) {
      auto const& mcParts = fmp.at(hit.key());

      const unsigned int nMCP = mcParts.size();

      if (nMCP > maxMCParticleMatches) maxMCParticleMatches = nMCP;

      if (nMCP == 0) { noMatches++; }
      else if (nMCP == 1) {
        hitsMatchedToOneMCP++;
      }
      else {
        bool atLeastOneDuplicate = false;
        int i = 0;
        std::cout << "Hit " << hit_count << '\n';
        for (auto const& p : mcParts) {
          auto const& d = fmp.data(hit.key()).at(i);
          auto const id = p->TrackId();
          auto const pdg = p->PdgCode();
          std::cout << "    TrackId " << id << '\n';
          std::cout << "    PDG " << pdg << '\n';
          std::cout << "    IDE Fraction " << d->ideFraction << '\n';
          std::cout << "    isMaxIDE " << d->isMaxIDE << '\n';
          std::cout << "    numElectrons " << d->numElectrons << '\n';
          std::cout << "    energy " << d->energy << "\n\n";
          i++;
          if (mcpIds.find(id) == mcpIds.end()) {
            mcpIds.insert(id);
            if (fPdgSetMap.find(pdg) == fPdgSetMap.end()) fPdgSetMap[pdg] = 0;
          }
          else {
            fPdgSetMap[pdg]++;
            atLeastOneDuplicate = true;
          }
          cosmicTags.push_back(IsProbablyACosmic(*p));
        }
        if (atLeastOneDuplicate) duplicates++;
        hitsMatchedToMultipleMCP++;
      }
      hit_count++;
    }
    const double unmatchedFraction = nHits > 0 ? noMatches / (double)nHits : -1.;
    const double matchedFraction = nHits > 0 ? hitsMatchedToOneMCP / (double)nHits : -1.;
    const double multimatchedFraction = nHits > 0 ? hitsMatchedToMultipleMCP / (double)nHits : -1.;
    const double multimatchedFractionContainingDupes =
      (nHits > 0) && (hitsMatchedToMultipleMCP > 0) ?
        duplicates / (double)hitsMatchedToMultipleMCP :
        -1.;

    const long unsigned int cosmicCount = std::count(cosmicTags.begin(), cosmicTags.end(), true);

    fIsCosmicSlice.push_back(cosmicCount > cosmicTags.size() / 2);
    fSliceNumbers.push_back(sliceID);
    fNHitsInSlice.push_back(nHits);

    fUnmatchedHits.push_back(noMatches);
    fMatchedHits.push_back(hitsMatchedToOneMCP);
    fMultiMatchedHits.push_back(hitsMatchedToMultipleMCP);

    fUnmatchedHitFraction.push_back(unmatchedFraction);
    fMatchedHitFraction.push_back(matchedFraction);
    fMultiMatchedHitFraction.push_back(multimatchedFraction);
    fMultiMatchedHitFractionContainingDupes.push_back(multimatchedFractionContainingDupes);

    fMaxMCParticleMatches.push_back(maxMCParticleMatches);
    fDuplicates.push_back(duplicates);
  }

  void SliceHitPurityAna::analyze(art::Event const& evt)
  {
    for (auto const& cryo : fCryostats) {
      std::string label = fSliceLabel + cryo.data();
      std::string mcLabel = fMCLabel + cryo.data();

      art::Handle<std::vector<recob::Slice>> slcHandle;
      evt.getByLabel(fSliceLabel + cryo.data(), slcHandle);
      if (!slcHandle.isValid()) {
        mf::LogWarning(fModuleName.data())
          << "Unabled to locate recob::Slice in " << fSliceLabel + cryo.data() << '\n';
      }

      art::Handle<std::vector<recob::Hit>> hitHandle;
      evt.getByLabel(fHitLabel + cryo.data(), hitHandle);
      if (!hitHandle.isValid()) {
        mf::LogWarning(fModuleName.data())
          << "Unabled to locate recob::Hit in " << fHitLabel + cryo.data() << '\n';
      }

      art::Handle<std::vector<recob::PFParticle>> pfpHandle;
      evt.getByLabel(fSliceLabel + cryo.data(), pfpHandle);
      if (!pfpHandle.isValid()) {
        mf::LogWarning(fModuleName.data())
          << "Unabled to locate recob::Hit in " << fHitLabel + cryo.data() << '\n';
      }

      art::FindMany<simb::MCParticle, anab::BackTrackerHitMatchingData> fmp(
        hitHandle, evt, mcLabel);

      art::FindManyP<recob::Hit> fmh(slcHandle, evt, label);
      art::FindManyP<recob::Vertex> fmvtx(pfpHandle, evt, label);
      art::FindManyP<recob::PFParticle> fmpfp(slcHandle, evt, label);

      std::map<int, std::set<int>> sliceToMCTrackId;
      std::map<int, std::map<int, int>> sliceToMCTrackIdAndHits;
      std::map<int, std::size_t> sliceToTotalHits;
      std::map<int, TVector3> sliceToVtx;
      // std::map<int,
      std::map<int, int> trackIdToPdg;
      std::map<int, std::map<int, std::set<int>>> trackIdSeenInViews;
      // make sure view 0
      // true start/stop z of MCParticles and see if they overlap
      for (auto const& slc : *slcHandle) {
        const int sliceId = slc.ID();
        auto const& hits = fmh.at(sliceId);
        auto const& pfps = fmpfp.at(sliceId);
        if (pfps.empty()) { continue ; }
        auto const& prim =
          std::find_if(pfps.begin(), pfps.end(), [](const art::Ptr<recob::PFParticle>& pfp) {
            return pfp->IsPrimary();
          });
        auto const& vertices = fmvtx.at(prim->key());

        if (vertices.size() == 1) {
          sliceToVtx[sliceId] = TVector3(
            vertices[0]->position().X(), vertices[0]->position().Y(), vertices[0]->position().Z());
        }

        // sliceToTotalHits[sliceId] = hits.size();
        // mf::LogInfo(fModuleName.data()) << "============= Slice " << sliceID << " (" << cryo << ")"
        // << " =============\n";
        std::map<int, std::set<int>> planeTPCSets;

        // DuplicateStudy(hits, fmp, sliceID);
        //
        // Make a pair of TPCID and nhits to calculate fraction of mixing cand?

        std::map<int, int> trackIdToNHits;
        // std::map<int, std::set<WirePlane>> seenInViews;
        std::set<int> uniqueTrackIds;
        sliceToTotalHits[sliceId] = 0;
        for (auto const& hit : hits) {
          auto const& wireID = hit->WireID();
          auto const& planeID = wireID.asPlaneID().deepestIndex();
          auto const& tpcID = wireID.asTPCID().deepestIndex();

          // if (planeID != WirePlane::IND1) { continue; }

          sliceToTotalHits[sliceId]++;

          planeTPCSets[planeID].insert(tpcID);

          // double matchedFraction = 0.;
          auto const& mcParts = fmp.at(hit.key());
          for (std::size_t i = 0; i < mcParts.size(); ++i) {
            // const auto& md = fmp.data(hit.key()).at(i);

            const int trackId = mcParts[i]->TrackId();
            if (uniqueTrackIds.find(trackId) == uniqueTrackIds.end()) {
              uniqueTrackIds.insert(trackId);
              sliceToMCTrackIdAndHits[sliceId][trackId] = 0;
              trackIdToPdg[trackId] = mcParts[i]->PdgCode();
            }
            trackIdSeenInViews[sliceId][trackId].insert(planeID);
          }
          for (auto const& id : uniqueTrackIds) {
            sliceToMCTrackIdAndHits[sliceId][id]++;
          }
          sliceToMCTrackId[sliceId] = uniqueTrackIds;
        }
      } // Slice

      for (auto it1 = sliceToMCTrackId.begin(); it1 != sliceToMCTrackId.end(); ++it1) {
        for (auto it2 = std::next(it1); it2 != sliceToMCTrackId.end(); ++it2) {
          const auto& slc1 = it1->first;
          const auto& slc2 = it2->first;

          std::set<int> commonIds;

          std::set_intersection(it1->second.begin(),
                                it1->second.end(),
                                it2->second.begin(),
                                it2->second.end(),
                                std::inserter(commonIds, commonIds.begin()));

          if (commonIds.empty()) { continue; }

          std::cout << "Slices " << slc1 << ", " << slc2 << '\n';
          std::cout << "  { ";
          for (auto const& id : it1->second) {
            std::cout << id << ' ';
          }
          std::cout << " }, { ";
          for (auto const& id : it2->second) {
            std::cout << id << ' ';
          }
          std::cout << " }\n  Common:\n";

          for (auto const& id : commonIds) {
            const auto frac1 = sliceToMCTrackIdAndHits[slc1][id] / (double)sliceToTotalHits[slc1];
            const auto frac2 = sliceToMCTrackIdAndHits[slc2][id] / (double)sliceToTotalHits[slc2];
            std::cout << "    trackId " << id << " (pdg = " << trackIdToPdg[id] << ")\n";
            std::cout << "    (Slice " << slc1 << ") Position (" << sliceToVtx[slc1].X() << ", "
                      << sliceToVtx[slc1].Y() << ", " << sliceToVtx[slc1].Z() << ")\n";
            std::cout << "    (Slice " << slc2 << ") Position (" << sliceToVtx[slc2].X() << ", "
                      << sliceToVtx[slc2].Y() << ", " << sliceToVtx[slc2].Z() << ")\n";
            std::cout << "    (Slice " << slc1 << ") hits fraction "
                      << sliceToMCTrackIdAndHits[slc1][id] << " / " << sliceToTotalHits[slc1]
                      << " = " << frac1 << '\n';
            std::cout << "    (Slice " << slc2 << ") hits fraction "
                      << sliceToMCTrackIdAndHits[slc2][id] << " / " << sliceToTotalHits[slc2]
                      << " = " << frac2 << '\n';
            const auto& seenViews1 = trackIdSeenInViews[slc1][id];
            const auto& seenViews2 = trackIdSeenInViews[slc2][id];
            std::cout << "    (Slice " << slc1 << ") This id seen in planes { ";
            for (auto const& view : seenViews1) {
              std::cout << view << ' ';
            }
            std::cout << "}\n";
            std::cout << "    (Slice " << slc2 << ") This id seen in planes { ";
            for (auto const& view : seenViews2) {
              std::cout << view << ' ';
            }
            std::cout << "}\n";
          }
        }
      }
    } // East/West Cryostat

    fRun = evt.id().run();
    fSubRun = evt.id().subRun();
    fEvent = evt.id().event();
    fTree->Fill();

    mf::LogInfo(fModuleName.data()) << "Done." << std::endl;
  } // end Analyze()

} // namespace SliceHitPurity

DEFINE_ART_MODULE(SliceHitPurity::SliceHitPurityAna)
