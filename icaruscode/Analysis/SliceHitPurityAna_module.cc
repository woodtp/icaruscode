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
#include <memory>
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
    int fMixedCandidateCount = 0;
    enum WirePlane { IND1 = 0, IND2 = 1, COLL = 2 };
    static constexpr std::string_view fModuleName = "SliceHitPurityAna";
    const std::array<std::string_view, 2> fCryostats{"E", "W"};
    const std::string fSliceLabel = "pandoraGausCryo";
    const std::string fHitLabel = "cluster3DCryo";
    const std::string fMCLabel = "mcassociationsGausCryo";
    constexpr static double fZGap = 200.; // for avoiding the z-Gap, e.g., std::abs(z) < fZGap

    art::Handle<std::vector<recob::Slice>> fSliceHandle;

    std::unique_ptr<art::FindManyP<recob::Hit>> fFindManyHits;
    std::unique_ptr<art::FindMany<simb::MCParticle, anab::BackTrackerHitMatchingData>>
      fFindManyMCParticles;
    std::unique_ptr<art::FindManyP<recob::PFParticle>> fFindManyPFParticles;
    std::unique_ptr<art::FindManyP<recob::Vertex>> fFindManyVertices;

    struct SliceMaps {
      std::map<int, std::set<int>> sliceToMCTrackId;
      std::map<int, std::map<int, int>> sliceToMCTrackIdAndHits;
      std::map<int, std::size_t> sliceToTotalHits;
      std::map<int, TVector3> sliceToVtx;
      // std::map<int, int> trackIdToPdg;
      std::map<int, std::map<int, std::set<int>>> trackIdSeenInViews;
      std::map<int, const simb::MCParticle*> trackIdToMCParticle;
    };

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
    void SetupFindManyPointers(art::Event const& evt, std::string_view cryo);
    void DuplicateStudy(
      const std::vector<art::Ptr<recob::Hit>>& hits,
      const art::FindMany<simb::MCParticle, anab::BackTrackerHitMatchingData>& fmp,
      const int sliceID);
    static bool IsProbablyACosmic(const simb::MCParticle& particle);
    [[nodiscard]] int SliceHitStudy(const SliceMaps& slcMaps);
    void DumpToTerminal();
  };

  SliceHitPurityAna::SliceHitPurityAna(fhicl::ParameterSet const& p)
    : EDAnalyzer{p}
    , fFindManyHits(nullptr)
    , fFindManyMCParticles(nullptr)
    , fFindManyPFParticles(nullptr)
    , fFindManyVertices(nullptr)
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

  void SliceHitPurityAna::SetupFindManyPointers(art::Event const& evt, std::string_view cryo)
  {
    std::string label = fSliceLabel + cryo.data();
    std::string mcLabel = fMCLabel + cryo.data();

    evt.getByLabel(fSliceLabel + cryo.data(), fSliceHandle);
    if (!fSliceHandle.isValid()) {
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

    fFindManyMCParticles =
      std::make_unique<art::FindMany<simb::MCParticle, anab::BackTrackerHitMatchingData>>(
        hitHandle, evt, mcLabel);

    fFindManyHits = std::make_unique<art::FindManyP<recob::Hit>>(fSliceHandle, evt, label);
    fFindManyVertices = std::make_unique<art::FindManyP<recob::Vertex>>(pfpHandle, evt, label);
    fFindManyPFParticles =
      std::make_unique<art::FindManyP<recob::PFParticle>>(fSliceHandle, evt, label);
  }

  bool SliceHitPurityAna::IsProbablyACosmic(const simb::MCParticle& particle)
  {
    const bool isMuon = std::abs(particle.PdgCode());
    const bool generatedOutsideVolume = std::fabs(particle.Vx()) > 360. ||
                                        std::fabs(particle.Vz()) > 900. || particle.Vy() > 135. ||
                                        particle.Vy() < -185.;

    return isMuon && generatedOutsideVolume;
  }

  void SliceHitPurityAna::analyze(art::Event const& evt)
  {
    for (auto const& cryo : fCryostats) {
      SetupFindManyPointers(evt, cryo);

      SliceMaps slcMaps;

      for (auto const& slc : *fSliceHandle) {
        const int sliceId = slc.ID();
        auto const& hits = fFindManyHits->at(sliceId);
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

        slcMaps.sliceToVtx[sliceId] = vtx;

        std::map<int, std::set<int>> planeTPCSets;

        // DuplicateStudy(hits, fmp, sliceID);

        std::map<int, int> trackIdToNHits;
        std::set<int> uniqueTrackIds;
        slcMaps.sliceToTotalHits[sliceId] = 0;
        std::vector<geo::WireID> wires;
        for (auto const& hit : hits) {
          auto const& wireID = hit->WireID();
          auto const& planeID = wireID.asPlaneID().deepestIndex();
          auto const& tpcID = wireID.asTPCID().deepestIndex();

          if (planeID == WirePlane::IND1) { wires.push_back(wireID); }

          slcMaps.sliceToTotalHits[sliceId]++;

          planeTPCSets[planeID].insert(tpcID);

          // double matchedFraction = 0.;
          auto const& mcParts = fFindManyMCParticles->at(hit.key());
          for (std::size_t i = 0; i < mcParts.size(); ++i) {
            auto const& mcp = mcParts[i];
            if (std::abs(mcp->Vz()) < fZGap || std::abs(mcp->EndZ()) < fZGap) { continue; }
            // auto const& md = fmp.data(hit.key()).at(i);

            const int trackId = mcp->TrackId();
            if (uniqueTrackIds.find(trackId) == uniqueTrackIds.end()) {
              uniqueTrackIds.insert(trackId);
              slcMaps.sliceToMCTrackIdAndHits[sliceId][trackId] = 0;
              // slcMaps.trackIdToPdg[trackId] = mcp->PdgCode();
              slcMaps.trackIdToMCParticle[trackId] = mcp;
            }
            slcMaps.trackIdSeenInViews[sliceId][trackId].insert(planeID);
          }
          for (auto const& id : uniqueTrackIds) {
            slcMaps.sliceToMCTrackIdAndHits[sliceId][id]++;
          }
          slcMaps.sliceToMCTrackId[sliceId] = uniqueTrackIds;
        } // Hits
      }   // Slices
      const int count = SliceHitStudy(slcMaps);
      fMixedCandidateCount += count;
      mf::LogInfo(fModuleName.data()) << "Final count for slice: " << count << '\n';
    } // East/West Cryostat

    fRun = evt.id().run();
    fSubRun = evt.id().subRun();
    fEvent = evt.id().event();
    fTree->Fill();

    mf::LogInfo(fModuleName.data()) << "Done.\n"
                                    << "FINAL COUNT " << fMixedCandidateCount << std::endl;
  } // end Analyze()

  int SliceHitPurityAna::SliceHitStudy(const SliceMaps& slcMaps)
  {
    int count = 0;
    for (auto it1 = slcMaps.sliceToMCTrackId.begin(); it1 != slcMaps.sliceToMCTrackId.end();
         ++it1) {
      for (auto it2 = std::next(it1); it2 != slcMaps.sliceToMCTrackId.end(); ++it2) {
        auto const& slc1 = it1->first;
        auto const& slc2 = it2->first;
        auto const& trackIds1 = it1->second;
        auto const& trackIds2 = it2->second;

        std::set<int> commonIds;

        std::set_intersection(trackIds1.begin(),
                              trackIds1.end(),
                              trackIds2.begin(),
                              trackIds2.end(),
                              std::inserter(commonIds, commonIds.begin()));

        if (commonIds.empty()) { continue; }

        std::set<int> otherIds;
        std::set_symmetric_difference(trackIds1.begin(),
                                      trackIds1.end(),
                                      trackIds2.begin(),
                                      trackIds2.end(),
                                      std::inserter(otherIds, otherIds.begin()));

        for (auto const& id : commonIds) {
          auto const noAssocHits1 = slcMaps.sliceToMCTrackIdAndHits.at(slc1).at(id);
          auto const noAssocHits2 = slcMaps.sliceToMCTrackIdAndHits.at(slc2).at(id);

          auto const noTotalHits1 = slcMaps.sliceToTotalHits.at(slc1);
          auto const noTotalHits2 = slcMaps.sliceToTotalHits.at(slc2);

          auto const assocHitsFrac1 = noAssocHits1 / (double)noTotalHits1;
          auto const assocHitsFrac2 = noAssocHits2 / (double)noTotalHits2;

          auto const& seenViews1 = slcMaps.trackIdSeenInViews.at(slc1).at(id);
          auto const& seenViews2 = slcMaps.trackIdSeenInViews.at(slc2).at(id);

          const bool trackIdSlc1SeenOnce = seenViews1.size() == 1;
          const bool trackIdSlc2SeenOnce = seenViews2.size() == 1;

          auto const trackIdSeenInduction1Only =
            (trackIdSlc1SeenOnce && *seenViews1.begin() == 0) ||
            (trackIdSlc2SeenOnce && *seenViews2.begin() == 0);

          auto const otherTrackIdSeenInMultipleViews =
            trackIdSeenInduction1Only && (seenViews1.size() > 1 || seenViews2.size() > 1);

          if (!otherTrackIdSeenInMultipleViews) { continue; }

          auto const& multiViews = seenViews1.size() > 1 ? seenViews1 : seenViews2;

          bool isAnyInd1 = false;
          for (auto const& view : multiViews) {
            if (view == 0) { isAnyInd1 = true; }
          }

          if (!isAnyInd1) { continue; }

          auto const& singleViewMCP = slcMaps.trackIdToMCParticle.at(id);

          const double startX = singleViewMCP->Vx();
          const double endX = singleViewMCP->EndX();

          const double startZ = singleViewMCP->Vz();
          const double endZ = singleViewMCP->EndZ();

          const double maxX = std::max(startX, endX);
          const double maxZ = std::max(startZ, endZ);
          const double minX = std::min(startX, endX);
          const double minZ = std::min(startZ, endZ);

          for (auto const& otherId : otherIds) {
            if (otherId == id) { continue; }
            auto const& otherMCP = slcMaps.trackIdToMCParticle.at(otherId);
            const double otherStartX = otherMCP->Vx();
            const double otherEndX = otherMCP->EndX();

            const double otherStartZ = otherMCP->Vz();
            const double otherEndZ = otherMCP->EndZ();

            const double otherMaxX = std::max(otherStartX, otherEndX);
            const double otherMaxZ = std::max(otherStartZ, otherEndZ);

            const double otherMinX = std::min(otherStartX, otherEndX);
            const double otherMinZ = std::min(otherStartZ, otherEndZ);

            const double overlapX = std::min(maxX, otherMaxX) - std::max(minX, otherMinX);
            const double overlapZ = std::min(maxZ, otherMaxZ) - std::max(minZ, otherMinZ);

            const bool noOverlapX = overlapX < std::numeric_limits<double>::epsilon();
            const bool noOverlapZ = overlapZ < std::numeric_limits<double>::epsilon();

            // if doesn't overlap in x or does overlap in z, continue.
            if (noOverlapX || !noOverlapZ) { continue; }

            count++;

            const char* const tab = "    ";
            std::cout << "Slices " << slc1 << ", " << slc2 << '\n';
            std::cout << tab << "{ ";
            for (auto const& id : trackIds1) {
              std::cout << id << ' ';
            }

            std::cout << "}, { ";
            for (auto const& id : trackIds2) {
              std::cout << id << ' ';
            }
            std::cout << " }\n";

            std::cout << tab << "COMMON TrackId " << id << " (pdg = " << singleViewMCP->PdgCode()
                      << ")\n";

            std::cout << tab << "Start Stop X: { " << singleViewMCP->Vx() << ", "
                      << singleViewMCP->EndX() << " }\n";
            std::cout << tab << "Start Stop Z: { " << singleViewMCP->Vz() << ", "
                      << singleViewMCP->EndZ() << " }\n\n";

            std::cout << tab << "NO Z OVERLAP FOUND\n";
            std::cout << tab << "trackId " << otherId << " (pdg = " << otherMCP->PdgCode() << ")\n";
            std::cout << tab << "Start Stop X: { " << otherMCP->Vx() << ", " << otherMCP->EndX()
                      << " }\n";
            std::cout << tab << "Start Stop Z: { " << otherMCP->Vz() << ", " << otherMCP->EndZ()
                      << " }\n";

            std::cout << tab << "Overlap X? FALSE\n";
            std::cout << tab << "Overlap Z? TRUE\n\n";

            std::cout << tab << "Slice " << slc1 << "\n";

            std::cout << tab << tab << "Vtx Position (" << slcMaps.sliceToVtx.at(slc1).X() << ", "
                      << slcMaps.sliceToVtx.at(slc1).Y() << ", " << slcMaps.sliceToVtx.at(slc1).Z()
                      << ")\n";

            std::cout << tab << tab << "Assoc Hits Fraction "
                      << slcMaps.sliceToMCTrackIdAndHits.at(slc1).at(id) << " / "
                      << slcMaps.sliceToTotalHits.at(slc1) << " = " << assocHitsFrac1 << '\n';

            std::cout << tab << tab << "This id seen in planes { ";

            for (auto const& view : seenViews1) {
              std::cout << view << ' ';
            }
            std::cout << "}\n";

            std::cout << tab << "Slice " << slc2 << "\n";

            std::cout << tab << tab << "Vtx Position (" << slcMaps.sliceToVtx.at(slc2).X() << ", "
                      << slcMaps.sliceToVtx.at(slc2).Y() << ", " << slcMaps.sliceToVtx.at(slc2).Z()
                      << ")\n";

            std::cout << tab << tab << "Assoc Hits Fraction "
                      << slcMaps.sliceToMCTrackIdAndHits.at(slc2).at(id) << " / "
                      << slcMaps.sliceToTotalHits.at(slc2) << " = " << assocHitsFrac2 << '\n';

            std::cout << tab << tab << "This id seen in planes { ";

            for (auto const& view : seenViews2) {
              std::cout << view << ' ';
            }
            std::cout << "}\n";
          } // other TrackIds
        }   // common TrackIds
      }     // slice-TrackId pair2
    }       // slice-TrackId pair1
    return count;
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

} // namespace SliceHitPurity

DEFINE_ART_MODULE(SliceHitPurity::SliceHitPurityAna)
