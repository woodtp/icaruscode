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
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include <iostream>

namespace SliceHitPurity
{
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

};


SliceHitPurityAna::SliceHitPurityAna(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}  // ,
  // More initializers here.
{
  // Call appropriate consumes<>() for any products to be retrieved by this module.
}

void SliceHitPurityAna::analyze(art::Event const& e)
{
  // Implementation of required member function here.
  std::cout << "Done." << std::endl;
}

} // namespace SliceHitPurity

DEFINE_ART_MODULE(SliceHitPurity::SliceHitPurityAna)

