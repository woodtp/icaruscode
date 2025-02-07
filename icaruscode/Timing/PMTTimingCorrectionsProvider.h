/**
 * @file   icaruscode/Timing/PMTTimingCorrectionsProvider.h
 * @brief  Service for the PMT timing corrections.
 * @author Andrea Scarpelli (ascarpell@bnl.gov), Matteo Vicenzi (mvicenzi@bnl.gov)
 */

#ifndef ICARUSCODE_TIMING_PMTIMINGCORRECTIONSPROVIDER_H
#define ICARUSCODE_TIMING_PMTIMINGCORRECTIONSPROVIDER_H

// Framework includes
#include "art/Framework/Services/Registry/ServiceDefinitionMacros.h"
#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Registry/ServiceDeclarationMacros.h"
#include "art/Framework/Principal/Run.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib_except/exception.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"
#include "cetlib_except/exception.h"

// Local
#include "icaruscode/Timing/PMTTimingCorrections.h"

// Database interface helpers
#include "larevt/CalibrationDBI/Providers/DBFolder.h"

// C/C++ standard libraries
#include <string>

namespace icarusDB::details {
    
  /// Structure for single channel corrections
  struct PMTTimeCorrectionsDB { 

    double triggerCableDelay=0;  ///< [&micro;s]
    double resetCableDelay=0;    ///< [&micro;s]
    double laserCableDelay=0;    ///< [&micro;s]
    double cosmicsCorrections=0; ///< [&micro;s]
    
  };
  
} // icarusDB::details

namespace icarusDB{ class PMTTimingCorrectionsProvider; }
/**
 * @brief 
 * 
 * This module reads the PMT timing corrections from the database.
 * Corrections are picked according to the run number being processed.  
 *
 * All time corrections are offsets (in microseconds) that need to be _added_ to the uncorrected time.
 * 
 * Configuration parameters
 * -------------------------
 * * `Tag` (default: `false`): Tag for database versioning
 * * `Verbose` (default: `false`): Print-out the corrections read from the database.
 * * `LogCategory` (default: `PMTTimingCorrection")
 *
 */
class icarusDB::PMTTimingCorrectionsProvider : public PMTTimingCorrections {

    public: 

        PMTTimingCorrectionsProvider(const fhicl::ParameterSet& pset);

	/// Read timing corrections from the database
        void readTimeCorrectionDatabase(const art::Run& run);

	/// Get time delay on the trigger line
        double getTriggerCableDelay( unsigned int channelID ) const override {
            return getChannelCorrOrDefault(channelID).triggerCableDelay;
        };

	/// Get time delay on the PPS reset line 
        double getResetCableDelay( unsigned int channelID ) const override {
            return getChannelCorrOrDefault(channelID).resetCableDelay;
        };

	/// Get timing corrections from laser data
        double getLaserCorrections( unsigned int channelID ) const override {
            return getChannelCorrOrDefault(channelID).laserCableDelay;
        };

	/// Get timing corrections from cosmics data
        double getCosmicsCorrections( unsigned int channelID ) const override {
            return getChannelCorrOrDefault(channelID).cosmicsCorrections;
        };

    private:
        
        using PMTTimeCorrectionsDB = details::PMTTimeCorrectionsDB;
        static constexpr PMTTimeCorrectionsDB CorrectionDefaults {}; ///< Default values

        bool fVerbose = false; ///< Whether to print the configuration we read.
        std::string fLogCategory; ///< Category tag for messages.
	fhicl::ParameterSet fTags; ///< List of database tags
	std::string fCablesTag;  ///< Tag for cable corrections database.	
	std::string fLaserTag;   ///< Tag for laser corrections database.
	std::string fCosmicsTag; ///< Tag for cosmics corrections database.	

	/// Map of corrections by channel
        std::map<unsigned int, PMTTimeCorrectionsDB> fDatabaseTimingCorrections;
        
        /// Internal access to the channel correction record; returns defaults if not present.
        PMTTimeCorrectionsDB const& getChannelCorrOrDefault
            (unsigned int channelID) const
            {
                auto const it = fDatabaseTimingCorrections.find(channelID);
                return (it == fDatabaseTimingCorrections.end())? CorrectionDefaults: it->second;
            }

	/// Convert run number to internal database
	uint64_t RunToDatabaseTimestamp(uint32_t run);

        void ReadPMTCablesCorrections(uint32_t run);

        void ReadLaserCorrections(uint32_t run);

        void ReadCosmicsCorrections(uint32_t run);

}; // services class

#endif 
