/* \class ConversionCountFilter
 *
 * Filters events if at least N conversions
 *
 * \author: Andre Stahl, LLR
 *
 */
#include "FWCore/Framework/interface/MakerMacros.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "DataFormats/EgammaCandidates/interface/ConversionFwd.h"
#include "CommonTools/UtilAlgos/interface/ObjectCountFilter.h"

 typedef ObjectCountFilter<
           reco::ConversionCollection
         >::type ConversionCountFilter;

DEFINE_FWK_MODULE( ConversionCountFilter );
