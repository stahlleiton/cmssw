/* \class ConversionSelector
 *
 * Selects conversions with a configurable string-based cut.
 * Saves clones of the selected conversions
 * Warning: this module can read anything that inherits from reco::Conversion, but it will
 *   only clone the reco::Conversion part of the object, the rest is lost.
 *
 * \author: Andre Stahl, LLR
 *
 * usage:
 *
 * module bestConversions = ConversionSelector {
 *   src = allConversions
 *   string cut = "(tracks().size() == 2) && conversionVertex().isValid() && refittedPair4Momentum().pt() < 10."
 * }
 *
 * for more details about the cut syntax, see the documentation
 * page below:
 *
 *   https://twiki.cern.ch/twiki/bin/view/CMS/SWGuidePhysicsCutParser
 *
 *
 */

#include "FWCore/Framework/interface/MakerMacros.h"
#include "CommonTools/UtilAlgos/interface/SingleObjectSelector.h"
#include "CommonTools/UtilAlgos/interface/StringCutObjectSelector.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "DataFormats/EgammaCandidates/interface/ConversionFwd.h"

 typedef SingleObjectSelector<
           edm::View<reco::Conversion>, 
           StringCutObjectSelector<reco::Conversion>,
           reco::ConversionCollection
         > ConversionSelector;

DEFINE_FWK_MODULE( ConversionSelector );
