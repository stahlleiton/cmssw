#ifndef HLTrigger_HLTfilters_TriggerExpressionPathReader_h
#define HLTrigger_HLTfilters_TriggerExpressionPathReader_h

#include <vector>
#include <string>

#include "DataFormats/Common/interface/PathStatus.h"
#include "FWCore/Utilities/interface/EDGetToken.h"

#include "HLTrigger/HLTcore/interface/TriggerExpressionEvaluator.h"

namespace triggerExpression {

  class PathReader : public Evaluator {
  public:
    PathReader(const std::string& pattern) : m_pattern(pattern), m_triggers() {}

    bool operator()(const Data& data) const override;

    void init(edm::ConsumesCollector&& iC) override;

    void init(const Data& data) override;

    void dump(std::ostream& out) const override;

  private:
    std::string m_pattern;
    std::vector<std::pair<std::string, unsigned int> > m_triggers;
    edm::EDGetTokenT<edm::PathStatus> m_token;
  };

}  // namespace triggerExpression

#endif  // HLTrigger_HLTfilters_TriggerExpressionPathReader_h
