#ifndef DEDXCALIBRATION_H
#define DEDXCALIBRATION_H
#include "CondFormats/Serialization/interface/Serializable.h"

#include <vector>
class DeDxCalibration {
public:
  DeDxCalibration();
  virtual ~DeDxCalibration() {}
  std::vector<double> thr;
  std::vector<double> alpha;
  std::vector<double> sigma;
  typedef std::pair<uint32_t, unsigned char> ChipId;
  std::map<ChipId, float> gain;

  COND_SERIALIZABLE;
};
#endif
