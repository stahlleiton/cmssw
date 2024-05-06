#include "CondFormats/SiStripObjects/interface/DeDxCalibration.h"
DeDxCalibration::DeDxCalibration() {
  thr.reserve(5);
  alpha.reserve(5);
  sigma.reserve(5);
}
