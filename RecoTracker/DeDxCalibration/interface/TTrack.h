#ifndef _TTrack_h_
#define _TTrack_h_

#include <utility>
#include "TObject.h"

#include "TPixelHit.h"
#include "TStripHit.h"

class TTrack : public TObject
{
 public:
  TTrack();
  virtual ~TTrack();

  short int charge;

  float eta;
  float pt;
  float phi;

  float chi2;
  short int ndf;

  std::vector<TPixelHit> pixelHits;
  std::vector<TStripHit> stripHits;

  bool isHighPurity;

  ClassDef(TTrack,1)
};

#endif
