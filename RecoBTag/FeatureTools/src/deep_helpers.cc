#include "RecoBTag/FeatureTools/interface/deep_helpers.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "FWCore/Utilities/interface/isFinite.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"
#include "TrackingTools/GeomPropagators/interface/AnalyticalTrajectoryExtrapolatorToLine.h"
#include "TrackingTools/GeomPropagators/interface/AnalyticalImpactPointExtrapolator.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/Candidate/interface/VertexCompositePtrCandidate.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "TLorentzVector.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"
#include "RecoVertex/VertexTools/interface/VertexDistanceXY.h"
#include "RecoVertex/VertexPrimitives/interface/ConvertToFromReco.h"
#include "RecoVertex/VertexPrimitives/interface/VertexState.h"
namespace btagbtvdeep {

  constexpr static int qualityMap[8] = {1, 0, 1, 1, 4, 4, 5, 6};

  enum qualityFlagsShiftsAndMasks {
    assignmentQualityMask = 0x7,
    assignmentQualityShift = 0,
    trackHighPurityMask = 0x8,
    trackHighPurityShift = 3,
    lostInnerHitsMask = 0x30,
    lostInnerHitsShift = 4,
    muonFlagsMask = 0x0600,
    muonFlagsShift = 9
  };

  // remove infs and NaNs with value  (adapted from DeepNTuples)
  const float catch_infs(const float in, const float replace_value) {
    if (in == in) {  // check if NaN
      if (std::isinf(in))
        return replace_value;
      else if (in < -1e32 || in > 1e32)
        return replace_value;
      return in;
    }
    return replace_value;
  }

  // remove infs/NaN and bound (adapted from DeepNTuples)
  const float catch_infs_and_bound(const float in,
                                   const float replace_value,
                                   const float lowerbound,
                                   const float upperbound,
                                   const float offset,
                                   const bool use_offsets) {
    float withoutinfs = catch_infs(in, replace_value);
    if (withoutinfs + offset < lowerbound)
      return lowerbound;
    if (withoutinfs + offset > upperbound)
      return upperbound;
    if (use_offsets)
      withoutinfs += offset;
    return withoutinfs;
  }

  // 2D distance between SV and PV (adapted from DeepNTuples)
  Measurement1D vertexDxy(const reco::VertexCompositePtrCandidate &svcand, const reco::Vertex &pv) {
    VertexDistanceXY dist;
    reco::Vertex::CovarianceMatrix csv;
    svcand.fillVertexCovariance(csv);
    reco::Vertex svtx(svcand.vertex(), csv);
    return dist.distance(svtx, pv);
  }

  //3D distance between SV and PV (adapted from DeepNTuples)
  Measurement1D vertexD3d(const reco::VertexCompositePtrCandidate &svcand, const reco::Vertex &pv) {
    VertexDistance3D dist;
    reco::Vertex::CovarianceMatrix csv;
    svcand.fillVertexCovariance(csv);
    reco::Vertex svtx(svcand.vertex(), csv);
    return dist.distance(svtx, pv);
  }

  // dot product between SV and PV (adapted from DeepNTuples)
  float vertexDdotP(const reco::VertexCompositePtrCandidate &sv, const reco::Vertex &pv) {
    reco::Candidate::Vector p = sv.momentum();
    reco::Candidate::Vector d(sv.vx() - pv.x(), sv.vy() - pv.y(), sv.vz() - pv.z());
    return p.Unit().Dot(d.Unit());
  }

  // compute minimum dr between SVs and a candidate (from DeepNTuples, now polymorphic)
  float mindrsvpfcand(const std::vector<reco::VertexCompositePtrCandidate> &svs,
                      const reco::Candidate *cand,
                      float mindr) {
    for (unsigned int i0 = 0; i0 < svs.size(); ++i0) {
      float tempdr = reco::deltaR(svs[i0], *cand);
      if (tempdr < mindr) {
        mindr = tempdr;
      }
    }
    return mindr;
  }

  // compute minimum distance between SVs and a candidate (from DeepNTuples, now polymorphic)
  float mindistsvpfcand(const std::vector<reco::VertexCompositePtrCandidate> &svs, const reco::TransientTrack track) {
    float mindist_ = 999.999;
    float out_dist = 0.0;
    for (unsigned int i = 0; i < svs.size(); ++i) {
      if (!track.isValid()) {
        continue;
      }
      reco::Vertex::CovarianceMatrix csv;
      svs[i].fillVertexCovariance(csv);
      reco::Vertex vertex(svs[i].vertex(), csv);
      if (!vertex.isValid()) {
        continue;
      }

      GlobalVector direction(svs[i].px(), svs[i].py(), svs[i].pz());

      AnalyticalImpactPointExtrapolator extrapolator(track.field());
      TrajectoryStateOnSurface tsos =
          extrapolator.extrapolate(track.impactPointState(), RecoVertex::convertPos(vertex.position()));

      VertexDistance3D dist;

      if (!tsos.isValid()) {
        continue;
      }
      GlobalPoint refPoint = tsos.globalPosition();
      GlobalError refPointErr = tsos.cartesianError().position();
      GlobalPoint vertexPosition = RecoVertex::convertPos(vertex.position());
      GlobalError vertexPositionErr = RecoVertex::convertError(vertex.error());

      std::pair<bool, Measurement1D> result(
          true, dist.distance(VertexState(vertexPosition, vertexPositionErr), VertexState(refPoint, refPointErr)));
      if (!result.first) {
        continue;
      }

      GlobalPoint impactPoint = tsos.globalPosition();
      GlobalVector IPVec(impactPoint.x() - vertex.x(), impactPoint.y() - vertex.y(), impactPoint.z() - vertex.z());
      double prod = IPVec.dot(direction);
      double sign = (prod >= 0) ? 1. : -1.;

      if (result.second.value() < mindist_) {
        out_dist = sign * result.second.value();
        mindist_ = result.second.value();
      }
    }
    return out_dist;
  }

  // instantiate template
  template bool sv_vertex_comparator<reco::VertexCompositePtrCandidate, reco::Vertex>(
      const reco::VertexCompositePtrCandidate &, const reco::VertexCompositePtrCandidate &, const reco::Vertex &);

  float vtx_ass_from_pfcand(const reco::PFCandidate &pfcand, int pv_ass_quality, const reco::VertexRef &pv) {
    float vtx_ass = pat::PackedCandidate::PVAssociationQuality(qualityMap[pv_ass_quality]);
    if (pfcand.trackRef().isNonnull() && pv->trackWeight(pfcand.trackRef()) > 0.5 && pv_ass_quality == 7) {
      vtx_ass = pat::PackedCandidate::UsedInFitTight;
    }
    return vtx_ass;
  }

  float quality_from_pfcand(const reco::PFCandidate &pfcand) {
    const auto &pseudo_track = (pfcand.bestTrack()) ? *pfcand.bestTrack() : reco::Track();
    // conditions from PackedCandidate producer
    bool highPurity = pfcand.trackRef().isNonnull() && pseudo_track.quality(reco::Track::highPurity);
    // do same bit operations than in PackedCandidate
    uint16_t qualityFlags = 0;
    qualityFlags = (qualityFlags & ~trackHighPurityMask) | ((highPurity << trackHighPurityShift) & trackHighPurityMask);
    bool isHighPurity = (qualityFlags & trackHighPurityMask) >> trackHighPurityShift;
    // to do as in TrackBase
    uint8_t quality = (1 << reco::TrackBase::loose);
    if (isHighPurity) {
      quality |= (1 << reco::TrackBase::highPurity);
    }
    return quality;
  }

  float lost_inner_hits_from_pfcand(const reco::PFCandidate &pfcand) {
    const auto &pseudo_track = (pfcand.bestTrack()) ? *pfcand.bestTrack() : reco::Track();
    // conditions from PackedCandidate producer
    bool highPurity = pfcand.trackRef().isNonnull() && pseudo_track.quality(reco::Track::highPurity);
    // do same bit operations than in PackedCandidate
    uint16_t qualityFlags = 0;
    qualityFlags = (qualityFlags & ~trackHighPurityMask) | ((highPurity << trackHighPurityShift) & trackHighPurityMask);
    return int16_t((qualityFlags & lostInnerHitsMask) >> lostInnerHitsShift) - 1;
  }

  std::pair<float, float> getDRSubjetFeatures(const reco::Jet &jet, const reco::Candidate *cand) {
    const auto *patJet = dynamic_cast<const pat::Jet *>(&jet);
    std::pair<float, float> features;
    // Do Subjets
    if (patJet) {
      if (patJet->nSubjetCollections() > 0) {
        auto subjets = patJet->subjets();
        std::nth_element(
            subjets.begin(),
            subjets.begin() + 1,
            subjets.end(),
            [](const edm::Ptr<pat::Jet> &p1, const edm::Ptr<pat::Jet> &p2) { return p1->pt() > p2->pt(); });
        features.first = !subjets.empty() ? reco::deltaR(*cand, *subjets[0]) : -1;
        features.second = subjets.size() > 1 ? reco::deltaR(*cand, *subjets[1]) : -1;
      } else {
        features.first = -1;
        features.second = -1;
      }
    } else {
      features.first = -1;
      features.second = -1;
    }
    return features;
  }
}  // namespace btagbtvdeep
