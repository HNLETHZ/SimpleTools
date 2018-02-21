#define G__DICTIONARY

#include "DataFormats/Common/interface/Wrapper.h"
#include "SimpleTools/HNLNtuple/interface/HNLKinematicVertexFitter.h"

#include "TLorentzVector.h"
#include "DataFormats/PatCandidates/interface/Muon.h"

namespace {
  struct HNL_SimpleTools {
    HNLKinematicVertexFitter hnlKinVtx_;

    std::pair<edm::Ptr<pat::Muon>,TLorentzVector> pppmtlv;
    edm::Wrapper<std::pair<edm::Ptr<pat::Muon>,TLorentzVector> > wpppmtlv;

    std::vector<std::pair<edm::Ptr<pat::Muon>,TLorentzVector> > vpppmtlv;
    edm::Wrapper<std::vector<std::pair<edm::Ptr<pat::Muon>,TLorentzVector> > > wvpppmtlv;

  };
}


