#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "TrackingTools/PatternTools/interface/ClosestApproachInRPhi.h"
#include "TrackingTools/PatternTools/interface/TSCBLBuilderNoMaterial.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "MagneticField/ParametrizedEngine/src/OAEParametrizedMagneticField.h"

#include "RecoVertex/KinematicFitPrimitives/interface/ParticleMass.h"
#include "RecoVertex/KinematicFitPrimitives/interface/MultiTrackKinematicConstraint.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "RecoVertex/KinematicFit/interface/KinematicConstrainedVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/TwoTrackMassKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleFitter.h"
#include "RecoVertex/KinematicFit/interface/MassKinematicConstraint.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicVertex.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParametersError.h"
#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"
#include "RecoVertex/VertexTools/interface/VertexDistanceXY.h"

class HNLKinematicVertexFitter {

  public:
    HNLKinematicVertexFitter() {};
    virtual ~HNLKinematicVertexFitter() {};

    reco::TransientTrack getTransientTrack(const reco::TrackRef& trackRef) {    
      reco::TransientTrack transientTrack(trackRef, paramField);
      return transientTrack;
    }

    RefCountedKinematicTree Fit(const pat::MuonCollection & muons){

      KinematicParticleFactoryFromTransientTrack pFactory;  
      std::vector<RefCountedKinematicParticle> XParticles_Tau3Mu;

      // do tau vertex fit
      std::vector<reco::TransientTrack> tau_tks;
      for (pat::MuonCollection::const_iterator imu = muons.begin(); imu != muons.end(); ++imu){
        XParticles_Tau3Mu.push_back(pFactory.particle(getTransientTrack(imu->track()), muon_mass, chi, ndf, muon_sigma));
      }

      // maybe there are other tracks to consider
      // http://cmslxr.fnal.gov/source/DataFormats/PatCandidates/interface/Muon.h#0070


      KinematicConstrainedVertexFitter kvFitterTau;
      RefCountedKinematicTree Tau3MuKinVtx = kvFitterTau.fit(XParticles_Tau3Mu); 
      
      return Tau3MuKinVtx;
        
    }

  private:
    OAEParametrizedMagneticField *paramField = new OAEParametrizedMagneticField("3_8T");
    // The mass of a muon and the insignificant mass sigma to avoid singularities in the covariance matrix.
    float muon_mass  = 0.10565837; 
    float muon_sigma = muon_mass*1.e-6;
    // initial chi2 and ndf before kinematic fits. The chi2 of the reconstruction is not considered 
    float chi        = 0.;
    float ndf        = 0.;

};

