// -*- C++ -*-
//
// Package:    SoftElectronInJetAnalyzer/ElectronID
// Class:      ElectronID
// 
/**\class ElectronID ElectronID.cc SoftElectronInJetAnalyzer/ElectronID/plugins/ElectronID.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Simon De Visscher
//         Created:  Thu, 08 Jan 2015 15:34:42 GMT
//
//

#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/ParticleFlowReco/interface/PFDisplacedTrackerVertex.h"
#include "DataFormats/ParticleFlowReco/interface/PFDisplacedVertexFwd.h"
#include "DataFormats/ParticleFlowReco/interface/PFDisplacedVertex.h"
#include "TrackingTools/PatternTools/interface/Trajectory.h"
#include "DataFormats/TrajectorySeed/interface/TrajectorySeed.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "DataFormats/ParticleFlowReco/interface/PFClusterFwd.h"
#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "RecoParticleFlow/PFTracking/interface/ElectronSeedMerger.h"
#include "DataFormats/TrajectorySeed/interface/PropagationDirection.h"
#include "TrackingTools/Records/interface/TransientRecHitRecord.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
// //#include "RecoParticleFlow/PFTracking/interface/GoodSeedProducer.h"
#include "DataFormats/EgammaReco/interface/ElectronSeed.h"
#include "DataFormats/EgammaReco/interface/ElectronSeedFwd.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/ParticleFlowReco/interface/PreId.h"
#include "DataFormats/ParticleFlowReco/interface/PreIdFwd.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"
#include "TrackingTools/PatternTools/interface/TSCBLBuilderNoMaterial.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DataFormats/Math/interface/deltaR.h"
//
//
// //add for the input variables
#include "RecoParticleFlow/PFTracking/interface/PFTrackTransformer.h"
#include "RecoParticleFlow/PFClusterTools/interface/PFResolutionMap.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/FileInPath.h"
#include "DataFormats/TrajectorySeed/interface/TrajectorySeedCollection.h"
#include "TrackingTools/TrackFitters/interface/TrajectoryFitter.h"
#include "TrackingTools/PatternTools/interface/TrajectorySmoother.h"
#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"
#include "FastSimulation/BaseParticlePropagator/interface/BaseParticlePropagator.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "RecoParticleFlow/PFClusterTools/interface/LinkByRecHit.h"
#include "RecoParticleFlow/PFClusterTools/interface/PFEnergyCalibration.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecTrack.h"
#include "TrackingTools/GsfTools/interface/MultiTrajectoryStateMode.h"
#include "TrackingTools/GsfTools/interface/MultiTrajectoryStateTransform.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticleFwd.h"
#include "DataFormats/HeavyIonEvent/interface/EvtPlane.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "SimDataFormats/Track/interface/SimTrack.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "SimDataFormats/Track/interface/CoreSimTrack.h"
#include "SimDataFormats/EncodedEventId/interface/EncodedEventId.h"
#include "SimDataFormats/Vertex/interface/SimVertex.h"
#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"
#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"
#include "SimDataFormats/TrackingHit/interface/UpdatablePSimHit.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticleFwd.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertexContainer.h"
#include "SimTracker/TrackAssociation/interface/TrackAssociatorBase.h"
#include "SimTracker/TrackAssociation/interface/TrackAssociatorByChi2.h"
#include "SimTracker/TrackAssociation/interface/QuickTrackAssociatorByHits.h"
#include "SimTracker/TrackerHitAssociation/interface/TrackerHitAssociator.h"
#include "SimTracker/Records/interface/TrackAssociatorRecord.h"
#include "SimTracker/TrackAssociation/plugins/ParametersDefinerForTPESProducer.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
//
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
//
#include "TrackingTools/PatternTools/interface/TwoTrackMinimumDistance.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
//
#include "RecoVertex/ConfigurableVertexReco/interface/ConfigurableVertexReconstructor.h"
//
#include "DataFormats/EgammaReco/interface/ElectronSeedFwd.h"
#include "DataFormats/EgammaReco/interface/ElectronSeed.h"
//
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrackFwd.h"
//
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/TauReco/interface/PFTau.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Common/interface/Association.h"
//
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
//
#include "RecoParticleFlow/PFProducer/interface/Utils.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
//
  #include "TrackingTools/TransientTrack/interface/TransientTrack.h"
  #include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
  #include "TrackingTools/Records/interface/TransientTrackRecord.h"
  #include "TrackingTools/IPTools/interface/IPTools.h"
  #include "DataFormats/GeometryVector/interface/GlobalVector.h"
//
  #include "RecoVertex/AdaptiveVertexFit/interface/AdaptiveVertexFitter.h"
  #include "RecoVertex/KalmanVertexFit/interface/KalmanVertexUpdator.h"
  #include "RecoVertex/KalmanVertexFit/interface/KalmanVertexTrackCompatibilityEstimator.h"
  #include "RecoVertex/KalmanVertexFit/interface/KalmanVertexSmoother.h"
  #include "RecoVertex/MultiVertexFit/interface/MultiVertexFitter.h"
//
  #include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
//
  #include "DataFormats/BTauReco/interface/SoftLeptonTagInfo.h"
//
  #include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
//
  #include "DataFormats/BTauReco/interface/JetTag.h"
  #include "DataFormats/JetReco/interface/Jet.h"
  #include "SimDataFormats/JetMatching/interface/JetFlavour.h"
  #include "SimDataFormats/JetMatching/interface/JetFlavourMatching.h"
//
//
//  // Muons
  #include "DataFormats/MuonReco/interface/Muon.h"
  #include "DataFormats/MuonReco/interface/MuonFwd.h"
  #include "DataFormats/MuonReco/interface/MuonSelectors.h"
//
  #include <TLorentzVector.h>
//
  #include "RecoToGenClass.h"
  #include <TFile.h>
  #include <TTree.h>
  #include <vector>
  #include <string>
  #include <sstream>
  #include <fstream>
  #include "TMath.h"
  #include "TVector2.h"
  #include "Math/VectorUtil.h"
  #include "TKey.h"
//

// class declaration
//

using namespace edm;
using namespace std;
using namespace reco;



class ElectronID : public edm::EDAnalyzer {
   public:
      explicit ElectronID(const edm::ParameterSet&);
      ~ElectronID();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      bool trackFilter(const reco::TrackRef &track) const;
      virtual bool isInnerMost(const reco::Track track1, const reco::Track track2, bool& sameLayer);
      bool findMatch(const reco::Candidate& genCand, Handle<TrackingParticleCollection> TPCollectionH,reco::SimToRecoCollection q,reco::Track& track,int& indexTrack,float& dR,float& sharedHits,RefToBase<Track> &trRefToBase);
      void GEDGsfElecFiller();
      bool  isaBhadron( int );
      bool  isaDhadron( int );
      bool  isaV( int );
      int GetOrigin(reco::GenParticle&);
      bool UseRECO, Verbose;
      InputTag PVerTag,gedgsfElectronTag_,genParticleTag_,TrackTag_,TrackingParticleTag_,allConversionsTag,offlineBeamSpotTag;
      float maxLIP;
      int minHits;
      RecoToGenFiller *rtgf;
      Handle<reco::VertexCollection> PrimaryVetices;
      ESHandle<TransientTrackBuilder> builder;
      Handle<edm::View<reco::Track> > theTrackColl;
      Handle<TrackingParticleCollection> theTrackingPartColl;
      ESHandle<TrackAssociatorBase> theHitsAssociator;
      Handle<GenParticleCollection> GPC;
      Handle<reco::GsfElectronCollection> gedgsfCandidates;
      Handle<reco::ConversionCollection> hConversions;
      edm::Handle<BeamSpot> beamSpot;
      const reco::Vertex* vertex;
      int nPV;
      const TransientTrackBuilder* transientTrackBuilder;
      reco::RecoToSimCollection theRecoToSimColl;
      GenParticleCollection gpc ;
      reco::VertexCollection pvc;
      bool goodvertex;
      reco::GsfPFRecTrackCollection  gsfpfTk;
      std::vector<TransientTrack> tts;

      // ----------member data ---------------------------
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
ElectronID::ElectronID(const edm::ParameterSet& iConfig):
        UseRECO(iConfig.getParameter<bool>         ("UseRECO") ),
        Verbose(iConfig.getParameter<bool>         ("Verbose") ),
        PVerTag(iConfig.getParameter<edm::InputTag>("PVerTag") ),
        gedgsfElectronTag_(iConfig.getParameter<edm::InputTag>("gedgsfElectronTag")),
        genParticleTag_(iConfig.getParameter<edm::InputTag>("genParticleTag")),
        TrackTag_(iConfig.getParameter<edm::InputTag>("tracksTag")),
        TrackingParticleTag_(iConfig.getParameter<edm::InputTag>("TrackingParticleTag")),
	allConversionsTag(iConfig.getParameter<edm::InputTag>("conversionTag")),
	offlineBeamSpotTag(iConfig.getParameter<edm::InputTag>("beamSpotTag")),
	maxLIP(iConfig.getParameter<double>("maximumLongitudinalImpactParameter")),
        minHits(iConfig.getParameter<unsigned int>("minHits"))
{
        std::cout<<"in constr"<<std::endl;
        rtgf= new RecoToGenFiller("GEDGSF");
        std::cout<<"in constr"<<std::endl;

}




ElectronID::~ElectronID()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
ElectronID::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
  std::cout<<"getbylabel 1"<<std::endl;
        iEvent.getByLabel(PVerTag, PrimaryVetices);
        iEvent.getByLabel(genParticleTag_, GPC);
        iEvent.getByLabel(gedgsfElectronTag_, gedgsfCandidates);
	iEvent.getByLabel(allConversionsTag, hConversions);
	iEvent.getByLabel(offlineBeamSpotTag, beamSpot);

        std::cout<<"test0"<<std::endl;
        if(!PrimaryVetices.isValid() || PrimaryVetices->empty()) return;
        vertex=&PrimaryVetices->front();
        nPV=PrimaryVetices->size();
        std::cout<<"test1"<<std::endl;
        iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", builder);
        transientTrackBuilder=builder.product();
        std::cout<<"test15"<<std::endl;
        const TrackAssociatorBase* theAssociatorByHits(NULL);
        std::cout<<"test16"<<std::endl;
        if(UseRECO){
                std::cout<<"UR0"<<std::endl;
                iEvent.getByLabel(TrackTag_, theTrackColl);
                iEvent.getByLabel(TrackingParticleTag_,theTrackingPartColl);
                iSetup.get<TrackAssociatorRecord>().get("quickTrackAssociatorByHits", theHitsAssociator);
                theAssociatorByHits = (TrackAssociatorBase*)theHitsAssociator.product();
                std::cout<<"UR05"<<std::endl;
                theRecoToSimColl=theAssociatorByHits->associateRecoToSim(theTrackColl, theTrackingPartColl, &iEvent, &iSetup);
                std::cout<<"UR1"<<std::endl;
        }
        gpc = *(GPC.product());
        pvc = *(PrimaryVetices);
        std::cout<<"test17"<<std::endl;
        goodvertex=!pvc.empty()?true:false;
        if(goodvertex)vertex=&pvc.front();
        const View<Track> tC = *(theTrackColl.product());
        std::cout<<"test2"<<std::endl;
        GEDGsfElecFiller();
}


void ElectronID::GEDGsfElecFiller(){
        int u=0;
        for (reco::GsfElectronCollection::const_iterator  gsf=gedgsfCandidates->begin(); gsf!=gedgsfCandidates->end(); ++gsf)
        {
                rtgf->initRecoToGenFillerObject();
                RefToBase<Track> tkRef  = RefToBase<Track>(gsf->gsfTrack() );
                if(theRecoToSimColl.find(tkRef) != theRecoToSimColl.end()){
                        TrackingParticleRef tpr = theRecoToSimColl[tkRef].begin()->first;
                        rtgf->tppdgId=tpr->pdgId();
                        TrackingParticle::genp_iterator j, b = tpr->genParticle_begin(), e = tpr->genParticle_end();
                        for( j = b; j != e; ++ j ) {
                                const reco::GenParticle * p = j->get();
                                reco::GenParticle *ncp=(reco::GenParticle*)p;
                                rtgf->Vtx=1;
                                rtgf->pdgId=ncp->pdgId();
                                rtgf->origin=GetOrigin(*ncp);
                                rtgf->ptGen=ncp->pt();
                                rtgf->pGen=ncp->p();
                                rtgf->etaGen=ncp->eta();
                        }
                }
                GsfElectronRef elref( gedgsfCandidates, u);
//                cout<<"electron pt="<<gsf->pt()<<" MVASoftOutput:"<< (mvasoft)[elref]<<endl;
//                rtgf->mva_e_pi=(mvasoft)[elref];
                rtgf->RecoPt=gsf->pt();
                rtgf->RecoEta=gsf->eta();
                rtgf->RecoPhi=gsf->phi();
                rtgf->EtotOvePin=gsf->eSuperClusterOverP();
                rtgf->EClusOverPout=gsf->eEleClusterOverPout();
                rtgf->fbrem=gsf->fbrem()>-1 ? gsf->fbrem():-1;
                float etot=gsf->eSuperClusterOverP()*gsf->trackMomentumAtVtx().R();
                float eEcal=gsf->eEleClusterOverPout()*gsf->trackMomentumAtEleClus().R();
                float dP=gsf->trackMomentumAtVtx().R()-gsf->trackMomentumAtEleClus().R();
                rtgf->RecoE=etot;
                rtgf->TotalEnergy=gsf->ecalEnergy();
                rtgf->EBremOverDeltaP=(etot-eEcal)/dP ;
                rtgf->logSigmaEtaEta=log(gsf->sigmaEtaEta());
                rtgf->DeltaEtaTrackEcalSeed=gsf->deltaEtaEleClusterTrackAtCalo();
                rtgf->HOverE=gsf->hadronicOverEm();
                rtgf->Chi2GSF=gsf->gsfTrack()->normalizedChi2()<200?gsf->gsfTrack()->normalizedChi2():200;
                bool validKF= false;
                reco::TrackRef myTrackRef = gsf->closestCtfTrackRef();
                validKF = (myTrackRef.isAvailable());
                validKF = (myTrackRef.isNonnull());
                float tempchi2kf=(validKF) ? myTrackRef->normalizedChi2() : 0 ;
                rtgf->CHi2KF=tempchi2kf<10?tempchi2kf:10 ;
                rtgf->nHits=(validKF) ? myTrackRef->hitPattern().trackerLayersWithMeasurement() : -1. ;
                rtgf->SigmaPtOverPt=gsf->gsfTrack().get()->ptModeError()/gsf->gsfTrack().get()->ptMode() ;

                rtgf->deta=gsf->deltaEtaSuperClusterTrackAtVtx()<0.06?gsf->deltaEtaSuperClusterTrackAtVtx():0.06;
                rtgf->dphi=gsf->deltaPhiSuperClusterTrackAtVtx()<0.6 ?gsf->deltaPhiSuperClusterTrackAtVtx():0.6;
                rtgf->detacalo=gsf->deltaEtaSeedClusterTrackAtCalo();
                rtgf->see=gsf->sigmaIetaIeta();
                rtgf->spp=gsf->sigmaIphiIphi();
                rtgf->etawidth=gsf->superCluster()->etaWidth();
                rtgf->phiwidth =  gsf->superCluster()->phiWidth();
                float temp=(gsf->e5x5()) !=0. ? 1.-(gsf->e1x5()/gsf->e5x5()) : -1. ;
                rtgf->e1x5e5x5=  temp>-1?temp:-1;
                rtgf->e1x5e5x5=  temp<2?temp:2;
                rtgf->R9              =  gsf->r9();
                rtgf->HoE             =  gsf->hcalOverEcalBc();
                rtgf->EoP             =  gsf->eSuperClusterOverP()<20? gsf->eSuperClusterOverP():20;
                rtgf->IoEmIoP         =  (1.0/gsf->ecalEnergy()) - (1.0 / gsf->p());
                rtgf->eleEoPout       =  gsf->eEleClusterOverPout()<20?gsf->eEleClusterOverPout():20;
                rtgf->PreShowerOverRaw = gsf->superCluster()->preshowerEnergy() / gsf->superCluster()->rawEnergy();
                bool fromConv=ConversionTools::hasMatchedConversion((*gsf),hConversions,beamSpot->position());
                rtgf->fromConversion = fromConv ? 1:0;
/*                std::cout<<"The weight is "<<GetWeight((*gsf))<<endl;


                rtgf->weight          = GetWeight((*gsf));*/
                if (gsf->gsfTrack().isNonnull()) {
                        rtgf->d0 = (-1.0)*gsf->gsfTrack()->dxy(pvc[0].position());
                } else if (gsf->closestCtfTrackRef().isNonnull()) {
                        rtgf->d0 = (-1.0)*gsf->closestCtfTrackRef()->dxy(pvc[0].position());
                } else {
                        rtgf->d0 = -9999.0;
                }



                if (gsf->gsfTrack().isNonnull()) {
                        const double gsfsign   = ( (-gsf->gsfTrack()->dxy(pvc[0].position()))   >=0 ) ? 1. : -1.;

                        const reco::TransientTrack tt = transientTrackBuilder->build(gsf->gsfTrack());
                        const std::pair<bool,Measurement1D> &ip3dpv =  IPTools::absoluteImpactParameter3D(tt,pvc[0]);
                        if (ip3dpv.first) {
                                double ip3d = gsfsign*ip3dpv.second.value();
                                double ip3derr = ip3dpv.second.error();
                                rtgf->ip3d = ip3d;
                                rtgf->ip3dSig = ip3d/ip3derr;
                        }
                }

                rtgf->nPV=pvc.size();

                rtgf->tree_purity->Fill();
                if(fabs(rtgf->pdgId)==11 && (rtgf->origin==4 || rtgf->origin==6)){
                        rtgf->weight=1;
                        rtgf->tree_purity_GenLepB->Fill();
                        cout<<"An gen elec from B"<<endl;
                }
                else if (fabs(rtgf->pdgId)==211 && (rtgf->origin>0)){
                        rtgf->tree_purity_GenPi->Fill();
                        cout<<"An gen pion from B"<<endl;
                }
                else if (rtgf->Vtx==1){
                        rtgf->tree_purity_GenX->Fill();
                        cout<<"An gen X"<<endl;
                }
                else if (fabs(rtgf->tppdgId)==11){
                        rtgf->tree_purity_NonGenElec->Fill();
                        cout<<"An non gen elec"<<endl;
                }
                else if (fabs(rtgf->tppdgId)==211){
                        rtgf->tree_purity_NonGenPi->Fill();
                        cout<<"An non gen pion"<<endl;
                }
                else {
                        rtgf->tree_purity_NonGenX->Fill();
                        cout<<"An non gen X"<<endl;
                }

                if(!(myTrackRef.isAvailable() && myTrackRef.isNonnull()))
                        continue;
                if( std::abs(myTrackRef->dz(pvc[0].position())) > maxLIP)
                        continue;
                TransientTrack tt = transientTrackBuilder->build(gsf->closestCtfTrackRef());
                tt.setBeamSpot(*beamSpot);
                tts.push_back(tt);
                u++;
        }
}



bool ElectronID::isInnerMost(const reco::Track track1, const reco::Track track2, bool& sameLayer) {
  reco::HitPattern hitPattern1 = track1.hitPattern();
  reco::HitPattern hitPattern2 = track2.hitPattern();

  int hitCounter1 = 0;
  trackingRecHit_iterator hitsIt1;
  for(hitsIt1 = track1.recHitsBegin(); hitsIt1 != track1.recHitsEnd(); ++hitsIt1, ++hitCounter1){
        if (((**hitsIt1).isValid())) break;
  }
  int hitCounter2 = 0;
  trackingRecHit_iterator hitsIt2;
  for(hitsIt2 = track2.recHitsBegin(); hitsIt2 != track2.recHitsEnd(); ++hitsIt2, ++hitCounter2){
        if (((**hitsIt2).isValid())) break;
  }
  uint32_t hit1 = hitPattern1.getHitPattern(HitPattern::TRACK_HITS,hitCounter1);
  uint32_t hit2 = hitPattern2.getHitPattern(HitPattern::TRACK_HITS,hitCounter2);
  if ( hitPattern1.getSubStructure(hit1) != hitPattern2.getSubStructure(hit2) )
    return hitPattern2.getSubStructure(hit2) < hitPattern1.getSubStructure(hit1);
  else if ( hitPattern1.getLayer(hit1) != hitPattern2.getLayer(hit2) )
    return hitPattern2.getLayer(hit2) < hitPattern1.getLayer(hit1);
  else {
    sameLayer = true;
    return false;
  }
}

int ElectronID::GetOrigin(reco::GenParticle& part){
        int isBanAncestor=0;
        int isDanAncestor=0;
        int isVanAncestor=0;
        if(part.numberOfMothers()!=0 ){
                const reco::Candidate *m1=part.mother(0);
                if(isaBhadron(fabs(m1->pdgId())))isBanAncestor=1;
                if(isaDhadron(fabs(m1->pdgId())))isDanAncestor=1;
                if(isaV(fabs(m1->pdgId())))isVanAncestor=1;
                if(m1->numberOfMothers()!=0 ){
                       const reco::Candidate * m2=m1->mother(0);
                        if(isaBhadron(fabs(m2->pdgId())))isBanAncestor=1;
                        if(isaDhadron(fabs(m2->pdgId())))isDanAncestor=1;
                        if(isaV(fabs(m2->pdgId())))isVanAncestor=1;
                        if(m2->numberOfMothers()!=0 ){
                                const reco::Candidate *m3=m2->mother(0);
                                if(isaBhadron(fabs(m3->pdgId())))isBanAncestor=1;
                                if(isaDhadron(fabs(m3->pdgId())))isDanAncestor=1;
                                if(isaV(fabs(m3->pdgId())))isVanAncestor=1;
                                if(m3->numberOfMothers()!=0 ){
                                        const reco::Candidate *m4=m3->mother(0);
                                        if(isaBhadron(fabs(m4->pdgId())))isBanAncestor=1;
                                        if(isaDhadron(fabs(m4->pdgId())))isDanAncestor=1;
                                        if(isaV(fabs(m4->pdgId())))isVanAncestor=1;
                                        if(m4->numberOfMothers()!=0 ){
                                                const reco::Candidate *m5=m4->mother(0);
                                                if(isaBhadron(fabs(m5->pdgId())))isBanAncestor=1;
                                                if(isaDhadron(fabs(m5->pdgId())))isDanAncestor=1;
                                                if(isaV(fabs(m5->pdgId())))isVanAncestor=1;
                                                if(m5->numberOfMothers()!=0 ){
                                                        const reco::Candidate *m6=m5->mother(0);
                                                        if(isaBhadron(fabs(m6->pdgId())))isBanAncestor=1;
                                                        if(isaDhadron(fabs(m6->pdgId())))isDanAncestor=1;
                                                        if(isaV(fabs(m6->pdgId())))isVanAncestor=1;
                                                }
                                        }

                                }
                        }

                }
        }
        int result=4*isBanAncestor+2*isDanAncestor+isVanAncestor;
        return result;
}

bool ElectronID::isaV(int pidAbs){
        int res=false;
        if(pidAbs==24 || pidAbs==23 || pidAbs==22)res=true;
        return res;
}


bool ElectronID::isaBhadron(int pidAbs){

        bool isB = false;
        if(pidAbs>500){

                pidAbs/= 100;

                if(pidAbs<60 && pidAbs>50)pidAbs/= 10;
                int mod10 = pidAbs % 5;

                if(mod10 == 0) {
                        isB = true;
                }
        }
        else  {
        }
        return isB;
}

bool ElectronID::isaDhadron(int pidAbs){

        bool isD = false;
        if(pidAbs>400){
                pidAbs/= 100;
                if(pidAbs<50 && pidAbs>40)pidAbs/= 10;
                int mod10 = pidAbs % 4;

                if(mod10 == 0 && pidAbs<10) {
                        isD = true;
                }
        }
        else{
        }
        return isD;
}

bool ElectronID::trackFilter(const reco::TrackRef &track) const
{
        if (track->hitPattern().numberOfValidHits() < (int)minHits)
                return false;
        if (track->pt() < 0.8 )
                return false;

        return true;
}


bool ElectronID::findMatch(const reco::Candidate& genCand,
                              Handle<TrackingParticleCollection> theTrackingPartColl,
                              reco::SimToRecoCollection q,
                              reco::Track& track,
                              int& indexTrack,
                              float& dR,
                              float& sharedHits,
                              RefToBase<Track> &trRefToBase) {

  const TrackingParticleCollection tPC = *(theTrackingPartColl.product());
  if ( tPC.size() == 0 ) return false;
  cout<<"in FM: playing with "<<tPC.size()<<" trackingParticles"<<endl;
  unsigned int nTPs = tPC.size();
  for(unsigned int iTP = 0; iTP < nTPs; ++iTP) {
    const TrackingParticle& lTP(tPC[iTP]);
    float deta = lTP.eta() - genCand.eta();
    float dphi = acos(cos(lTP.phi() - genCand.phi()));
    dR = sqrt(deta*deta + dphi*dphi);
    if ( lTP.pdgId() != genCand.pdgId() || dR > 0.05 ) continue;
        cout<<"in FM3: found a TP close enough to the gen particle "<<lTP.pdgId()<<" "<<genCand.pdgId()<<" "<<dR<<endl;
    edm::Ref<TrackingParticleCollection> tp(theTrackingPartColl, iTP);
    try {
        std::vector<std::pair<RefToBase<Track>, double> > trackV = q[tp];
        cout<<"in FM3: "<<trackV.size()<<endl;
      if ( trackV.size() == 1 ) {
        std::vector<std::pair<RefToBase<Track>, double> >::const_iterator it = trackV.begin();
        RefToBase<Track> tr = it->first;
        indexTrack = tr.key();
        track = *tr;
        trRefToBase=tr;
        sharedHits = it->second;
        return true;
      }
      if ( trackV.size() > 1 ) {
        std::vector<std::pair<RefToBase<Track>, double> >::const_iterator it = trackV.begin();
        RefToBase<Track> tr = it->first;
        Track track1 = *tr;
        int indexTrack1 = tr.key();
        float nHitsRatio1 = it->second;
        ++it; tr = it->first;
        int indexTrack2 = tr.key();
        Track track2 = *tr;
        float nHitsRatio2 = it->second;
        bool sameLayer = false;
        bool result = isInnerMost(track1, track2, sameLayer);
        Track finalTrack;
        if ( !sameLayer ) {
          if ( result ) {
            finalTrack = track2;
            sharedHits = nHitsRatio2;
            indexTrack = indexTrack2;
          }
          else {
            finalTrack = track1;
            sharedHits = nHitsRatio1;
            indexTrack = indexTrack1;
          }
        } else {
          if ( nHitsRatio2 > nHitsRatio1 ) {
            finalTrack = track2;
            sharedHits = nHitsRatio2;
            indexTrack = indexTrack2;
          }
          else {
            finalTrack = track1;
            sharedHits = nHitsRatio1;
            indexTrack = indexTrack1;
          }
        }
        track = finalTrack;
        return true;
      }
    } catch (Exception event) {
      return false;
    }
  }
  return false;
}



// ------------ method called once each job just before starting event loop  ------------
void 
ElectronID::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
ElectronID::endJob() 
{
	rtgf->WriteInFileAndCloseIt();
}

// ------------ method called when starting to processes a run  ------------
/*
void 
ElectronID::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void 
ElectronID::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void 
ElectronID::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void 
ElectronID::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
ElectronID::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(ElectronID);
