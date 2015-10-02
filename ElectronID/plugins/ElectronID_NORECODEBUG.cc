// -*- C++ -*-
//
// Package:    SoftElectronInJetAnalyzer/ElectronID_NORECODEBUG
// Class:      ElectronID_NORECODEBUG
// 
/**\class ElectronID_NORECODEBUG ElectronID_NORECODEBUG.cc SoftElectronInJetAnalyzer/ElectronID_NORECODEBUG/plugins/ElectronID_NORECODEBUG.cc

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

#include "DataFormats/BTauReco/interface/JetTag.h"
#include "DataFormats/JetReco/interface/Jet.h"
#include "SimDataFormats/JetMatching/interface/JetFlavour.h"
#include "SimDataFormats/JetMatching/interface/JetFlavourMatching.h"

  #include <TLorentzVector.h>
//
  #include "GenToRecoClass.h"
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



class ElectronID_NORECODEBUG : public edm::EDAnalyzer {
   public:
      explicit ElectronID_NORECODEBUG(const edm::ParameterSet&);
      ~ElectronID_NORECODEBUG();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;
      void PFElecFiller();
      bool  isaBhadron( int );
      bool  isaDhadron( int );
      bool  isaV( int );
      int GetOrigin(reco::GenParticle&);
      bool Verbose;
      InputTag PVerTag,PFCandidateTag_,GEDGsfElecTag_,genParticleTag_,TrackTag_,allConversionsTag,offlineBeamSpotTag;
      float maxLIP;
      int minHits;
      GenToRecoFiller *rtgf;
      Handle<reco::VertexCollection> PrimaryVetices;
      Handle<GenParticleCollection> GPC;
      Handle<reco::PFCandidateCollection> PFCandidates;
 edm::Handle<edm::View<reco::GsfElectron> > GEDGsfElecs;//     Handle<reco::GsfElectronCollection> GEDGsfElecs;	
      Handle<reco::ConversionCollection> hConversions;
      edm::Handle<BeamSpot> beamSpot;

  edm::Handle<edm::ValueMap<bool> > veto_id_decisions;
  edm::Handle<edm::ValueMap<bool> > loose_id_decisions;
  edm::Handle<edm::ValueMap<bool> > medium_id_decisions;
  edm::Handle<edm::ValueMap<bool> > tight_id_decisions; 
  edm::Handle<edm::ValueMap<float> > mvaValues;
      const reco::Vertex* vertex;
      int nPV;
      GenParticleCollection gpc ;
      reco::VertexCollection pvc;
      bool goodvertex;


      edm::EDGetTokenT<edm::ValueMap<bool> > eleVetoIdMapToken_;
      edm::EDGetTokenT<edm::ValueMap<bool> > eleLooseIdMapToken_;
      edm::EDGetTokenT<edm::ValueMap<bool> > eleMediumIdMapToken_;
      edm::EDGetTokenT<edm::ValueMap<bool> > eleTightIdMapToken_;
      edm::EDGetTokenT<edm::ValueMap<float> > mvaValuesMapToken_;
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
ElectronID_NORECODEBUG::ElectronID_NORECODEBUG(const edm::ParameterSet& iConfig):
        Verbose(iConfig.getParameter<bool>         ("Verbose") ),
        PVerTag(iConfig.getParameter<edm::InputTag>("PVerTag") ),
        PFCandidateTag_(iConfig.getParameter<edm::InputTag>("PFCandidateTag")),
	GEDGsfElecTag_(iConfig.getParameter<edm::InputTag>("GEDGsfElecTag")),
        genParticleTag_(iConfig.getParameter<edm::InputTag>("genParticleTag")),
	allConversionsTag(iConfig.getParameter<edm::InputTag>("conversionTag")),
	offlineBeamSpotTag(iConfig.getParameter<edm::InputTag>("beamSpotTag")),
	maxLIP(iConfig.getParameter<double>("maximumLongitudinalImpactParameter")),
        minHits(iConfig.getParameter<unsigned int>("minHits")),
  eleVetoIdMapToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleVetoIdMap"))),
  eleLooseIdMapToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleLooseIdMap"))),
  eleMediumIdMapToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleMediumIdMap"))),
  eleTightIdMapToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleTightIdMap"))),
  mvaValuesMapToken_(consumes<edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>("mvaValuesMap")))
{
	cout<<"testcons"<<endl;
        rtgf= new GenToRecoFiller("PF");
	cout<<"testcons ok"<<endl;

}




ElectronID_NORECODEBUG::~ElectronID_NORECODEBUG()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
ElectronID_NORECODEBUG::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
        iEvent.getByLabel(PVerTag, PrimaryVetices);
        iEvent.getByLabel(genParticleTag_, GPC);
        iEvent.getByLabel(PFCandidateTag_, PFCandidates);
	iEvent.getByLabel(GEDGsfElecTag_,GEDGsfElecs);
	iEvent.getByLabel(allConversionsTag, hConversions);
	iEvent.getByLabel(offlineBeamSpotTag, beamSpot);

  iEvent.getByToken(eleVetoIdMapToken_ ,veto_id_decisions);
  iEvent.getByToken(eleLooseIdMapToken_ ,loose_id_decisions);
  iEvent.getByToken(eleMediumIdMapToken_,medium_id_decisions);
  iEvent.getByToken(eleTightIdMapToken_,tight_id_decisions);
  iEvent.getByToken(mvaValuesMapToken_,mvaValues);

        if(!PrimaryVetices.isValid() || PrimaryVetices->empty()) return;
        vertex=&PrimaryVetices->front();
        nPV=PrimaryVetices->size();
        gpc = *(GPC.product());
        pvc = *(PrimaryVetices);
        goodvertex=!pvc.empty()?true:false;
        if(goodvertex)vertex=&pvc.front();
        PFElecFiller();

}


void ElectronID_NORECODEBUG::PFElecFiller(){
	for (int j = 0 ; j < (int)gpc.size(); ++j)
        {
		cout<<"ok"<<endl;
		rtgf->initGenToRecoFillerObject();
                bool isFinal=gpc[j].status()==1;
                bool isElecOrPionOrKaon=fabs(gpc[j].pdgId())==11 || fabs(gpc[j].pdgId())==211 || fabs(gpc[j].pdgId())==321;
                bool isKineOk=gpc[j].pt()>1.5 && fabs(gpc[j].eta())<2.4;
                if(isFinal && isElecOrPionOrKaon && isKineOk){
			cout<<"---->okok"<<endl;
			rtgf->pdgId=fabs(gpc[j].pdgId());
                        rtgf->origin=GetOrigin(gpc[j]);
                        rtgf->ptGen=gpc[j].pt();
                        rtgf->pGen=gpc[j].p();
                        rtgf->etaGen=gpc[j].eta();
			for (reco::PFCandidateCollection::const_iterator  pfe=PFCandidates->begin(); pfe!=PFCandidates->end(); ++pfe){
				if(pfe->particleId()==2 && deltaR(gpc[j].eta(), gpc[j].phi(), pfe->eta(),pfe->phi())<0.01){
					cout<<"--------> okokok"<<endl;
					GsfElectronRef elec=pfe->gsfElectronRef();
					rtgf->ecalseed=elec->core()->ecalDrivenSeed();
					rtgf->trkseed=elec->core()->trackerDrivenSeed();
					rtgf->pf_pt=elec->pt();
					float pt_surr=0;
					for (reco::PFCandidateCollection::const_iterator  pfe_surr=PFCandidates->begin(); pfe_surr!=PFCandidates->end(); ++pfe_surr){
						float dr=deltaR(pfe_surr->eta(), pfe_surr->phi(), pfe->eta(),pfe->phi());
						if(dr!=0 && dr<0.3)pt_surr=+pfe_surr->pt();
					}
					cout<<"--------> okokok a"<<endl;
					rtgf->pf_isol=pt_surr/elec->pt();
					rtgf->id_veto= (*veto_id_decisions)[elec];
					cout<<"--------> okokok b"<<endl;
    					rtgf->id_loose  = (*loose_id_decisions)[elec];
    					rtgf->id_medium  = (*medium_id_decisions)[elec];
					rtgf->id_tight = (*tight_id_decisions)[elec];
					cout<<"--------> okokok c"<<endl;
					rtgf->mva = (*mvaValues)[elec];
					cout<<"--------> okokok d"<<endl;
				}
			}
			for (size_t i = 0; i < GEDGsfElecs->size(); ++i){
//			for (reco::GsfElectronCollection::const_iterator  gsf=GEDGsfElecs->begin(); gsf!=GEDGsfElecs->end(); ++gsf){
				const auto el = GEDGsfElecs->ptrAt(i);
                                if(deltaR(gpc[j].eta(), gpc[j].phi(), el->eta(),el->phi())<0.01){
                                        rtgf->mva_gsf=(*mvaValues)[el];
                                        rtgf->ecalseed_gsf=el->core()->ecalDrivenSeed();
                                        rtgf->trkseed_gsf=el->core()->trackerDrivenSeed();
                                }
                       // }
			}
			rtgf->tree_efficiency->Fill();
		}
	}
	
}


int ElectronID_NORECODEBUG::GetOrigin(reco::GenParticle& part){
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

bool ElectronID_NORECODEBUG::isaV(int pidAbs){
        int res=false;
        if(pidAbs==24 || pidAbs==23 || pidAbs==22)res=true;
        return res;
}


bool ElectronID_NORECODEBUG::isaBhadron(int pidAbs){

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

bool ElectronID_NORECODEBUG::isaDhadron(int pidAbs){

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



void
ElectronID_NORECODEBUG::beginJob()
{
}

void
ElectronID_NORECODEBUG::endJob()
{
        rtgf->WriteInFileAndCloseIt();
}

void
ElectronID_NORECODEBUG::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
DEFINE_FWK_MODULE(ElectronID_NORECODEBUG);
