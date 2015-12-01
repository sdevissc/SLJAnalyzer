// -*- C++ -*-
//
// Package:    SoftElectronInJetAnalyzer/ElecIsolation
// Class:      ElecIsolation
// 
/**\class ElecIsolation ElecIsolation.cc SoftElectronInJetAnalyzer/ElecIsolation/plugins/ElecIsolation.cc

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

#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/PatCandidates/interface/Electron.h"



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



class ElecIsolation : public edm::EDAnalyzer {
   public:
      explicit ElecIsolation(const edm::ParameterSet&);
      ~ElecIsolation();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;
      void PFElecFiller(int);
      bool  isaBhadron( int );
      bool  isaDhadron( int );
      bool  isaV( int );
      bool isAncestor(const reco::Candidate* , const reco::Candidate *);
      int GetOrigin(const pat::PackedGenParticle);
      edm::EDGetTokenT<edm::View<pat::Electron> > elecToken;
      edm::EDGetTokenT<edm::View<pat::PackedGenParticle> > packedGenToken;
      edm::EDGetTokenT<edm::View<reco::GenParticle> > prunedToken;
      edm::EDGetTokenT<edm::View<pat::PackedCandidate> > PFToken;
      GenToRecoFiller *rtgf;

      Handle<edm::View<pat::Electron> > elec;
      Handle<edm::View<reco::GenParticle> > pruned;
      Handle<edm::View<pat::PackedGenParticle> > packed;
      Handle<edm::View<pat::PackedCandidate > > PF;

      edm::Handle<edm::ValueMap<float> > mvaValues;
      edm::EDGetTokenT<edm::ValueMap<float> > mvaValuesMapToken_;
      int type; 
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
ElecIsolation::ElecIsolation(const edm::ParameterSet& iConfig):
	elecToken(consumes<edm::View<pat::Electron> >(iConfig.getParameter<edm::InputTag>("GEDGsfElecTag"))),
	packedGenToken(consumes<edm::View<pat::PackedGenParticle> >(iConfig.getParameter<edm::InputTag>("packed"))),
        prunedToken(consumes<edm::View<reco::GenParticle> >(iConfig.getParameter<edm::InputTag>("genParticleTag"))),
	PFToken(consumes<edm::View<pat::PackedCandidate> >(iConfig.getParameter<edm::InputTag>("packedPFCandidatesTag"))),
	mvaValuesMapToken_(consumes<edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>("mvaValuesMap"))),
	type(iConfig.getParameter<unsigned int>("Type"))
{
        rtgf= new GenToRecoFiller("PF");
}




ElecIsolation::~ElecIsolation()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
ElecIsolation::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
        iEvent.getByToken(packedGenToken,packed);
        iEvent.getByToken(prunedToken,pruned);
	iEvent.getByToken(elecToken,elec);
	iEvent.getByToken(PFToken,PF);
	iEvent.getByToken(mvaValuesMapToken_,mvaValues);
 	
        PFElecFiller(type);

}

bool ElecIsolation::isAncestor(const reco::Candidate* ancestor, const reco::Candidate * particle)
{    
	if(ancestor == particle ) return true;
 	for(size_t i=0;i< particle->numberOfMothers();i++){
		cout<<"ancestor "<<particle->pdgId();
                if(isAncestor(ancestor,particle->mother(i))) return true;
        }
 return false;
}



void ElecIsolation::PFElecFiller(int type){
	cout<<"A new event, first initializing"<<endl;
	for(int lep=0;lep<(int)elec->size();lep++){
		rtgf->initGenToRecoFillerObject();
		rtgf->isPF=((*elec)[lep].isPF()==true)?1:0;
		rtgf->origin=0;
		cout		<<" FOUND a RECO electron at reco level with pt="<<(*elec)[lep].pt()
				<<" and eta= "<<(*elec)[lep].eta()<<" phi="<<(*elec)[lep].phi()
				<<", will now run on the "<<packed->size()<<" gen particles"<<endl;
		for (int j = 0 ; j < (int)packed->size(); ++j){
			float dR=deltaR((*elec)[lep].eta(),(*elec)[lep].phi(),(*packed)[j].eta(),(*packed)[j].phi());
			if(dR>0.01)continue;
			cout<<"Tried with pdgId="<<(*packed)[j].pdgId()
				<<" Pt="<<(*packed)[j].pt()
				<<" eta="<<(*packed)[j].eta()
				<<" phi="<<(*packed)[j].phi()
				<<"===> dr="<<dR<<endl;
			for(size_t i = 0; i < pruned->size(); i++){
				cout<<"Prune "<<i<<endl;
				float Zancestor=abs((*pruned)[i].pdgId()) == 23;
				float Bancestor=abs((*pruned)[i].pdgId())>500 && abs((*pruned)[i].pdgId()) <600;
				cout<<"Prune "<<i<<" "<<Zancestor<<" "<<Bancestor<<endl;
                               	if( Zancestor || Bancestor){
					cout<<"Now checking the ancestorship"<<endl;
                                       	const Candidate * MassiveGuy = &(*pruned)[i];
                                       	const Candidate * motherInPrunedCollection = (*packed)[j].mother(0) ;
                                       	if(motherInPrunedCollection != nullptr && isAncestor( MassiveGuy , motherInPrunedCollection)){
                                       	        if(Zancestor){
							rtgf->origin=1;cout<<"--->---> matched to a Z"<<endl;
							rtgf->pdgId=abs((*packed)[j].pdgId());
						}
                                       	        else if(Bancestor){
							rtgf->origin=4;cout<<"--->---> matched to a B"<<endl;
							rtgf->pdgId=abs((*packed)[j].pdgId());
						}
                                       	}
                               	}
                        }
		}
		float sumpt=0;
		//float ptPF=0;
	
		cout<<"Setting charged, neutral and PU to zero"<<endl;
		double charged = 0, neutral = 0, pileup  = 0;
		std::vector<reco::CandidatePtr> footprint;
		cout<<"Now building the footprint"<<endl;
        	for (unsigned int i = 0, n = (*elec)[lep].numberOfSourceCandidatePtrs(); i < n; ++i) {
	            cout<<"-->find a PF, pushing back to the footprint vector"<<endl;
        	    footprint.push_back((*elec)[lep].sourceCandidatePtr(i));
		 }
		cout<<"Now running on the PF vector, looking at elements sufficiently close (0.3)"<<endl;
		for (unsigned int i = 0, n = PF->size(); i < n; ++i) {
            		const pat::PackedCandidate &pf = (*PF)[i];
            		if (deltaR(pf,(*elec)[lep]) < 0.3) {
				cout<<"--->dR ok--> is that an element of the footprint? (pt="<<pf.pt()<<")"<<endl;
				if (std::find(footprint.begin(), footprint.end(), reco::CandidatePtr(PF,i)) != footprint.end()) {
                   			continue;
                		}
				cout<<"--->no--> continue"<<endl;
                		if (pf.charge() == 0) {
					cout<<"--->---> neutral PF, adding to neutral is pt>0.5"<<endl;
                    			if (pf.pt() > 0.5) neutral += pf.pt();
					cout<<"--->---> neutral is now "<<neutral<<endl;
                			} 
				else if (pf.fromPV() >= 2) {
					cout<<"--->---> charged and from PV, adding to charged"<<endl;
                    			charged += pf.pt();
					cout<<"--->---> charged is now "<<charged<<endl;
                		} else {
					 cout<<"--->---> charged and from PU, adding to PU"<<endl;
                    			if (pf.pt() > 0.5) pileup += pf.pt();
					cout<<"--->---> PU  is now "<<pileup<<endl;
                		}
			}
		}
		double iso = charged + neutral;
		cout<<"new computation of PFiso is "<<iso/(*elec)[lep].pt()<<endl;
		for (int k = 0 ; k < (int)PF->size(); ++k){
			float dR=deltaR((*elec)[lep].eta(),(*elec)[lep].phi(),(*PF)[k].eta(),(*PF)[k].phi());	
			if(dR<0.3)sumpt+=(*PF)[k].pt();
		}
		cout<<"Summary: origin="<<rtgf->origin<<" pdg="<<rtgf->pdgId<<endl;
                rtgf->ptGen=(*elec)[lep].pt();
                rtgf->etaGen=abs((*elec)[lep].eta());
		float dr03TkSumPt=(*elec)[lep].dr03TkSumPt();
                float dr03EcalRecHitSumEt=(*elec)[lep].dr03EcalRecHitSumEt();
                float dr03HcalTowerSumEt=(*elec)[lep].dr03HcalTowerSumEt();
		rtgf->detIso=( dr03TkSumPt + dr03EcalRecHitSumEt + dr03HcalTowerSumEt) / (*elec)[lep].pt();
                rtgf->dr03TkSumPt=(*elec)[lep].dr03TkSumPt();
                rtgf->dr03EcalRecHitSumEt=(*elec)[lep].dr03EcalRecHitSumEt();
                rtgf->dr03HcalTowerSumEt=(*elec)[lep].dr03HcalTowerSumEt();
                rtgf->ecalPFClusterIso=(*elec)[lep].ecalPFClusterIso();
                rtgf->hcalPFClusterIso=(*elec)[lep].hcalPFClusterIso();
                rtgf->mva = (*mvaValues)[elec->ptrAt(lep)];
		GsfElectron::PflowIsolationVariables pfIso = (*elec)[lep].pfIsolationVariables();
                float isoChargedHadrons = pfIso.sumChargedHadronPt;
                float isoNeutralHadrons = pfIso.sumNeutralHadronEt;
                float isoPhotons = pfIso.sumPhotonEt;
                float isoChargedHadFromPileup = pfIso.sumPUPt;
		rtgf->isoChargedHadrons = pfIso.sumChargedHadronPt;
                rtgf->isoNeutralHadrons = pfIso.sumNeutralHadronEt;
                rtgf->isoPhotons = pfIso.sumPhotonEt;
                rtgf->isoChargedHadFromPileup = pfIso.sumPUPt;
		rtgf->ManualIsoPF = iso/(*elec)[lep].pt();
                rtgf->pfIso=(isoChargedHadrons+max(0.0,isoNeutralHadrons+isoPhotons-0.5*isoChargedHadFromPileup))/(*elec)[lep].pt();
               	cout<<"pfIso="<<rtgf->pfIso  <<" isoChargedHadrons="<<isoChargedHadrons
                	<<" isoNeutralHadrons="<<isoNeutralHadrons
                        <<" isoPhotons="<<isoPhotons
                        <<" isoChargedHadFromPileup="<<isoChargedHadFromPileup
                        <<" detIso="<<rtgf->detIso<<endl;
                rtgf->tree_efficiency->Fill();	
		cout<<"***********FINISH*****************"<<endl;
		
	}


}


int ElecIsolation::GetOrigin(const pat::PackedGenParticle part){
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

bool ElecIsolation::isaV(int pidAbs){
        int res=false;
        if(pidAbs==24 || pidAbs==23 || pidAbs==22)res=true;
        return res;
}


bool ElecIsolation::isaBhadron(int pidAbs){

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

bool ElecIsolation::isaDhadron(int pidAbs){

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
ElecIsolation::beginJob()
{
}

void
ElecIsolation::endJob()
{
        rtgf->WriteInFileAndCloseIt();
}

void
ElecIsolation::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
DEFINE_FWK_MODULE(ElecIsolation);
