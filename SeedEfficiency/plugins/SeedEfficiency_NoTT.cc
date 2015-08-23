// -*- C++ -*-
//
// Package:    SoftElectronInJetAnalyzer/SeedEfficiency_noTT
// Class:      SeedEfficiency_noTT
// 
/**\class SeedEfficiency_noTT SeedEfficiency_noTT.cc SoftElectronInJetAnalyzer/SeedEfficiency_noTT/plugins/SeedEfficiency_noTT.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Simon De Visscher
//         Created:  Wed, 17 Dec 2014 12:38:17 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticleFwd.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "SimTracker/TrackAssociation/interface/TrackAssociatorBase.h"
#include "SimTracker/Records/interface/TrackAssociatorRecord.h"
#include "GenToRecoClass.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrackFwd.h"
#include "DataFormats/ParticleFlowReco/interface/GsfPFRecTrack.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "SimTracker/Records/interface/TrackAssociatorRecord.h"
#include "SimTracker/TrackAssociation/interface/QuickTrackAssociatorByHits.h"
#include "SimTracker/TrackAssociation/interface/TrackAssociatorByChi2.h"
#include "SimTracker/TrackAssociation/interface/QuickTrackAssociatorByHits.h"
#include "DataFormats/EgammaReco/interface/ElectronSeedFwd.h"
#include "DataFormats/EgammaReco/interface/ElectronSeed.h"
#include "DataFormats/JetReco/interface/GenJet.h"
//
// class declaration
//
using namespace edm;
using namespace std;
using namespace reco;

class SeedEfficiency_noTT : public edm::EDAnalyzer {
   public:
      explicit SeedEfficiency_noTT(const edm::ParameterSet&);
      ~SeedEfficiency_noTT();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      bool trackFilter(const reco::TrackRef &track) const;
      virtual bool isInnerMost(const reco::Track track1, const reco::Track track2, bool& sameLayer);
      bool findMatch(const reco::Candidate& genCand, Handle<TrackingParticleCollection> TPCollectionH,reco::SimToRecoCollection q,reco::Track& track,int& indexTrack,float& dR,float& sharedHits,RefToBase<Track> &trRefToBase);
      void EffEstimator();
      bool  isaBhadron( int );
      bool  isaDhadron( int );
      bool  isaV( int );
      int GetOrigin(reco::GenParticle&);
      bool UseRECO, Verbose;
	InputTag PVerTag,gedgsfElectronTag_,genParticleTag_,genJetTag_,TrackTag_,TrackingParticleTag_;
      int minHits;
      GenToRecoFiller *gtrf;
      Handle<reco::VertexCollection> PrimaryVetices;
      ESHandle<TransientTrackBuilder> builder;
      Handle<edm::View<reco::Track> > theTrackColl;
      Handle<TrackingParticleCollection> theTrackingPartColl;
      ESHandle<TrackAssociatorBase> theHitsAssociator;
      Handle<GenParticleCollection> GPC;
      Handle<reco::GsfElectronCollection> gedgsfCandidates;
      Handle<reco::GsfPFRecTrackCollection> gsfPfRecTracksCollection ;
      Handle<reco::GsfTrackCollection> gTCol;
	Handle<reco::GenJetCollection> JetMC;

      const reco::Vertex* vertex;	
      int nPV;
      const TransientTrackBuilder* transientTrackBuilder;
      reco::SimToRecoCollection theSimToRecoColl;
      reco::RecoToSimCollection theRecoToSimColl;
      GsfElectronCollection gedgsfCollection ;
      GenParticleCollection gpc ;
	GenJetCollection gjc;
      reco::VertexCollection pvc;
      bool goodvertex;
      reco::GsfPFRecTrackCollection  gsfpfTk;
      reco::GsfTrackCollection gsftc;
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
SeedEfficiency_noTT::SeedEfficiency_noTT(const edm::ParameterSet& iConfig):
        UseRECO(iConfig.getParameter<bool>         ("UseRECO") ),
        Verbose(iConfig.getParameter<bool>         ("Verbose") ),
	PVerTag(iConfig.getParameter<edm::InputTag>("PVerTag") ),
        gedgsfElectronTag_(iConfig.getParameter<edm::InputTag>("gedgsfElectronTag")),
        genParticleTag_(iConfig.getParameter<edm::InputTag>("genParticleTag")),
	genJetTag_(iConfig.getParameter<edm::InputTag>("genJetTag")),
	TrackTag_(iConfig.getParameter<edm::InputTag>("tracksTag")),
        TrackingParticleTag_(iConfig.getParameter<edm::InputTag>("TrackingParticleTag")),
        minHits(iConfig.getParameter<unsigned int>("minHits"))
{
        std::cout<<"in constr"<<std::endl;
        gtrf= new GenToRecoFiller("GTRC");
        std::cout<<"in constr"<<std::endl;
	
}



SeedEfficiency_noTT::~SeedEfficiency_noTT()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
SeedEfficiency_noTT::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
	cout<<"test"<<endl;
   using namespace edm;
        iEvent.getByLabel(PVerTag, PrimaryVetices);
        iEvent.getByLabel(genParticleTag_, GPC);
        iEvent.getByLabel(gedgsfElectronTag_, gedgsfCandidates);
        iEvent.getByLabel("pfTrackElec",gsfPfRecTracksCollection);
	iEvent.getByLabel( "electronGsfTracks", gTCol);
	iEvent.getByLabel(genJetTag_, JetMC);
	iEvent.getByLabel(TrackTag_, theTrackColl);
	cout<<"test1"<<endl;
        if(!PrimaryVetices.isValid() || PrimaryVetices->empty()) return;
        vertex=&PrimaryVetices->front();
        nPV=PrimaryVetices->size();
	gedgsfCollection = *(gedgsfCandidates.product());
        //Loading the genParticles
        cout<<"test2"<<endl;
        gpc = *(GPC.product());
	gjc = *(JetMC.product());
        pvc = *(PrimaryVetices);
	gsftc=*(gTCol.product());
        goodvertex=!pvc.empty()?true:false;
        if(goodvertex)vertex=&pvc.front();
	cout<<"test3"<<endl;
        const View<Track> tC = *(theTrackColl.product());
        gsfpfTk=*(gsfPfRecTracksCollection.product());
	cout<<"test4"<<endl;
	EffEstimator();
}

void SeedEfficiency_noTT::EffEstimator(){
        for (int j = 0 ; j < (int)gpc.size(); ++j)
        {
                bool isFinal=gpc[j].status()==1;
                bool isElecOrPionOrKaon=fabs(gpc[j].pdgId())==11 || fabs(gpc[j].pdgId())==211 || fabs(gpc[j].pdgId())==321;
                bool isKineOk=gpc[j].pt()>1.5 && fabs(gpc[j].eta())<2.4;
                if(isFinal && isElecOrPionOrKaon && isKineOk){
			gtrf->initGenToRecoFillerObject();
                        gtrf->pdgId=gpc[j].pdgId();
                        gtrf->origin=GetOrigin(gpc[j]);
			if(fabs(gpc[j].pdgId())==11 && isKineOk)std::cout<<"We have an electron or pion or kaon with pdgId="<<gtrf->pdgId<<" origin="<<gtrf->origin<<" pt="<<gpc[j].pt()<<std::endl;
                        gtrf->ptGen=gpc[j].pt();
                        gtrf->etaGen=gpc[j].eta();
                        gtrf->phiGen=gpc[j].phi();
                        gtrf->pGen=gpc[j].p();
                     //   bool inRecoJet=false;
                        Track track;
                     //   float sharedHits = 0;
                     //   float dR = 100;
                     //   int indexTrack;
                     	cout<<"Checking the jet (size="<<gjc.size()<<")"<<endl;
                        for (int k = 0 ; k < (int)gjc.size(); ++k){
				if(gjc[k].pt()>25 ){
                                	cout<<"jet "<<k<<" Pt= "<<gjc[k].pt()<<endl;
                                	if(deltaR(gpc[j].eta(), gpc[j].phi(), gjc[k].eta(),gjc[k].phi())<0.5){
                                	        cout<<"Distance to a GenJet smaller than 0.5---> inGenJet=1"<<endl;
                                	        gtrf->inGenJet=1;
                                	}
				}
                        }

                        RefToBase<Track> tr;
			if(fabs(gpc[j].pdgId())==11 && isKineOk)std::cout<<"Checking result from findMatch"<<std::endl;
/*                        bool result = findMatch(gpc[j], theTrackingPartColl, theSimToRecoColl,track, indexTrack, dR, sharedHits,tr);
                        if ( result ) {
                                gtrf->gen_match_ = 1;
                                gtrf->sharedHits=sharedHits;
				if(fabs(gpc[j].pdgId())==11 && isKineOk)std::cout<<"-->matched with TT: deltaR="<<cout<<"GsfPFRecTrack: deltaR="
					<<deltaR(gpc[j].eta(), gpc[j].phi(), tr.get()->eta(), tr.get()->phi())
					<<"  pt diff rel"<<fabs(gpc[j].pt()-tr.get()->pt())/gpc[j].pt()<<endl;
                        }
                        else          gtrf->gen_match_ = 0;
                        gtrf->inRecoJet=inRecoJet;
			if(fabs(gpc[j].pdgId())==11 && isKineOk)std::cout<<"Now checking the matching with gsfpfrectracks("<<gsfpfTk.size()<<" in total)"<<std::endl;
*/
                        for (int u = 0 ; u < (int)gsfpfTk.size(); ++u)
                        {
				if(fabs(gpc[j].pdgId())==11 && isKineOk){
					cout<<"GsfPFRecTrack: deltaR="<<deltaR(gpc[j].eta(), gpc[j].phi(), gsfpfTk[u].gsfTrackRef().get()->eta(), gsfpfTk[u].gsfTrackRef().get()->phi())
						<<"  pt diff rel"<<fabs(gpc[j].pt()-gsfpfTk[u].gsfTrackRef().get()->pt())/gpc[j].pt()<<endl;
				}
				if(deltaR(gpc[j].eta(), gpc[j].phi(), gsfpfTk[u].gsfTrackRef().get()->eta(), gsfpfTk[u].gsfTrackRef().get()->phi())<0.01 ){
                                        if(fabs(gpc[j].pdgId())==11 && isKineOk)std::cout<<"-->TRACKER DRIVEN SEEDED DeltaRMatched with seed number "<<u<<std::endl;
                                        gtrf->isDeltaRMatchedWithASeed=1;
					
                                }
				
				if(gtrf->isDeltaRMatchedWithASeed==0 && gtrf->isMatchedWithASeed==1)cout<<"************************* DeltaR Gen Track is "<<deltaR(gpc[j].eta(), gpc[j].phi(), gsfpfTk[u].gsfTrackRef().get()->eta(), gsfpfTk[u].gsfTrackRef().get()->phi())<<endl;
                        }
			if(fabs(gpc[j].pdgId())==11 && isKineOk)cout<<"Now checking the ecalSeeding("<<gsftc.size()<<" in total)"<<endl;
			for (int u = 0 ; u < (int)gsftc.size(); ++u)
                        {
				if(deltaR(gpc[j].eta(), gpc[j].phi(), gsftc[u].eta(), gsftc[u].phi())<0.01){
                                        const GsfTrackRef trackRef = edm::Ref<GsfTrackCollection>(&gsftc, u);
                                        if(fabs(gpc[j].pdgId())==11 && isKineOk)cout<<"-->Got a GEN-RECO deltaR of 0 -->Found it"<<endl;
                                        reco::ElectronSeedRef seedRef=  trackRef->seedRef().castTo<reco::ElectronSeedRef>();
                                        if(fabs(gpc[j].pdgId())==11 && isKineOk)std::cout<<"----> DeltaR match: seedRef.isAvailable() ?"<<seedRef.isAvailable()<<std::endl;
                                        if(seedRef.isAvailable() && seedRef->isEcalDriven()) {
                                                gtrf->isDeltaREcalDrivenSeeded=1;
                                                if(fabs(gpc[j].pdgId())==11 && isKineOk)std::cout<<"----> ECAL-DRIVEN SEEDED"<<std::endl;
                                        }
                                }
/*				if(!tr) continue;
				if(deltaR(tr.get()->eta(), tr.get()->phi(), gsftc[u].eta(), gsftc[u].phi())==0){
					gtrf->DistanceBetweenGenAndRECOTrack=deltaR(gpc[j].eta(), gpc[j].phi(), gsftc[u].eta(), gsftc[u].phi());
					const GsfTrackRef trackRef = edm::Ref<GsfTrackCollection>(&gsftc, u);
                                        if(fabs(gpc[j].pdgId())==11 && isKineOk)cout<<"-->Got a TT-RECO deltaR of 0 -->Found it"<<endl;
					reco::ElectronSeedRef seedRef=  trackRef->seedRef().castTo<reco::ElectronSeedRef>();
					if(fabs(gpc[j].pdgId())==11 && isKineOk)std::cout<<"----> seedRef.isAvailable() ?"<<seedRef.isAvailable()<<std::endl;
                                	if(seedRef.isAvailable() && seedRef->isEcalDriven()) {
						gtrf->isEcalDrivenSeeded=1;
						if(fabs(gpc[j].pdgId())==11 && isKineOk)std::cout<<"----> ECAL-DRIVEN SEEDED"<<std::endl;
					}
                                }
*/
				
                        }
			// Trial with no TrackingTruth Info, just DeltaR
			
						


                        gtrf->nPV=pvc.size();
                        gtrf->tree_efficiency->Fill();
                        if(fabs(gpc[j].pdgId())==11 && isKineOk)cout<<"#####################################################################"<<endl;
                }
         }
}


bool SeedEfficiency_noTT::isInnerMost(const reco::Track track1, const reco::Track track2, bool& sameLayer) {
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

int SeedEfficiency_noTT::GetOrigin(reco::GenParticle& part){
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
	cout<<"returning ============> "<<result<<endl;
        return result;
}

bool SeedEfficiency_noTT::isaV(int pidAbs){
        int res=false;
        if(pidAbs==24 || pidAbs==23 || pidAbs==22)res=true;
        return res;
}


bool SeedEfficiency_noTT::isaBhadron(int pidAbs){
        cout<<"is a B hadron? pdgId="<<pidAbs<<endl;
        bool isB = false;
        if(pidAbs>500){
                cout<<" pdgId="<<pidAbs<<endl;
                pidAbs/= 100;

                if(pidAbs<60 && pidAbs>50)pidAbs/= 10;
                cout<<"   after division by 100 or 1000: "<<pidAbs<<endl;
                int mod10 = pidAbs % 5;
                cout<<"   now the modulo is : "<<mod10<<endl;
                if(mod10 == 0) {
                        isB = true;
                        cout<<"   ===> is a B!"<<mod10<<endl;
                }
        }
        else  {
		cout<<"    NO"<<pidAbs<<endl;
        }
        return isB;
}


bool SeedEfficiency_noTT::isaDhadron(int pidAbs){
	cout<<"is a D hadron? pdgId="<<pidAbs<<endl;
        bool isD = false;
        if(pidAbs>400){
		cout<<" pdgId="<<pidAbs<<endl;
                pidAbs/= 100;
                if(pidAbs<50 && pidAbs>40)pidAbs/= 10;
		cout<<"   after division by 100 or 1000: "<<pidAbs<<endl;
                int mod10 = pidAbs % 4;
		cout<<"   now the modulo is : "<<mod10<<endl;

                if(mod10 == 0 && pidAbs<10) {
                        isD = true;
			cout<<"   ===> is a D!"<<mod10<<endl;
                }
        }
        else{
		cout<<"    NO"<<pidAbs<<endl;

        }
        return isD;
}

bool SeedEfficiency_noTT::trackFilter(const reco::TrackRef &track) const
{
        if (track->hitPattern().numberOfValidHits() < (int)minHits)
                return false;
        if (track->pt() < 0.8 )
                return false;

        return true;
}

bool SeedEfficiency_noTT::findMatch(const reco::Candidate& genCand,
                              Handle<TrackingParticleCollection> theTrackingPartColl,
                              reco::SimToRecoCollection q,
                              reco::Track& track,
                              int& indexTrack,
                              float& dR,
                              float& sharedHits,
                              RefToBase<Track> &trRefToBase) {

  const TrackingParticleCollection tPC = *(theTrackingPartColl.product());
  if ( tPC.size() == 0 ) return false;
//  cout<<"in FM: playing with "<<tPC.size()<<" trackingParticles"<<endl;
  unsigned int nTPs = tPC.size();
  for(unsigned int iTP = 0; iTP < nTPs; ++iTP) {
    const TrackingParticle& lTP(tPC[iTP]);
    float deta = lTP.eta() - genCand.eta();
    float dphi = acos(cos(lTP.phi() - genCand.phi()));
    dR = sqrt(deta*deta + dphi*dphi);
    if ( lTP.pdgId() != genCand.pdgId() || dR > 0.05 ) continue;
//	cout<<"in FM3: found a TP close enough to the gen particle "<<lTP.pdgId()<<" "<<genCand.pdgId()<<" "<<dR<<endl;
    edm::Ref<TrackingParticleCollection> tp(theTrackingPartColl, iTP);
    try {
        std::vector<std::pair<RefToBase<Track>, double> > trackV = q[tp];
//	cout<<"in FM3: "<<trackV.size()<<endl;
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
SeedEfficiency_noTT::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
SeedEfficiency_noTT::endJob() 
{
	gtrf->WriteInFileAndCloseIt();
}

// ------------ method called when starting to processes a run  ------------
/*
void 
SeedEfficiency_noTT::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void 
SeedEfficiency_noTT::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void 
SeedEfficiency_noTT::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void 
SeedEfficiency_noTT::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
SeedEfficiency_noTT::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(SeedEfficiency_noTT);
