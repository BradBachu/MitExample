#include "TString.h"
#include "TH1F.h"
#include "TFile.h"
#include "TTree.h"
#include "TROOT.h"
#include "TMath.h"
#include "MitExample/DataFormats/interface/TnPEvent.h"
#include "TLorentzVector.h"
#include "RooDataSet.h"
#include "RooRealVar.h"
#include "RooWorkspace.h"

using namespace mithep;
using namespace RooFit;

void makeFlatNtuples(TString fileName = "/home/snarayan/cms/hist/egamma_v3/t2mit/filefi/032/SingleElectron+Run2012A-22Jan2013-v1+AOD/ntuple.root",
                     TString fileOutName = "/home/snarayan/cms/hist/egamma_v3/t2mit/filefi/032/SingleElectron+Run2012A-22Jan2013-v1+AOD/flattened.root") {
  TFile * fIn = new TFile(fileName);
  TTree * events = (TTree*)fIn->FindObjectAny("events");
  int nEvents = events->GetEntries();
  Float_t tagpx[64];
  Float_t tagpy[64];
  Float_t tagpz[64];
  Float_t tagenergy[64];
  Float_t probepx[64];
  Float_t probepy[64];
  Float_t probepz[64];
  Float_t probeenergy[64];
  UInt_t nPairs;
  UInt_t nPairsYesVeto;
  UInt_t nPairsNoVeto;
  UInt_t nVertices;
  
  events->SetBranchAddress("tag.px",&tagpx);
  events->SetBranchAddress("tag.py",&tagpy);
  events->SetBranchAddress("tag.pz",&tagpz);
  events->SetBranchAddress("tag.energy",&tagenergy);
  events->SetBranchAddress("probe.px",&probepx);
  events->SetBranchAddress("probe.py",&probepy);
  events->SetBranchAddress("probe.pz",&probepz);
  events->SetBranchAddress("probe.energy",&probeenergy);
  events->SetBranchAddress("nPairs",&nPairs);
  events->SetBranchAddress("nPairsNoVeto",&nPairsNoVeto);
  events->SetBranchAddress("nPairsYesVeto",&nPairsYesVeto);
  events->SetBranchAddress("nVertices",&nVertices);

  Float_t mass; 
  Float_t eta;
  Float_t p_T;
  Float_t phi;
  Int_t electronVetoApplied;
  TFile* fOut = new TFile(fileOutName,"RECREATE");
  TTree * tOut = new TTree("eventsSkimmed","eventsSkimmed");
  tOut->Branch("mass",&mass);
  tOut->Branch("eta",&eta);
  tOut->Branch("p_T",&p_T);
  tOut->Branch("phi",&phi);
  tOut->Branch("nVertices",&nVertices);
  tOut->Branch("electronVetoApplied",&electronVetoApplied);
  for (unsigned i=0;i <nEvents; i++) {
    events->GetEntry(i);
    for (unsigned j=0;j<(nPairsNoVeto+nPairsYesVeto); j++) {
      TLorentzVector vTag(tagpx[j],tagpy[j],tagpz[j],tagenergy[j]);
      TLorentzVector vProbe(probepx[j],probepy[j],probepz[j],probeenergy[j]);
      mass = (vTag+vProbe).M();
      if (mass<0. || mass > pow(10,10)) continue; // probably some stupid shit
      eta = TMath::Log((probeenergy[j]+probepz[j])/(probeenergy[j]-probepz[j]))/2.;
      p_T = TMath::Sqrt(probepx[j]*probepx[j] + probepy[j]*probepy[j]);
      phi = TMath::ATan(probepy[j]/probepx[j]);
      if (j >= nPairsNoVeto) {
        electronVetoApplied = 1;
      } else {
        electronVetoApplied = 0;
      }
      tOut->Fill();
    }
  }
  
  fOut->cd();
  tOut->Write("events");
  fOut->Close(); 

}
