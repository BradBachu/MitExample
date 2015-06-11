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

void makeFlatNtuples(TString fileName = "/home/snarayan/cms/cmssw/040/CMSSW_7_4_0/src/MitExample/macros/ntuples.root",
                     TString fileOutName = "skimmed.root") {
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
  UInt_t nVerticesIn;
  UInt_t nVerticesOut;
  events->SetBranchAddress("tag.px",&tagpx);
  events->SetBranchAddress("tag.py",&tagpy);
  events->SetBranchAddress("tag.pz",&tagpz);
  events->SetBranchAddress("tag.energy",&tagenergy);
  events->SetBranchAddress("probe.px",&probepx);
  events->SetBranchAddress("probe.py",&probepy);
  events->SetBranchAddress("probe.pz",&probepz);
  events->SetBranchAddress("probe.energy",&probeenergy);
  events->SetBranchAddress("nPairs",&nPairs);
  events->SetBranchAddress("nVertices",&nVerticesIn);
  Float_t mass; 
  Float_t eta;
  Float_t p_T;
  TFile* fOut = new TFile(fileOutName,"RECREATE");
  TTree * tOut = new TTree("eventsSkimmed","eventsSkimmed");
  tOut->Branch("mass",&mass);
  tOut->Branch("eta",&eta);
  tOut->Branch("p_T",&p_T);
  tOut->Branch("nVertices",&nVerticesOut);
  for (unsigned i=0;i <nEvents; i++) {
    events->GetEntry(i);
    for (unsigned j=0;j<nPairs; j++) {
      TLorentzVector vTag(tagpx[j],tagpy[j],tagpz[j],tagenergy[j]);
      TLorentzVector vProbe(probepx[j],probepy[j],probepz[j],probeenergy[j]);
      mass = (vTag+vProbe).M();
      eta = TMath::Log((probeenergy[j]+probepz[j])/(probeenergy[j]-probepz[j]))/2.;
      p_T = TMath::Sqrt(probepx[j]*probepx[j] + probepy[j]*probepy[j]);
      nVerticesOut = nVerticesIn;
      tOut->Fill();
    }
  }
  
  fOut->cd();
  tOut->Write("events");
  fOut->Close(); 

}
