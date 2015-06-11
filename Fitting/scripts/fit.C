#include "TString.h"
#include "TH1F.h"
#include "TFile.h"
#include "TTree.h"
#include "TROOT.h"
#include "MitExample/DataFormats/interface/TnPEvent.h"
#include "MitExample/Fitting/interface/FittingUtils.h"
#include "TLorentzVector.h"
#include "RooDataSet.h"
#include "RooRealVar.h"
#include "RooWorkspace.h"
#include "RooAbsPdf.h"
#include "TCanvas.h"

using namespace mithep;
using namespace RooFit;
// using namespace RooStats;

void fit() {
  TFile * fIn = new TFile("mass.root");
  TTree * events = (TTree*)fIn->Get("events");
  int nEvents = events->GetEntries();

  RooWorkspace * workspace = new RooWorkspace("workspace");  
  RooArgSet * argset = (RooArgSet*)workspace->factory("{mass[0., 350.]}");
  workspace->factory("BreitWigner::bw(mass,m0[91.2,80.,100.],gamma[4.,0.,10.])");
  workspace->factory("RooCBShape::cb(mass,m1[0.,-20.,20.],s1[1.6,-100.,100.],a1[1000.,0.,10000.],n1[1.,0,30000])");
  workspace->factory("Exponential::bg(mass,alpha[0,-10.,.001])");
  workspace->factory("FCONV::sig(mass,bw,cb)");
  workspace->factory("SUM::model(f[0,1]*sig,bg)");

  RooDataSet * dataset = FittingUtils::createDataSet(events, workspace->var("mass"));
  workspace->pdf("model")->fitTo(*dataset,Minimizer("Minuit2","Migrad"));
  TCanvas * c1 = new TCanvas("c1","c1",600,400);
  c1->cd();
  RooPlot * frame = workspace->var("mass")->frame();
  dataset->plotOn(frame);
  workspace->pdf("model")->plotOn(frame);
  frame->Draw();
  c1->SaveAs("ds.pdf");

}
