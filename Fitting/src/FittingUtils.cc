#include "MitExample/Fitting/interface/FittingUtils.h"
#include "MitExample/DataFormats/interface/TnPEvent.h"

#include "RooRealVar.h"
#include "RooArgSet.h"
#include "RooDataSet.h"
#include "TTree.h"
#include "TLorentzVector.h"

RooDataSet*
mithep::FittingUtils::createDataSet(TTree* _source, RooRealVar* x, char const* _name/* = "dataset"*/, char const* _title/* = "mass"*/)
{
  Float_t massVal;
  _source->SetBranchAddress("mass",&massVal);

  RooArgSet argset(*x);
  RooDataSet* dataset(new RooDataSet(_name, _title, argset));

  long iEntry(0);
  while (_source->GetEntry(iEntry++) > 0) {
    x->setVal((Double_t)massVal);
    dataset->add(argset);
  }

  return dataset;
}
