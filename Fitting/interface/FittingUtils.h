#ifndef MITEXAMPLE_FITTING_FITTINGUTILS_H
#define MITEXAMPLE_FITTING_FITTINGUTILS_H

class RooDataSet;
class RooArgSet;
class RooRealVar;
class TTree;

namespace mithep {
  class FittingUtils {
    public:
      static RooDataSet* createDataSet(TTree*, RooRealVar*, char const* name = "dataset", char const* title ="data with mass");
  };

}

#endif
