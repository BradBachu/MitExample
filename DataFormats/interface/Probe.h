#ifndef MITEXAMPLE_DATAFORMATS_PROBE_H
#define MITEXAMPLE_DATAFORMATS_PROBE_H

#include "Rtypes.h"

namespace mithep {

  class TnPEvent;

  class Probe {
  public:
    Probe(TnPEvent&, UInt_t);

    Float_t& pt;
    Float_t& eta;
    Float_t& phi;
    Float_t& px;
    Float_t& py;
    Float_t& pz;
    Float_t& energy;
  private:
    Probe(); // hide default constructor
  };

}

#endif
