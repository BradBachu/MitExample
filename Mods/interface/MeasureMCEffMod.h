#ifndef MITEXAMPLE_MODS_MCEFFMOD_H
#define MITEXAMPLE_MODS_MCEFFMOD_H

#include "MitAna/TreeMod/interface/BaseMod.h"
#include "MitAna/DataTree/interface/ElectronCol.h"
#include "MitAna/DataTree/interface/PhotonCol.h"
#include "MitAna/DataTree/interface/MCParticleCol.h"
#include "MitAna/DataTree/interface/MCParticle.h"
#include "MitAna/DataTree/interface/VertexFwd.h"

#include "TTree.h"
#include "TString.h"

namespace mithep {
  
  class MeasureMCEffMod : public BaseMod {
  public:
    MeasureMCEffMod(char const* name = "MeasureMCEffMod", char const* title = "Measures electron veto efficiency in MC");
    void SetProbePhotonsNoVetoName(char const* _name) { fProbePhotonsNoVetoName = _name; }
    void SetTagElectronsName(char const* _name) { fTagElectronsName = _name; }
    void SetProbePhotonsYesVetoName(char const* _name) { fProbePhotonsYesVetoName = _name; }

  protected:
    void Process() override;
    void SlaveBegin() override;
    void SlaveTerminate() override;

    TString fTagElectronsName;
    TString fProbePhotonsNoVetoName;
    TString fProbePhotonsYesVetoName;
    TString                        fMCParticlesName;

    ElectronCol const* fTagElectrons;
    PhotonCol const* fProbePhotonsNoVeto;
    PhotonCol const* fProbePhotonsYesVeto;
    const MCParticleCol           *fMCParticles;

    Float_t fEta;
    Int_t fVetoed;
    TTree* fNtuples;

    ClassDef(MeasureMCEffMod, 1);
  };

}

#endif
