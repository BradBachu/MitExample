#ifndef MITEXAMPLE_MODS_NTUPLESMOD_H
#define MITEXAMPLE_MODS_NTUPLESMOD_H

#include "MitAna/TreeMod/interface/BaseMod.h"
#include "MitAna/DataTree/interface/ElectronCol.h"
#include "MitAna/DataTree/interface/PhotonCol.h"
#include "MitAna/DataTree/interface/TriggerObjectCol.h"
#include "MitExample/DataFormats/interface/TnPEvent.h"
#include "MitAna/DataTree/interface/VertexFwd.h"

#include "TTree.h"
#include "TString.h"

namespace mithep {
  
  class NtuplesMod : public BaseMod {
  public:
    NtuplesMod(char const* name = "NtuplesMod", char const* title = "Flat-tree ntuples producer");
    void SetTagElectronsName(char const* _name) { fTagElectronsName = _name; }
    void SetProbePhotonsNoVetoName(char const* _name) { fProbePhotonsNoVetoName = _name; }
    void SetProbePhotonsYesVetoName(char const* _name) { fProbePhotonsYesVetoName = _name; }
    void SetTriggerObjectsName(char const* _name) { fTriggerObjectsName = _name; }
    void SetTriggerMatchName(char const* _name) { fTriggerMatchName = _name; }
    void                SetPVName(const char *n)               { fPVName = n;                }

  protected:
    void Process() override;
    void SlaveBegin() override;
    void SlaveTerminate() override;

    TString fTagElectronsName;
    TString fProbePhotonsNoVetoName;
    TString fProbePhotonsYesVetoName;
    TString fTriggerObjectsName;
    TString fTriggerMatchName;

    VertexCol const* fVertices;
    ElectronCol const* fTagElectrons;
    PhotonCol const* fProbePhotonsNoVeto;
    PhotonCol const* fProbePhotonsYesVeto;
    TriggerObjectCol const* fTriggerObjects;

    TnPEvent fEvent;
    TTree* fNtuples;
    TString fPVName;

    ClassDef(NtuplesMod, 0)
  };

}

#endif
