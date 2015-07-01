#include "MitExample/Mods/interface/MeasureMCEffMod.h"
#include "TVector2.h"
#include "TTree.h"

ClassImp(mithep::MeasureMCEffMod);

  mithep::MeasureMCEffMod::MeasureMCEffMod(char const* name, char const* title) :
    BaseMod(name, title),
    fTagElectronsName("TagElectrons"),
    fProbePhotonsNoVetoName("ProbePhotons"),
    fProbePhotonsYesVetoName("ProbePhotons"),
    fMCParticlesName("MCParticles"),
    fTagElectrons(0),
    fProbePhotonsNoVeto(0),
    fProbePhotonsYesVeto(0),
    fMCParticles(0),
    fEta(0),
    fVetoed(0),
    fNtuples(0)
  {
  }

  void
  mithep::MeasureMCEffMod::Process()
  {
    LoadEventObject(fTagElectronsName, fTagElectrons);
    LoadEventObject(fProbePhotonsNoVetoName, fProbePhotonsNoVeto);
    LoadEventObject(fProbePhotonsYesVetoName, fProbePhotonsYesVeto);
    LoadBranch(fMCParticlesName);

    std::vector<MCParticle const*> MCElectrons;
    for (unsigned i=0; i < fMCParticles->GetEntries(); i++) {
      MCParticle const& eleCand = *fMCParticles->At(i);
      if (eleCand.Status()!=1) 
        continue;
      if (eleCand.PdgId()==11 || eleCand.PdgId()==-11)
        MCElectrons.push_back(&eleCand);
    }

  if (fTagElectrons->GetEntries()==0) return; // no 'tag' electrons
  std::vector<Electron const*> tags;
  for (unsigned iE(0); iE != fTagElectrons->GetEntries(); ++iE) {
    Electron const& inEle(*fTagElectrons->At(iE));
    tags.push_back(&inEle);
  }

  std::vector<Photon const*> probesNoVeto;
  for (unsigned iP(0); iP != fProbePhotonsNoVeto->GetEntries(); ++iP) {
    Photon const& inPh(*fProbePhotonsNoVeto->At(iP));
    bool pass = 1;
    for(Electron const*tag : tags) {
      if (tag->SCluster() == inPh.SCluster()) {
        pass=0;
        break;
      }
    }
    if (pass) 
      probesNoVeto.push_back(&inPh);
  }

  std::vector<Photon const*> probesYesVeto;
  for (unsigned iP(0); iP != fProbePhotonsYesVeto->GetEntries(); ++iP) {
    Photon const& inPh(*fProbePhotonsYesVeto->At(iP));
    bool pass = 1;
    for(Electron const*tag : tags) {
      if (tag->SCluster() == inPh.SCluster()) {
        pass=0;
        break;
      }
    }
    if (pass) 
      probesYesVeto.push_back(&inPh);
  }

  for (Photon const* probe : probesNoVeto) {
    for (MCParticle const* ele : MCElectrons) {
      // check that they're close
      double dEta(ele->Eta() - probe->Eta());
      double dPhi(TVector2::Phi_mpi_pi(ele->Phi() - probe->Phi()));
        if (dEta * dEta + dPhi * dPhi < 0.15 * 0.15) {
          fEta = ele->Eta();
          fVetoed = 0;
          fNtuples->Fill();
          break;
        }
    }
  }
  
  for (Photon const* probe : probesYesVeto) {
    for (MCParticle const* ele : MCElectrons) {
      // check that they're close
      double dEta(ele->Eta() - probe->Eta());
      double dPhi(TVector2::Phi_mpi_pi(ele->Phi() - probe->Phi()));
        if (dEta * dEta + dPhi * dPhi < 0.15 * 0.15) {
          fEta = ele->Eta();
          fVetoed = 1;
          fNtuples->Fill();
          break;
        }
    }
  }
}
  void
  mithep::MeasureMCEffMod::SlaveBegin()
  {
    ReqEventObject(fMCParticlesName, fMCParticles,   true);
    fNtuples = new TTree("events","electrons faking photons");
    fNtuples->Branch("eta",&fEta,"eta/F");
    fNtuples->Branch("vetoed",&fVetoed,"vetoed/i");

    AddOutput(fNtuples);
  }

  void
  mithep::MeasureMCEffMod::SlaveTerminate()
  {
  }

