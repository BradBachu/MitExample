import ROOT
from sys import argv

try:
	fileName = argv[1]
	nBatch = int(argv[2])
	nSkip = int(argv[3])
	outFileName = '~/cms/hist/egammaLocal_v4/'+fileName.split('/')[-1].split('.')[0]+'_'+str(nSkip)+'.root'
except IndexError:
  #	fileName='/local/snarayan/3EE27797-3773-E211-B2A6-00266CF27430.root'
  fileName='/mnt/hadoop/cms/store/user/paus/filefi/032/SingleElectron+Run2012A-22Jan2013-v1+AOD/06BA13DA-4573-E211-8F2B-003048D3C880.root'
	nBatch=1000
	nSkip=0
	outFileName = 'ntuples.root'

ROOT.gROOT.SetBatch(True)

ROOT.gSystem.Load('libMitAnaDataTree.so')
ROOT.gSystem.Load('libMitPhysicsMods.so')
ROOT.gSystem.Load('libMitExampleMods.so')

mithep = ROOT.mithep

print "Setting up analysis"
analysis = mithep.Analysis()
analysis.SetOutputName(outFileName)

analysis.AddFile(fileName)
analysis.SetProcessNEvents(nBatch)
analysis.SetSkipFirstNEvents(nSkip)

hltMod = mithep.HLTMod()
hltMod.SetBitsName('HLTBits')
hltMod.SetTrigObjsName('SingleElectronTriggerObjects')
hltMod.AddTrigger('HLT_Ele27_WP80_v*')

goodPVMod = mithep.GoodPVFilterMod()
goodPVMod.SetMinVertexNTracks(0)
goodPVMod.SetMinNDof(4)
goodPVMod.SetMaxAbsZ(24.0)
goodPVMod.SetMaxRho(2.0)
goodPVMod.SetIsMC(False)
goodPVMod.SetVertexesName('PrimaryVertexes')
goodPVMod.SetOutputName('GoodVertexes')

eleIdMod = mithep.ElectronIDMod()
eleIdMod.SetPtMin(30.)
eleIdMod.SetEtaMax(2.5)
eleIdMod.SetApplyEcalFiducial(True)
eleIdMod.SetIDType('CustomTight')
eleIdMod.SetIsoType('PFIso')
eleIdMod.SetApplyConversionFilterType1(False)
eleIdMod.SetApplyConversionFilterType2(False)
eleIdMod.SetChargeFilter(False)
eleIdMod.SetApplyD0Cut(True)
eleIdMod.SetApplyDZCut(True)
eleIdMod.SetWhichVertex(0)
eleIdMod.SetNExpectedHitsInnerCut(2)
eleIdMod.SetElectronsFromBranch(True)
eleIdMod.SetInputName('Electrons')
eleIdMod.SetGoodElectronsName('TightElectrons')
eleIdMod.SetRhoType(mithep.RhoUtilities.CMS_RHO_RHOKT6PFJETS)

phoIdMod = mithep.PhotonIDMod()
phoIdMod.SetPtMin(10.0)
phoIdMod.SetOutputName('MediumPhotonsNoEVeto')
phoIdMod.SetIDType('EgammaMedium')
phoIdMod.SetIsoType('MITPUCorrected')
phoIdMod.SetApplyElectronVeto(False)
phoIdMod.SetApplyElectronVetoConvRecovery(False)
phoIdMod.SetApplyPixelSeed(False)
phoIdMod.SetApplyConversionId(False)
phoIdMod.SetApplyFiduciality(True)
phoIdMod.SetIsData(True)
phoIdMod.SetPhotonsFromBranch(True)

phoIdMod2 = mithep.PhotonIDMod()
phoIdMod2.SetPtMin(10.0)
phoIdMod2.SetOutputName('MediumPhotonsYesEVeto')
phoIdMod2.SetIDType('EgammaMedium')
phoIdMod2.SetIsoType('MITPUCorrected')
phoIdMod2.SetApplyElectronVeto(False)
phoIdMod2.SetApplyElectronVetoConvRecovery(False)
phoIdMod2.SetApplyPixelSeed(True)
phoIdMod2.SetApplyConversionId(False)
phoIdMod2.SetApplyFiduciality(True)
phoIdMod2.SetIsData(True)
phoIdMod2.SetPhotonsFromBranch(True)

ntuplesMod = mithep.NtuplesMod('NtuplesMod', 'Flat ntuples producer')
ntuplesMod.SetPVName('GoodVertexes')
ntuplesMod.SetTagElectronsName('TightElectrons')
ntuplesMod.SetProbePhotonsNoVetoName('MediumPhotonsNoEVeto')
ntuplesMod.SetProbePhotonsYesVetoName('MediumPhotonsYesEVeto')
ntuplesMod.SetTriggerObjectsName('SingleElectronTriggerObjects')
ntuplesMod.SetTriggerMatchName('hltEle27WP80TrackIsoFilter')

analysis.AddSuperModule(hltMod)
hltMod.Add(goodPVMod)
goodPVMod.Add(eleIdMod)
eleIdMod.Add(phoIdMod)
phoIdMod.Add(phoIdMod2)
phoIdMod2.Add(ntuplesMod)

analysis.Run(False)
