#!/usr/bin/env python

import ROOT
inFile = open('SingleElectron_Run2012A_FileList.txt','r') # list of all relevant files, separated by newline
outFile = open('SingleElectron_Run2012A_FileList_nEntries.txt','w') # file to be read by condor submission script
for l in inFile:
	if not l:
		break
	f0 = ROOT.TFile(l[:-1])
	eventTree = f0.Get("Events")
	nEntries = eventTree.GetEntries()
	outFile.write("%s\t%d\n"%(l[:-1],nEntries))
