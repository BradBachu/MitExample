#!/usr/bin/env python

import ROOT
inFile = open('fileList.txt','r') # list of all relevant files, separated by newline
outFile = open('fileList_nEntries.txt','w') # file to be read by condor submission script
for l in inFile:
  print l
  if not l:
    break
  f0 = ROOT.TFile(l[:-1])
  eventTree = f0.Get("Events")
  nEntries = eventTree.GetEntries()
  outFile.write("%s\t%d\n"%(l[:-1],nEntries))
