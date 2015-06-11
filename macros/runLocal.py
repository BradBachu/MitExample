#!/usr/bin/python

from os import system
from sys import stderr

fIn = open("fileList_nEntries.txt")

totalEvts=0

for l in fIn:
    fileName = l[:-1].split()[0]
    nEvtsPerFile = int(l[:-1].split()[1])
    stderr.write("Getting %d events from %s\n"%(nEvtsPerFile,fileName))
    system('python $CMSSW_BASE/src/MitExample/macros/makeNtuples.py %s %d 0'%(fileName,nEvtsPerFile))
    totalEvts+=nEvtsPerFile
