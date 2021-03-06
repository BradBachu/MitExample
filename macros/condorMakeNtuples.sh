#!/bin/bash

fileName=$1
scramDir='/home/snarayan/cms/cmssw/040/CMSSW_7_4_0'
outDir='/scratch5/snarayan/data/SingleElectron+Run2012A-22Jan2013-v1+AOD'
rootDir='/home/snarayan/cms/root'
nBatch=${2}
let nSkip=${3}*nBatch

cd $scramDir
eval `scramv1 runtime -sh`
cd $outDir
pwd
# cp ${rootDir}/.rootlogon.C .
python ${scramDir}/src/MitExample/macros/makeNtuples.py $fileName $nBatch $nSkip
