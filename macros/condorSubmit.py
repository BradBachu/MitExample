#!/usr/bin/env python

from sys import argv
from os import system
nEvents = int(argv[1])
fileList = open(argv[2],'r')
nBatch=1000
nEventsProcessed=0

for l in fileList:
	fileName = l[:-1].split()[0]
	nEvtsPerFile = int(l[:-1].split()[1])
	nToDo=min(nEvents-nEventsProcessed,nEvtsPerFile)
  	nBatchTmp = min(nBatch,nToDo)
	nEventsProcessed+=nToDo
	fileShortName=fileName.split('/')[-1].split('.')[0]
	print "Submitting %d events from %s"%(nToDo,fileName)

	condorSubmitFile = open('submit.cfg','w')
	condorSubmitFile.write('Executable  = condorMakeNtuples.sh \nUniverse  = vanilla\n')
	condorSubmitFile.write('Requirements = OpSysAndVer == "SL6" && ARCH == "X86_64"\n')
	condorSubmitFile.write('Error = /scratch5/snarayan/logs/%s_$(Process).err\n'%(fileShortName))
	condorSubmitFile.write('Log = /scratch5/snarayan/logs/%s_$(Process).log\n'%(fileShortName))
	condorSubmitFile.write('Output = /scratch5/snarayan/logs/%s_$(Process).out\n'%(fileShortName))
	condorSubmitFile.write('GetEnv=True\n')
	condorSubmitFile.write('Arguments = "%s %d $(Process)"\n'%(fileName,nBatchTmp))
	condorSubmitFile.write('+AccountingGroup= "group_cmsuser.snarayan"\n')
	if nToDo%nBatch:
		condorSubmitFile.write('Queue %d\n'%(int(nToDo/nBatch)+1))
	else:
		condorSubmitFile.write('Queue %d\n'%(int(nToDo/nBatch)))
	condorSubmitFile.close()
	system('condor_submit submit.cfg')
	if nEventsProcessed >= nEvents:
		break
