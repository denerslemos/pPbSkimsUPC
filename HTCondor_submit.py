#!/usr/bin/env python

'''Add important stuff'''
import os.path
import optparse
import subprocess

''' Inputs for the skim code '''
usage = 'usage: %prog [options]'
parser = optparse.OptionParser(usage)
parser.add_option('-i', '--infiles', dest='infiles', help='input list of files', default='', type='string')
parser.add_option('-o', '--outfiles', dest='outfiles', help='output files', default='', type='string')
parser.add_option('-m', '--ismontecarlo', dest='ismontecarlo', help='0 for data and >0 for MC', default='0', type='int')
parser.add_option('-n', '--multiplicity', dest='multiplicity', help='Here only tracks bellow 40', default='0', type='int')
parser.add_option('-s', '--subfiles', dest='subfiles', help='HTCondor submission file', default='', type='string')
parser.add_option('-x', '--side', dest='psides', help='side', default='0', type='int')
parser.add_option('-z', '--zdcfile', dest='zdc', help='zdc input file', default='', type='string')

(opt, args) = parser.parse_args()
inFiles = opt.infiles
outFiles = opt.outfiles
isMC = opt.ismontecarlo
mult = opt.multiplicity
subFiles = opt.subfiles
protonSide = opt.psides
ZDCfiles = opt.zdc

''' Read list of files '''
listOfFiles = open(inFiles+'.txt', 'r')
Lines = listOfFiles.readlines()
print ("Number of files/jobs: "+str(len(Lines)))


''' Start the write submission file '''
fsubfile = open(subFiles+".sub", "w")
command_lines = '''universe   = vanilla
getenv     = True
executable = sub_skim.sh
+JobFlavour           = "tomorrow"
requirements = ((OpSysAndVer =?= "AlmaLinux9") && (CERNEnvironment =?= "qa"))
RequestCpus = 1
transfer_input_files  = voms_proxy.txt
environment = "X509_USER_PROXY=voms_proxy.txt"
'''

''' Loop over files '''
i=0
for line in Lines:
    outtempfiles = open(inFiles+"_part"+str(i)+".txt", "w")
    outtempfiles.write(line)
    outtempfiles.close()
    temp = '''
log        = cond/'''+subFiles+'''_part_'''+str(i)+'''.log
output     = cond/'''+subFiles+'''_part_'''+str(i)+'''.out
error      = cond/'''+subFiles+'''_part_'''+str(i)+'''.err
arguments = '''+inFiles+'''_part'''+str(i)+'''.txt '''+outFiles+'''_'''+str(i)+'''.root '''+str(isMC)+'''  '''+str(mult)+'''  '''+str(protonSide)+'''  '''+ZDCfiles+'''.txt
queue
'''
    command_lines += temp
    i=i+1

fsubfile.write(command_lines)
fsubfile.close()
subprocess.call(["condor_submit", subFiles+".sub"])
