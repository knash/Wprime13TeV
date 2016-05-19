#! /usr/bin/env python
import re
import os
import subprocess
from os import listdir
from os.path import isfile, join
import glob
import math
import ROOT
from ROOT import *
import sys
from DataFormats.FWLite import Events, Handle
from optparse import OptionParser
parser = OptionParser()
parser.add_option('-c', '--cuts', metavar='F', type='string', action='store',
                  default	=	'rate_default',
                  dest		=	'cuts',
                  help		=	'Cuts type (ie default, rate, etc)')
(options, args) = parser.parse_args()

cuts = options.cuts

import Wprime_Functions	
from Wprime_Functions import *

#Load up cut values based on what selection we want to run 
Cons = LoadConstants()
lumi = Cons['lumi']
kfac = Cons['kfac']
xsec_wpr = Cons['xsec_wpr']
xsec_wpl = Cons['xsec_wpl']
xsec_wplr = Cons['xsec_wplr']
xsec_ttbar = Cons['xsec_ttbar']
xsec_qcd = Cons['xsec_qcd']
xsec_st = Cons['xsec_st']
nev_wpr = Cons['nev_wpr']
nev_wpl = Cons['nev_wpl']
nev_wplr = Cons['nev_wplr']
nev_ttbar = Cons['nev_ttbar']
nev_qcd = Cons['nev_qcd']
nev_st = Cons['nev_st']


files = sorted(glob.glob("*job*of*.root"))

j = []
print files
for f in files:
	j.append(f.replace('_jo'+re.search('_jo(.+?)_PSET', f).group(1),""))

files_to_sum = list(set(j))

commands = []
commands.append('rm *.log') 
commands.append('rm temprootfiles/*.root')
commands.append('rm -rf notneeded')
for f in files_to_sum:
	commands.append('rm '+f) 
	commands.append('hadd ' + f + " " + f.replace('_PSET','_job*_PSET') )
	commands.append('mv ' +  f.replace('_PSET','_job*_PSET') + ' temprootfiles/')
commands.append('mv TBratefiledata_PSET_*.root rootfiles/')


for s in commands :
    print 'executing ' + s
    subprocess.call( [s], shell=True )

commands = []

# ttbar weighting and organization
commands.append('rm rootfiles/TBratefilettbar_PSET_'+cuts+'weighted.root') #removes old weighted file with same name in /rootfiles/
commands.append('python HistoWeight.py -i TBratefilettbar_PSET_'+cuts+'.root -o rootfiles/TBratefilettbar_PSET_'+cuts+'weighted.root -w ' + str(lumi*xsec_ttbar['PH']/nev_ttbar['PH']))
commands.append('mv TBratefilettbar_PSET_'+cuts+'.root temprootfiles/')
commands.append('mv TBratefilettbar_PSET_'+cuts+'weighted.root rootfiles/')

# QCD weighting and organization
#ptarray = [300, 470, 600, 800, 1000, 1400]

#commands.append('rm ' + 'TBratefileQCD_PSET_'+cuts+'weighted.root')
#commands.append('rm ' + 'TBratefileQCD_PSET_'+cuts+'.root')
#commands.append('hadd ' + 'TBratefileQCD_PSET_'+cuts+'.root' + " " +'TBratefileQCDPT*_PSET_'+cuts+'.root')	#adds the separated pt files into one
#for pti in ptarray:
#	pt = str(pti)
#	commands.append('rm ' + 'TBratefileQCDPT'+pt+'_PSET_'+cuts+'weighted.root')
#	commands.append('python HistoWeight.py -i TBratefileQCDPT'+pt+'_PSET_'+cuts+'.root -o TBratefileQCDPT'+pt+'_PSET_'+cuts+'weighted.root -w ' + str(lumi*xsec_qcd[pt]/nev_qcd[pt])) #weights individual pt files by their appropriate weight

#commands.append('hadd ' + 'TBratefileQCD_PSET_'+cuts+'weighted.root' + " " + 'TBratefileQCDPT*_PSET_'+cuts+'weighted.root') #adds the separated weighted files together	
#commands.append('mv ' + 'TBratefileQCDPT*_PSET_'+cuts+'.root' + " " + 'temprootfiles/')		#moves the individual pt files to temp
#commands.append('mv ' + 'TBratefileQCDPT*_PSET_'+cuts+'weighted.root' + " " + 'temprootfiles/')	#moves the invididual weighted pt files to temp

#commands.append('rm temprootfiles/TBratefileQCD_PSET_'+cuts+'.root')		#removes prexisting ratefile
#commands.append('rm rootfiles/TBratefileQCD_PSET_'+cuts+'weighted.root')	#removes prexisting weighted ratefile
#commands.append('mv TBratefileQCD_PSET_'+cuts+'.root temprootfiles/')		#moves new ratefile to /temprootfiles
#commands.append('mv TBratefileQCD_PSET_'+cuts+'weighted.root rootfiles/')	#moves new weighted ratefile to /rootfiles

STfiles = sorted(glob.glob('TBratefileS*_PSET_'+cuts+'.root'))
for f in STfiles:
	
	chan = f.replace('TBratefileST','').replace('_PSET_'+cuts+'.root','')	
	if chan=='':
		continue 
	print "Single top " + chan + " channel"

	xsec_ST = xsec_st[chan]
	nev_ST = nev_st[chan]
	commands.append('rm ' + f.replace('TBratefileST','TBratefileweightedST'))	 
	commands.append('python HistoWeight.py -i '+f+' -o '+f.replace('TBratefileST','TBratefileweightedST')+' -w ' + str(lumi*xsec_ST/nev_ST))
	commands.append('mv '+f+' temprootfiles/')
commands.append('rm TBratefileST_PSET_'+cuts+'.root')
commands.append('hadd TBratefileST_PSET_'+cuts+'.root TBratefileweightedST*'+cuts+'.root')
commands.append('mv TBratefileST_PSET_'+cuts+'.root rootfiles/')
commands.append('mv TBratefileweightedST*'+cuts+'.root temprootfiles/')



qcdfiles = sorted(glob.glob('TBratefileQCDHT*_PSET_'+cuts+'.root'))
for f in qcdfiles:
	pt = f.lstrip('TBratefileQCDHT').rstrip('_PSET_'+cuts+'.root')
	print "QCD HT " + pt
	xsec_QCD = xsec_qcd['HT'+pt]
	nev_QCD = nev_qcd['HT'+pt]
	commands.append('rm ' + f.replace('TBratefileQCDHT','TBratefileweightedQCDHT'))	 
	commands.append('python HistoWeight.py -i '+f+' -o '+f.replace('TBratefileQCDHT','TBratefileweightedQCDHT')+' -w ' + str(lumi*xsec_QCD/nev_QCD))
	commands.append('mv '+f+' temprootfiles/')
commands.append('hadd TBratefileQCD_PSET_'+cuts+'.root TBratefileweightedQCDHT*'+cuts+'.root')
commands.append('mv TBratefileQCD_PSET_'+cuts+'.root rootfiles/')
commands.append('mv TBratefileweightedQCDHT*'+cuts+'.root temprootfiles/')



for coup in ['right']:
	sigfiles = sorted(glob.glob('TBratefilesignal'+coup+'*_PSET_'+cuts+'.root'))
	for f in sigfiles:
		mass = f.lstrip('TBratefilesignal'+coup).rstrip('_PSET_'+cuts+'.root')
		print "Signal " + coup + " at " + mass 
		xsec_sig = xsec_wpr[mass]*kfac
		nev_sig = nev_wpr[mass]
		commands.append('rm ' + f.replace('TBratefilesignal'+coup,'TBratefileweightedsignal'+coup))	 
		commands.append('python HistoWeight.py -i '+f+' -o '+f.replace('TBratefilesignal'+coup,'TBratefileweightedsignal'+coup)+' -w ' + str(lumi*xsec_sig/nev_sig))
		commands.append('mv '+f+' temprootfiles/')
		commands.append('mv '+f.replace('TBratefilesignal'+coup,'TBratefileweightedsignal'+coup)+' rootfiles/')
for s in commands :
    print 'executing ' + s
    subprocess.call( [s], shell=True )







