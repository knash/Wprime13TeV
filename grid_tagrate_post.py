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
ttagsf = Cons['ttagsf']
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
commands.append('rm TBratefilettbar_PSET_'+cuts+'.root')
commands.append('python HistoWeight.py -i TBratefilettbar700_PSET_'+cuts+'.root -o temprootfiles/TBratefilettbar700_PSET_'+cuts+'weighted.root -w ' + str(lumi*xsec_ttbar['700']*ttagsf/nev_ttbar['700']))
commands.append('python HistoWeight.py -i TBratefilettbar1000_PSET_'+cuts+'.root -o temprootfiles/TBratefilettbar1000_PSET_'+cuts+'weighted.root -w ' + str(lumi*xsec_ttbar['1000']*ttagsf/nev_ttbar['1000']))
commands.append('hadd TBratefilettbar_PSET_'+cuts+'.root temprootfiles/TBratefilettbar700_PSET_'+cuts+'weighted.root temprootfiles/TBratefilettbar1000_PSET_'+cuts+'weighted.root')
commands.append('mv TBratefilettbar700_PSET_'+cuts+'.root TBratefilettbar1000_PSET_'+cuts+'.root temprootfiles/')
commands.append('mv TBratefilettbar_PSET_'+cuts+'.root rootfiles/')
for coup in ['right','left','mixed']:
	sigfiles = sorted(glob.glob('TBratefilesignal'+coup+'*_PSET_'+cuts+'.root'))
	for f in sigfiles:
		mass = f.lstrip('TBratefilesignal'+coup).rstrip('_PSET_'+cuts+'.root')
		xsec_sig = xsec_wpr[mass]
		nev_sig = nev_wpr[mass]
		commands.append('rm ' + f.replace('TBratefilesignal'+coup,'TBratefileweightedsignal'+coup))	 
		commands.append('python HistoWeight.py -i '+f+' -o '+f.replace('TBratefilesignal'+coup,'TBratefileweightedsignal'+coup)+' -w ' + str(lumi*xsec_sig*ttagsf/nev_sig))
		commands.append('mv '+f+' temprootfiles/')
		commands.append('mv '+f.replace('TBratefilesignal'+coup,'TBratefileweightedsignal'+coup)+' rootfiles/')
for s in commands :
    print 'executing ' + s
    subprocess.call( [s], shell=True )







