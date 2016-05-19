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
                  default	=	'default',
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

filestrs = ['none_','JES_up_','JES_down_','JER_up_','JER_down_','_pileup_up_','_pileup_down_','none_pdf_NNPDF_up_','none_pdf_NNPDF_down_']

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
commands.append('mv TBanalyzerdata_*.root rootfiles/')


for s in commands :
    print 'executing ' + s
    subprocess.call( [s], shell=True )

commands = []

for f in filestrs:
	f = f.replace('NNPDF','')
	commands.append('rm rootfiles/TBanalyzerttbar_'+f+'PSET_'+cuts+'weighted.root') #removes old file with same name in /rootfiles/
	commands.append('python HistoWeight.py -i TBanalyzerttbar_Trigger_nominal_'+f+'PSET_'+cuts+'.root -o rootfiles/TBanalyzerttbar_'+f+'PSET_'+cuts+'weighted.root -w ' + str(lumi*xsec_ttbar['PH']/nev_ttbar['PH']))
	commands.append('mv TBanalyzerttbar_Trigger_nominal_'+f+'PSET_'+cuts+'.root temprootfiles/')


for scale in ['scaleup','scaledown']:
	commands.append('rm rootfiles/TBanalyzerttbar'+scale+'_Trigger_nominal_none_PSET_'+cuts+'weighted.root') #removes old file with same name in /rootfiles/
	commands.append('python HistoWeight.py -i TBanalyzerttbar'+scale+'_Trigger_nominal_none_PSET_'+cuts+'.root -o rootfiles/TBanalyzerttbar'+scale+'_Trigger_nominal_none_PSET_'+cuts+'weighted.root -w ' + str(lumi*xsec_ttbar['PH'+scale]/nev_ttbar['PH'+scale]))
	commands.append('mv TBanalyzerttbar'+scale+'_Trigger_nominal_none_PSET_'+cuts+'.root temprootfiles/')

STfiles = sorted(glob.glob('TBanalyzerS*_Trigger_nominal_none_PSET_'+cuts+'.root'))
for f in STfiles:
	
	chan = f.replace('TBanalyzerST','').replace('_Trigger_nominal_none_PSET_'+cuts+'.root','')	
	if chan=='':
		continue 
	print "Single top " + chan + " channel"

	xsec_ST = xsec_st[chan]
	nev_ST = nev_st[chan]
	commands.append('rm ' + f.replace('TBanalyzerST','TBanalyzerweightedST'))	 
	commands.append('python HistoWeight.py -i '+f+' -o '+f.replace('TBanalyzerST','TBanalyzerweightedST')+' -w ' + str(lumi*xsec_ST/nev_ST))
	commands.append('mv '+f+' temprootfiles/')
commands.append('rm TBanalyzerST_Trigger_nominal_none_PSET_'+cuts+'.root')
commands.append('hadd TBanalyzerST_Trigger_nominal_none_PSET_'+cuts+'.root TBanalyzerweightedST*_Trigger_nominal_none_PSET_'+cuts+'.root')
commands.append('mv TBanalyzerST_Trigger_nominal_none_PSET_'+cuts+'.root rootfiles/')
commands.append('mv TBanalyzerweightedST*_Trigger_nominal_none_PSET_'+cuts+'.root temprootfiles/')



for mod in ['','_modm_down','_modm_up']:	
	qcdfiles = sorted(glob.glob('TBanalyzerQCDHT*_Trigger_nominal_none'+mod+'_PSET_'+cuts+'.root'))
	for f in qcdfiles:
		pt = f.lstrip('TBanalyzerQCDHT').rstrip('_Trigger_nominal_none'+mod+'_PSET_'+cuts+'.root')
		print "QCD HT " + pt
		xsec_QCD = xsec_qcd['HT'+pt]
		nev_QCD = nev_qcd['HT'+pt]
		commands.append('rm ' + f.replace('TBanalyzerQCDHT','TBanalyzerweightedQCDHT'))	 
		commands.append('python HistoWeight.py -i '+f+' -o '+f.replace('TBanalyzerQCDHT','TBanalyzerweightedQCDHT')+' -w ' + str(lumi*xsec_QCD/nev_QCD))
		commands.append('mv '+f+' temprootfiles/')
	commands.append('hadd TBanalyzerQCD_Trigger_nominal_none'+mod+'_PSET_'+cuts+'.root TBanalyzerweightedQCDHT*_Trigger_nominal_none'+mod+'_PSET_'+cuts+'.root')
	commands.append('mv TBanalyzerQCD_Trigger_nominal_none'+mod+'_PSET_'+cuts+'.root rootfiles/')
	commands.append('mv TBanalyzerweightedQCDHT*_Trigger_nominal_none'+mod+'_PSET_'+cuts+'.root temprootfiles/')



for g in filestrs:
    for coup in ['right','left','mixed']:
	sigfiles = sorted(glob.glob('TBanalyzersignal'+coup+'*_'+g+'PSET_'+cuts+'.root'))
	for f in sigfiles:
		mass = f.replace('TBanalyzersignal'+coup,'')[:4]
		
		if coup =='right':
			xsec_sig = xsec_wpr[mass]*kfac
			nev_sig = nev_wpr[mass]
		if coup =='left':
			xsec_sig = xsec_wpl[mass]*kfac
			nev_sig = nev_wpl[mass]
		if coup =='mixed':
			xsec_sig = xsec_wplr[mass]*kfac
			nev_sig = nev_wplr[mass]
		commands.append('rm ' + f.replace('TBanalyzersignal'+coup,'TBanalyzerweightedsignal'+coup))	 
		commands.append('python HistoWeight.py -i '+f+' -o '+f.replace('TBanalyzersignal'+coup,'TBanalyzerweightedsignal'+coup)+' -w ' + str(lumi*xsec_sig/nev_sig))
		commands.append('mv '+f+' temprootfiles/')
		commands.append('mv '+f.replace('TBanalyzersignal'+coup,'TBanalyzerweightedsignal'+coup)+' rootfiles/')




for s in commands :
    print 'executing ' + s
    subprocess.call( [s], shell=True )







