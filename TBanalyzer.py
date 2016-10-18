#! /usr/bin/env python

###################################################################
##								 ##
## Name: TBanalyzer.py	   			                 ##
## Author: Kevin Nash 						 ##
## Date: 6/5/2012						 ##
## Purpose: This program performs the main analysis.  		 ##
##	    It takes the tagrates created by  	 		 ##
##          TBrate_Maker.py stored in fitdata, and uses 	 ##
##          them to weigh pre b tagged samples to create a 	 ##
##	    QCD background estimate along with the full event    ##
##	    selection to product Mtb inputs to Theta		 ##
##								 ##
###################################################################

import os
import glob
import math
from math import sqrt
#import quickroot
#from quickroot import *
import ROOT 
from ROOT import *

import sys
from DataFormats.FWLite import Events, Handle
from optparse import OptionParser

parser = OptionParser()

parser.add_option('-s', '--set', metavar='F', type='string', action='store',
                  default	=	'data',
                  dest		=	'set',
                  help		=	'data or ttbar')
parser.add_option('-x', '--pileup', metavar='F', type='string', action='store',
                  default	=	'on',
                  dest		=	'pileup',
                  help		=	'If not data do pileup reweighting?')
parser.add_option('-n', '--num', metavar='F', type='string', action='store',
                  default	=	'all',
                  dest		=	'num',
                  help		=	'job number')

parser.add_option('-B', '--bkg', metavar='F', type='string', action='store',
                  default	=	'nominal',
                  dest		=	'bkg',
                  help		=	'cuts to use for rate file')

parser.add_option('-j', '--jobs', metavar='F', type='string', action='store',
                  default	=	'1',
                  dest		=	'jobs',
                  help		=	'number of jobs')


parser.add_option('-y', '--modmass', metavar='F', type='string', action='store',
                  default	=	'nominal',
                  dest		=	'modmass',
                  help		=	'nominal up or down')

parser.add_option('-Y', '--modtb', metavar='F', type='string', action='store',
                  default	=	'nominal',
                  dest		=	'modtb',
                  help		=	'nominal up or down')



parser.add_option('-t', '--tname', metavar='F', type='string', action='store',
                  default	=	'HLT_PFHT800_v2,HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV0p45_v3',
                  dest		=	'tname',
                  help		=	'trigger name')
parser.add_option('-S', '--split', metavar='F', type='string', action='store',
                  default	=	'file',
                  dest		=	'split',
                  help		=	'split by event of file')

#parser.add_option('-t', '--trigger', metavar='F', type='string', action='store',
#                  default	=	'none',
#                  dest		=	'trigger',
#                  help		=	'none, nominal, up, or down')

parser.add_option('-m', '--modulesuffix', metavar='F', type='string', action='store',
                  default	=	'none',
                  dest		=	'modulesuffix',
                  help		=	'ex. PtSmearUp')

parser.add_option('-g', '--grid', metavar='F', type='string', action='store',
                  default	=	'off',
                  dest		=	'grid',
                  help		=	'running on grid off or on')
parser.add_option('-u', '--ptreweight', metavar='F', type='string', action='store',
                  default	=	'none',
                  dest		=	'ptreweight',
                  help		=	'on or off')

parser.add_option('-p', '--pdfweights', metavar='F', type='string', action='store',
                  default	=	'nominal',
                  dest		=	'pdfweights',
                  help		=	'nominal, up, or down')
parser.add_option('-J', '--JES', metavar='F', type='string', action='store',
                  default	=	'nominal',
                  dest		=	'JES',
                  help		=	'nominal, up, or down')
parser.add_option('-R', '--JER', metavar='F', type='string', action='store',
                  default	=	'nominal',
                  dest		=	'JER',
                  help		=	'nominal, up, or down')
parser.add_option('-z', '--pdfset', metavar='F', type='string', action='store',
                  default	=	'',
                  dest		=	'pdfset',
                  help		=	'pdf set')
parser.add_option('--printEvents', metavar='F', action='store_true',
                  default=False,
                  dest='printEvents',
                  help='Print events that pass selection (run:lumi:event)')
parser.add_option('-c', '--cuts', metavar='F', type='string', action='store',
                  default	=	'default',
                  dest		=	'cuts',
                  help		=	'Cuts type (ie default, rate, etc)')
parser.add_option('-b', '--bx', metavar='F', type='string', action='store',
                  default	=	'25ns',
                  dest		=	'bx',
                  help		=	'bunch crossing 50ns or 25ns')
parser.add_option('-T', '--tempop', metavar='F', type='string', action='store',
                  default	=	'norm',
                  dest		=	'tempop',
                  help		=	'temp option')
parser.add_option('--twojet', metavar='F', action='store_true',
                  default=False,
                  dest='twojet',
                  help='twojet')

(options, args) = parser.parse_args()

tname = options.tname.split(',')
tnamestr = ''
for iname in range(0,len(tname)):
	tnamestr+=tname[iname]
	if iname!=len(tname)-1:
		tnamestr+='OR'

if tnamestr=='HLT_PFHT800_v2ORHLT_AK8DiPFJet280_200_TrimMass30_BTagCSV0p45_v3':
	tnameformat='nominal'
elif tnamestr=='':
	tnameformat='none'
else:
	tnameformat=tnamestr


print "Options summary"
print "=================="
for  opt,value in options.__dict__.items():
	#print str(option)+ ": " + str(options[option]) 
	print str(opt) +': '+ str(value)
print "=================="
print ""
di = ""
if options.grid == 'on':
	di = "tardir/"
	sys.path.insert(0, 'tardir/')

gROOT.Macro(di+"rootlogon.C")

import Wprime_Functions	
from Wprime_Functions import *

#Load up cut values based on what selection we want to run 
Cuts = LoadCuts(options.cuts)
bpt = Cuts['bpt']
tpt = Cuts['tpt']
ptmincut = Cuts['ptmincut']
dy = Cuts['dy']
tmass = Cuts['tmass']
nsubjets = Cuts['nsubjets']
tau32 = Cuts['tau32']
minmass = Cuts['minmass']
sjbtag = Cuts['sjbtag']
bmass = Cuts['bmass']
btag = Cuts['btag']
eta1 = Cuts['eta1']
eta2 = Cuts['eta2']
eta3 = Cuts['eta3']



#For large datasets we need to parallelize the processing
jobs=int(options.jobs)
if jobs != 1:
	num=int(options.num)
	jobs=int(options.jobs)
	print "Running over " +str(jobs)+ " jobs"
	print "This will process job " +str(num)
else:
	print "Running over all events"

#This section defines some strings that are used in naming the optput files

mmstr = ""
if options.modmass!="nominal":
	print "using modm uncertainty"
	mmstr = "_modm_"+options.modmass

mmtbstr = ""
if options.modtb!="nominal":
	print "using modtb uncertainty"
	mmtbstr = "_modtb_"+options.modtb


TJ = ''
if options.twojet:
	TJ='_twojetmod'

mod = ''
post = ''
if options.JES!='nominal':
	mod = mod + 'JES_' + options.JES
	post='jes'+options.JES
if options.JER!='nominal':
	mod = mod + 'JER_' + options.JER
	post='jer'+options.JER

TO = ""
if options.tempop!='norm':
	TO=options.tempop

pstr = ""
if options.pdfweights!="nominal":
	print "using pdf uncertainty"
	pstr = "_pdf_"+options.pdfset+"_"+options.pdfweights

pustr = ""
if options.pileup=='off':
	pustr = "_pileup_unweighted"
if options.pileup=='up':
	pustr = "_pileup_up"
if options.pileup=='down':
	pustr = "_pileup_down"
mod = mod+pustr
if mod == '':
	mod = options.modulesuffix


ModFile = ROOT.TFile(di+"ModMassFile_PSET_rate_"+options.cuts +".root")
ModPlot = ModFile.Get("rtmass")


ModtbFile = ROOT.TFile(di+"bkgwQCD"+options.cuts+".root")
ModtbPlot = ModtbFile.Get("bkgw")

run_b_SF = True
#Based on what set we want to analyze, we find all Ntuple root files 

files = Load_Ntuples(options.set, options.bx) 
bkgset='data'
if (options.set.find('ttbar') != -1) or (options.set.find('ST') != -1):
	settype = 'ttbar'
elif (options.set.find('QCD') != -1):
	settype ='ttbar'
	bkgset='QCD'
	run_b_SF = False
else :
	
	settype = options.set.replace('right','').replace('left','').replace('mixed','')



print 'The type of set is ' + settype

if options.set != 'data':
	#Load up scale factors (to be used for MC only)
	TrigFile = TFile(di+"Triggerweight_databtagsmass.root")
	TrigPlot = TrigFile.Get("TriggerWeight_"+tnamestr+"_pre_HLT_PFHT475_v2")


	#UNCOMMENT LATER
	if options.set=="WJETS":
		PileFile = TFile(di+"PileUp_Ratio_ttbar.root")
	else:
		PileFile = TFile(di+"PileUp_Ratio_"+settype+".root")
	if options.pileup == 'up':
		PilePlot = PileFile.Get("Pileup_Ratio_up")
	elif options.pileup == 'down':
		PilePlot = PileFile.Get("Pileup_Ratio_down")
	else:
		PilePlot = PileFile.Get("Pileup_Ratio")
jobiter = 0
splitfiles = []

if jobs != 1 and options.split=="file":
    for ifile in range(1,len(files)+1):
    	if (ifile-1) % jobs == 0:
		jobiter+=1
	count_index = ifile  - (jobiter-1)*jobs
	if count_index==num:
		splitfiles.append(files[ifile-1])

    events = Events(splitfiles)
if options.split=="event" or jobs == 1:	  
	events = Events(files)

#Here we load up handles and labels.
#These are used to grab entries from the Ntuples.
#To see all the current types in an Ntuple use edmDumpEventContent /PathtoNtuple/Ntuple.root

AK8HL = Initlv("jetsAK8",post)


PtnomHandle 	= 	Handle (  "vector<float> "  )
PtnomLabel  	= 	( "jetsAK8" , "jetAK8Pt")

BDiscHandle 	= 	Handle (  "vector<float> "  )
BDiscLabel  	= 	( "jetsAK8" , "jetAK8CSV")

GeneratorHandle 	= 	Handle (  "GenEventInfoProduct")
GeneratorLabel  	= 	( "generator" , "")

puHandle    	= 	Handle("int")
puLabel     	= 	( "eventUserData", "puNtrueInt" )

minmassHandle 	= 	Handle (  "vector<float> "  )
minmassLabel  	= 	( "jetsAK8" , "jetAK8minmass")

nSubjetsHandle 	= 	Handle (  "vector<float> "  )
nSubjetsLabel  	= 	( "jetsAK8" , "jetAK8nSubJets")


pdfHandle 	= 	Handle (  "vector<float> "  )
pdfLabel  	= 	( "weights" , "pdfWeights"+options.pdfset)

partonFlavourHandle 	= 	Handle (  "vector<float> "  )
partonFlavourLabel  	= 	( "jetsAK8" , "jetAK8PartonFlavour")

softDropMassHandle 	= 	Handle (  "vector<float> "  )
softDropMassLabel  	= 	( "jetsAK8" , "jetAK8softDropMass")

tau1Handle 	= 	Handle (  "vector<float> "  )
tau1Label  	= 	( "jetsAK8" , "jetAK8tau1")

tau2Handle 	= 	Handle (  "vector<float> "  )
tau2Label  	= 	( "jetsAK8" , "jetAK8tau2")

tau3Handle 	= 	Handle (  "vector<float> "  )
tau3Label  	= 	( "jetsAK8" , "jetAK8tau3")

topMassHandle 	= 	Handle (  "vector<float> "  )
topMassLabel  	= 	( "jetsAK8" , "jetAK8topMass")

subjetsCSVHandle 	= 	Handle (  "vector<float> "  )
subjetsCSVLabel  	= 	( "subjetsCmsTopTag" , "subjetCmsTopTagCSV")

subjets0indexHandle 	= 	Handle (  "vector<float> "  )
subjets0indexLabel  	= 	( "jetsAK8" , "jetAK8topSubjetIndex0")

subjets1indexHandle 	= 	Handle (  "vector<float> "  )
subjets1indexLabel  	= 	( "jetsAK8" , "jetAK8topSubjetIndex1")

subjets2indexHandle 	= 	Handle (  "vector<float> "  )
subjets2indexLabel  	= 	( "jetsAK8" , "jetAK8topSubjetIndex2")

subjets3indexHandle 	= 	Handle (  "vector<float> "  )
subjets3indexLabel  	= 	( "jetsAK8" , "jetAK8topSubjetIndex3")





softDropMassuncorrHandle 	= 	Handle (  "vector<float> "  )
softDropMassuncorrLabel  	= 	( "jetsAK8" , "jetAK8softDropMassuncorr")

vsubjets0indexHandle 	= 	Handle (  "vector<float> "  )
vsubjets0indexLabel  	= 	( "jetsAK8" , "jetAK8vSubjetIndex0")

vsubjets1indexHandle 	= 	Handle (  "vector<float> "  )
vsubjets1indexLabel  	= 	( "jetsAK8" , "jetAK8vSubjetIndex1")

subjetsAK8CSVHandle 	= 	Handle (  "vector<float> "  )
subjetsAK8CSVLabel  	= 	( "subjetsAK8" , "subjetAK8CSV")



#---------------------------------------------------------------------------------------------------------------------#

if jobs != 1:
	f = TFile( "TBanalyzer"+options.set+"_Trigger_"+tnameformat+"_"+mod+pstr+mmstr+mmtbstr+TO+"_job"+options.num+"of"+options.jobs+TJ+"_PSET_"+options.cuts+".root", "recreate" )
else:
	f = TFile( "TBanalyzer"+options.set+"_Trigger_"+tnameformat+"_"+mod +pstr+mmstr+mmtbstr+TO+TJ+"_PSET_"+options.cuts+".root", "recreate" )


if options.bkg=='nominal':
	ratestr = options.cuts
else:
	ratestr = options.bkg

#Load up the average b-tagging rates -- Takes parameters from text file and makes a function
BTR = BTR_Init('Bifpoly','rate_'+ratestr,di,bkgset)
BTR_err = BTR_Init('Bifpoly_err','rate_'+ratestr,di,bkgset)
fittitles = ["pol0","pol2","pol3","FIT","Bifpoly","expofit"]
fits = []
for fittitle in fittitles:
	fits.append(BTR_Init(fittitle,'rate_'+ratestr,di,bkgset))

print "Creating histograms"

#Define Histograms

TagFile1 = TFile(di+"Tagrate2D"+bkgset+'rate_'+ratestr+".root")
TagPlot2de1= TagFile1.Get("tagrateeta1")
TagPlot2de2= TagFile1.Get("tagrateeta2")
TagPlot2de3= TagFile1.Get("tagrateeta3")

f.cd()
#---------------------------------------------------------------------------------------------------------------------#
Mtb		= TH1F("Mtb",		"mass of tb",     	  	      	140, 500, 4000 )
Mb		= TH1F("Mb",		"mass of b",     	  	      	200, 0, 400 )

Nevents		= TH1F("Nevents",	"mass of tb",     	  	          5, 0., 5. )

Nsdsj		= TH1F("Nsdsj",	"mass of tb",     	  	          5, 0., 5. )


Mtbptup		= TH1F("Mtbptup",	"mass of tb ttbar pt reweighting up",	140, 500, 4000 )
Mtbptdown	= TH1F("Mtbptdown",     "mass of tb ttbar pt reweighting up",   140, 500, 4000 )

MtbBup		= TH1F("MtbBup",	"mass of tb B tag SF up",     	  	140, 500, 4000 )
MtbBdown	= TH1F("MtbBdown",	"mass of tb B tag SF up",     	  	140, 500, 4000 )

MtbTup		= TH1F("MtbTup",	"mass of tb B tag SF up",     	  	140, 500, 4000 )
MtbTdown	= TH1F("MtbTdown",	"mass of tb B tag SF up",     	  	140, 500, 4000 )


Mtbtrigup		= TH1F("Mtbtrigup",	"mass of tb B tag SF up",     	  	140, 500, 4000 )
Mtbtrigdown	= TH1F("Mtbtrigdown",	"mass of tb B tag SF up",     	  	140, 500, 4000 )


QCDbkg		= TH1F("QCDbkg",	"QCD background estimate",		140, 500, 4000 )
QCDbkgh		= TH1F("QCDbkgh",	"QCD background estimate up error",	140, 500, 4000 )
QCDbkgl		= TH1F("QCDbkgl",	"QCD background estimate down error",	140, 500, 4000 )
QCDbkg2D	= TH1F("QCDbkg2D",	"QCD background estimate 2d error",	140, 500, 4000 )

pflavpre = TH1F("pflavpre",	"pflavpre",		30, -0.5, 29.5 )
pflavpost = TH1F("pflavpost",	"pflavpost",		30, -0.5, 29.5 )

pflavpre.Sumw2()
pflavpost.Sumw2()

Mtb.Sumw2()

Mtbtrigup.Sumw2()
Mtbtrigdown.Sumw2()


Mtbptup.Sumw2()
Mtbptdown.Sumw2()

MtbBup.Sumw2()
MtbBdown.Sumw2()

QCDbkg.Sumw2()
QCDbkgh.Sumw2()
QCDbkgl.Sumw2()

QCDbkg_ARR = []



for ihist in range(0,len(fittitles)):
	QCDbkg_ARR.append(TH1F("QCDbkg"+str(fittitles[ihist]),     "mass W' in b+1 pt est etabin",     	  	      140, 500, 4000 ))
	QCDbkg_ARR[ihist].Sumw2()

#---------------------------------------------------------------------------------------------------------------------#

# loop over events
#---------------------------------------------------------------------------------------------------------------------#

count = 0

print "Start looping"
#initialize the ttree variables
tree_vars = {"bpt":array('d',[0.]),"bmass":array('d',[0.]),"btag":array('d',[0.]),"tpt":array('d',[0.]),"tmass":array('d',[0.]),"nsubjets":array('d',[0.]),"sjbtag":array('d',[0.]),"run":array('d',[0.]),"lumi":array('d',[0.]),"event":array('d',[0.]),"weight":array('d',[0.])}
Tree = Make_Trees(tree_vars)

goodEvents = []
#totevents = events.size()
#print str(totevents)  +  ' Events total'

#usegenweight = False
#if options.set == "QCDFLAT7000":
#	usegenweight = True
#	print "Using gen weight"
nfail = 0 
for event in events:
    weightSFb = 1.0
    errorSFb = 0.0
    count	= 	count + 1
    
    m = 0
    t = 0
    if count > 100000:
	break

    if count % 100000 == 0 :
      print  '--------- Processing Event ' + str(count) #+'   -- percent complete ' + str(100*count/totevents) + '% -- '

    #Here we split up event processing based on number of jobs 
    #This is set up to have jobs range from 1 to the total number of jobs (ie dont start at job 0)

    if jobs != 1 and options.split=="event":
    	if (count - 1) % jobs == 0:
		jobiter+=1
	count_index = count - (jobiter-1)*jobs
	if count_index!=num:
		continue 

 #   if usegenweight:
#
#		try:
#			event.getByLabel (GeneratorLabel, GeneratorHandle)
 #   			gen 		= 	GeneratorHandle.product()
#			Nevents.Fill(0.,gen.weightProduct())
#		except:
#			continue 

    try:
   	AK8LV = Makelv(AK8HL,event)
    except:
	nfail+=1
	continue 
    # We load up the relevant handles and labels and create collections

    if len(AK8LV)==0:
	continue

    tindex,bindex = Hemispherize(AK8LV,AK8LV)


    bJetsh1 = []
    bJetsh0  =  []
    topJetsh1 = []
    topJetsh0  = []
    for i in range(0,len(bindex[1])):
   	bJetsh1.append(AK8LV[bindex[1][i]])
    for i in range(0,len(bindex[0])):
    	bJetsh0.append(AK8LV[bindex[0][i]])
    for i in range(0,len(tindex[1])):
    	topJetsh1.append(AK8LV[tindex[1][i]])
    for i in range(0,len(tindex[0])):
    	topJetsh0.append(AK8LV[tindex[0][i]])
    
    bjh0 = 0
    bjh1 = 0

    tjh0 = 0
    tjh1 = 0

    #Require 1 pt>150 jet in each hemisphere (top jets already have the 150GeV requirement) 
    for bjet in bJetsh0:
	if ptmincut[0] < bjet.Perp() < ptmincut[1]:
		bjh0+=1
    for tjet in topJetsh0:
	if ptmincut[0] < tjet.Perp()  < ptmincut[1]:
		tjh0+=1

    for bjet in bJetsh1:
	if ptmincut[0] < bjet.Perp()  < ptmincut[1]:
		bjh1+=1

    for tjet in topJetsh1:
	if ptmincut[0] < tjet.Perp()  < ptmincut[1]:
		tjh1+=1
    if options.twojet:
   	 njets11b0 	= 	((tjh1 == 1) and (bjh0 == 1))
   	 njets11b1 	= 	((tjh0 == 1) and (bjh1 == 1))
    else:
   	 njets11b0 	= 	((tjh1 >= 1) and (bjh0 >= 1))
   	 njets11b1 	= 	((tjh0 >= 1) and (bjh1 >= 1))
    for hemis in ['hemis0','hemis1']:
    	if hemis == 'hemis0'   :
		if not njets11b0 :
			continue 
		#The Ntuple entries are ordered in pt, so [0] is the highest pt entry
		#We are calling a candidate b jet (highest pt jet in hemisphere0)  

		tindexval = tindex[1][0]
		bindexval = bindex[0][0]

		bjet = bJetsh0[0]
		tjet = topJetsh1[0]


    	if hemis == 'hemis1'  :
		if not njets11b1:
			continue 

		tindexval = tindex[0][0]
		bindexval = bindex[1][0]

		bjet = bJetsh1[0]
		tjet = topJetsh0[0]
	
	if abs(bjet.Eta()) > 2.4 or abs(tjet.Eta()) > 2.4:
		continue


	
    	weight=1.0
     	weightSFptup=1.0
     	weightSFptdown=1.0

    	bpt_cut = bpt[0]<bjet.Perp()<bpt[1]
    	tpt_cut = tpt[0]<tjet.Perp()<tpt[1]



    	dy_cut = dy[0]<=abs(tjet.Rapidity()-bjet.Rapidity())<dy[1]
    	if bpt_cut and tpt_cut and dy_cut: 

    		if options.pdfweights != "nominal" :

            		event.getByLabel( pdfLabel, pdfHandle )
            		pdfs = pdfHandle.product()
			if options.tempop=="norm":
				iweight = PDF_Lookup( pdfs , options.pdfweights)
			if options.tempop=="ave":
				iweight = PDF_LookupAVE( pdfs , options.pdfweights)
			if options.tempop=="max":
				iweight = PDF_LookupMAX( pdfs , options.pdfweights)
			if options.tempop=="std":
				iweight = PDF_LookupSTD( pdfs , options.pdfweights)

            		weight *= iweight


		weightSFb = 1.0
		weightSFbdown = 1.0
		weightSFbup = 1.0

		weightSFt = 1.0
		weightSFtdown = 1.0
		weightSFtup = 1.0


		if options.set!="data":

		
			event.getByLabel (puLabel, puHandle)
    			PileUp 		= 	puHandle.product()
               		bin1 = PilePlot.FindBin(float(PileUp[0])) 

			if options.pileup != 'off':
				weight *= PilePlot.GetBinContent(bin1)


			if run_b_SF :
				#btagging scale factor reweighting done here
				SFB = SFB_Lookup( bjet.Perp() )
				weightSFb = SFB[0]
				weightSFbdown = SFB[1]
				weightSFbup = SFB[2]
		
			if options.cuts=="default" or options.cuts=="sideband3" :
				#btagging scale factor reweighting done here
				SFT = SFT_Lookup( tjet.Perp() )
				weightSFt = SFT[0]
				weightSFtdown = SFT[1]
				weightSFtup = SFT[2]

        	event.getByLabel (softDropMassLabel, softDropMassHandle)
        	topJetMass 	= 	softDropMassHandle.product()


        	event.getByLabel (softDropMassuncorrLabel, softDropMassuncorrHandle)
        	topJetMassuncorr 	= 	softDropMassuncorrHandle.product()

        	event.getByLabel (PtnomLabel, PtnomHandle)
        	Ptnom 	= 	PtnomHandle.product()
		
		tcorr = 1.0
		bcorr = 1.0
		if options.JER!='nominal':
			tcorr=tjet.Perp()/Ptnom[tindexval]
			bcorr=bjet.Perp()/Ptnom[bindexval]
		if options.JES!='nominal':
			bcorr=bjet.Perp()/Ptnom[bindexval]

		#####TEMPORARY#####
		[tmassc,bmassc] = [topJetMassuncorr[tindexval]*tcorr,topJetMass[bindexval]*bcorr]
	

		tmass_cut = tmass[0]<tmassc<tmass[1]

		if tmass_cut :

        		event.getByLabel ( nSubjetsLabel , nSubjetsHandle )
    			nSubjets 		= 	nSubjetsHandle.product()
        		event.getByLabel (minmassLabel, minmassHandle)
    			topJetminmass 	= 	minmassHandle.product()

			minmass_cut = minmass[0]<=topJetminmass[tindexval]<minmass[1]
			nsubjets_cut = nsubjets[0]<=nSubjets[tindexval]<nsubjets[1]
		
	 		if minmass_cut and nsubjets_cut:

   				event.getByLabel (BDiscLabel, BDiscHandle)
    				bJetBDisc 	= 	BDiscHandle.product()

                		btag_cut = btag[0]<bJetBDisc[bindexval]<=btag[1]

				ht = tjet.Perp() + bjet.Perp()



				weighttrigup=1.0
				weighttrigdown=1.0
				if options.tname != 'none' and options.set!='data':


					TRW = Trigger_Lookup( ht , TrigPlot )[0]
					TRWup = Trigger_Lookup( ht , TrigPlot )[1]
					TRWdown = Trigger_Lookup( ht , TrigPlot )[2]

					weighttrigup=weight*TRWup
					weighttrigdown=weight*TRWdown
					weight*=TRW

				if options.ptreweight == "on":
					#ttbar pt reweighting done here
					event.getByLabel( GenLabel, GenHandle )
					GenParticles = GenHandle.product()
					PTW = PTW_Lookup( GenParticles )
					weight*=PTW
     					weightSFptup=max(0.0,weight*(2*PTW-1))
     					weightSFptdown=weight
				#Done differently than TBkinematics
				weightSFtup=weight*weightSFb*weightSFtup
				weightSFtdown=weight*weightSFb*weightSFtdown
				weight*=weightSFt

				weightb = weight*weightSFb
				weightSFptup*=weightSFb
				weightSFptdown*=weightSFb
				weightSFbup = weight*(weightSFbup)
				weightSFbdown = weight*(weightSFbdown)
				weighttrigup*=weightSFb*weightSFt
				weighttrigdown*=weightSFb*weightSFt

    				#event.getByLabel (subjets0indexLabel, subjets0indexHandle)
    				#subjets0index 		= 	subjets0indexHandle.product() 

    				#event.getByLabel (subjets1indexLabel, subjets1indexHandle)
    				#subjets1index 		= 	subjets1indexHandle.product() 

    				#event.getByLabel (subjets2indexLabel, subjets2indexHandle)
    				#subjets2index 		= 	subjets2indexHandle.product() 

    				#event.getByLabel (subjets3indexLabel, subjets3indexHandle)
    				#subjets3index 		= 	subjets3indexHandle.product()

    				#event.getByLabel (subjetsCSVLabel, subjetsCSVHandle)
    				#subjetsCSV 		= 	subjetsCSVHandle.product()  



				#SJ_csvs = [subjets0index,subjets1index,subjets2index,subjets3index]
			

    				event.getByLabel (vsubjets0indexLabel,vsubjets0indexHandle )
    				vsubjets0index 		= 	vsubjets0indexHandle.product() 

    				event.getByLabel (vsubjets1indexLabel,vsubjets1indexHandle )
    				vsubjets1index 		= 	vsubjets1indexHandle.product() 

    				event.getByLabel (subjetsAK8CSVLabel,subjetsAK8CSVHandle )
    				subjetsAK8CSV		= 	subjetsAK8CSVHandle.product() 

				Nsdsj.Fill(len(subjetsAK8CSV))

				if len(subjetsAK8CSV)==0:
					continue
				if len(subjetsAK8CSV)<2:
					SJ_csvvals = [subjetsAK8CSV[int(vsubjets0index[tindexval])]]
				else:
    					SJ_csvvals = [subjetsAK8CSV[int(vsubjets0index[tindexval])],subjetsAK8CSV[int(vsubjets1index[tindexval])]]
				SJ_csvmax = max(SJ_csvvals)





				#SJ_csvvals = []
				#for icsv in range(0,int(nSubjets[tindexval])):
			#		if int(SJ_csvs[icsv][tindexval])!=-1:
			#			SJ_csvvals.append(subjetsCSV[int(SJ_csvs[icsv][tindexval])])
			#		else:
			#			SJ_csvvals.append(0.)
			#	SJ_csvmax = max(SJ_csvvals)

				sjbtag_cut = sjbtag[0]<SJ_csvmax<=sjbtag[1]



    				event.getByLabel (tau1Label, tau1Handle)
    				tau1 		= 	tau1Handle.product()  


    				event.getByLabel (tau2Label, tau2Handle)
    				tau2 		= 	tau2Handle.product()  


    				event.getByLabel (tau3Label, tau3Handle)
    				tau3 		= 	tau3Handle.product()  


				if tau2[tindexval]>0.:
					tau32_cut =  tau32[0]<=tau3[tindexval]/tau2[tindexval]<tau32[1]
				else:
					tau32_cut =  tau32[0]<=1.0<tau32[1]


				if sjbtag_cut :

					if tau32_cut:

						bmass_cut = bmass[0]<=bmassc<bmass[1]
						if bmass_cut:

							eta_regions = [eta1,eta2,eta3]
							BTRweight = bkg_weight(bjet,BTR,eta_regions)
							BTRweightsigsq = bkg_weight(bjet,BTR_err,eta_regions)

							BTRweighterrup = BTRweight+sqrt(BTRweightsigsq)
							BTRweighterrdown = BTRweight-sqrt(BTRweightsigsq)

							modm = bmassc
							if options.modmass=='nominal':
                						massw = ModPlot.Interpolate(modm)
							if options.modmass=='up':
                						massw = 1 + 0.5*(ModPlot.Interpolate(modm)-1)
							if options.modmass=='down':
                						massw = max(0.0,1 + 1.5*(ModPlot.Interpolate(modm)-1))


							modmtb = (tjet+bjet).M()
							tbbin = ModtbPlot.FindBin(modmtb)
							tberr = ModtbPlot.GetBinError(tbbin)
							if options.modtb=='nominal':
                						tbw = ModtbPlot.GetBinContent(tbbin)
							if options.modtb=='up':
                						tbw = 1.0

							if options.modtb=='down':
                						tbw = max(0.0,1 + 2.0*(ModtbPlot.GetBinContent(tbbin)-1))
							if bkgset!='QCD':
								massw*=tbw
					
							eta1_cut = eta1[0]<=abs(bjet.Eta())<eta1[1]
							eta2_cut = eta2[0]<=abs(bjet.Eta())<eta2[1]
							eta3_cut = eta3[0]<=abs(bjet.Eta())<eta3[1]
							if options.set!="data":
    								event.getByLabel (partonFlavourLabel, partonFlavourHandle)
    								partonFlavour 		= 	partonFlavourHandle.product()  
				
							if not btag_cut:
								if options.set!="data":
									pflavpre.Fill(abs(partonFlavour[bindexval]))
								if (eta1_cut) :
									xbin = TagPlot2de1.GetXaxis().FindBin(bjet.Perp())
									ybin = TagPlot2de1.GetYaxis().FindBin((tjet+bjet).M())
									tagrate2d = TagPlot2de1.GetBinContent(xbin,ybin)
									QCDbkg2D.Fill((tjet+bjet).M(),tagrate2d*weight*massw)
			
								if (eta2_cut):
									xbin = TagPlot2de2.GetXaxis().FindBin(bjet.Perp())
									ybin = TagPlot2de2.GetYaxis().FindBin((tjet+bjet).M())
									tagrate2d = TagPlot2de2.GetBinContent(xbin,ybin)
									QCDbkg2D.Fill((tjet+bjet).M(),tagrate2d*weight*massw)

								if (eta3_cut):
	
									xbin = TagPlot2de3.GetXaxis().FindBin(bjet.Perp())
									ybin = TagPlot2de3.GetYaxis().FindBin((tjet+bjet).M())
									tagrate2d = TagPlot2de3.GetBinContent(xbin,ybin)
									QCDbkg2D.Fill((tjet+bjet).M(),tagrate2d*weight*massw)

							
								for ifit in range(0,len(fittitles)):
									tempweight = bkg_weight(bjet,fits[ifit],eta_regions)
									QCDbkg_ARR[ifit].Fill((tjet+bjet).M(),tempweight*weight*massw) 
								
								QCDbkg.Fill((tjet+bjet).M(),BTRweight*weight*massw)
								QCDbkgh.Fill((tjet+bjet).M(),BTRweighterrup*weight*massw)
								QCDbkgl.Fill((tjet+bjet).M(),BTRweighterrdown*weight*massw) 
	
        				        	if btag_cut:

								if options.set!="data":
								
									pflavpost.Fill(abs(partonFlavour[bindexval]))


                                      				goodEvents.append( [ event.object().id().run(), event.object().id().luminosityBlock(), event.object().id().event() ] )
								Mtb.Fill((tjet+bjet).M(),weightb) 



								MtbBup.Fill((tjet+bjet).M(),weightSFbup) 
								MtbBdown.Fill((tjet+bjet).M(),weightSFbdown) 


								MtbTup.Fill((tjet+bjet).M(),weightSFtup) 
								MtbTdown.Fill((tjet+bjet).M(),weightSFtdown) 


								Mtbptup.Fill((tjet+bjet).M(),weightSFptup) 
								Mtbptdown.Fill((tjet+bjet).M(),weightSFptdown) 

								Mtbtrigup.Fill((tjet+bjet).M(),weighttrigup)
								Mtbtrigdown.Fill((tjet+bjet).M(),weighttrigdown)
								#print weightb
								#print weightSFb
								#print weightSFt
								#print PilePlot.GetBinContent(bin1)
								#print float(PileUp[0])
								#print 
								temp_variables = {"bpt":bjet.Perp(),"bmass":bmassc,"btag":bJetBDisc[bindexval],"tpt":tjet.Perp(),"tmass":tmassc,"nsubjets":nSubjets[tindexval],"sjbtag":SJ_csvmax,"run":event.object().id().run(),"lumi":event.object().id().luminosityBlock(),"event":event.object().id().event(),"weight":weightb }		


								for tv in tree_vars.keys():
									tree_vars[tv][0] = temp_variables[tv]
								Tree.Fill()


f.cd()
f.Write()
f.Close()

print "number of events: " + str(count)


if options.printEvents:
    Outf1   =   open("DataEvents"+options.num+".txt", "w")
    sys.stdout = Outf1
    for goodEvent in goodEvents :
        print '{0:12.0f}:{1:12.0f}:{2:12.0f}'.format(
            goodEvent[0], goodEvent[1], goodEvent[2]
        )
