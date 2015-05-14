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
parser.add_option('-j', '--jobs', metavar='F', type='string', action='store',
                  default	=	'1',
                  dest		=	'jobs',
                  help		=	'number of jobs')
parser.add_option('-t', '--trigger', metavar='F', type='string', action='store',
                  default	=	'none',
                  dest		=	'trigger',
                  help		=	'none, nominal, up, or down')

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
parser.add_option('-z', '--pdfset', metavar='F', type='string', action='store',
                  default	=	'cteq66',
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




(options, args) = parser.parse_args()

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



import Wprime_Functions	
from Wprime_Functions import *

#Load up cut values based on what selection we want to run 
Cuts = LoadCuts(options.cuts)
bpt = Cuts['bpt']
tpt = Cuts['tpt']
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
mod = "ttbsmAna"
if options.modulesuffix != "none" :
	mod = mod + options.modulesuffix

pstr = ""
if options.pdfweights!="nominal":
	print "using pdf uncertainty"
	pstr = "_pdf_"+options.pdfset+"_"+options.pdfweights

pustr = ""
if options.pileup=='off':
	pustr = "_pileup_unweighted"


run_b_SF = True
#Based on what set we want to analyze, we find all Ntuple root files 

files = Load_Ntuples(options.set)

if (options.set.find('ttbar') != -1) or (options.set.find('singletop') != -1):
	settype = 'ttbar'
elif (options.set.find('QCD') != -1):
	settype ='QCD'
	run_b_SF = False
else :
	
	settype = options.set.replace('right','').replace('left','').replace('mixed','')

print 'The type of set is ' + settype

if options.set != 'data':
	#Load up scale factors (to be used for MC only)

	TrigFile = TFile(di+"TRIG_EFFICWPHTdata_dijet8TeV.root")
	TrigPlot = TrigFile.Get("r11")

	PileFile = TFile(di+"PileUp_Ratio_"+settype+".root")
	PilePlot = PileFile.Get("Pileup_Ratio")

# We select all the events:    
events = Events (files)

#Here we load up handles and labels.
#These are used to grab entries from the Ntuples.
#To see all the current types in an Ntuple use edmDumpEventContent /PathtoNtuple/Ntuple.root

pdfHandle 		= 	Handle("std::vector<double>")
pdfLabel 		= 	( "pdfWeights",options.pdfset)

puHandle    	= 	Handle("int")
puLabel     	= 	( mod, "npvRealTrue" )

npvHandle    	= 	Handle("unsigned int")
npvLabel     	= 	( mod, "npv" )

CA8Handle     	= 	Handle( "vector<ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > >")
CA8Label      	= 	( "nsub", "CA8P4" )

TopTau2Handle       = 	Handle( "std::vector<double>" )
TopTau2Label    	= 	( "nsub" , "Tau2")

TopTau3Handle       = 	Handle( "std::vector<double>" )
TopTau3Label    	= 	( "nsub" , "Tau3")

TopRcnHandle       = 	Handle( "std::vector<double>" )
TopRcnLabel    	= 	( "nsub" , "Rcn")



#Load all hemisphere 0 objects
#---------------------------------------------------------------------------------------------------------------------#

hemis0topHandle     	= 	Handle( "vector<ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > >")
hemis0topLabel      	= 	( mod, "topTagP4Hemis0" )

hemis0NSubJetsHandle 	= 	Handle (  "vector<double> "  )
hemis0NSubJetsLabel  	= 	( mod , "topTagNSubjetsHemis0")

hemis0TopBDiscsj0CSVHandle       = 	Handle( "std::vector<double>" )
hemis0TopBDiscsj0CSVLabel    	= 	( mod, "topTagsj0BDiscCSVHemis0")

hemis0TopBDiscsj1CSVHandle       = 	Handle( "std::vector<double>" )
hemis0TopBDiscsj1CSVLabel    	= 	( mod, "topTagsj1BDiscCSVHemis0")

hemis0TopBDiscsj2CSVHandle       = 	Handle( "std::vector<double>" )
hemis0TopBDiscsj2CSVLabel    	= 	( mod, "topTagsj2BDiscCSVHemis0")

hemis0TopBDiscsj3CSVHandle       = 	Handle( "std::vector<double>" )
hemis0TopBDiscsj3CSVLabel    	= 	( mod, "topTagsj3BDiscCSVHemis0")

hemis0MinMassHandle     = 	Handle( "std::vector<double>" )
hemis0MinMassLabel  	= 	( mod, "topTagMinMassHemis0" )

hemis0TopMassHandle     = 	Handle( "std::vector<double>" )
hemis0TopMassLabel  	= 	( mod, "topTagTopMassHemis0" )

hemis0prCA8Handle  	= 	Handle ("vector<ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > >")
hemis0prCA8Label    	=	( mod, "wTagP4Hemis0" )

hemis0BDiscHandle       = 	Handle( "std::vector<double>" )
hemis0BDiscLabel    	= 	( mod , "wTagBDiscHemis0CSV")


#---------------------------------------------------------------------------------------------------------------------#


# Load all Hemisphere 1 objects
#---------------------------------------------------------------------------------------------------------------------#

hemis1topHandle     	= 	Handle( "vector<ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > >")
hemis1topLabel      	= 	( mod, "topTagP4Hemis1" )

hemis1NSubJetsHandle 	= 	Handle (  "vector<double> "  )
hemis1NSubJetsLabel  	= 	( mod , "topTagNSubjetsHemis1")

hemis1MinMassHandle     = 	Handle( "std::vector<double>" )
hemis1MinMassLabel  	= 	( mod, "topTagMinMassHemis1" )

hemis1TopMassHandle     = 	Handle( "std::vector<double>" )
hemis1TopMassLabel  	= 	( mod, "topTagTopMassHemis1" )

hemis1prCA8Handle     	= 	Handle( "vector<ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > >")
hemis1prCA8Label      	= 	( mod, "wTagP4Hemis1" )

hemis1TopBDiscsj0CSVHandle       = 	Handle( "std::vector<double>" )
hemis1TopBDiscsj0CSVLabel    	= 	( mod, "topTagsj0BDiscCSVHemis1")

hemis1TopBDiscsj1CSVHandle       = 	Handle( "std::vector<double>" )
hemis1TopBDiscsj1CSVLabel    	= 	( mod, "topTagsj1BDiscCSVHemis1")

hemis1TopBDiscsj2CSVHandle       = 	Handle( "std::vector<double>" )
hemis1TopBDiscsj2CSVLabel    	= 	( mod, "topTagsj2BDiscCSVHemis1")

hemis1TopBDiscsj3CSVHandle       = 	Handle( "std::vector<double>" )
hemis1TopBDiscsj3CSVLabel    	= 	( mod, "topTagsj3BDiscCSVHemis1")

hemis1BDiscHandle       = 	Handle( "std::vector<double>" )
hemis1BDiscLabel    	= 	( mod , "wTagBDiscHemis1CSV")

GenHandle = Handle( "vector<reco::GenParticle>" )
GenLabel = ( "prunedGenParticles", "" )

#---------------------------------------------------------------------------------------------------------------------#

if jobs != 1:
	f = TFile( "TBanalyzer"+options.set+"_Trigger_"+options.trigger+"_"+options.modulesuffix +pustr+pstr+"_job"+options.num+"of"+options.jobs+"_PSET_"+options.cuts+".root", "recreate" )
else:
	f = TFile( "TBanalyzer"+options.set+"_Trigger_"+options.trigger+"_"+options.modulesuffix +pustr+pstr+"_PSET_"+options.cuts+".root", "recreate" )

#Load up the average b-tagging rates -- Takes parameters from text file and makes a function
BTR = BTR_Init('Bifpoly','rate_'+options.cuts,di)
BTR_err = BTR_Init('Bifpoly_err','rate_'+options.cuts,di)
fittitles = ["pol0","pol2","pol3","FIT","Bifpoly","expofit"]
fits = []
for fittitle in fittitles:
	fits.append(BTR_Init(fittitle,'rate_'+options.cuts,di))

print "Creating histograms"

#Define Histograms

TagFile1 = TFile(di+"Tagrate2D.root")
TagPlot2de1= TagFile1.Get("tagrateeta1")
TagPlot2de2= TagFile1.Get("tagrateeta2")
TagPlot2de3= TagFile1.Get("tagrateeta3")

f.cd()
#---------------------------------------------------------------------------------------------------------------------#
Mtb	    = TH1F("Mtb",     "mass of tb",     	  	      140, 500, 4000 )

Mtbptup	    = TH1F("Mtbptup",     "mass of tb ttbar pt reweighting up",     	  	      140, 500, 4000 )
Mtbptdown	    = TH1F("Mtbptdown",     "mass of tb ttbar pt reweighting up",     	  	      140, 500, 4000 )

MtbBup	    = TH1F("MtbBup",     "mass of tb B tag SF up",     	  	      140, 500, 4000 )
MtbBDown	    = TH1F("MtbBDown",     "mass of tb B tag SF up",     	  	      140, 500, 4000 )


QCDbkg= TH1F("QCDbkg",     "QCD background estimate",     	  	      140, 500, 4000 )
QCDbkgh= TH1F("QCDbkgh",     "QCD background estimate up error",     	  	      140, 500, 4000 )
QCDbkgl= TH1F("QCDbkgl",     "QCD background estimate down error",     	  	      140, 500, 4000 )
QCDbkg2D= TH1F("QCDbkg2D",     "QCD background estimate 2d error",     	  	      140, 500, 4000 )


Mtb.Sumw2()

Mtbptup.Sumw2()
Mtbptdown.Sumw2()

MtbBup.Sumw2()
MtbBDown.Sumw2()

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
jobiter = 0
print "Start looping"
#initialize the ttree variables
tree_vars = {"bpt":array('d',[0.]),"bmass":array('d',[0.]),"btag":array('d',[0.]),"tpt":array('d',[0.]),"tmass":array('d',[0.]),"nsubjets":array('d',[0.]),"sjbtag":array('d',[0.])}
Tree = Make_Trees(tree_vars)

goodEvents = []
totevents = events.size()
print str(totevents)  +  ' Events total'
for event in events:
    weightSFb = 1.0
    errorSFb = 0.0
    count	= 	count + 1
    m = 0
    t = 0
    #if count > 1000:
	#break

    if count % 100000 == 0 :
      print  '--------- Processing Event ' + str(count) +'   -- percent complete ' + str(100*count/totevents) + '% -- '

    #Here we split up event processing based on number of jobs 
    #This is set up to have jobs range from 1 to the total number of jobs (ie dont start at job 0)

    if jobs != 1:
    	if (count - 1) % jobs == 0:
		jobiter+=1
	count_index = count - (jobiter-1)*jobs
	if count_index!=num:
		continue 
    # We load up the relevant handles and labels and create collections

    event.getByLabel (hemis1topLabel, hemis1topHandle)
    topJetsh1 		= 	hemis1topHandle.product()

    event.getByLabel (hemis0topLabel, hemis0topHandle)
    topJetsh0 		= 	hemis0topHandle.product() 

    event.getByLabel (hemis1prCA8Label, hemis1prCA8Handle)
    bJetsh1 		= 	hemis1prCA8Handle.product()

    event.getByLabel (hemis0prCA8Label, hemis0prCA8Handle)
    bJetsh0 		= 	hemis0prCA8Handle.product()

    bjh0 = 0
    bjh1 = 0
    #Require 1 pt>150 jet in each hemisphere (top jets already have the 150GeV requirement) 
    for bjet in bJetsh0:
	if bjet.pt() > 150.0:
		bjh0+=1
    for bjet in bJetsh1:
	if bjet.pt() > 150.0:
		bjh1+=1

    njets11b0 	= 	((len(topJetsh1) == 1) and (bjh0 == 1))
    njets11b1 	= 	((len(topJetsh0) == 1) and (bjh1 == 1))

    for hemis in ['hemis0','hemis1']:
    	if hemis == 'hemis0'   :
		if not njets11b0:
			continue 
		#The Ntuple entries are ordered in pt, so [0] is the highest pt entry
		#We are calling a candidate b jet (highest pt jet in hemisphere0)  
		bjet = bJetsh0[0]
		BDiscLabel = hemis0BDiscLabel
		BDiscHandle = hemis0BDiscHandle

		tjet = topJetsh1[0]
        	TopMassLabel = hemis1TopMassLabel
		TopMassHandle = hemis1TopMassHandle
        	NSubJetsLabel = hemis1NSubJetsLabel
		NSubJetsHandle = hemis1NSubJetsHandle 
        	MinMassLabel = hemis1MinMassLabel
		MinMassHandle = hemis1MinMassHandle
    		TopBDiscsj0CSVLabel = hemis1TopBDiscsj0CSVLabel
		TopBDiscsj0CSVHandle = hemis1TopBDiscsj0CSVHandle
    		TopBDiscsj1CSVLabel = hemis1TopBDiscsj1CSVLabel
		TopBDiscsj1CSVHandle = hemis1TopBDiscsj1CSVHandle
    		TopBDiscsj2CSVLabel = hemis1TopBDiscsj2CSVLabel
		TopBDiscsj2CSVHandle = hemis1TopBDiscsj2CSVHandle
    		TopBDiscsj3CSVLabel = hemis1TopBDiscsj3CSVLabel
		TopBDiscsj3CSVHandle = hemis1TopBDiscsj3CSVHandle


    	if hemis == 'hemis1'  :
		if not njets11b1:
			continue 
		bjet = bJetsh1[0]
		BDiscLabel = hemis1BDiscLabel
		BDiscHandle = hemis1BDiscHandle

		tjet = topJetsh0[0]
        	TopMassLabel = hemis0TopMassLabel
		TopMassHandle = hemis0TopMassHandle
        	NSubJetsLabel = hemis0NSubJetsLabel
		NSubJetsHandle = hemis0NSubJetsHandle 
        	MinMassLabel = hemis0MinMassLabel
		MinMassHandle = hemis0MinMassHandle
    		TopBDiscsj0CSVLabel = hemis0TopBDiscsj0CSVLabel
		TopBDiscsj0CSVHandle = hemis0TopBDiscsj0CSVHandle
    		TopBDiscsj1CSVLabel = hemis0TopBDiscsj1CSVLabel
		TopBDiscsj1CSVHandle = hemis0TopBDiscsj1CSVHandle
    		TopBDiscsj2CSVLabel = hemis0TopBDiscsj2CSVLabel
		TopBDiscsj2CSVHandle = hemis0TopBDiscsj2CSVHandle
    		TopBDiscsj3CSVLabel = hemis0TopBDiscsj3CSVLabel
		TopBDiscsj3CSVHandle = hemis0TopBDiscsj3CSVHandle



    	weight=1.0
     	weightSFptup=1.0
     	weightSFptdown=1.0

    	bpt_cut = bpt[0]<bjet.pt()<bpt[1]
    	tpt_cut = tpt[0]<tjet.pt()<tpt[1]
    	dy_cut = dy[0]<=abs(tjet.Rapidity()-bjet.Rapidity())<dy[1]
    	if bpt_cut and tpt_cut and dy_cut: 

    		if options.pdfweights != "nominal" :
            		event.getByLabel( pdfLabel, pdfHandle )
            		pdfs = pdfHandle.product()
			iweight = PDF_Lookup( pdfs , options.pdfweights )
            		weight *= iweight


		if options.set!="data":

			event.getByLabel (puLabel, puHandle)
    			PileUp 		= 	puHandle.product()
               		bin1 = PilePlot.FindBin(PileUp[0]) 

			if options.pileup != 'off':
				weight *= PilePlot.GetBinContent(bin1)


			if run_b_SF :
				#btagging scale factor reweighting done here
				SFB = SFB_Lookup( bjet.pt() )
				weightSFb = SFB[0]
				errorSFb = SFB[1]

        	event.getByLabel (TopMassLabel, TopMassHandle)
        	topJetMass 	= 	TopMassHandle.product()


		tmass_cut = tmass[0]<topJetMass[0]<tmass[1]

		if tmass_cut :
	 		event.getByLabel (MinMassLabel, MinMassHandle)
         		topJetMinMass 	= 	MinMassHandle.product()
         		event.getByLabel ( NSubJetsLabel , NSubJetsHandle )
    	 		NSubJets 		= 	NSubJetsHandle.product()

			minmass_cut = minmass[0]<=topJetMinMass[0]<minmass[1]
			nsubjets_cut = nsubjets[0]<=NSubJets[0]<nsubjets[1]
	 		if minmass_cut and nsubjets_cut:

   				event.getByLabel (BDiscLabel, BDiscHandle)
    				bJetBDisc 	= 	BDiscHandle.product()
		
                		btag_cut = btag[0]<bJetBDisc[0]<=btag[1]

				ht = tjet.pt() + bjet.pt()
				if options.trigger != "none" :
					#Trigger reweighting done here
					TRW = Trigger_Lookup( ht , TrigPlot , options.trigger )
					weight*=TRW
		
				if options.ptreweight == "on":
					#ttbar pt reweighting done here
					event.getByLabel( GenLabel, GenHandle )
					GenParticles = GenHandle.product()
					PTW = PTW_Lookup( GenParticles )
					weight*=PTW
     					weightSFptup=max(0.0,weight*(2*PTW-1))
     					weightSFptdown=weight

				weightb = weight*weightSFb
				weightSFptup*=weightSFb
				weightSFptdown*=weightSFb
				weightSFbup = weight*(weightSFb+errorSFb)
				weightSFbdown = weight*(weightSFb-errorSFb)




    				event.getByLabel (CA8Label, CA8Handle)
    				CA8Jets 		= 	CA8Handle.product() 
	
    				event.getByLabel (TopBDiscsj0CSVLabel, TopBDiscsj0CSVHandle)
    				Topsj0BDiscCSV 		= 	TopBDiscsj0CSVHandle.product() 

    				event.getByLabel (TopBDiscsj1CSVLabel, TopBDiscsj1CSVHandle)
    				Topsj1BDiscCSV 		= 	TopBDiscsj1CSVHandle.product() 

    				event.getByLabel (TopBDiscsj2CSVLabel, TopBDiscsj2CSVHandle)
    				Topsj2BDiscCSV 		= 	TopBDiscsj2CSVHandle.product() 

    				event.getByLabel (TopBDiscsj3CSVLabel, TopBDiscsj3CSVHandle)
    				Topsj3BDiscCSV 		= 	TopBDiscsj3CSVHandle.product() 

				SJ_csvmax = max(Topsj0BDiscCSV[0],Topsj1BDiscCSV[0],Topsj2BDiscCSV[0],Topsj3BDiscCSV[0])
				sjbtag_cut = sjbtag[0]<SJ_csvmax<=sjbtag[1]

    				event.getByLabel (TopTau2Label, TopTau2Handle)
    				Tau2		= 	TopTau2Handle.product() 
		
    				event.getByLabel (TopTau3Label, TopTau3Handle)
    				Tau3		= 	TopTau3Handle.product() 

				index = -1

				for ijet in range(0,len(CA8Jets)):
					if (abs(ROOT.Math.VectorUtil.DeltaR(CA8Jets[ijet],tjet))<0.5):
						index = ijet
						break

				tau32_cut =  tau32[0]<=Tau3[index]/Tau2[index]<tau32[1]
				if sjbtag_cut :
					if tau32_cut:
						bmass_cut = bmass[0]<=bjet.mass()<bmass[1]
						if bmass_cut:
							eta_regions = [eta1,eta2,eta3]
							BTRweight = bkg_weight(bjet,BTR,eta_regions)
							BTRweightsigsq = bkg_weight(bjet,BTR_err,eta_regions)

							BTRweighterrup = BTRweight+sqrt(BTRweightsigsq)
							BTRweighterrdown = BTRweight-sqrt(BTRweightsigsq)


							eta1_cut = eta1[0]<=abs(bjet.eta())<eta1[1]
							eta2_cut = eta2[0]<=abs(bjet.eta())<eta2[1]
							eta3_cut = eta3[0]<=abs(bjet.eta())<eta3[1]

							if (eta1_cut) :
								xbin = TagPlot2de1.GetXaxis().FindBin(bjet.pt())
								ybin = TagPlot2de1.GetYaxis().FindBin((tjet+bjet).mass())
								tagrate2d = TagPlot2de1.GetBinContent(xbin,ybin)
								QCDbkg2D.Fill((tjet+bjet).mass(),tagrate2d*weight)
			
							if (eta2_cut):
								xbin = TagPlot2de2.GetXaxis().FindBin(bjet.pt())
								ybin = TagPlot2de2.GetYaxis().FindBin((tjet+bjet).mass())
								tagrate2d = TagPlot2de2.GetBinContent(xbin,ybin)
								QCDbkg2D.Fill((tjet+bjet).mass(),tagrate2d*weight)

							if (eta3_cut):

								xbin = TagPlot2de3.GetXaxis().FindBin(bjet.pt())
								ybin = TagPlot2de3.GetYaxis().FindBin((tjet+bjet).mass())
								tagrate2d = TagPlot2de3.GetBinContent(xbin,ybin)
								QCDbkg2D.Fill((tjet+bjet).mass(),tagrate2d*weight)


							for ifit in range(0,len(fittitles)):
									tempweight = bkg_weight(bjet,fits[ifit],eta_regions)
									QCDbkg_ARR[ifit].Fill((tjet+bjet).mass(),tempweight*weight) 

							QCDbkg.Fill((tjet+bjet).mass(),BTRweight*weight)
							QCDbkgh.Fill((tjet+bjet).mass(),BTRweighterrup*weight)
							QCDbkgl.Fill((tjet+bjet).mass(),BTRweighterrdown*weight)  
        				        	if btag_cut:
                                      				goodEvents.append( [ event.object().id().run(), event.object().id().luminosityBlock(), event.object().id().event() ] )
								Mtb.Fill((tjet+bjet).mass(),weightb) 
								MtbBup.Fill((tjet+bjet).mass(),weightSFbup) 
								MtbBDown.Fill((tjet+bjet).mass(),weightSFbdown) 
								Mtbptup.Fill((tjet+bjet).mass(),weightSFptup) 
								Mtbptdown.Fill((tjet+bjet).mass(),weightSFptdown) 
								temp_variables = {"bpt":bjet.pt(),"bmass":bjet.mass(),"btag":bJetBDisc[0],"tpt":tjet.pt(),"tmass":topJetMass[0],"nsubjets":NSubJets[0],"sjbtag":SJ_csvmax}

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
