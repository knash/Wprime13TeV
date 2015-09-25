
###################################################################
##								 ##
## Name: WPrime_Functions.py	   			         ##
## Author: Kevin Nash 						 ##
## Date: 5/13/2015						 ##
## Purpose: This contains all functions used by the              ##
##	    analysis.  A method is generally placed here if 	 ##
##	    it is called more than once in reproducing all	 ##
##	    analysis results.  The functions contained here 	 ##
##	    Are capable of tuning the analysis - such as changing##
##	    cross sections, updating lumi, changing file	 ##
##	    locations, etc. with all changes propegating 	 ##
##	    to all relevant files automatically.  		 ##
##								 ##
###################################################################

import cppyy
import os
import array
import glob
import math
import ROOT
import sys
import math
from math import sqrt
from array import *
from ROOT import *

from DataFormats.FWLite import Events, Handle
#This is the most impostant Function.  Correct information here is essential to obtaining valid results.
#In order we have Luminosity, top tagging scale factor, cross sections for wprime right,left,mixed,ttbar,qcd, and singletop and their corresponding event numbers
#If I wanted to access the left handed W' cross section at 1900 GeV I could do Xsecl1900 = LoadConstants()['xsec_wpl']['1900']
def LoadConstants():
	 return  {
		'lumi':20.38,
		'ttagsf':1.0,
		'xsec_wpr':{'1300':2.9615,'2000':0.21351,'2700':0.041294},
		'xsec_wpl':{'1400': 0.},
		'xsec_wplr':{'1400': 0.},
		'xsec_ttbar':{'MG':831.76 },
		'xsec_qcd':{'300':7823,'470':648.2,'600':186.9,'800':32.293,'1000':9.4183,'1400':0.84265,'800_BROKEN':32.293,'FLAT7000':2022100000},
		'xsec_st':{'s':3.79,'sB':1.76,'t':56.4,'tB':30.7,'tW':11.1,'tWB':11.1},
		'nev_wpr':{'1200':200000,'1400':199200,'1600':198200,'1800':200000,'2000':198400,'2400':197600,'2600':200000,'2800':200000,'3000':199200},
		'nev_wpl':{'2000':474565,},
		'nev_wplr':{'2000':468663,},
		'nev_ttbar':{'MG':11339232},
		'nev_qcd':{'300':2933611,'470':1936832,'600':1964128,'800':1937216,'1000':1487136,'1400':197959},
		'nev_st':{'s':259176,'sB':139604,'t':3748155,'tB':1930185,'tW':495559,'tWB':491463},
		}



#This is also a very impostant Function.  The analysis runs on "PSETS", which correspond to the TYPE variable here.
#These each load a cut profile.  For instance 'default' is the standard selection used to set limits
def LoadCuts(TYPE):
	if TYPE=='default':
 		return  {
			'bpt':[350.0,float("inf")],
			'tpt':[350.0,float("inf")],
			'dy':[0.0,1.6],
			'tmass':[130.0,200.0],
			'nsubjets':[3,10],
			'tau32':[0.0,0.61],
			'minmass':[50.0,float("inf")],
			'sjbtag':[0.890,1.0],
			'bmass':[0.0,70.0],
			'btag':[0.890,1.0],
			'eta1':[0.0,0.5],
			'eta2':[0.5,1.15],	
			'eta3':[1.15,2.4]
			}
	if TYPE=='rate_default':
 		return  {
			'bpt':[350.0,float("inf")],
			'tpt':[350.0,float("inf")],
			'dy':[0.0,1.6],
			'tmass':[130.0,200.0],
			'nsubjets':[1,3],
			'tau32':[0.0,1.0],
			'minmass':[0.0,float("inf")],
			'sjbtag':[0.890,1.0],
			'bmass':[0.0,70.0],
			'btag':[0.890,1.0],
			'eta1':[0.0,0.5],
			'eta2':[0.5,1.15],	
			'eta3':[1.15,2.4]
			}

	if TYPE=='loose':
 		return  {
			'bpt':[350.0,float("inf")],
			'tpt':[350.0,float("inf")],
			'dy':[0.0,1.6],
			'tmass':[130.0,200.0],
			'nsubjets':[3,10],
			'tau32':[0.0,1.0],
			'minmass':[50.0,float("inf")],
			'sjbtag':[0.0,1.0],
			'bmass':[0.0,70.0],
			'btag':[0.890,1.0],
			'eta1':[0.0,0.5],
			'eta2':[0.5,1.15],	
			'eta3':[1.15,2.4]
			}
	if TYPE=='rate_loose':
 		return  {
			'bpt':[350.0,float("inf")],
			'tpt':[350.0,float("inf")],
			'dy':[0.0,1.6],
			'tmass':[130.0,200.0],
			'nsubjets':[1,3],
			'tau32':[0.0,1.0],
			'minmass':[0.0,float("inf")],
			'sjbtag':[0.0,1.0],
			'bmass':[0.0,70.0],
			'btag':[0.890,1.0],
			'eta1':[0.0,0.5],
			'eta2':[0.5,1.15],	
			'eta3':[1.15,2.4]
			}

	if TYPE=='sideband1':	
 		return  {
			'bpt':[350.0,float("inf")],
			'tpt':[350.0,float("inf")],
			'dy':[0.0,1.6],
			'tmass':[130.0,200.0],
			'nsubjets':[3,10],
			'tau32':[0.61,1.0],	#
			'minmass':[0.0,50.0],
			'sjbtag':[0.890,1.0],
			'bmass':[0.0,70.0],
			'btag':[0.890,1.0],
			'eta1':[0.0,0.5],
			'eta2':[0.5,1.15],	
			'eta3':[1.15,2.4]
			}

	if TYPE=='rate_sideband1':	#minmass inv
 		return  {
			'bpt':[350.0,float("inf")],
			'tpt':[350.0,float("inf")],
			'dy':[0.0,1.6],
			'tmass':[130.0,200.0],
			'nsubjets':[0,3],
			'tau32':[0.0,1.0],	#
			'minmass':[0.0,float("inf")],
			'sjbtag':[0.890,1.0],
			'bmass':[0.0,70.0],
			'btag':[0.890,1.0],
			'eta1':[0.0,0.5],
			'eta2':[0.5,1.15],	
			'eta3':[1.15,2.4]
			}

	if TYPE=='sideband2':
 		return  {
			'bpt':[350.0,float("inf")],
			'tpt':[350.0,float("inf")],
			'dy':[0.0,1.6],
			'tmass':[130.0,200.0],
			'nsubjets':[3,10],
			'tau32':[0.0,0.61],	#
			'minmass':[50.0,float("inf")],
			'sjbtag':[0.0,0.890],
			'bmass':[0.0,70.0],
			'btag':[0.890,1.0],
			'eta1':[0.0,0.5],
			'eta2':[0.5,1.15],	
			'eta3':[1.15,2.4]
			}

	if TYPE=='rate_sideband2':	#Sjbtag inv
 		return  {
			'bpt':[350.0,float("inf")],
			'tpt':[350.0,float("inf")],
			'dy':[0.0,1.6],
			'tmass':[130.0,200.0],
			'nsubjets':[0,3],
			'tau32':[0.0,1.0],	#
			'minmass':[0.0,float("inf")],
			'sjbtag':[0.0,0.890],
			'bmass':[0.0,70.0],
			'btag':[0.890,1.0],
			'eta1':[0.0,0.5],
			'eta2':[0.5,1.15],	
			'eta3':[1.15,2.4]
			}

	if TYPE=='sideband3':
 		return  {
			'bpt':[350.0,float("inf")],
			'tpt':[350.0,float("inf")],
			'dy':[0.0,1.6],
			'tmass':[130.0,200.0],
			'nsubjets':[3,10],
			'tau32':[0.0,0.61],
			'minmass':[50.0,float("inf")],
			'sjbtag':[0.890,1.0],
			'bmass':[70.0,float("inf")],	#inverted for lots of ttbar
			'btag':[0.890,1.0],
			'eta1':[0.0,0.5],
			'eta2':[0.5,1.15],	
			'eta3':[1.15,2.4]
			}
	if TYPE=='rate_sideband3':
 		return  {
			'bpt':[350.0,float("inf")],
			'tpt':[350.0,float("inf")],
			'dy':[0.0,1.6],
			'tmass':[130.0,200.0],
			'nsubjets':[1,3],
			'tau32':[0.0,1.0],
			'minmass':[0.0,float("inf")],
			'sjbtag':[0.890,1.0],
			'bmass':[0.0,70.0],	
			'btag':[0.890,1.0],
			'eta1':[0.0,0.5],
			'eta2':[0.5,1.15],	
			'eta3':[1.15,2.4]
			}


#This function loads up Ntuples based on what type of set you want to analyze.  
#This needs to be updated whenever new Ntuples are produced (unless the file locations are the same).
def Load_Ntuples(string,bx):
	print 'running on ' + string 
	if string == 'data':
		files = glob.glob("/eos/uscms/store/user/knash/JetHT/crab_JetHT_Run2015D_25ns_Sep24_7411/150924_220444/0000/*.root")
		#files += glob.glob("/uscms_data/d3/knash/WPrime8TeV/data/CMSSW_5_3_18/src/Analysis/TTBSMPatTuples/test/Run2012B-22Jan2013/res/*.root")
	
 	if string == 'ttbar':
 		files = glob.glob("/eos/uscms/store/user/knash/TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_b2ganafw741_TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v2/150608_221852/0000/*.root")
 	if string == 'QCDPT300':
 		files = glob.glob("/eos/uscms/store/user/lcorcodi/QCD_Pt_300to470_TuneCUETP8M1_13TeV_pythia8/crab_b2ganafw74xV2_QCD_Pt_300to470_TuneCUETP8M1_13TeV_pythia8_RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/150824_212922/0000/*.root")
 	if string == 'QCDPT470':
 		files = glob.glob("/eos/uscms/store/user/lcorcodi/QCD_Pt_470to600_TuneCUETP8M1_13TeV_pythia8/crab_b2ganafw74xV2_QCD_Pt_470to600_TuneCUETP8M1_13TeV_pythia8_RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v2/150825_210455/0000/*.root")
 	if string == 'QCDPT600':
 		files = glob.glob("/eos/uscms/store/user/lcorcodi/QCD_Pt_600to800_TuneCUETP8M1_13TeV_pythia8/crab_b2ganafw74xV2_QCD_Pt_600to800_TuneCUETP8M1_13TeV_pythia8_RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v3/150824_205008/0000/*.root")
 	if string == 'QCDPT800':
 		files = glob.glob("/eos/uscms/store/user/lcorcodi/QCD_Pt_800to1000_TuneCUETP8M1_13TeV_pythia8/crab_b2ganafw74xV2_QCD_Pt_800to1000_TuneCUETP8M1_13TeV_pythia8_RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v2/150826_123934/0000/*.root")
 	if string == 'QCDPT1000':
 		files = glob.glob("/eos/uscms/store/user/lcorcodi/QCD_Pt_1000to1400_TuneCUETP8M1_13TeV_pythia8/crab_b2ganafw74xV2_QCD_Pt_1000to1400_TuneCUETP8M1_13TeV_pythia8_RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/150824_213158/0000/*.root")
 	if string == 'QCDPT1400':
 		files = glob.glob("/eos/uscms/store/user/lcorcodi/QCD_Pt_1400to1800_TuneCUETP8M1_13TeV_pythia8/crab_b2ganafw74xV2_QCD_Pt_1400to1800_TuneCUETP8M1_13TeV_pythia8_RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/150824_213219/0000/*.root")



	if string == 'signalright1200':
		files = glob.glob("/eos/uscms/store/user/knash/WprimeToTB_TToHad_M-1200_RH_TuneCUETP8M1_13TeV-comphep-pythia8/crab_WPrime13TeV_B2GAnaFW_V5p2_M1200_RH_25ns/150923_161540/0000/*.root")
	if string == 'signalright1400':
		files = glob.glob("/eos/uscms/store/user/knash/WprimeToTB_TToHad_M-1400_RH_TuneCUETP8M1_13TeV-comphep-pythia8/crab_WPrime13TeV_B2GAnaFW_V5p2_M1400_RH_25ns/150923_154748/0000/*.root")
	if string == 'signalright1600':
		files = glob.glob("/eos/uscms/store/user/knash/WprimeToTB_TToHad_M-1600_RH_TuneCUETP8M1_13TeV-comphep-pythia8/crab_WPrime13TeV_B2GAnaFW_V5p2_M1600_RH_25ns/150918_190442/0000/*.root")
	if string == 'signalright1800':
		files = glob.glob("/eos/uscms/store/user/knash/WprimeToTB_TToHad_M-1800_RH_TuneCUETP8M1_13TeV-comphep-pythia8/crab_WPrime13TeV_B2GAnaFW_V5p2_M1800_RH_25ns/150918_190502/0000/*.root")
	if string == 'signalright2000':
		files = glob.glob("/eos/uscms/store/user/knash/WprimeToTB_TToHad_M-2000_RH_TuneCUETP8M1_13TeV-comphep-pythia8/crab_WPrime13TeV_B2GAnaFW_V5p2_M2000_RH_25ns/150924_210409/0000/*.root")
	if string == 'signalright2400':
		files = glob.glob("/eos/uscms/store/user/knash/WprimeToTB_TToHad_M-2400_RH_TuneCUETP8M1_13TeV-comphep-pythia8/crab_WPrime13TeV_B2GAnaFW_V5p2_M2400_RH_25ns/150918_190557/0000/*.root")
	if string == 'signalright2600':
		files = glob.glob("/eos/uscms/store/user/knash/WprimeToTB_TToHad_M-2600_RH_TuneCUETP8M1_13TeV-comphep-pythia8/crab_WPrime13TeV_B2GAnaFW_V5p2_M2600_RH_25ns/150918_190616/0000/*.root")
	if string == 'signalright2800':
		files = glob.glob("/eos/uscms/store/user/knash/WprimeToTB_TToHad_M-2800_RH_TuneCUETP8M1_13TeV-comphep-pythia8/crab_WPrime13TeV_B2GAnaFW_V5p2_M2800_RH_25ns/150923_160258/0000/*.root")
	if string == 'signalright3000':
		files = glob.glob("/eos/uscms/store/user/knash/WprimeToTB_TToHad_M-3000_RH_TuneCUETP8M1_13TeV-comphep-pythia8/crab_WPrime13TeV_B2GAnaFW_V5p2_M3000_RH_25ns/150918_190657/0000/*.root")

	try:
		print 'A total of ' + str(len(files)) + ' files'
	except:
		print 'Bad files option'
	return files




#This function initializes the average b tagging rates used for QCD determination
#It tages the type of functional form as an argument.  The default fit is Bifpoly

#This is a poorly written function, but I cant think of a better way to do this 
#It works, but you should be able to just have one input

def BTR_Init(ST,CUT,di,setval):

	if setval == "ttbar":
		setval = "data"


	if ST == 'Bifpoly':
		TRBPE1 = open(di+"fitdata/bpinput"+setval+"eta1_PSET_"+CUT+".txt")
		TRBPE1.seek(0)
		TRBPE2 = open(di+"fitdata/bpinput"+setval+"eta2_PSET_"+CUT+".txt")
		TRBPE2.seek(0)
		TRBPE3 = open(di+"fitdata/bpinput"+setval+"eta3_PSET_"+CUT+".txt")
		TRBPE3.seek(0)
		eta1fit = TF1("eta1fit",BifPoly,0,1400,5)
		eta2fit = TF1("eta2fit",BifPoly,0,1400,5)
		eta3fit = TF1("eta3fit",BifPoly,0,1400,5)
		Params = 5
	if ST == 'Bifpoly_err':
		TRBPE1 = open(di+"fitdata/bperrorinput"+setval+"eta1_PSET_"+CUT+".txt")
		TRBPE1.seek(0)
		TRBPE2 = open(di+"fitdata/bperrorinput"+setval+"eta2_PSET_"+CUT+".txt")
		TRBPE2.seek(0)
		TRBPE3 = open(di+"fitdata/bperrorinput"+setval+"eta3_PSET_"+CUT+".txt")
		TRBPE3.seek(0)
		eta1fit=TF1("eta1fit",BifPolyErr,0,1400,10)
		eta2fit=TF1("eta2fit",BifPolyErr,0,1400,10)
		eta3fit=TF1("eta3fit",BifPolyErr,0,1400,10)
		Params = 10

	if ST == 'pol0':
		TRBPE1 = open(di+"fitdata/pol0input"+setval+"eta1_PSET_"+CUT+".txt")
		TRBPE1.seek(0)
		TRBPE2 = open(di+"fitdata/pol0input"+setval+"eta2_PSET_"+CUT+".txt")
		TRBPE2.seek(0)
		TRBPE3 = open(di+"fitdata/pol0input"+setval+"eta3_PSET_"+CUT+".txt")
		TRBPE3.seek(0)
		eta1fit = TF1("eta1fit",'pol0',0,1400)
		eta2fit = TF1("eta2fit",'pol0',0,1400)
		eta3fit = TF1("eta3fit",'pol0',0,1400)
		Params = 1

	if ST == 'pol2':
		TRBPE1 = open(di+"fitdata/pol2input"+setval+"eta1_PSET_"+CUT+".txt")
		TRBPE1.seek(0)
		TRBPE2 = open(di+"fitdata/pol2input"+setval+"eta2_PSET_"+CUT+".txt")
		TRBPE2.seek(0)
		TRBPE3 = open(di+"fitdata/pol2input"+setval+"eta3_PSET_"+CUT+".txt")
		TRBPE3.seek(0)
		eta1fit = TF1("eta1fit",'pol2',0,1400)
		eta2fit = TF1("eta2fit",'pol2',0,1400)
		eta3fit = TF1("eta3fit",'pol2',0,1400)
		Params = 3

	if ST == 'pol3':
		TRBPE1 = open(di+"fitdata/pol3input"+setval+"eta1_PSET_"+CUT+".txt")
		TRBPE1.seek(0)
		TRBPE2 = open(di+"fitdata/pol3input"+setval+"eta2_PSET_"+CUT+".txt")
		TRBPE2.seek(0)
		TRBPE3 = open(di+"fitdata/pol3input"+setval+"eta3_PSET_"+CUT+".txt")
		TRBPE3.seek(0)
		eta1fit = TF1("eta1fit",'pol3',0,1400)
		eta2fit = TF1("eta2fit",'pol3',0,1400)
		eta3fit = TF1("eta3fit",'pol3',0,1400)
		Params = 4
	if ST == 'FIT':
		TRBPE1 = open(di+"fitdata/newfitinput"+setval+"eta1_PSET_"+CUT+".txt")
		TRBPE1.seek(0)
		TRBPE2 = open(di+"fitdata/newfitinput"+setval+"eta2_PSET_"+CUT+".txt")
		TRBPE2.seek(0)
		TRBPE3 = open(di+"fitdata/newfitinput"+setval+"eta3_PSET_"+CUT+".txt")
		TRBPE3.seek(0)
		eta1fit = TF1("eta1fit",'[0]*([1]+x)/([2]+x)+[3]*x',0,1400)
		eta2fit = TF1("eta2fit",'[0]*([1]+x)/([2]+x)+[3]*x',0,1400)
		eta3fit = TF1("eta3fit",'[0]*([1]+x)/([2]+x)+[3]*x',0,1400)
		Params = 4
	if ST == 'expofit':
		TRBPE1 = open(di+"fitdata/expoconinput"+setval+"eta1_PSET_"+CUT+".txt")
		TRBPE1.seek(0)
		TRBPE2 = open(di+"fitdata/expoconinput"+setval+"eta2_PSET_"+CUT+".txt")
		TRBPE2.seek(0)
		TRBPE3 = open(di+"fitdata/expoconinput"+setval+"eta3_PSET_"+CUT+".txt")
		TRBPE3.seek(0)
		eta1fit = TF1("eta1fit",'expo(0) + pol0(2)',0,1400)
		eta2fit = TF1("eta2fit",'expo(0) + pol0(2)',0,1400)
		eta3fit = TF1("eta3fit",'expo(0) + pol0(2)',0,1400)
		Params = 3

	TBP1 = TRBPE1.read()
	TBP2 = TRBPE2.read()
	TBP3 = TRBPE3.read()
	
	for i in range(0,Params):

		eta1fit.SetParameter(i,float(TBP1.split('\n')[i]) )
		eta2fit.SetParameter(i,float(TBP2.split('\n')[i]) )
		eta3fit.SetParameter(i,float(TBP3.split('\n')[i]) )

	return [eta1fit.Clone(),eta2fit.Clone(),eta3fit.Clone()] 

#This takes the average b tagging rates that are initialized in the above function and produces 
#A QCD background estimate based on them 
def bkg_weight(blv, funcs, etabins):
	for ibin in range(0,len(etabins)):
		if (etabins[ibin][0] <= abs(blv.Eta()) < etabins[ibin][1]) :
			tagratept = funcs[ibin].Eval(blv.Perp())	

	try :
		return tagratept
	except :
		print blv.Eta()
		tagratept=0.0	

#This is the bifurcated polynomial function and its associated uncertainty 
def BifPoly( x, p ):
	xx=x[0]
	if xx<p[4]:
      		return p[0]+p[1]*xx+p[2]*(xx-p[4])**2
        else:
		return p[0]+p[1]*xx+p[3]*(xx-p[4])**2
def BifPolyErr( x, p ):
	xx=x[0]
	if xx<p[9]:
      		return p[0]+p[1]*xx**2+p[2]*(xx-p[9])**4+p[3]*xx+p[4]*(xx-p[9])**2+p[5]*xx*(xx-p[9])**2
        else:
		return p[0]+p[1]*xx**2+p[6]*(xx-p[9])**4+p[3]*xx+p[7]*(xx-p[9])**2+p[8]*xx*(xx-p[9])**2

#This is the first in a series of functions used to extract Monte Carlo to data scale factors and their uncertainty 
#This looks up the b tagging scale factor 
def SFB_Lookup( Y ):
	ptminsfb = [320, 400, 500, 600]
	ptmaxsfb = [400, 500, 600, 800]	
	SFb_error = [0.0313175,0.0415417,0.0740446,0.0596716]
	SFb = TFormula("SFb","(0.938887+(0.00017124*x))+(-2.76366e-07*(x*x))")

	if Y <= 800.0:
		weightSFb  = SFb.Eval(Y)
		if ptminsfb[0] < Y <= ptmaxsfb[0]:
			errorSFb = SFb_error[0]
		if ptminsfb[1] < Y <= ptmaxsfb[1]:
			errorSFb = SFb_error[1]
		if ptminsfb[2] < Y <= ptmaxsfb[2]:
			errorSFb = SFb_error[2]
		if ptminsfb[3] < Y <= ptmaxsfb[3]:
			errorSFb = SFb_error[3]
	else: 
		weightSFb  = SFb.Eval(800.0)
		errorSFb = 2*SFb_error[3]
	return [weightSFb,errorSFb]
#This looks up the PDF uncertainty
def PDF_Lookup( pdfs , pdfOP ):
	iweight = 0.0
        if pdfOP == "up" :
       		for pdf in pdfs[1::2] :
              		iweight = iweight + pdf
        else :
        	for pdf in pdfs[2::2] :
        		iweight = iweight + pdf
        return (iweight/pdfs[0]) / (len(pdfs)-1) * 2.0
#This looks up the b tagging scale factor or uncertainty
def Trigger_Lookup( H , TRP):
        Weight = 1.0
	
        if H < 1300.0:
                bin0 = TRP.FindBin(H) 
                jetTriggerWeight = TRP.GetBinContent(bin0)
                Weight = jetTriggerWeight
	return Weight

def Trigger_Pass(tnamestr,trigs,bits):
	###TAKE OUT!
	tnamestr = ['HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV0p45_v2','HLT_PFHT800_v1']

	for t in range(0,len(trigs)):
		for tname in tnamestr:	
			if trigs[t]==tname and bits[t] :
				return True
	return False

#This looks up the ttbar pt reweighting scale factor 
def PTW_Lookup( GP ):
		genTpt = -100.
		genTBpt = -100	
		for ig in GP :
			isT = ig.pdgId() == 6 and ig.status() == 3
			isTB = ig.pdgId() == -6 and ig.status() == 3
			if isT:
				genTpt = ig.Perp()
			if isTB:
				genTBpt = ig.Perp()	
		if (genTpt<0) or (genTBpt<0):
			print "ERROR"

      		wTPt = exp(0.156-0.00137*genTpt)
      		wTbarPt = exp(0.156-0.00137*genTBpt)
      		return sqrt(wTPt*wTbarPt)


#This is just a quick function to automatically make a tree
#This is used right now to automatically output branches used to validate the cuts used in a run
def Make_Trees(Floats):
        t = TTree("Tree", "Tree");
	print "Booking trees"
	for F in Floats.keys():
    		t.Branch(F, Floats[F], F+"/D")
	return t

#This takes all of the alternative fit forms for the average b tagging rate and 
#Compares them to the chosen nominal fit (bifpoly).  It outputs the mean squared error uncertainty from this comparison 
def Fit_Uncertainty(List):
	sigmah	    = List[0]
	fits=len(List)-1
	for ihist in range(0,len(List)):
		print List[ihist].GetName()
		if List[ihist].GetName().find('Bifpoly') != -1:
			nominalhist = List[ihist]
	for ibin in range(0,nominalhist.GetXaxis().GetNbins()+1):

		mse=0.0
		sigma=0.0
		sumsqdiff = 0.0
		for ihist in range(0,len(List)):
			if List[ihist].GetName().find('Bifpoly') == -1:
				sumsqdiff+=(List[ihist].GetBinContent(ibin)-nominalhist.GetBinContent(ibin))*(List[ihist].GetBinContent(ibin)-nominalhist.GetBinContent(ibin))
		mse = sumsqdiff/fits
		sigma = sqrt(mse)
		sigmah.SetBinContent(ibin,sigma)
	
	return sigmah

#Makes the blue pull plots
def Make_Pull_plot( DATA,BKG,BKGUP,BKGDOWN ):
	pull = DATA.Clone("pull")
	pull.Add(BKG,-1)
	sigma = 0.0
	FScont = 0.0
	BKGcont = 0.0
	for ibin in range(1,pull.GetNbinsX()+1):
		FScont = DATA.GetBinContent(ibin)
		BKGcont = BKG.GetBinContent(ibin)
		if FScont>=BKGcont:
			FSerr = DATA.GetBinErrorLow(ibin)
			BKGerr = abs(BKGUP.GetBinContent(ibin)-BKG.GetBinContent(ibin))
		if FScont<BKGcont:
			FSerr = DATA.GetBinErrorUp(ibin)
			BKGerr = abs(BKGDOWN.GetBinContent(ibin)-BKG.GetBinContent(ibin))
		sigma = sqrt(FSerr*FSerr + BKGerr*BKGerr)
		if FScont == 0.0:
			pull.SetBinContent(ibin, 0.0 )	
		else:
			if sigma != 0 :
				pullcont = (pull.GetBinContent(ibin))/sigma
				pull.SetBinContent(ibin, pullcont)
			else :
				pull.SetBinContent(ibin, 0.0 )
	return pull

def Initlv(string):
	PtHandle 	= 	Handle (  "vector<float> "  )
	PtLabel  	= 	( "jetsAK8" , "jetAK8Pt")

	EtaHandle 	= 	Handle (  "vector<float> "  )
	EtaLabel  	= 	( "jetsAK8" , "jetAK8Eta")

	PhiHandle 	= 	Handle (  "vector<float> "  )
	PhiLabel  	= 	( "jetsAK8" , "jetAK8Phi")

	MassHandle 	= 	Handle (  "vector<float> "  )
	MassLabel  	= 	( "jetsAK8" , "jetAK8Mass")

	return [[PtHandle,PtLabel],[EtaHandle,EtaLabel],[PhiHandle,PhiLabel],[MassHandle,MassLabel]]

def Makelv(vector,event):

    	event.getByLabel (vector[0][1], vector[0][0])
    	Pt 		= 	vector[0][0].product()

    	event.getByLabel (vector[1][1], vector[1][0])
    	Eta 		= 	vector[1][0].product()

    	event.getByLabel (vector[2][1], vector[2][0])
    	Phi 		= 	vector[2][0].product()

    	event.getByLabel (vector[3][1], vector[3][0])
    	Mass 		= 	vector[3][0].product()

	lvs = []
	for i in range(0,len(Pt)):

		#lvs.append(ROOT.Math.LorentzVector('ROOT::Math::PtEtaPhiM4D<double>')(Pt[i],Eta[i],Phi[i],Mass[i]))

		lvs.append(TLorentzVector())
		lvs[i].SetPtEtaPhiM(Pt[i],Eta[i],Phi[i],Mass[i])
	return lvs


def Hemispherize(LV1,LV2):
	tjets = [[],[]]
	bjets = [[],[]]
	for iLV1 in range(0,len(LV1)):
		if Math.VectorUtil.DeltaPhi(LV1[0],LV1[iLV1])> TMath.Pi()/2:
			tjets[1].append(iLV1)
		else:
			tjets[0].append(iLV1)
	for iLV2 in range(0,len(LV2)):
		if Math.VectorUtil.DeltaPhi(LV1[0],LV2[iLV2])> TMath.Pi()/2:
			bjets[1].append(iLV2)
		else:
			bjets[0].append(iLV2)
	return tjets,bjets


#Some lazy string formatting functions 
def strf( x ):
	return '%.2f' % x

def strf1( x ):
	return '%.0f' % x

