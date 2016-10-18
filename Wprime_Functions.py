
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
		#Oct: 925 PR4: 1628
		
		'lumi':2553.0,
		'kfac':1.2,
		'xsec_wpr':{'1000':4.52,'1100':3.09,'1200':2.17,'1300':1.55,'1400':1.129,'1500':0.834,'1600':0.623,'1700':0.471,'1800':0.359,'1900':0.276,'2000':0.214,'2100':0.166,'2200':0.130,'2300':0.103,'2400':0.0813 ,'2500':0.0646,'2600':0.0516,'2700':0.0413,'2800':0.0318,'2900':0.0267,'3000':0.0216},
		'xsec_wpl':{'1400': 0.},
		'xsec_wplr':{'1400': 0.},
		'xsec_ttbar':{'MG':831.76,'PH':831.76,'PHscaleup':831.76,'PHscaledown':831.76 },
		'xsec_qcd':{'300':7823,'470':648.2,'600':186.9,'800':32.293,'1000':9.4183,'1400':0.84265,'1800':0.1091,'HT500':31630,'HT700':6802,'HT1000':1206,'HT1500':120.4,'HT2000':25.25},
		'xsec_st':{'S':11.36,'T':216.97,'TW':35.85,'TWB':35.85},
		'xsec_wjets':{'HT600':95.14},
		'nev_wpr':{'1000':198400,'1100':200000,'1200':200000,'1300':197600,'1400':157200,'1500':104800,'1600':198200,'1700':200000,'1800':200000,'1900':200000,'2000':198400,'2100':200000,'2200':200000,'2300':200000,'2400':197600,'2500':200000,'2600':200000,'2700':194200,'2800':200000,'2900':199200,'3000':199200},
		'nev_wpl':{'2000':474565,},
		'nev_wplr':{'2000':468663,},
		'nev_wjets':{'HT600':1008034},
		'nev_ttbar':{'MG':11339232,'PH':19757190,'PHscaleup':9921174,'PHscaledown':9860774},
		'nev_qcd':{'300':2933611,'470':1936832,'600':1964128,'800':1937216,'1000':1487136,'1400':197959,'1800':193608,'HT500':19664159,'HT700':15077752,'HT1000':4963895,'HT1500':3868886,'HT2000':1912529},
		'nev_st':{'S':984400,'T':49858384,'TW':995600 ,'TWB':988500},
		}



#This is also a very impostant Function.  The analysis runs on "PSETS", which correspond to the TYPE variable here.
#These each load a cut profile.  For instance 'default' is the standard selection used to set limits
def LoadCuts(TYPE):
	if TYPE=='default':
 		return  {
			'bpt':[350.0,float("inf")],
			'tpt':[350.0,float("inf")],
			'ptmincut':[200.,float("inf")],
			'dy':[0.0,1.3],
			'tmass':[110.0,210.0],
			#'nsubjets':[3,10],
			'nsubjets':[0,10],
			'tau32':[0.0,0.61],
			#'minmass':[50.0,float("inf")],
			'minmass':[-float("inf"),float("inf")],
			'sjbtag':[0.760,1.0],
			'bmass':[0.0,70.0],
			#'btag':[0.605,1.0],
			'btag':[0.605,1.0],
			'eta1':[0.0,0.5],
			'eta2':[0.5,1.15],	
			'eta3':[1.15,2.4]
			}
	if TYPE=='rate_default':
 		return  {
			'bpt':[350.0,float("inf")],
			'tpt':[350.0,float("inf")],
			'ptmincut':[200.,float("inf")],
			'dy':[0.0,1.3],
			'tmass':[50.0,170.0],
			#'nsubjets':[1,3],
			'nsubjets':[0,10],
			'tau32':[0.75,1.0],
			#'minmass':[50.0,float("inf")],
			'minmass':[-float("inf"),float("inf")],
			'sjbtag':[0.760,1.0],
			'bmass':[0.0,70.0],
			#'btag':[0.605,1.0],
			'btag':[0.605,1.0],
			'eta1':[0.0,0.5],
			'eta2':[0.5,1.15],	
			'eta3':[1.15,2.4]
			}

	if TYPE=='highmass':
 		return  {
			'bpt':[350.0,float("inf")],
			'tpt':[350.0,float("inf")],
			'ptmincut':[200.,float("inf")],
			'dy':[0.0,1.6],
			'tmass':[110.0,210.0],
			#'nsubjets':[3,10],
			'nsubjets':[0,10],
			'tau32':[0.0,0.61],
			#'minmass':[50.0,float("inf")],
			'minmass':[-float("inf"),float("inf")],
			'sjbtag':[0.760,1.0],
			'bmass':[0.0,70.0],
			#'btag':[0.605,1.0],
			'btag':[0.89,1.0],
			'eta1':[0.0,0.5],
			'eta2':[0.5,1.15],	
			'eta3':[1.15,2.4]
			}
	if TYPE=='rate_highmass':
 		return  {
			'bpt':[350.0,float("inf")],
			'tpt':[350.0,float("inf")],
			'ptmincut':[200.,float("inf")],
			'dy':[0.0,1.6],
			'tmass':[50.0,170.0],
			#'nsubjets':[1,3],
			'nsubjets':[0,10],
			'tau32':[0.75,1.0],
			#'minmass':[50.0,float("inf")],
			'minmass':[-float("inf"),float("inf")],
			'sjbtag':[0.760,1.0],
			'bmass':[0.0,70.0],
			#'btag':[0.605,1.0],
			'btag':[0.89,1.0],
			'eta1':[0.0,0.5],
			'eta2':[0.5,1.15],	
			'eta3':[1.15,2.4]
			}


	if TYPE=='loose':
 		return  {
			'bpt':[350.0,float("inf")],
			'tpt':[350.0,float("inf")],
			'ptmincut':[200.,float("inf")],
			'dy':[0.0,1.3],
			'tmass':[110.0,210.0],
			#'nsubjets':[3,10],
			'nsubjets':[0,10],
			'tau32':[0.0,1.0],
			#'minmass':[50.0,float("inf")],
			'minmass':[-float("inf"),float("inf")],
			'sjbtag':[0.0,1.0],
			'bmass':[0.0,70.0],
		#	'btag':[0.605,1.0],
			'btag':[0.605,1.0],
			'eta1':[0.0,0.5],
			'eta2':[0.5,1.15],	
			'eta3':[1.15,2.4]
			}
	if TYPE=='rate_loose':
 		return  {
			'bpt':[350.0,float("inf")],
			'tpt':[350.0,float("inf")],
			'ptmincut':[200.,float("inf")],
			'dy':[0.0,1.3],
			'tmass':[50.0,170.0],
			#'nsubjets':[1,3],
			'nsubjets':[0,10],
			'tau32':[0.75,1.0],
			'minmass':[-float("inf"),float("inf")],
			'sjbtag':[0.0,1.0],
			'bmass':[0.0,70.0],
			#'btag':[0.605,1.0],
			'btag':[0.605,1.0],
			'eta1':[0.0,0.5],
			'eta2':[0.5,1.15],	
			'eta3':[1.15,2.4]
			}

	if TYPE=='sideband2':
 		return  {
			'bpt':[350.0,float("inf")],
			'tpt':[350.0,float("inf")],
			'ptmincut':[200.,float("inf")],
			'dy':[0.0,1.3],
			'tmass':[110.0,210.0],
			#'nsubjets':[3,10],
			'nsubjets':[0,10],
			'tau32':[0.0,0.61],	#
			#'minmass':[50.0,float("inf")],
			'minmass':[-float("inf"),float("inf")],
			'sjbtag':[0.0,0.760],
			'bmass':[0.0,70.0],
			#'btag':[0.605,1.0],
			'btag':[0.605,1.0],
			'eta1':[0.0,0.5],
			'eta2':[0.5,1.15],	
			'eta3':[1.15,2.4]
			}

	if TYPE=='rate_sideband2':	#Sjbtag inv
 		return  {
			'bpt':[350.0,float("inf")],
			'tpt':[350.0,float("inf")],
			'ptmincut':[200.,float("inf")],
			'dy':[0.0,1.3],
			'tmass':[50.0,170.0],
			#'nsubjets':[1,3],
			'nsubjets':[0,10],
			'tau32':[0.75,1.0],	#
			'minmass':[-float("inf"),float("inf")],
			'sjbtag':[0.0,0.760],
			'bmass':[0.0,70.0],
			#'btag':[0.605,1.0],
			'btag':[0.605,1.0],
			'eta1':[0.0,0.5],
			'eta2':[0.5,1.15],	
			'eta3':[1.15,2.4]
			}

	if TYPE=='sideband3':
 		return  {
			'bpt':[350.0,float("inf")],
			'tpt':[350.0,float("inf")],
			'ptmincut':[200.,float("inf")],
			'dy':[0.0,1.3],
			'tmass':[110.0,210.0],
		#	'nsubjets':[3,10],
			'nsubjets':[0,10],
			'tau32':[0.0,0.61],
			#'minmass':[50.0,float("inf")],
			'minmass':[-float("inf"),float("inf")],
			'sjbtag':[0.760,1.0],
			'bmass':[70.0,float("inf")],	#inverted for lots of ttbar
			#'btag':[0.605,1.0],
			'btag':[0.605,1.0],
			'eta1':[0.0,0.5],
			'eta2':[0.5,1.15],	
			'eta3':[1.15,2.4]
			}
	if TYPE=='rate_sideband3':
 		return  {
			'bpt':[350.0,float("inf")],
			'tpt':[350.0,float("inf")],
			'ptmincut':[200.,float("inf")],
			'dy':[0.0,1.3],
			'tmass':[50.0,170.0],
			#'nsubjets':[1,3],
			'nsubjets':[0,10],
			'tau32':[0.75,1.0],
			'minmass':[-float("inf"),float("inf")],
			'sjbtag':[0.760,1.0],
			'bmass':[70.0,float("inf")],	#inverted for lots of ttbar
			#'btag':[0.605,1.0],
			'btag':[0.605,1.0],
			'eta1':[0.0,0.5],
			'eta2':[0.5,1.15],	
			'eta3':[1.15,2.4]
			}


	if TYPE=='antibtag':
 		return  {
			'bpt':[350.0,float("inf")],
			'tpt':[350.0,float("inf")],
			'ptmincut':[200.,float("inf")],
			'dy':[0.0,1.3],
			'tmass':[110.0,210.0],
		#	'nsubjets':[3,10],
			'nsubjets':[0,10],
			'tau32':[0.0,0.61],
			#'minmass':[50.0,float("inf")],
			'minmass':[-float("inf"),float("inf")],
			'sjbtag':[0.760,1.0],
			'bmass':[70.0,float("inf")],	#inverted for lots of ttbar
			#'btag':[0.605,1.0],
			'btag':[0.0,0.650],
			'eta1':[0.0,0.5],
			'eta2':[0.5,1.15],	
			'eta3':[1.15,2.4]
			}
	if TYPE=='rate_antibtag':
 		return  {
			'bpt':[350.0,float("inf")],
			'tpt':[350.0,float("inf")],
			'ptmincut':[200.,float("inf")],
			'dy':[0.0,1.3],
			'tmass':[50.0,170.0],
			#'nsubjets':[1,3],
			'nsubjets':[0,10],
			'tau32':[0.75,1.0],
			'minmass':[-float("inf"),float("inf")],
			'sjbtag':[0.760,1.0],
			'bmass':[70.0,float("inf")],	#inverted for lots of ttbar
			#'btag':[0.605,1.0],
			'btag':[0.0,0.650],
			'eta1':[0.0,0.5],
			'eta2':[0.5,1.15],	
			'eta3':[1.15,2.4]
			}




	if TYPE=='cutmod':
 		return  {
			'bpt':[350.0,float("inf")],
			'tpt':[350.0,float("inf")],
			'ptmincut':[200.,float("inf")],
			'dy':[0.0,1.3],
			'tmass':[110.0,210.0],
			#'nsubjets':[3,10],
			'nsubjets':[0,10],
			'tau32':[0.0,0.61],
			#'minmass':[50.0,float("inf")],
			'minmass':[-float("inf"),float("inf")],
			'sjbtag':[0.760,1.0],
			'bmass':[0.0,70.0],
			#'btag':[0.605,1.0],
			'btag':[0.605,1.0],
			'eta1':[0.0,0.5],
			'eta2':[0.5,1.15],	
			'eta3':[1.15,2.4]
			}
	if TYPE=='cutmod_default':
 		return  {
			'bpt':[350.0,float("inf")],
			'tpt':[350.0,float("inf")],
			'ptmincut':[200.,float("inf")],
			'dy':[0.0,1.3],
			'tmass':[50.0,170.0],
			#'nsubjets':[1,3],
			'nsubjets':[0,10],
			'tau32':[0.75,1.0],
			'minmass':[-float("inf"),float("inf")],
			'sjbtag':[0.760,1.0],
			'bmass':[0.0,70.0],
			#'btag':[0.605,1.0],
			'btag':[0.605,1.0],
			'eta1':[0.0,0.5],
			'eta2':[0.5,1.15],	
			'eta3':[1.15,2.4]
			}



#This function loads up Ntuples based on what type of set you want to analyze.  
#This needs to be updated whenever new Ntuples are produced (unless the file locations are the same).
def Load_Ntuples(string,bx):
	print 'running on ' + string 
	if string == 'data':
		files = glob.glob("/uscms_data/d3/knash/SlimNtuples/JetHT/crab_JetHT_Run2015D-05Oct2015-v1_B2GAnaFW_V8p4_Slim_V10/*/0000/*.root")
		files += glob.glob("/uscms_data/d3/knash/SlimNtuples/JetHT/crab_JetHT_Run2015D-PromptReco-v4_B2GAnaFW_V8p4_Slim_V10/*/0000/*.root")
		#files = glob.glob("/eos/uscms/store/user/algomez/JetHT/Run2015D-05Oct2015-v1_RUNA_v09/151116_122733/*/*.root")
		#files += glob.glob("/eos/uscms/store/user/algomez/JetHT/Run2015D-PromptReco-v4_RUNA_v09/151117_100001/*/*.root")

 		#files = glob.glob("/uscms_data/d3/knash/WPrime13TeV/B2GAnaFW/SlimNtuples_test/CMSSW_7_4_15/src/Analysis/B2GAnaFW/test/B2GEDMNtuple_slimtest.root")
 		#files = glob.glob("/uscms_data/d3/knash/WPrime13TeV/B2GAnaFW/SlimNtuples_test/CMSSW_7_4_15/src/Analysis/B2GAnaFW/test/B2GEDMNtuple_slimregular.root")


	if string == 'B2GsignalLH1200':
		files = glob.glob("/eos/uscms/store/user/lcorcodi/BstarToTW_M-1200_LH_TuneCUETP8M1_13TeV-madgraph-pythia8/RunIISpring16MiniAODv2_B2GAnaFW_80x_V1p0/160624*/0000/*.root")
	if string == 'B2GsignalLH1400':
		files = glob.glob("/eos/uscms/store/user/lcorcodi/BstarToTW_M-1400_LH_TuneCUETP8M1_13TeV-madgraph-pythia8/RunIISpring16MiniAODv2_B2GAnaFW_80x_V1p0/160624*/0000/*.root")
	if string == 'B2GsignalLH1600':
		files = glob.glob("/eos/uscms/store/user/lcorcodi/BstarToTW_M-1600_LH_TuneCUETP8M1_13TeV-madgraph-pythia8/RunIISpring16MiniAODv2_B2GAnaFW_80x_V1p0/160624*/0000/*.root")
	if string == 'B2GsignalLH1800':
		files = glob.glob("/eos/uscms/store/user/lcorcodi/BstarToTW_M-1800_LH_TuneCUETP8M1_13TeV-madgraph-pythia8/RunIISpring16MiniAODv2_B2GAnaFW_80x_V1p0/160624*/0000/*.root")
	if string == 'B2GsignalLH2000':
		files = glob.glob("/eos/uscms/store/user/lcorcodi/BstarToTW_M-2000_LH_TuneCUETP8M1_13TeV-madgraph-pythia8/RunIISpring16MiniAODv2_B2GAnaFW_80x_V1p0/160624*/0000/*.root")
	if string == 'B2GsignalLH2200':
		files = glob.glob("/eos/uscms/store/user/lcorcodi/BstarToTW_M-2200_LH_TuneCUETP8M1_13TeV-madgraph-pythia8/RunIISpring16MiniAODv2_B2GAnaFW_80x_V1p0/160624*/0000/*.root")
	if string == 'B2GsignalLH2400':
		files = glob.glob("/eos/uscms/store/user/lcorcodi/BstarToTW_M-2400_LH_TuneCUETP8M1_13TeV-madgraph-pythia8/RunIISpring16MiniAODv2_B2GAnaFW_80x_V1p0/160624*/0000/*.root")
	if string == 'B2GsignalLH2600':
		files = glob.glob("/eos/uscms/store/user/lcorcodi/BstarToTW_M-2600_LH_TuneCUETP8M1_13TeV-madgraph-pythia8/RunIISpring16MiniAODv2_B2GAnaFW_80x_V1p0/160624*/0000/*.root")
	if string == 'B2GsignalLH2800':
		files = glob.glob("/eos/uscms/store/user/lcorcodi/BstarToTW_M-2800_LH_TuneCUETP8M1_13TeV-madgraph-pythia8/RunIISpring16MiniAODv2_B2GAnaFW_80x_V1p0/160624*/0000/*.root")
	if string == 'B2GsignalLH3000':
		files = glob.glob("/eos/uscms/store/user/lcorcodi/BstarToTW_M-3000_LH_TuneCUETP8M1_13TeV-madgraph-pythia8/RunIISpring16MiniAODv2_B2GAnaFW_80x_V1p0/160624*/0000/*.root")
	if string == 'ttbar80X':
 		files = glob.glob("/uscms_data/d3/lcorcodi/BStar13TeV/SlimTuples/crab_TT_TuneCUETP8M1_13TeV-powheg-pythia8_B2GAnaFW_V1p1_80x_Slim_V2/*/0000/*.root")

	if string == 'ttbar':
 		files = glob.glob("/uscms_data/d3/knash/SlimNtuples/TT_TuneCUETP8M1_13TeV-powheg-pythia8/crab_TT_TuneCUETP8M1_13TeV-powheg-pythia8_B2GAnaFW_V8p4_Slim_V12/*/0000/*.root")
	if string == 'ttbarscaleup':
 		files = glob.glob("/uscms_data/d3/knash/SlimNtuples/TT_TuneCUETP8M1_13TeV-powheg-scaleup-pythia8/crab_TT_TuneCUETP8M1_13TeV-powheg-scaleup-pythia8_B2GAnaFW_V8p4_Slim_V12/*/0000/*.root")
	if string == 'ttbarscaledown':
 		files = glob.glob("/uscms_data/d3/knash/SlimNtuples/TT_TuneCUETP8M1_13TeV-powheg-scaledown-pythia8/crab_TT_TuneCUETP8M1_13TeV-powheg-scaledown-pythia8_B2GAnaFW_V8p4_Slim_V12/*/0000/*.root")

 #	if string == 'QCDPT300':
 #		files = glob.glob("/uscms_data/d3/knash/SlimNtuples/QCD_Pt_300to470_TuneCUETP8M1_13TeV_pythia8/crab_QCD_Pt_300to470_TuneCUETP8M1_13TeV_pythia8_v74x_V6_25ns/150928_144140/0000/*.root")
 #	if string == 'QCDPT470':
 #		files = glob.glob("/uscms_data/d3/knash/SlimNtuples/QCD_Pt_470to600_TuneCUETP8M1_13TeV_pythia8/crab_QCD_Pt_470to600_TuneCUETP8M1_13TeV_pythia8_v74x_V6_25ns/150928_144250/0000/*.root")
 #	if string == 'QCDPT600':
 #		files = glob.glob("/uscms_data/d3/knash/SlimNtuples/QCD_Pt_600to800_TuneCUETP8M1_13TeV_pythia8/crab_QCD_Pt_600to800_TuneCUETP8M1_13TeV_pythia8_v74x_V6_25ns/150928_144306/0000/*.root")
 #	if string == 'QCDPT800':
 #		files = glob.glob("/uscms_data/d3/knash/SlimNtuples/QCD_Pt_800to1000_TuneCUETP8M1_13TeV_pythia8/crab_QCD_Pt_800to1000_TuneCUETP8M1_13TeV_pythia8_v74x_V6_25ns/150928_200957/0000/*.root")
#	if string == 'QCDPT1000':
 #		files = glob.glob("/uscms_data/d3/knash/SlimNtuples/QCD_Pt_1000to1400_TuneCUETP8M1_13TeV_pythia8/crab_QCD_Pt_1000to1400_TuneCUETP8M1_13TeV_pythia8_v74x_V6_25ns/150928_144156/0000/*.root")
 #	if string == 'QCDPT1400':
 #		files = glob.glob("/uscms_data/d3/knash/SlimNtuples/QCD_Pt_1400to1800_TuneCUETP8M1_13TeV_pythia8/crab_QCD_Pt_1400to1800_TuneCUETP8M1_13TeV_pythia8_v74x_V6_25ns/150928_144216/0000/*.root")
#	if string == 'QCDPT1800':
 #		files = glob.glob("/uscms_data/d3/knash/SlimNtuples/QCD_Pt_1800to2400_TuneCUETP8M1_13TeV_pythia8/crab_QCD_Pt_1800to2400_TuneCUETP8M1_13TeV_pythia8_v74x_V6_25ns/150928_144234/0000/*.root")


	if string == 'WJETS':
 		files = glob.glob("/uscms_data/d3/knash/SlimNtuples/WJetsToQQ_HT-600ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_WJetsToQQ_HT-600ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_V8p4_Slim_V12/160629_214546/0000/*.root")
 

	if string == 'QCDHT300':
 		files = glob.glob("/uscms_data/d3/knash/SlimNtuples/QCD_HT300to500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_QCD_HT300to500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_V8p4_Slim_V12/*/0000/*.root")
 
	if string == 'QCDHT500':
 		files = glob.glob("/uscms_data/d3/knash/SlimNtuples/QCD_HT500to700_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_QCD_HT500to700_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_V8p4_Slim_V12/*/0000/*.root")
 
	if string == 'QCDHT700':
 		files = glob.glob("/uscms_data/d3/knash/SlimNtuples/QCD_HT700to1000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_QCD_HT700to1000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_V8p4_Slim_V11/*/0000/*.root")
 	if string == 'QCDHT1000':
 		files = glob.glob("/uscms_data/d3/knash/SlimNtuples/QCD_HT1000to1500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_QCD_HT1000to1500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_V8p4_Slim_V12/*/0000/*.root")
 	if string == 'QCDHT1500':
 		files = glob.glob("/uscms_data/d3/knash/SlimNtuples/QCD_HT1500to2000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_QCD_HT1500to2000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_V8p4_Slim_V11/*/0000/*.root")
	if string == 'QCDHT2000':
 		files = glob.glob("/uscms_data/d3/knash/SlimNtuples/QCD_HT2000toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_QCD_HT2000toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_V8p4_Slim_V12/*/0000/*.root")
 
  


	#if string == 'STS':
 	#	files = glob.glob("/uscms_data/d3/knash/SlimNtuples/ST_s-channel_4f_leptonDecays_13TeV-amcatnlo-pythia8_TuneCUETP8M1/crab_ST_s-channel_4f_leptonDecays_13TeV-amcatnlo-pythia8_TuneCUETP8M1_V8p4_Slim_V12/160113_205221/0000/*.root")
	#if string == 'STT':
 	#	files = glob.glob("/uscms_data/d3/knash/SlimNtuples/ST_t-channel_4f_leptonDecays_13TeV-amcatnlo-pythia8_TuneCUETP8M1/crab_ST_t-channel_4f_leptonDecays_13TeV-amcatnlo-pythia8_TuneCUETP8M1_V8p4_Slim_V12/160113_205147/0000/*.root")
 	#	files += glob.glob("/uscms_data/d3/knash/SlimNtuples/ST_t-channel_4f_leptonDecays_13TeV-amcatnlo-pythia8_TuneCUETP8M1/crab_ST_t-channel_4f_leptonDecays_13TeV-amcatnlo-pythia8_TuneCUETP8M1_ext1_V8p4_Slim_V12/160114_161134/0000/*.root")
	if string == 'STTW':
 		files = glob.glob("/uscms_data/d3/knash/SlimNtuples/ST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/crab_ST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1_V8p4_Slim_V12/*/0000/*.root")
	if string == 'STTWB':
 		files = glob.glob("/uscms_data/d3/knash/SlimNtuples/ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/crab_ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1_V8p4_Slim_V12/*/0000/*.root")

	if string.find('signalright')!=-1:
		sigmass = string.replace('signalright','')
		if string == 'signalright1000' or string == 'signalright1100': 
			files = glob.glob("/uscms_data/d3/knash/SlimNtuples/WprimeToTB_TToHad_M-"+sigmass+"_RH_TuneCUETP8M1_13TeV-comphep-pythia8/crab_WPrime13TeV_RH_"+sigmass+"_B2GAnaFW_V8p4_Slim_V12/*/*/*.root")
		else:
			files = glob.glob("/uscms_data/d3/knash/SlimNtuples/WprimeToTB_TToHad_M-"+sigmass+"_RH_TuneCUETP8M1_13TeV-comphep-pythia8/crab_WPrime13TeV_RH_"+sigmass+"_B2GAnaFW_V8p4_Slim_V12/*/*/*.root")


	#if string == 'signalleft2800':
	#	files = glob.glob("/eos/uscms/store/user/knash/Wprimeleft280013TeV_pythia8/crab_WPrime13TeV_LeftCut_Fullsim_B2GAnaFW_V6_2800_LH_25ns/151123_180958/0000/*.root")
	for i in range(0,len(files)):
		files[i] = files[i].replace('/eos/uscms','root://cmsxrootd.fnal.gov//')

	try:
		print 'A total of ' + str(len(files)) + ' files'
	except:
		print 'Bad files option'
	return files




#This function initializes the average b tagging rates used for QCD determination
#It tages the type of functional form as an argument.  The default fit is Bifpoly





def BTR_Init(ST,CUT,di,setval):

	if setval.find("QCD")==-1:
		setval = "data"


	if ST == 'Bifpoly':
		TRBPE1 = open("./"+di+"fitdata/bpinput"+setval+"eta1_PSET_"+CUT+".txt")
		TRBPE1.seek(0)
		TRBPE2 = open("./"+di+"fitdata/bpinput"+setval+"eta2_PSET_"+CUT+".txt")
		TRBPE2.seek(0)
		TRBPE3 = open("./"+di+"fitdata/bpinput"+setval+"eta3_PSET_"+CUT+".txt")
		TRBPE3.seek(0)
		eta1fit = TF1("eta1fit",BifPoly,0,2000,5)
		eta2fit = TF1("eta2fit",BifPoly,0,2000,5)
		eta3fit = TF1("eta3fit",BifPoly,0,2000,5)
		Params = 5
	if ST == 'Bifpoly_err':
		TRBPE1 = open("./"+di+"fitdata/bperrorinput"+setval+"eta1_PSET_"+CUT+".txt")
		TRBPE1.seek(0)
		TRBPE2 = open("./"+di+"fitdata/bperrorinput"+setval+"eta2_PSET_"+CUT+".txt")
		TRBPE2.seek(0)
		TRBPE3 = open("./"+di+"fitdata/bperrorinput"+setval+"eta3_PSET_"+CUT+".txt")
		TRBPE3.seek(0)
		eta1fit=TF1("eta1fit",BifPolyErr,0,2000,10)
		eta2fit=TF1("eta2fit",BifPolyErr,0,2000,10)
		eta3fit=TF1("eta3fit",BifPolyErr,0,2000,10)
		Params = 10

	if ST == 'pol0':
		TRBPE1 = open("./"+di+"fitdata/pol0input"+setval+"eta1_PSET_"+CUT+".txt")
		TRBPE1.seek(0)
		TRBPE2 = open("./"+di+"fitdata/pol0input"+setval+"eta2_PSET_"+CUT+".txt")
		TRBPE2.seek(0)
		TRBPE3 = open("./"+di+"fitdata/pol0input"+setval+"eta3_PSET_"+CUT+".txt")
		TRBPE3.seek(0)
		eta1fit = TF1("eta1fit",'pol0',0,2000)
		eta2fit = TF1("eta2fit",'pol0',0,2000)
		eta3fit = TF1("eta3fit",'pol0',0,2000)
		Params = 1

	if ST == 'pol2':
		TRBPE1 = open("./"+di+"fitdata/pol2input"+setval+"eta1_PSET_"+CUT+".txt")
		TRBPE1.seek(0)
		TRBPE2 = open("./"+di+"fitdata/pol2input"+setval+"eta2_PSET_"+CUT+".txt")
		TRBPE2.seek(0)
		TRBPE3 = open("./"+di+"fitdata/pol2input"+setval+"eta3_PSET_"+CUT+".txt")
		TRBPE3.seek(0)
		eta1fit = TF1("eta1fit",'pol2',0,2000)
		eta2fit = TF1("eta2fit",'pol2',0,2000)
		eta3fit = TF1("eta3fit",'pol2',0,2000)
		Params = 3

	if ST == 'pol3':
		TRBPE1 = open("./"+di+"fitdata/pol3input"+setval+"eta1_PSET_"+CUT+".txt")
		TRBPE1.seek(0)
		TRBPE2 = open("./"+di+"fitdata/pol3input"+setval+"eta2_PSET_"+CUT+".txt")
		TRBPE2.seek(0)
		TRBPE3 = open("./"+di+"fitdata/pol3input"+setval+"eta3_PSET_"+CUT+".txt")
		TRBPE3.seek(0)
		eta1fit = TF1("eta1fit",'pol3',0,2000)
		eta2fit = TF1("eta2fit",'pol3',0,2000)
		eta3fit = TF1("eta3fit",'pol3',0,2000)
		Params = 4
	if ST == 'FIT':
		TRBPE1 = open("./"+di+"fitdata/newfitinput"+setval+"eta1_PSET_"+CUT+".txt")
		TRBPE1.seek(0)
		TRBPE2 = open("./"+di+"fitdata/newfitinput"+setval+"eta2_PSET_"+CUT+".txt")
		TRBPE2.seek(0)
		TRBPE3 = open("./"+di+"fitdata/newfitinput"+setval+"eta3_PSET_"+CUT+".txt")
		TRBPE3.seek(0)
		eta1fit = TF1("eta1fit",'[0]*([1]+x)/([2]+x)+[3]*x',0,2000)
		eta2fit = TF1("eta2fit",'[0]*([1]+x)/([2]+x)+[3]*x',0,2000)
		eta3fit = TF1("eta3fit",'[0]*([1]+x)/([2]+x)+[3]*x',0,2000)
		Params = 4
	if ST == 'expofit':
		TRBPE1 = open("./"+di+"fitdata/expoconinput"+setval+"eta1_PSET_"+CUT+".txt")
		TRBPE1.seek(0)
		TRBPE2 = open("./"+di+"fitdata/expoconinput"+setval+"eta2_PSET_"+CUT+".txt")
		TRBPE2.seek(0)
		TRBPE3 = open("./"+di+"fitdata/expoconinput"+setval+"eta3_PSET_"+CUT+".txt")
		TRBPE3.seek(0)
		eta1fit = TF1("eta1fit",'expo(0) + pol0(2)',0,2000)
		eta2fit = TF1("eta2fit",'expo(0) + pol0(2)',0,2000)
		eta3fit = TF1("eta3fit",'expo(0) + pol0(2)',0,2000)
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

	ptminsfb = [300]
	ptmaxsfb = [670]	
	SFb = TFormula("SFb","0.908299+(2.70877e-06*(log(x+370.144)*(log(x+370.144)*(3-(-(104.614*log(x+370.144)))))))")
	#SFb = TFormula("SFb","0.901434*((1.+(0.0852659*x))/(1.+(0.0769021*x)))" )
	SFb_down =	[
			TFormula("SFb_down_1","0.908299+((2.70877e-06*(log(x+370.144)*(log(x+370.144)*(3-(-(104.614*log(x+370.144)))))))-0.047329759448766708)")
		#	TFormula("SFb_down_1", "(0.901434*((1.+(0.0852659*x))/(1.+(0.0769021*x))))-0.11035951972007751") 
			#TFormula("SFb_down_1","-(0.0443172)+((0.00496634*(log(x+1267.85)*(log(x+1267.85)*(3-(-(0.110428*log(x+1267.85)))))))-0.064764265418052673)") 
			]
	SFb_up =	[
			TFormula("SFb_up_1","(0.908299+(2.70877e-06*(log(x+370.144)*(log(x+370.144)*(3-(-(104.614*log(x+370.144))))))))+0.047329759448766708")			#TFormula("SFb_up_1","(-(0.0443172)+(0.00496634*(log(x+1267.85)*(log(x+1267.85)*(3-(-(0.110428*log(x+1267.85))))))))+0.064764265418052673")
		#	TFormula("SFb_up_1","(0.901434*((1.+(0.0852659*x))/(1.+(0.0769021*x))))+0.11035951972007751")
			]

	if Y <= 670.0:
		weightSFb  = SFb.Eval(Y)
		if Y <= ptmaxsfb[0]:
			weightSFb_down = SFb_down[0].Eval(Y)
			weightSFb_up = SFb_up[0].Eval(Y)
	else: 
		weightSFb  = SFb.Eval(670.0)

		deltaSFb_down = SFb_down[0].Eval(670.0) - weightSFb
		deltaSFb_up = SFb_up[0].Eval(670.0) - weightSFb

		weightSFb_down = SFb.Eval(670.0) + 2*deltaSFb_down
		weightSFb_up = SFb.Eval(670.0) + 2*deltaSFb_up

	return [weightSFb,weightSFb_down,weightSFb_up]
#This looks up the PDF uncertainty
def SFT_Lookup( pttop ):

	ttagsf = [[0.88,0.11],[1.00,0.23]]
	ttagsfregions = [[0,550],[550,float("inf")]]

	for ipttop in range(0,len(ttagsfregions)):
		if ttagsfregions[ipttop][0]<pttop<=ttagsfregions[ipttop][1]:
			return [ttagsf[ipttop][0],ttagsf[ipttop][0]-ttagsf[ipttop][1],ttagsf[ipttop][0]+ttagsf[ipttop][1]]

#def PDF_Lookup( pdfs , pdfOP ):
	#print pdfs
	#for pd in pdfs:
#		print pd
#	iweight = 0.0

 #       if pdfOP == "up" :
  #     		for pdf in pdfs[1::2] :
   #           		iweight = iweight + pdf*pdf
    #    else :
     #   	for pdf in pdfs[2::2] :
      #  		iweight = iweight + pdf*pdf
	##print sqrt((iweight) / (len(pdfs)-1) * 2.0)
        #return sqrt((iweight) / (len(pdfs)-1) * 2.0)
def PDF_LookupMAX( pdfs , pdfOP ):
	#print pdfs
	#for pd in pdfs:
#		print pd
	iweight=1.0
        if pdfOP == "up" :
       		for pdf in pdfs :
              		iweight = max(iweight,pdf)

        else :
        	for pdf in pdfs :
              		iweight = min(iweight,pdf)

        return iweight

def PDF_LookupAVE( pdfs , pdfOP ):
	iweight = 0.0
        if pdfOP == "up" :
       		for pdf in pdfs[1::2] :
              		iweight = iweight + pdf
        else :
        	for pdf in pdfs[2::2] :
        		iweight = iweight + pdf
        return (iweight) / (len(pdfs)-1) * 2.0
def PDF_Lookup( pdfs , pdfOP ):
	iweight = 0.0
	#print "LEN"
	#print len(pdfs)
	ave =  pdfs
	ave =  reduce(lambda x, y: x + y, pdfs) / len(pdfs)
	#print ave
       	for pdf in pdfs :
             	iweight = iweight + (pdf-ave)*(pdf-ave)

        if pdfOP == "up" :
        	return 1+sqrt((iweight) / (len(pdfs)))
        else :
          	return 1-sqrt((iweight) / (len(pdfs)))
#This looks up the b tagging scale factor or uncertainty
def Trigger_Lookup( H , TRP):
        Weight = 1.0
	Weightup = 1.0
	Weightdown = 1.0
        if H < 1200.0:
                bin0 = TRP.FindBin(H) 
                jetTriggerWeight = TRP.GetBinContent(bin0)
                Weight = jetTriggerWeight

		deltaTriggerEff  = 0.5*(1.0-jetTriggerWeight)
                Weightup  =   min(1.0,jetTriggerWeight + deltaTriggerEff)
                Weightdown  =   max(0.0,jetTriggerWeight - deltaTriggerEff)
		
	return [Weight,Weightup,Weightdown]

def Trigger_Pass(tnamestr,trigs,bits):
	###TAKE OUT!
	#tnamestr = ['HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV0p45_v3','HLT_PFHT800_v2']
	foundname=False
	for t in range(0,len(trigs)):
		for tname in tnamestr:	
			if trigs[t]==tname :
				foundname=True
			if trigs[t]==tname and bits[t] :
				return True
	if foundname==False:
		print "TRIGGER NOT IN FILE"
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

def Initlv(string,post=''):
	PtHandle 	= 	Handle (  "vector<float> "  )
	PtLabel  	= 	( string+post , string.replace("jets","jet")+"Pt")

	EtaHandle 	= 	Handle (  "vector<float> "  )
	EtaLabel  	= 	( string+post , string.replace("jets","jet")+"Eta")

	PhiHandle 	= 	Handle (  "vector<float> "  )
	PhiLabel  	= 	( string+post , string.replace("jets","jet")+"Phi")

	MassHandle 	= 	Handle (  "vector<float> "  )
	MassLabel  	= 	( string+post , string.replace("jets","jet")+"Mass")

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

		if abs(Math.VectorUtil.DeltaPhi(LV1[0],LV1[iLV1]))> TMath.Pi()/2:
			tjets[1].append(iLV1)
		else:
			tjets[0].append(iLV1)

	for iLV2 in range(0,len(LV2)):
		if abs(Math.VectorUtil.DeltaPhi(LV1[0],LV2[iLV2]))> TMath.Pi()/2:
			bjets[1].append(iLV2)
		else:
			bjets[0].append(iLV2)
	return tjets,bjets


#Some lazy string formatting functions 
def strf( x ):
	return '%.2f' % x

def strf1( x ):
	return '%.0f' % x
def strf2( x ):
	return '%.1f' % x
