
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


import os
import array
import glob
import math
import ROOT
import sys
from array import *
from ROOT import *

#This is the most impostant Function.  Correct information here is essential to obtaining valid results.
#In order we have Luminosity, top tagging scale factor, cross sections for wprime right,left,mixed,ttbar,qcd, and singletop and their corresponding event numbers
#If I wanted to access the left handed W' cross section at 1900 GeV I could do Xsecl1900 = LoadConstants()['xsec_wpl']['1900']
def LoadConstants():
	 return  {
		'lumi':19757,
		'ttagsf':1.036,
		'xsec_wpr':{'1700': 0.12460496399999998, '2200': 0.021867502800000001, '2100': 0.030537407999999995, '1300': 0.58233251999999991, '1800': 0.086716607999999987, '1600': 0.18060609599999999, '2000': 0.042946991999999996, '2700': 0.0046708991999999993, '1500': 0.26379143999999999, '2900': 0.0027131543999999994, '3000': 0.0021096041999999998, '2400': 0.011475855599999999, '2300': 0.015779834399999998, '2600': 0.0062373959999999992, '2800': 0.0035367815999999995, '1400': 0.38991875999999992, '1900': 0.060816491999999993, '2500': 0.0084188544000000001},
		'xsec_wpl':{'1700': 0.20028183251999998, '2200': 0.12110620368, '2100': 0.12566885759999999, '1300': 0.50971714199999996, '1800': 0.1740819102, '1600': 0.25585124993999997, '2000': 0.13701113712000001, '2700': 0.11197266959999999, '1500': 0.31458947519999991, '2900': 0.11225782920000001, '3000': 0.11240040900000001, '2400': 0.11543941529999999, '2300': 0.11660361635999998, '2600': 0.11313865133999999, '2800': 0.11213426004, '1400': 0.40551726599999999, '1900': 0.14849230967999999, '2500': 0.11430079584},
		'xsec_wplr':{'1700': 0.32403986999999995, '2200': 0.14765102615999998, '2100': 0.16013479679999998, '1300': 1.1192727479999998, '1800': 0.26698274459999999, '1600': 0.43940786999999992, '2000': 0.18624247487999995, '2700': 0.11840769599999998, '1500': 0.56915749439999985, '2900': 0.11684897399999998, '3000': 0.11516685839999997, '2400': 0.13107319631999997, '2300': 0.13539302568, '2600': 0.12258562272, '2800': 0.11761920587999998, '1400': 0.81357839639999996, '1900': 0.21316229351999999, '2500': 0.12686314211999999},
		'xsec_ttbar':{'700':245.8*1.23*0.074,'1000':245.8*1.23*0.014},
		'xsec_qcd':{'300':1759.6,'470':113.9,'600':27.0,'800':3.57,'1000':0.734,'1400':0.03352235},
		'xsec_st':{'s':3.79,'sB':1.76,'t':56.4,'tB':30.7,'tW':11.1,'tWB':11.1},
		'nev_wpr':{'1300':489658,'1500':583622,'1700':582074,'1900':579561,'2100':581051,'2300':580410,'2700':581035,'3100':486557},
		'nev_wpl':{'1300':474565,'1500':483975,'1700':464047,'1900':447300,'2100':474159,'2300':478510,'2700':469821,'3100':476969},
		'nev_wplr':{'1300':468663,'1500':465828,'1700':456763,'1900':476482,'2100':480806,'2300':486335,'2700':480485,'3100':462627},
		'nev_ttbar':{'700':3058076,'1000':1233739,'700scaleup':2225727,'1000scaleup':1225662,'700scaledown':2153111,'1000scaledown':1292980},
		'nev_qcd':{'300':5908205,'470':3919113,'600':3902030,'800':3881338,'1000':1895936,'1400':1912782},
		'nev_st':{'s':259176,'sB':139604,'t':3748155,'tB':1930185,'tW':495559,'tWB':491463},
		}

#This is also a very impostant Function.  The analysis runs on "PSETS", which correspond to the TYPE variable here.
#These each load a cut profile.  For instance 'default' is the standard selection used to set limits
def LoadCuts(TYPE):
	if TYPE=='default':
 		return  {
			'bpt':[370.0,float("inf")],
			'tpt':[450.0,float("inf")],
			'dy':[0.0,1.6],
			'tmass':[140.0,250.0],
			'nsubjets':[3,10],
			'tau32':[0.0,0.55],
			'minmass':[50.0,float("inf")],
			'sjbtag':[0.679,1.0],
			'bmass':[0.0,70.0],
			'btag':[0.679,1.0],
			'eta1':[0.0,0.5],
			'eta2':[0.5,1.15],	
			'eta3':[1.15,2.4]
			}
	if TYPE=='rate_default':
 		return  {
			'bpt':[370.0,float("inf")],
			'tpt':[450.0,float("inf")],
			'dy':[0.0,1.6],
			'tmass':[140.0,250.0],
			'nsubjets':[0,3],
			'tau32':[0.0,1.0],
			'minmass':[0.0,float("inf")],
			'sjbtag':[0.679,1.0],
			'bmass':[0.0,70.0],
			'btag':[0.679,1.0],
			'eta1':[0.0,0.5],
			'eta2':[0.5,1.15],	
			'eta3':[1.15,2.4]
			}

	if TYPE=='sideband1':
 		return  {
			'bpt':[370.0,float("inf")],
			'tpt':[450.0,float("inf")],
			'dy':[0.0,1.6],
			'tmass':[140.0,250.0],
			'nsubjets':[3,10],
			'tau32':[0.55,1.0],
			'minmass':[0.0,50.0],
			'sjbtag':[0.679,1.0],
			'bmass':[0.0,70.0],
			'btag':[0.679,1.0],
			'eta1':[0.0,0.5],
			'eta2':[0.5,1.15],	
			'eta3':[1.15,2.4]
			}
	if TYPE=='sideband2':
 		return  {
			'bpt':[370.0,float("inf")],
			'tpt':[450.0,float("inf")],
			'dy':[0.0,1.6],
			'tmass':[140.0,250.0],
			'nsubjets':[3,10],
			'tau32':[0.0,0.55],
			'minmass':[50.0,float("inf")],
			'sjbtag':[0.0,0.679],
			'bmass':[0.0,70.0],
			'btag':[0.679,1.0],
			'eta1':[0.0,0.5],
			'eta2':[0.5,1.15],	
			'eta3':[1.15,2.4]
			}

	if TYPE=='sideband3':
 		return  {
			'bpt':[370.0,float("inf")],
			'tpt':[450.0,float("inf")],
			'dy':[0.0,1.6],
			'tmass':[140.0,250.0],
			'nsubjets':[3,10],
			'tau32':[0.0,0.55],
			'minmass':[50.0,float("inf")],
			'sjbtag':[0.679,1.0],
			'bmass':[70.0,float("inf")],
			'btag':[0.679,1.0],
			'eta1':[0.0,0.5],
			'eta2':[0.5,1.15],	
			'eta3':[1.15,2.4]
			}

#This function loads up Ntuples based on what type of set you want to analyze.  
#This needs to be updated whenever new Ntuples are produced (unless the file locations are the same).
def Load_Ntuples(string):
	print 'running on ' + string 
	if string == 'data':
		files = glob.glob("/uscms_data/d3/knash/WPrime8TeV/data/CMSSW_5_3_18/src/Analysis/TTBSMPatTuples/test/Run2012A-22Jan2013/res/*.root")
		files += glob.glob("/uscms_data/d3/knash/WPrime8TeV/data/CMSSW_5_3_18/src/Analysis/TTBSMPatTuples/test/Run2012B-22Jan2013/res/*.root")
		files += glob.glob("/uscms_data/d3/knash/WPrime8TeV/data/CMSSW_5_3_18/src/Analysis/TTBSMPatTuples/test/Run2012C-22Jan2013/res/*.root")
		files += glob.glob("/uscms_data/d3/knash/WPrime8TeV/data/CMSSW_5_3_18/src/Analysis/TTBSMPatTuples/test/Run2012D-22Jan2013/res/*.root")
	if string == 'QCDHT1000':
		files = glob.glob("/uscms_data/d3/knash/WPrime8TeV/CMSSW_5_3_18/src/Analysis/TTBSMPatTuples/test/QCD_HT1000_to_HTinf/res/*.root" )
	if string == 'ttbar700':
		files = glob.glob("/uscms_data/d3/knash/WPrime8TeV/CMSSW_5_3_18/src/Analysis/TTBSMPatTuples/test/ttbar_Mtt-700to1000/res/*.root" )
	if string == 'ttbar1000':
		files = glob.glob("/uscms_data/d3/knash/WPrime8TeV/CMSSW_5_3_18/src/Analysis/TTBSMPatTuples/test/ttbar_Mtt-1000toinf/res/*.root" )
	if string == 'ttbar700scaleup':
		files = glob.glob("/uscms_data/d3/knash/WPrime8TeV/CMSSW_5_3_18/src/Analysis/TTBSMPatTuples/test/ttbar_Mtt-700to1000scaleup/res/*.root" )
	if string == 'ttbar700scaledown':
		files = glob.glob("/uscms_data/d3/knash/WPrime8TeV/CMSSW_5_3_18/src/Analysis/TTBSMPatTuples/test/ttbar_Mtt-700to1000scaledown/res/*.root" )
	if string == 'ttbar1000scaleup':
		files = glob.glob("/uscms_data/d3/knash/WPrime8TeV/CMSSW_5_3_18/src/Analysis/TTBSMPatTuples/test/ttbar_Mtt-1000toInfscaleup/res/*.root")
	if string == 'ttbar1000scaledown':
		files = glob.glob("/uscms_data/d3/knash/WPrime8TeV/CMSSW_5_3_18/src/Analysis/TTBSMPatTuples/test/ttbar_Mtt-1000toInfscaledown/res/*.root" )
	if string == 'singletop_s':
		files = glob.glob("/uscms_data/d3/knash/WPrime8TeV/CMSSW_5_3_18/src/Analysis/TTBSMPatTuples/test/singletop_s/res/*.root" )
	if string == 'singletop_sB':
		files = glob.glob("/uscms_data/d3/knash/WPrime8TeV/CMSSW_5_3_18/src/Analysis/TTBSMPatTuples/test/singletop_sB/res/*.root" )
	if string == 'singletop_t':
		files = glob.glob("/uscms_data/d3/knash/WPrime8TeV/CMSSW_5_3_18/src/Analysis/TTBSMPatTuples/test/singletop_t/res/*.root" )
	if string == 'singletop_tB':
		files = glob.glob("/uscms_data/d3/knash/WPrime8TeV/CMSSW_5_3_18/src/Analysis/TTBSMPatTuples/test/singletop_tB/res/*.root" )
	if string == 'singletop_tW':
		files = glob.glob("/uscms_data/d3/knash/WPrime8TeV/CMSSW_5_3_18/src/Analysis/TTBSMPatTuples/test/singletop_tW/res/*.root" )
	if string == 'singletop_tWB':
		files = glob.glob("/uscms_data/d3/knash/WPrime8TeV/CMSSW_5_3_18/src/Analysis/TTBSMPatTuples/test/singletop_tWB/res/*.root" )

	if string == 'signalright1300':
		files = glob.glob("/uscms_data/d3/knash/WPrime8TeV/CMSSW_5_3_18/src/Analysis/TTBSMPatTuples/test/SingletopWprimeTToHad_M-1300/res/*.root" )
	if string == 'signalright1500':
		files = glob.glob("/uscms_data/d3/knash/WPrime8TeV/CMSSW_5_3_18/src/Analysis/TTBSMPatTuples/test/SingletopWprimeTToHad_M-1500/res/*.root" )
	if string == 'signalright1700':
		files = glob.glob("/uscms_data/d3/knash/WPrime8TeV/CMSSW_5_3_18/src/Analysis/TTBSMPatTuples/test/SingletopWprimeTToHad_M-1700/res/*.root" )
	if string == 'signalright1900':
		files = glob.glob("/uscms_data/d3/knash/WPrime8TeV/CMSSW_5_3_18/src/Analysis/TTBSMPatTuples/test/SingletopWprimeTToHad_M-1900/res/*.root" )
	if string == 'signalright2100':
		files = glob.glob("/uscms_data/d3/knash/WPrime8TeV/CMSSW_5_3_18/src/Analysis/TTBSMPatTuples/test/SingletopWprimeTToHad_M-2100/res/*.root" )
	if string == 'signalright2300':
		files = glob.glob("/uscms_data/d3/knash/WPrime8TeV/CMSSW_5_3_18/src/Analysis/TTBSMPatTuples/test/SingletopWprimeTToHad_M-2300/res/*.root" )
	if string == 'signalright2700':
		files = glob.glob("/uscms_data/d3/knash/WPrime8TeV/CMSSW_5_3_18/src/Analysis/TTBSMPatTuples/test/SingletopWprimeTToHad_M-2700/res/*.root" )
	if string == 'signalright3100':
		files = glob.glob("/uscms_data/d3/knash/WPrime8TeV/CMSSW_5_3_18/src/Analysis/TTBSMPatTuples/test/SingletopWprimeTToHad_M-3100/res/*.root" )

	if string == 'signalleft1300':
		files = glob.glob("/uscms_data/d3/knash/WPrime8TeV/CMSSW_5_3_18/src/Analysis/TTBSMPatTuples/test/SingletopWprimeTToHad_M-1300_left/res/*.root" )
	if string == 'signalleft1500':
		files = glob.glob("/uscms_data/d3/knash/WPrime8TeV/CMSSW_5_3_18/src/Analysis/TTBSMPatTuples/test/SingletopWprimeTToHad_M-1500_left/res/*.root" )
	if string == 'signalleft1700':
		files = glob.glob("/uscms_data/d3/knash/WPrime8TeV/CMSSW_5_3_18/src/Analysis/TTBSMPatTuples/test/SingletopWprimeTToHad_M-1700_left/res/*.root" )
	if string == 'signalleft1900':
		files = glob.glob("/uscms_data/d3/knash/WPrime8TeV/CMSSW_5_3_18/src/Analysis/TTBSMPatTuples/test/SingletopWprimeTToHad_M-1900_left/res/*.root" )
	if string == 'signalleft2100':
		files = glob.glob("/uscms_data/d3/knash/WPrime8TeV/CMSSW_5_3_18/src/Analysis/TTBSMPatTuples/test/SingletopWprimeTToHad_M-2100_left/res/*.root" )
	if string == 'signalleft2300':
		files = glob.glob("/uscms_data/d3/knash/WPrime8TeV/CMSSW_5_3_18/src/Analysis/TTBSMPatTuples/test/SingletopWprimeTToHad_M-2300_left/res/*.root" )
	if string == 'signalleft2700':
		files = glob.glob("/uscms_data/d3/knash/WPrime8TeV/CMSSW_5_3_18/src/Analysis/TTBSMPatTuples/test/SingletopWprimeTToHad_M-2700_left/res/*.root" )
	if string == 'signalleft3100':
		files = glob.glob("/uscms_data/d3/knash/WPrime8TeV/CMSSW_5_3_18/src/Analysis/TTBSMPatTuples/test/SingletopWprimeTToHad_M-3100_left/res/*.root" )
	

	if string == 'signalmixed1300':
		files = glob.glob("/uscms_data/d3/knash/WPrime8TeV/CMSSW_5_3_18/src/Analysis/TTBSMPatTuples/test/SingletopWprimeTToHad_M-1300_mixed/res/*.root" )
	if string == 'signalmixed1500':
		files = glob.glob("/uscms_data/d3/knash/WPrime8TeV/CMSSW_5_3_18/src/Analysis/TTBSMPatTuples/test/SingletopWprimeTToHad_M-1500_mixed/res/*.root" )
	if string == 'signalmixed1700':
		files = glob.glob("/uscms_data/d3/knash/WPrime8TeV/CMSSW_5_3_18/src/Analysis/TTBSMPatTuples/test/SingletopWprimeTToHad_M-1700_mixed/res/*.root" )
	if string == 'signalmixed1900':
		files = glob.glob("/uscms_data/d3/knash/WPrime8TeV/CMSSW_5_3_18/src/Analysis/TTBSMPatTuples/test/SingletopWprimeTToHad_M-1900_mixed/res/*.root" )
	if string == 'signalmixed2100':
		files = glob.glob("/uscms_data/d3/knash/WPrime8TeV/CMSSW_5_3_18/src/Analysis/TTBSMPatTuples/test/SingletopWprimeTToHad_M-2100_mixed/res/*.root" )
	if string == 'signalmixed2300':
		files = glob.glob("/uscms_data/d3/knash/WPrime8TeV/CMSSW_5_3_18/src/Analysis/TTBSMPatTuples/test/SingletopWprimeTToHad_M-2300_mixed/res/*.root" )
	if string == 'signalmixed2700':
		files = glob.glob("/uscms_data/d3/knash/WPrime8TeV/CMSSW_5_3_18/src/Analysis/TTBSMPatTuples/test/SingletopWprimeTToHad_M-2700_mixed/res/*.root" )
	if string == 'signalmixed3100':
		files = glob.glob("/uscms_data/d3/knash/WPrime8TeV/CMSSW_5_3_18/src/Analysis/TTBSMPatTuples/test/SingletopWprimeTToHad_M-3100_mixed/res/*.root" )
	
	try:
		print 'A total of ' + str(len(files)) + ' files'
	except:
		print 'Bad files option'
	return files




#This function initializes the average b tagging rates used for QCD determination
#It tages the type of functional form as an argument.  The default fit is Bifpoly

#This is a poorly written function, but I cant think of a better way to do this 
#It works, but you should be able to just have one input
def BTR_Init(ST,CUT,di):
	if ST == 'Bifpoly':
		TRBPE1 = open(di+"fitdata/bpinputeta1_PSET_"+CUT+".txt")
		TRBPE1.seek(0)
		TRBPE2 = open(di+"fitdata/bpinputeta2_PSET_"+CUT+".txt")
		TRBPE2.seek(0)
		TRBPE3 = open(di+"fitdata/bpinputeta3_PSET_"+CUT+".txt")
		TRBPE3.seek(0)
		eta1fit = TF1("eta1fit",BifPoly,0,1400,5)
		eta2fit = TF1("eta2fit",BifPoly,0,1400,5)
		eta3fit = TF1("eta3fit",BifPoly,0,1400,5)
		Params = 5
	if ST == 'Bifpoly_err':
		TRBPE1 = open(di+"fitdata/bperrorinputeta1_PSET_"+CUT+".txt")
		TRBPE1.seek(0)
		TRBPE2 = open(di+"fitdata/bperrorinputeta2_PSET_"+CUT+".txt")
		TRBPE2.seek(0)
		TRBPE3 = open(di+"fitdata/bperrorinputeta3_PSET_"+CUT+".txt")
		TRBPE3.seek(0)
		eta1fit=TF1("eta1fit",BifPolyErr,0,1400,10)
		eta2fit=TF1("eta2fit",BifPolyErr,0,1400,10)
		eta3fit=TF1("eta3fit",BifPolyErr,0,1400,10)
		Params = 10

	if ST == 'pol0':
		TRBPE1 = open(di+"fitdata/pol0inputeta1_PSET_"+CUT+".txt")
		TRBPE1.seek(0)
		TRBPE2 = open(di+"fitdata/pol0inputeta2_PSET_"+CUT+".txt")
		TRBPE2.seek(0)
		TRBPE3 = open(di+"fitdata/pol0inputeta3_PSET_"+CUT+".txt")
		TRBPE3.seek(0)
		eta1fit = TF1("eta1fit",'pol0',0,1400)
		eta2fit = TF1("eta2fit",'pol0',0,1400)
		eta3fit = TF1("eta3fit",'pol0',0,1400)
		Params = 1

	if ST == 'pol2':
		TRBPE1 = open(di+"fitdata/pol2inputeta1_PSET_"+CUT+".txt")
		TRBPE1.seek(0)
		TRBPE2 = open(di+"fitdata/pol2inputeta2_PSET_"+CUT+".txt")
		TRBPE2.seek(0)
		TRBPE3 = open(di+"fitdata/pol2inputeta3_PSET_"+CUT+".txt")
		TRBPE3.seek(0)
		eta1fit = TF1("eta1fit",'pol2',0,1400)
		eta2fit = TF1("eta2fit",'pol2',0,1400)
		eta3fit = TF1("eta3fit",'pol2',0,1400)
		Params = 3

	if ST == 'pol3':
		TRBPE1 = open(di+"fitdata/pol3inputeta1_PSET_"+CUT+".txt")
		TRBPE1.seek(0)
		TRBPE2 = open(di+"fitdata/pol3inputeta2_PSET_"+CUT+".txt")
		TRBPE2.seek(0)
		TRBPE3 = open(di+"fitdata/pol3inputeta3_PSET_"+CUT+".txt")
		TRBPE3.seek(0)
		eta1fit = TF1("eta1fit",'pol3',0,1400)
		eta2fit = TF1("eta2fit",'pol3',0,1400)
		eta3fit = TF1("eta3fit",'pol3',0,1400)
		Params = 4
	if ST == 'FIT':
		TRBPE1 = open(di+"fitdata/newfitinputeta1_PSET_"+CUT+".txt")
		TRBPE1.seek(0)
		TRBPE2 = open(di+"fitdata/newfitinputeta2_PSET_"+CUT+".txt")
		TRBPE2.seek(0)
		TRBPE3 = open(di+"fitdata/newfitinputeta3_PSET_"+CUT+".txt")
		TRBPE3.seek(0)
		eta1fit = TF1("eta1fit",'[0]*([1]+x)/([2]+x)+[3]*x',0,1400)
		eta2fit = TF1("eta2fit",'[0]*([1]+x)/([2]+x)+[3]*x',0,1400)
		eta3fit = TF1("eta3fit",'[0]*([1]+x)/([2]+x)+[3]*x',0,1400)
		Params = 4
	if ST == 'expofit':
		TRBPE1 = open(di+"fitdata/expoconinputeta1_PSET_"+CUT+".txt")
		TRBPE1.seek(0)
		TRBPE2 = open(di+"fitdata/expoconinputeta2_PSET_"+CUT+".txt")
		TRBPE2.seek(0)
		TRBPE3 = open(di+"fitdata/expoconinputeta3_PSET_"+CUT+".txt")
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
		eta2fit.SetParameter(i,float(TBP1.split('\n')[i]) )
		eta3fit.SetParameter(i,float(TBP1.split('\n')[i]) )

	return [eta1fit.Clone(),eta2fit.Clone(),eta3fit.Clone()] 

#This takes the average b tagging rates that are initialized in the above function and produces 
#A QCD background estimate based on them 
def bkg_weight(blv, funcs, etabins):
	for ibin in range(0,len(etabins)):
		if (etabins[ibin][0] <= abs(blv.eta()) < etabins[ibin][1]) :
			tagratept = funcs[ibin].Eval(blv.pt())		
	return tagratept

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
def Trigger_Lookup( H , TRP , TROP ):
        Weight = 1.0
        if H < 1300.0:
                bin0 = TRP.FindBin(H) 
                jetTriggerWeight = TRP.GetBinContent(bin0)
		deltaTriggerEff  = 0.5*(1.0-jetTriggerWeight)
                jetTriggerWeightUp  =   jetTriggerWeight + deltaTriggerEff
                jetTriggerWeightDown  = jetTriggerWeight - deltaTriggerEff
                jetTriggerWeightUp  = min(1.0,jetTriggerWeightUp)
                jetTriggerWeightDown  = max(0.0,jetTriggerWeightDown)
                if TROP  == "nominal" :
                	   Weight = jetTriggerWeight
                if TROP  == "up" :
                	   Weight = jetTriggerWeightUp
                if TROP == "down" :
                	   Weight = jetTriggerWeightDown
	return Weight
#This looks up the ttbar pt reweighting scale factor 
def PTW_Lookup( GP ):
		genTpt = -100.
		genTBpt = -100	
		for ig in GP :
			isT = ig.pdgId() == 6 and ig.status() == 3
			isTB = ig.pdgId() == -6 and ig.status() == 3
			if isT:
				genTpt = ig.pt()
			if isTB:
				genTBpt = ig.pt()	
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
		if List[ihist].GetName() == 'QCDbkgBifpoly':
			nominalhist = List[ihist]
	for ibin in range(0,nominalhist.GetXaxis().GetNbins()+1):

		mse=0.0
		sigma=0.0
		sumsqdiff = 0.0
		for ihist in range(0,len(List)):
			if List[ihist].GetName() != 'QCDbkgBifpoly':
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
		if FScont < 0.99:
			pull.SetBinContent(ibin, 0.0 )	
		else:
			if sigma != 0 :
				pullcont = (pull.GetBinContent(ibin))/sigma
				pull.SetBinContent(ibin, pullcont)
			else :
				pull.SetBinContent(ibin, 0.0 )
	return pull
#Some lazy string formatting functions 
def strf( x ):
	return '%.2f' % x

def strf1( x ):
	return '%.0f' % x

