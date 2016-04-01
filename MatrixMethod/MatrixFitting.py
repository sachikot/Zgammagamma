from ROOT import TH1, TH1F, TH2, TH2D, TFile, TCanvas, THStack, TLegend, TMath, TList
from ROOT import gROOT, gStyle
from ROOT import *
import array, glob, ROOT
import sys, os, re, uuid, math, copy, imp, random, collections, pickle, time
from array import array
from uncertainties import ufloat
from uncertainties import unumpy

ROOT.gROOT.SetBatch()
ROOT.gStyle.SetOptStat(0)

_sieie_cuts  = { 'EB' : (0.011,0.029), 'EE' : (0.033,0.087) }
_chIso_cuts  = { 'EB' : (1.5,15.0),    'EE' : (1.2,15.6) }

def get_default_cuts(var) :
    
    if var == 'Sieie' :
        return { 'EB' : { 'tight' : ( 0, _sieie_cuts['EB'][0]-0.0001 ), 'loose' : ( _sieie_cuts['EB'][0]+0.0001,_sieie_cuts['EB'][1]-0.0001 ) },
                 'EE' : { 'tight' : ( 0, _sieie_cuts['EE'][0]-0.0001 ), 'loose' : ( _sieie_cuts['EE'][0]+0.0001,_sieie_cuts['EE'][1]-0.0001 ) } 
               }

    if var == 'ChIso' :
        return { 'EB' : { 'tight' : ( 0.0, _chIso_cuts['EB'][0]-0.01 ), 'loose' : ( _chIso_cuts['EB'][0]+0.01,_chIso_cuts['EB'][1]-0.01 ) },
                 'EE' : { 'tight' : ( 0.0, _chIso_cuts['EE'][0]-0.01 ), 'loose' : ( _chIso_cuts['EE'][0]+0.01,_chIso_cuts['EE'][1]-0.01 ) } 
               }

def get_syst_uncertainty(sample) :
    if sample == 'bakg' :
        # use a flat 20% uncertainty for now
        return 0.20

    if sample == 'data' :
        return 0.0

def get2PhoHist(filename, histname) :
    file = TFile.Open(filename)
    phoHist = file.Get(histname)
    phoHist.SetDirectory(0)
    phoHist.SetFillColor(0)

    return phoHist

def getTemp1DHist(filename, histname) :
    file = TFile.Open(filename)
    temp1DHist = file.Get(histname)
    temp1DHist.SetDirectory(0)
    temp1DHist.SetFillColor(0)

    return temp1DHist

def getTemp2DHist(filename, histname, varMin, varMax, projAxis) :
    file = TFile.Open(filename)
    #    print "histoname=%s" % histname
    #    file.Print()
    Hist2D = file.Get(histname)
    #    print type (Hist2D)
    Hist2D.SetDirectory(0)

    if projAxis == 'Y' :
        firstBin = Hist2D.GetYaxis().FindBin(varMin)
        endBin = Hist2D.GetYaxis().FindBin(varMax)
    elif projAxis == 'X' :
        firstBin = Hist2D.GetXaxis().FindBin(varMin)
        endBin = Hist2D.GetXaxis().FindBin(varMax)
    else :
        print 'Error!!! I need to put X or Y on projAxis!!'
        return

    # project a 2D histogram into a 1D histogram
    if projAxis == 'Y' :
        temp2DHist = Hist2D.ProjectionX('projX_'+histname+'_'+str(varMin)+'_'+str(varMax), firstBin, endBin)
    else :
        temp2DHist = Hist2D.ProjectionY('projY_'+histname+'_'+str(varMin)+'_'+str(varMax), firstBin, endBin)

    temp2DHist.SetDirectory(0)
    temp2DHist.SetFillColor(0)
    
    return temp2DHist

def getWeightedHist(filename, histName, leftCut=None, rightCut=None):

    if leftCut==None and rightCut==None:
        # template is 1D histogram
        sumHist = getTemp1DHist(filename, histName)
        return sumHist
    else :
        # template is projected 2D histogram
        sumHistp = getTemp2DHist(filename, histName, leftCut, rightCut, projAxis='Y')
        return sumHistp

def run_diphoton_fit(templates, gg_hist, lead_eta, subl_eta, templates_corr, ndim, var, outputDir=None ) :

    accept_reg = ['EB', 'EE']
    if lead_eta not in accept_reg :
        print 'Lead region does not make sense'
        return
    if subl_eta not in accept_reg :
        print 'Subl region does not make sense'
        return

    # cut point
    cuts = get_default_cuts(var)

    if isinstance(gg_hist, dict) :
        bins_subl_tight = (gg_hist['leadPass'].GetYaxis().FindBin( cuts[subl_eta]['tight'][0] ), gg_hist['leadPass'].GetYaxis().FindBin( cuts[subl_eta]['tight'][1] ))
        bins_subl_loose = (gg_hist['leadPass'].GetYaxis().FindBin( cuts[subl_eta]['loose'][0] ), gg_hist['leadPass'].GetYaxis().FindBin( cuts[subl_eta]['loose'][1] ))

        print 'bins_subl_tight ', bins_subl_tight
        print 'bins_subl_loose ', bins_subl_loose

        # Integrate to the the data in the four regions
        Ndata_TT = gg_hist['leadPass'].Integral( bins_subl_tight[0], bins_subl_tight[1] )
        Ndata_TL = gg_hist['leadPass'].Integral( bins_subl_loose[0], bins_subl_loose[1] )
        Ndata_LT = gg_hist['leadFail'].Integral( bins_subl_tight[0], bins_subl_tight[1] )
        Ndata_LL = gg_hist['leadFail'].Integral( bins_subl_loose[0], bins_subl_loose[1] )

    else :

        # gg_hist is 2d histogram between lead and sublead photons
        bins_lead_tight = (gg_hist.GetXaxis().FindBin(cuts[lead_eta]['tight'][0]), gg_hist.GetXaxis().FindBin(cuts[lead_eta]['tight'][1]))
        bins_lead_loose = (gg_hist.GetXaxis().FindBin(cuts[lead_eta]['loose'][0]), gg_hist.GetXaxis().FindBin(cuts[lead_eta]['loose'][1]))
        bins_subl_tight = (gg_hist.GetYaxis().FindBin(cuts[subl_eta]['tight'][0]), gg_hist.GetYaxis().FindBin(cuts[subl_eta]['tight'][1]))
        bins_subl_loose = (gg_hist.GetYaxis().FindBin(cuts[subl_eta]['loose'][0]), gg_hist.GetYaxis().FindBin(cuts[subl_eta]['loose'][1]))
        
        print 'bins_lead_tight ', bins_lead_tight
        print 'bins_lead_loose ', bins_lead_loose
        print 'bins_subl_tight ', bins_subl_tight
        print 'bins_subl_loose ', bins_subl_loose
        
        Ndata_TT = gg_hist.Integral(bins_lead_tight[0], bins_lead_tight[1], bins_subl_tight[0], bins_subl_tight[1])
        Ndata_TL = gg_hist.Integral(bins_lead_tight[0], bins_lead_tight[1], bins_subl_loose[0], bins_subl_loose[1])
        Ndata_LT = gg_hist.Integral(bins_lead_loose[0], bins_lead_loose[1], bins_subl_tight[0], bins_subl_tight[1])
        Ndata_LL = gg_hist.Integral(bins_lead_loose[0], bins_lead_loose[1], bins_subl_loose[0], bins_subl_loose[1])
        
    # ufloat it!
    Ndata = {}
    Ndata['TT'] = ufloat( Ndata_TT, math.sqrt(Ndata_TT ), 'Ndata_TT' )
    Ndata['TL'] = ufloat( Ndata_TL, math.sqrt(Ndata_TL ), 'Ndata_TL' )
    Ndata['LT'] = ufloat( Ndata_LT, math.sqrt(Ndata_LT ), 'Ndata_LT' )
    Ndata['LL'] = ufloat( Ndata_LL, math.sqrt(Ndata_LL ), 'Ndata_LL' )

    print 'N data TT = ', Ndata['TT']
    print 'N data TL = ', Ndata['TL']
    print 'N data LT = ', Ndata['LT']
    print 'N data LL = ', Ndata['LL']

    eff_cuts = {}
    eff_cuts['lead'] = {}
    eff_cuts['subl'] = {}
    eff_cuts['lead']['tight'] = cuts[lead_eta]['tight']
    eff_cuts['lead']['loose'] = cuts[lead_eta]['loose']
    eff_cuts['subl']['tight'] = cuts[subl_eta]['tight']
    eff_cuts['subl']['loose'] = cuts[subl_eta]['loose']

    #get 2D efficiencies from 1D inputs
    eff_results = generate_2d_efficiencies( templates, eff_cuts, lead_eta, subl_eta, var=var )
    (eff_1d_stat_tight, eff_1d_stat_loose, eff_1d_syst_bkg, eff_1d_syst_temp) =generate_1d_efficiencies( templates, eff_cuts, lead_eta, subl_eta, var=var )

    if templates_corr is not None :
        eff_ff_2d_stat, eff_ff_2d_syst = generate_2d_corr_efficiencies(templates_corr, eff_cuts, lead_eta, subl_eta, var=var)

        # for correlated FF templates
        # the uncertainties from 1-d templates
        # don't apply.  Set those to zero
        eff_results['stat_tight']['eff_FF_TT'] = ufloat( eff_ff_2d_stat['eff_FF_TT'].n, 0.0 )
        eff_results['stat_tight']['eff_FF_TL'] = ufloat( eff_ff_2d_stat['eff_FF_TL'].n, 0.0 )
        eff_results['stat_tight']['eff_FF_LT'] = ufloat( eff_ff_2d_stat['eff_FF_LT'].n, 0.0 )
        eff_results['stat_tight']['eff_FF_LL'] = ufloat( eff_ff_2d_stat['eff_FF_LL'].n, 0.0 )
        eff_results['stat_loose']['eff_FF_TT'] = ufloat( eff_ff_2d_stat['eff_FF_TT'].n, 0.0 )
        eff_results['stat_loose']['eff_FF_TL'] = ufloat( eff_ff_2d_stat['eff_FF_TL'].n, 0.0 )
        eff_results['stat_loose']['eff_FF_LT'] = ufloat( eff_ff_2d_stat['eff_FF_LT'].n, 0.0 )
        eff_results['stat_loose']['eff_FF_LL'] = ufloat( eff_ff_2d_stat['eff_FF_LL'].n, 0.0 )
        eff_results['syst_bkg'  ]['eff_FF_TT'] = ufloat( eff_ff_2d_stat['eff_FF_TT'].n, 0.0 )
        eff_results['syst_bkg'  ]['eff_FF_TL'] = ufloat( eff_ff_2d_stat['eff_FF_TL'].n, 0.0 )
        eff_results['syst_bkg'  ]['eff_FF_LT'] = ufloat( eff_ff_2d_stat['eff_FF_LT'].n, 0.0 )
        eff_results['syst_bkg'  ]['eff_FF_LL'] = ufloat( eff_ff_2d_stat['eff_FF_LL'].n, 0.0 )
        eff_results['syst_temp' ]['eff_FF_TT'] = ufloat( eff_ff_2d_stat['eff_FF_TT'].n, 0.0 )
        eff_results['syst_temp' ]['eff_FF_TL'] = ufloat( eff_ff_2d_stat['eff_FF_TL'].n, 0.0 )
        eff_results['syst_temp' ]['eff_FF_LT'] = ufloat( eff_ff_2d_stat['eff_FF_LT'].n, 0.0 )
        eff_results['syst_temp' ]['eff_FF_LL'] = ufloat( eff_ff_2d_stat['eff_FF_LL'].n, 0.0 )

        # make a new uncertainty
        # with fake-fake uncertainties
        eff_results['stat_ff'] = {}
        eff_results['stat_ff']['eff_FF_TT'] = eff_ff_2d_stat['eff_FF_TT']
        eff_results['stat_ff']['eff_FF_TL'] = eff_ff_2d_stat['eff_FF_TL']
        eff_results['stat_ff']['eff_FF_LT'] = eff_ff_2d_stat['eff_FF_LT']
        eff_results['stat_ff']['eff_FF_LL'] = eff_ff_2d_stat['eff_FF_LL']
        # put in all of the other values, use stat_tight as template
        eff_results['stat_ff']['eff_RR_TT'] = ufloat( eff_results['stat_tight']['eff_RR_TT'].n, 0.0 )
        eff_results['stat_ff']['eff_RR_TL'] = ufloat( eff_results['stat_tight']['eff_RR_TL'].n, 0.0 )
        eff_results['stat_ff']['eff_RR_LT'] = ufloat( eff_results['stat_tight']['eff_RR_LT'].n, 0.0 )
        eff_results['stat_ff']['eff_RR_LL'] = ufloat( eff_results['stat_tight']['eff_RR_LL'].n, 0.0 )
        
        eff_results['stat_ff']['eff_RF_TT'] = ufloat( eff_results['stat_tight']['eff_RF_TT'].n, 0.0 )
        eff_results['stat_ff']['eff_RF_TL'] = ufloat( eff_results['stat_tight']['eff_RF_TL'].n, 0.0 )
        eff_results['stat_ff']['eff_RF_LT'] = ufloat( eff_results['stat_tight']['eff_RF_LT'].n, 0.0 )
        eff_results['stat_ff']['eff_RF_LL'] = ufloat( eff_results['stat_tight']['eff_RF_LL'].n, 0.0 )
        
        eff_results['stat_ff']['eff_FR_TT'] = ufloat( eff_results['stat_tight']['eff_FR_TT'].n, 0.0 )
        eff_results['stat_ff']['eff_FR_TL'] = ufloat( eff_results['stat_tight']['eff_FR_TL'].n, 0.0 )
        eff_results['stat_ff']['eff_FR_LT'] = ufloat( eff_results['stat_tight']['eff_FR_LT'].n, 0.0 )
        eff_results['stat_ff']['eff_FR_LL'] = ufloat( eff_results['stat_tight']['eff_FR_LL'].n, 0.0 )
    else :
        # put in all of the other values, use stat_tight as template
        eff_results['stat_ff'] = {}
        eff_results['stat_ff']['eff_FF_TT'] = ufloat( eff_results['stat_tight']['eff_FF_TT'].n, 0.0 )
        eff_results['stat_ff']['eff_FF_TL'] = ufloat( eff_results['stat_tight']['eff_FF_TL'].n, 0.0 )
        eff_results['stat_ff']['eff_FF_LT'] = ufloat( eff_results['stat_tight']['eff_FF_LT'].n, 0.0 )
        eff_results['stat_ff']['eff_FF_LL'] = ufloat( eff_results['stat_tight']['eff_FF_LL'].n, 0.0 )

        eff_results['stat_ff']['eff_RR_TT'] = ufloat( eff_results['stat_tight']['eff_RR_TT'].n, 0.0 )
        eff_results['stat_ff']['eff_RR_TL'] = ufloat( eff_results['stat_tight']['eff_RR_TL'].n, 0.0 )
        eff_results['stat_ff']['eff_RR_LT'] = ufloat( eff_results['stat_tight']['eff_RR_LT'].n, 0.0 )
        eff_results['stat_ff']['eff_RR_LL'] = ufloat( eff_results['stat_tight']['eff_RR_LL'].n, 0.0 )
        
        eff_results['stat_ff']['eff_RF_TT'] = ufloat( eff_results['stat_tight']['eff_RF_TT'].n, 0.0 )
        eff_results['stat_ff']['eff_RF_TL'] = ufloat( eff_results['stat_tight']['eff_RF_TL'].n, 0.0 )
        eff_results['stat_ff']['eff_RF_LT'] = ufloat( eff_results['stat_tight']['eff_RF_LT'].n, 0.0 )
        eff_results['stat_ff']['eff_RF_LL'] = ufloat( eff_results['stat_tight']['eff_RF_LL'].n, 0.0 )
        
        eff_results['stat_ff']['eff_FR_TT'] = ufloat( eff_results['stat_tight']['eff_FR_TT'].n, 0.0 )
        eff_results['stat_ff']['eff_FR_TL'] = ufloat( eff_results['stat_tight']['eff_FR_TL'].n, 0.0 )
        eff_results['stat_ff']['eff_FR_LT'] = ufloat( eff_results['stat_tight']['eff_FR_LT'].n, 0.0 )
        eff_results['stat_ff']['eff_FR_LL'] = ufloat( eff_results['stat_tight']['eff_FR_LL'].n, 0.0 )


    eff_2d_nouncert = {}
    for key, val in eff_results['stat_tight'].iteritems() :
        if type( val ) == type( ufloat(0,0) ) :
            eff_2d_nouncert[key] = ufloat( val.n, 0.0 )
        else :
            eff_2d_nouncert[key] = val

    if ndim == 3 :

        data = {'TL': Ndata['TL'], 'LT' : Ndata['LT'], 'LL' : Ndata['LL']}
        data_nostat = { }
        data_nostat['TL'] = ufloat( Ndata['TL'].n, 0.0 )
        data_nostat['LT'] = ufloat( Ndata['LT'].n, 0.0 )
        data_nostat['LL'] = ufloat( Ndata['LL'].n, 0.0 )

        results_stat_data       = run_fit( data       , eff_2d_nouncert)
        results_stat_temp_tight = run_fit( data_nostat, eff_results['stat_tight'] )
        results_stat_temp_loose = run_fit( data_nostat, eff_results['stat_loose'] )
        results_stat_ff         = run_fit( data_nostat, eff_results['stat_ff'] )
        results_syst_bkg        = run_fit( data_nostat, eff_results['syst_bkg'] )
        results_syst_temp       = run_fit( data_nostat, eff_results['syst_temp'] )

    if ndim == 4 :
        print 'RUNNING FIT WITH 4 DIM'
        # consider only SB uncertainties
        data_SB = { }
        data_SB['TT'] = ufloat( Ndata['TT'].n, 0.0 )
        data_SB['TL'] = Ndata['TL']
        data_SB['LT'] = Ndata['LT']
        data_SB['LL'] = Ndata['LL']

        # consider only SR uncertainties
        data_SR = { }
        data_SR['TT'] =  Ndata['TT']
        data_SR['TL'] = ufloat( Ndata['TL'].n, 0.0 )
        data_SR['LT'] = ufloat( Ndata['LT'].n, 0.0 )
        data_SR['LL'] = ufloat( Ndata['LL'].n, 0.0 )

        # no data uncertainties
        data_nostat = { }
        data_nostat['TT'] = ufloat( Ndata['TT'].n, 0.0 )
        data_nostat['TL'] = ufloat( Ndata['TL'].n, 0.0 )
        data_nostat['LT'] = ufloat( Ndata['LT'].n, 0.0 )
        data_nostat['LL'] = ufloat( Ndata['LL'].n, 0.0 )

        results_stat_dataSB     = run_fit( data_SB, eff_2d_nouncert)
        results_stat_dataSR     = run_fit( data_SR, eff_2d_nouncert)
        results_stat_temp_tight = run_fit( data_nostat, eff_results['stat_tight'] )
        results_stat_temp_loose = run_fit( data_nostat, eff_results['stat_loose'] )
        results_stat_ff         = run_fit( data_nostat, eff_results['stat_ff'] )
        results_syst_bkg        = run_fit( data_nostat, eff_results['syst_bkg'] )
        results_syst_temp       = run_fit( data_nostat, eff_results['syst_temp'] )

    if ndim == 3 :
        idxrf = 0
        idxfr = 1
        idxff = 2

    if ndim == 4 :
        idxrf = 1
        idxfr = 2
        idxff = 3

    #get fitted predictions
    if ndim == 3 :
         p_RR_TT = ufloat(1.0, 0.0 )
         p_RR_TL = ufloat(0.0, 0.0 )
         p_RR_LT = ufloat(0.0, 0.0 )
         p_RR_LL = ufloat(0.0, 0.0 )
         p_RF_TT = results_stat_dataSB.item(0)*eff_results['stat_tight']['eff_RF_TT']
         p_RF_TL = results_stat_dataSB.item(0)*eff_results['stat_tight']['eff_RF_TL']
         p_RF_LT = results_stat_dataSB.item(0)*eff_results['stat_tight']['eff_RF_LT']
         p_RF_LL = results_stat_dataSB.item(0)*eff_results['stat_tight']['eff_RF_LL']
         p_FR_TT = results_stat_dataSB.item(1)*eff_results['stat_tight']['eff_FR_TT']
         p_FR_TL = results_stat_dataSB.item(1)*eff_results['stat_tight']['eff_FR_TL']
         p_FR_LT = results_stat_dataSB.item(1)*eff_results['stat_tight']['eff_FR_LT']
         p_FR_LL = results_stat_dataSB.item(1)*eff_results['stat_tight']['eff_FR_LL']
         p_FF_TT = results_stat_dataSB.item(2)*eff_results['stat_tight']['eff_FF_TT']
         p_FF_TL = results_stat_dataSB.item(2)*eff_results['stat_tight']['eff_FF_TL']
         p_FF_LT = results_stat_dataSB.item(2)*eff_results['stat_tight']['eff_FF_LT']
         p_FF_LL = results_stat_dataSB.item(2)*eff_results['stat_tight']['eff_FF_LL']

    if ndim == 4 : #everybody moves down 1 if 4 dim
         p_RR_TT = results_stat_dataSB.item(0)*eff_results['stat_tight']['eff_RR_TT']
         p_RR_TL = results_stat_dataSB.item(0)*eff_results['stat_tight']['eff_RR_TL']
         p_RR_LT = results_stat_dataSB.item(0)*eff_results['stat_tight']['eff_RR_LT']
         p_RR_LL = results_stat_dataSB.item(0)*eff_results['stat_tight']['eff_RR_LL']
         p_RF_TT = results_stat_dataSB.item(1)*eff_results['stat_tight']['eff_RF_TT']
         p_RF_TL = results_stat_dataSB.item(1)*eff_results['stat_tight']['eff_RF_TL']
         p_RF_LT = results_stat_dataSB.item(1)*eff_results['stat_tight']['eff_RF_LT']
         p_RF_LL = results_stat_dataSB.item(1)*eff_results['stat_tight']['eff_RF_LL']
         p_FR_TT = results_stat_dataSB.item(2)*eff_results['stat_tight']['eff_FR_TT']
         p_FR_TL = results_stat_dataSB.item(2)*eff_results['stat_tight']['eff_FR_TL']
         p_FR_LT = results_stat_dataSB.item(2)*eff_results['stat_tight']['eff_FR_LT']
         p_FR_LL = results_stat_dataSB.item(2)*eff_results['stat_tight']['eff_FR_LL']
         p_FF_TT = results_stat_dataSB.item(3)*eff_results['stat_tight']['eff_FF_TT']
         p_FF_TL = results_stat_dataSB.item(3)*eff_results['stat_tight']['eff_FF_TL']
         p_FF_LT = results_stat_dataSB.item(3)*eff_results['stat_tight']['eff_FF_LT']
         p_FF_LL = results_stat_dataSB.item(3)*eff_results['stat_tight']['eff_FF_LL']

    print 'Npred FF_LL'
    print p_FF_LL
    print 'Npred FF_LT'
    print p_FF_LT
    print 'Npred FF_TL'
    print p_FF_TL
    print 'Npred FF_TT'
    print p_FF_TT

    print 'Npred RF_LL'
    print p_RF_LL
    print 'Npred RF_LT'
    print p_RF_LT
    print 'Npred RF_TL'
    print p_RF_TL
    print 'Npred RF_TT'
    print p_RF_TT

    print 'Npred FR_LL'
    print p_FR_LL
    print 'Npred FR_LT'
    print p_FR_LT
    print 'Npred FR_TL'
    print p_FR_TL
    print 'Npred FR_TT'
    print p_FR_TT

    print 'Npred RR_LL'
    print p_RR_LL
    print 'Npred RR_LT'
    print p_RR_LT
    print 'Npred RR_TL'
    print p_RR_TL
    print 'Npred RR_TT'
    print p_RR_TT

    print 'Npred Bkg Total' 
    print (p_FF_TT+p_RF_TT+p_FR_TT)

    create_2dEff_hist(eff_2d_nouncert, reg[0], reg[1], outputDir, ndim)

    text_results_stat_dataSB       = collect_results( results_stat_dataSB      , Ndata, eff_2d_nouncert , templates, eff_cuts, ndim)
    text_results_stat_dataSR       = collect_results( results_stat_dataSR      , Ndata, eff_2d_nouncert , templates, eff_cuts, ndim)
    text_results_stat_temp_tight = collect_results( results_stat_temp_tight, Ndata, eff_results['stat_tight'] , templates, eff_cuts, ndim)
    text_results_stat_temp_loose = collect_results( results_stat_temp_loose, Ndata, eff_results['stat_loose'] , templates, eff_cuts, ndim)
    text_results_stat_ff         = collect_results( results_stat_ff        , Ndata, eff_results['stat_ff'] , templates, eff_cuts, ndim)
    text_results_syst_bkg        = collect_results( results_syst_bkg       , Ndata, eff_results['syst_bkg']   , templates, eff_cuts, ndim)
    text_results_syst_temp       = collect_results( results_syst_temp      , Ndata, eff_results['syst_temp']  , templates, eff_cuts, ndim)

    return text_results_stat_dataSR,text_results_stat_dataSB, text_results_stat_temp_tight, text_results_stat_temp_loose, text_results_stat_ff, text_results_syst_bkg, text_results_syst_temp

def collect_results(results, data, efficiencies, templates, cuts, ndim) :

    text_results = collections.OrderedDict()

    for key, val in efficiencies.iteritems() :
        text_results[key] = val

    if ndim == 4 :
        
        text_results['Ndata_TT'] = data['TT']
        text_results['Ndata_TL'] = data['TL']
        text_results['Ndata_LT'] = data['LT']
        text_results['Ndata_LL'] = data['LL']

        text_results['alpha_RR'] = results.item(0)
        text_results['alpha_RF'] = results.item(1)
        text_results['alpha_FR'] = results.item(2)
        text_results['alpha_FF'] = results.item(3)

        text_results['Npred_RR_TT'] = text_results['alpha_RR']*text_results['eff_RR_TT']
        text_results['Npred_RR_TL'] = text_results['alpha_RR']*text_results['eff_RR_TL']
        text_results['Npred_RR_LT'] = text_results['alpha_RR']*text_results['eff_RR_LT']
        text_results['Npred_RR_LL'] = text_results['alpha_RR']*text_results['eff_RR_LL']

    else :
        text_results['Ndata_TT'] = ufloat(0, 0)
        text_results['Ndata_TL'] = data['TL']
        text_results['Ndata_LT'] = data['LT']
        text_results['Ndata_LL'] = data['LL']


        text_results['alpha_RF'] = results.item(0)
        text_results['alpha_FR'] = results.item(1)
        text_results['alpha_FF'] = results.item(2)


    text_results['Npred_RF_TT'] = text_results['alpha_RF']*text_results['eff_RF_TT']
    text_results['Npred_RF_TL'] = text_results['alpha_RF']*text_results['eff_RF_TL']
    text_results['Npred_RF_LT'] = text_results['alpha_RF']*text_results['eff_RF_LT']
    text_results['Npred_RF_LL'] = text_results['alpha_RF']*text_results['eff_RF_LL']

    text_results['Npred_FR_TT'] = text_results['alpha_FR']*text_results['eff_FR_TT']
    text_results['Npred_FR_TL'] = text_results['alpha_FR']*text_results['eff_FR_TL']
    text_results['Npred_FR_LT'] = text_results['alpha_FR']*text_results['eff_FR_LT']
    text_results['Npred_FR_LL'] = text_results['alpha_FR']*text_results['eff_FR_LL']

    text_results['Npred_FF_TT'] = text_results['alpha_FF']*text_results['eff_FF_TT']
    text_results['Npred_FF_TL'] = text_results['alpha_FF']*text_results['eff_FF_TL']
    text_results['Npred_FF_LT'] = text_results['alpha_FF']*text_results['eff_FF_LT']
    text_results['Npred_FF_LL'] = text_results['alpha_FF']*text_results['eff_FF_LL']

    # add the template integrals to results

    bins_lead_tight = ( templates['lead']['real'].GetXaxis().FindBin(cuts['lead']['tight'][0]), templates['lead']['real'].GetXaxis().FindBin(cuts['lead']['tight'][1]) )
    bins_lead_loose = ( templates['lead']['real'].GetXaxis().FindBin(cuts['lead']['loose'][0]), templates['lead']['real'].GetXaxis().FindBin(cuts['lead']['loose'][1]) )
    bins_subl_tight = ( templates['subl']['real'].GetXaxis().FindBin(cuts['subl']['tight'][0]), templates['subl']['real'].GetXaxis().FindBin(cuts['subl']['tight'][1]) )
    bins_subl_loose = ( templates['subl']['real'].GetXaxis().FindBin(cuts['subl']['loose'][0]), templates['subl']['real'].GetXaxis().FindBin(cuts['subl']['loose'][1]) )

    int_lead_real_loose = get_integral_and_error(templates['lead']['real'], bins_lead_loose )
    int_lead_real_tight = get_integral_and_error(templates['lead']['real'], bins_lead_tight )
    int_lead_fake_loose = get_integral_and_error(templates['lead']['fake'], bins_lead_loose )
    int_lead_fake_tight = get_integral_and_error(templates['lead']['fake'], bins_lead_tight )

    int_subl_real_loose = get_integral_and_error(templates['subl']['real'], bins_subl_loose )
    int_subl_real_tight = get_integral_and_error(templates['subl']['real'], bins_subl_tight )
    int_subl_fake_loose = get_integral_and_error(templates['subl']['fake'], bins_subl_loose )
    int_subl_fake_tight = get_integral_and_error(templates['subl']['fake'], bins_subl_tight )


    text_results['template_int_lead_real_loose'] = int_lead_real_loose
    text_results['template_int_lead_real_tight'] = int_lead_real_tight
    text_results['template_int_lead_fake_loose'] = int_lead_fake_loose
    text_results['template_int_lead_fake_tight'] = int_lead_fake_tight
    text_results['template_int_subl_real_loose'] = int_subl_real_loose
    text_results['template_int_subl_real_tight'] = int_subl_real_tight
    text_results['template_int_subl_fake_loose'] = int_subl_fake_loose
    text_results['template_int_subl_fake_tight'] = int_subl_fake_tight

    return text_results

def generate_1d_efficiencies( templates, eff_cuts, lead_eta, subl_eta, var=None) :

    (int_stat, int_syst) = get_template_integrals( templates, eff_cuts, lead_eta, subl_eta, var=var)
    (eff_1d_stat_tight, eff_1d_stat_loose, eff_1d_syst_bkg, eff_1d_syst_temp) = get_1d_loose_efficiencies( int_stat, int_syst, lead_eta, subl_eta, var=var)

    return eff_1d_stat_tight, eff_1d_stat_loose, eff_1d_syst_bkg, eff_1d_syst_temp


def generate_2d_efficiencies( templates, eff_cuts, lead_eta, subl_eta, var=None) :

    (int_stat, int_syst) = get_template_integrals( templates, eff_cuts, lead_eta, subl_eta, var=var)

    # get efficiencies with three sources of uncertainty
    # 1. statistical uncertainty on templates
    # 2. systematics uncertainty onbackground subtraction
    # 3. systematic uncertainty on template shapes
    # The integrals from get_template_integrals
    # give 1 and 2

    eff_results = {}

    (eff_1d_stat_tight, eff_1d_stat_loose, eff_1d_syst_bkg, eff_1d_syst_temp) = get_1d_loose_efficiencies( int_stat, int_syst, lead_eta, subl_eta, var=var )

    eff_stat_tight = collections.OrderedDict()
    eff_stat_loose = collections.OrderedDict()
    eff_syst_bkg   = collections.OrderedDict()
    eff_syst_temp  = collections.OrderedDict()

    stat_tight_eff_R_L_lead = eff_1d_stat_tight['eff_R_L_lead']
    stat_tight_eff_F_L_lead = eff_1d_stat_tight['eff_F_L_lead']
    stat_tight_eff_R_L_subl = eff_1d_stat_tight['eff_R_L_subl']
    stat_tight_eff_F_L_subl = eff_1d_stat_tight['eff_F_L_subl']

    stat_loose_eff_R_L_lead = eff_1d_stat_loose['eff_R_L_lead']
    stat_loose_eff_F_L_lead = eff_1d_stat_loose['eff_F_L_lead']
    stat_loose_eff_R_L_subl = eff_1d_stat_loose['eff_R_L_subl']
    stat_loose_eff_F_L_subl = eff_1d_stat_loose['eff_F_L_subl']

    syst_bkg_eff_R_L_lead = eff_1d_syst_bkg['eff_R_L_lead']
    syst_bkg_eff_F_L_lead = eff_1d_syst_bkg['eff_F_L_lead']
    syst_bkg_eff_R_L_subl = eff_1d_syst_bkg['eff_R_L_subl']
    syst_bkg_eff_F_L_subl = eff_1d_syst_bkg['eff_F_L_subl']

    syst_temp_eff_R_L_lead = eff_1d_syst_temp['eff_R_L_lead']
    syst_temp_eff_F_L_lead = eff_1d_syst_temp['eff_F_L_lead']
    syst_temp_eff_R_L_subl = eff_1d_syst_temp['eff_R_L_subl']
    syst_temp_eff_F_L_subl = eff_1d_syst_temp['eff_F_L_subl']

    # tight efficiencies are just 1-loose efficiencies
    stat_tight_eff_R_T_lead = ufloat(1.0, 0.0) - stat_tight_eff_R_L_lead
    stat_tight_eff_F_T_lead = ufloat(1.0, 0.0) - stat_tight_eff_F_L_lead
    stat_tight_eff_R_T_subl = ufloat(1.0, 0.0) - stat_tight_eff_R_L_subl
    stat_tight_eff_F_T_subl = ufloat(1.0, 0.0) - stat_tight_eff_F_L_subl

    stat_loose_eff_R_T_lead = ufloat(1.0, 0.0) - stat_loose_eff_R_L_lead
    stat_loose_eff_F_T_lead = ufloat(1.0, 0.0) - stat_loose_eff_F_L_lead
    stat_loose_eff_R_T_subl = ufloat(1.0, 0.0) - stat_loose_eff_R_L_subl
    stat_loose_eff_F_T_subl = ufloat(1.0, 0.0) - stat_loose_eff_F_L_subl

    syst_bkg_eff_R_T_lead = ufloat(1.0, 0.0)  - syst_bkg_eff_R_L_lead
    syst_bkg_eff_F_T_lead = ufloat(1.0, 0.0)  - syst_bkg_eff_F_L_lead
    syst_bkg_eff_R_T_subl = ufloat(1.0, 0.0)  - syst_bkg_eff_R_L_subl
    syst_bkg_eff_F_T_subl = ufloat(1.0, 0.0)  - syst_bkg_eff_F_L_subl

    syst_temp_eff_R_T_lead = ufloat(1.0, 0.0)  - syst_temp_eff_R_L_lead
    syst_temp_eff_F_T_lead = ufloat(1.0, 0.0)  - syst_temp_eff_F_L_lead
    syst_temp_eff_R_T_subl = ufloat(1.0, 0.0)  - syst_temp_eff_R_L_subl
    syst_temp_eff_F_T_subl = ufloat(1.0, 0.0)  - syst_temp_eff_F_L_subl

    print 'syst_bkg_eff_R_L_lead ',syst_bkg_eff_R_L_lead  
    print 'syst_bkg_eff_F_L_lead ',syst_bkg_eff_F_L_lead  
    print 'syst_bkg_eff_R_L_subl ',syst_bkg_eff_R_L_subl  
    print 'syst_bkg_eff_F_L_subl ',syst_bkg_eff_F_L_subl  

    print 'syst_bkg_eff_R_T_lead ',syst_bkg_eff_R_T_lead
    print 'syst_bkg_eff_F_T_lead ',syst_bkg_eff_F_T_lead
    print 'syst_bkg_eff_R_T_subl ',syst_bkg_eff_R_T_subl
    print 'syst_bkg_eff_F_T_subl ',syst_bkg_eff_F_T_subl

    print 'syst_temp_eff_R_L_lead ',syst_temp_eff_R_L_lead  
    print 'syst_temp_eff_F_L_lead ',syst_temp_eff_F_L_lead  
    print 'syst_temp_eff_R_L_subl ',syst_temp_eff_R_L_subl  
    print 'syst_temp_eff_F_L_subl ',syst_temp_eff_F_L_subl  

    print 'syst_temp_eff_R_T_lead ',syst_temp_eff_R_T_lead
    print 'syst_temp_eff_F_T_lead ',syst_temp_eff_F_T_lead
    print 'syst_temp_eff_R_T_subl ',syst_temp_eff_R_T_subl
    print 'syst_temp_eff_F_T_subl ',syst_temp_eff_F_T_subl

    # store the 1-d efficiencies in the
    eff_stat_tight['eff_1d'] = eff_1d_stat_tight
    eff_stat_loose['eff_1d'] = eff_1d_stat_loose
    eff_syst_bkg  ['eff_1d'] = eff_1d_syst_bkg
    eff_syst_temp ['eff_1d'] = eff_1d_syst_temp

    eff_stat_tight['eff_RR_TT'] = stat_tight_eff_R_T_lead*stat_tight_eff_R_T_subl 
    eff_stat_tight['eff_RR_TL'] = stat_tight_eff_R_T_lead*stat_tight_eff_R_L_subl 
    eff_stat_tight['eff_RR_LT'] = stat_tight_eff_R_L_lead*stat_tight_eff_R_T_subl 
    eff_stat_tight['eff_RR_LL'] = stat_tight_eff_R_L_lead*stat_tight_eff_R_L_subl 

    eff_stat_tight['eff_RF_TT'] = stat_tight_eff_R_T_lead*stat_tight_eff_F_T_subl 
    eff_stat_tight['eff_RF_TL'] = stat_tight_eff_R_T_lead*stat_tight_eff_F_L_subl 
    eff_stat_tight['eff_RF_LT'] = stat_tight_eff_R_L_lead*stat_tight_eff_F_T_subl 
    eff_stat_tight['eff_RF_LL'] = stat_tight_eff_R_L_lead*stat_tight_eff_F_L_subl 

    eff_stat_tight['eff_FR_TT'] = stat_tight_eff_F_T_lead*stat_tight_eff_R_T_subl 
    eff_stat_tight['eff_FR_TL'] = stat_tight_eff_F_T_lead*stat_tight_eff_R_L_subl 
    eff_stat_tight['eff_FR_LT'] = stat_tight_eff_F_L_lead*stat_tight_eff_R_T_subl 
    eff_stat_tight['eff_FR_LL'] = stat_tight_eff_F_L_lead*stat_tight_eff_R_L_subl 

    eff_stat_tight['eff_FF_TT'] = stat_tight_eff_F_T_lead*stat_tight_eff_F_T_subl 
    eff_stat_tight['eff_FF_TL'] = stat_tight_eff_F_T_lead*stat_tight_eff_F_L_subl 
    eff_stat_tight['eff_FF_LT'] = stat_tight_eff_F_L_lead*stat_tight_eff_F_T_subl 
    eff_stat_tight['eff_FF_LL'] = stat_tight_eff_F_L_lead*stat_tight_eff_F_L_subl 

    eff_stat_loose['eff_RR_TT'] = stat_loose_eff_R_T_lead*stat_loose_eff_R_T_subl 
    eff_stat_loose['eff_RR_TL'] = stat_loose_eff_R_T_lead*stat_loose_eff_R_L_subl 
    eff_stat_loose['eff_RR_LT'] = stat_loose_eff_R_L_lead*stat_loose_eff_R_T_subl 
    eff_stat_loose['eff_RR_LL'] = stat_loose_eff_R_L_lead*stat_loose_eff_R_L_subl 

    eff_stat_loose['eff_RF_TT'] = stat_loose_eff_R_T_lead*stat_loose_eff_F_T_subl 
    eff_stat_loose['eff_RF_TL'] = stat_loose_eff_R_T_lead*stat_loose_eff_F_L_subl 
    eff_stat_loose['eff_RF_LT'] = stat_loose_eff_R_L_lead*stat_loose_eff_F_T_subl 
    eff_stat_loose['eff_RF_LL'] = stat_loose_eff_R_L_lead*stat_loose_eff_F_L_subl 

    eff_stat_loose['eff_FR_TT'] = stat_loose_eff_F_T_lead*stat_loose_eff_R_T_subl 
    eff_stat_loose['eff_FR_TL'] = stat_loose_eff_F_T_lead*stat_loose_eff_R_L_subl 
    eff_stat_loose['eff_FR_LT'] = stat_loose_eff_F_L_lead*stat_loose_eff_R_T_subl 
    eff_stat_loose['eff_FR_LL'] = stat_loose_eff_F_L_lead*stat_loose_eff_R_L_subl 

    eff_stat_loose['eff_FF_TT'] = stat_loose_eff_F_T_lead*stat_loose_eff_F_T_subl 
    eff_stat_loose['eff_FF_TL'] = stat_loose_eff_F_T_lead*stat_loose_eff_F_L_subl 
    eff_stat_loose['eff_FF_LT'] = stat_loose_eff_F_L_lead*stat_loose_eff_F_T_subl 
    eff_stat_loose['eff_FF_LL'] = stat_loose_eff_F_L_lead*stat_loose_eff_F_L_subl 

    eff_syst_bkg['eff_RR_TT'] = syst_bkg_eff_R_T_lead*syst_bkg_eff_R_T_subl 
    eff_syst_bkg['eff_RR_TL'] = syst_bkg_eff_R_T_lead*syst_bkg_eff_R_L_subl 
    eff_syst_bkg['eff_RR_LT'] = syst_bkg_eff_R_L_lead*syst_bkg_eff_R_T_subl 
    eff_syst_bkg['eff_RR_LL'] = syst_bkg_eff_R_L_lead*syst_bkg_eff_R_L_subl 

    eff_syst_bkg['eff_RF_TT'] = syst_bkg_eff_R_T_lead*syst_bkg_eff_F_T_subl 
    eff_syst_bkg['eff_RF_TL'] = syst_bkg_eff_R_T_lead*syst_bkg_eff_F_L_subl 
    eff_syst_bkg['eff_RF_LT'] = syst_bkg_eff_R_L_lead*syst_bkg_eff_F_T_subl 
    eff_syst_bkg['eff_RF_LL'] = syst_bkg_eff_R_L_lead*syst_bkg_eff_F_L_subl 

    eff_syst_bkg['eff_FR_TT'] = syst_bkg_eff_F_T_lead*syst_bkg_eff_R_T_subl 
    eff_syst_bkg['eff_FR_TL'] = syst_bkg_eff_F_T_lead*syst_bkg_eff_R_L_subl 
    eff_syst_bkg['eff_FR_LT'] = syst_bkg_eff_F_L_lead*syst_bkg_eff_R_T_subl 
    eff_syst_bkg['eff_FR_LL'] = syst_bkg_eff_F_L_lead*syst_bkg_eff_R_L_subl 

    eff_syst_bkg['eff_FF_TT'] = syst_bkg_eff_F_T_lead*syst_bkg_eff_F_T_subl 
    eff_syst_bkg['eff_FF_TL'] = syst_bkg_eff_F_T_lead*syst_bkg_eff_F_L_subl 
    eff_syst_bkg['eff_FF_LT'] = syst_bkg_eff_F_L_lead*syst_bkg_eff_F_T_subl 
    eff_syst_bkg['eff_FF_LL'] = syst_bkg_eff_F_L_lead*syst_bkg_eff_F_L_subl 

    eff_syst_temp['eff_RR_TT'] = syst_temp_eff_R_T_lead*syst_temp_eff_R_T_subl 
    eff_syst_temp['eff_RR_TL'] = syst_temp_eff_R_T_lead*syst_temp_eff_R_L_subl 
    eff_syst_temp['eff_RR_LT'] = syst_temp_eff_R_L_lead*syst_temp_eff_R_T_subl 
    eff_syst_temp['eff_RR_LL'] = syst_temp_eff_R_L_lead*syst_temp_eff_R_L_subl 

    eff_syst_temp['eff_RF_TT'] = syst_temp_eff_R_T_lead*syst_temp_eff_F_T_subl 
    eff_syst_temp['eff_RF_TL'] = syst_temp_eff_R_T_lead*syst_temp_eff_F_L_subl 
    eff_syst_temp['eff_RF_LT'] = syst_temp_eff_R_L_lead*syst_temp_eff_F_T_subl 
    eff_syst_temp['eff_RF_LL'] = syst_temp_eff_R_L_lead*syst_temp_eff_F_L_subl 

    eff_syst_temp['eff_FR_TT'] = syst_temp_eff_F_T_lead*syst_temp_eff_R_T_subl 
    eff_syst_temp['eff_FR_TL'] = syst_temp_eff_F_T_lead*syst_temp_eff_R_L_subl 
    eff_syst_temp['eff_FR_LT'] = syst_temp_eff_F_L_lead*syst_temp_eff_R_T_subl 
    eff_syst_temp['eff_FR_LL'] = syst_temp_eff_F_L_lead*syst_temp_eff_R_L_subl 

    eff_syst_temp['eff_FF_TT'] = syst_temp_eff_F_T_lead*syst_temp_eff_F_T_subl 
    eff_syst_temp['eff_FF_TL'] = syst_temp_eff_F_T_lead*syst_temp_eff_F_L_subl 
    eff_syst_temp['eff_FF_LT'] = syst_temp_eff_F_L_lead*syst_temp_eff_F_T_subl 
    eff_syst_temp['eff_FF_LL'] = syst_temp_eff_F_L_lead*syst_temp_eff_F_L_subl 

    eff_results['stat_tight'] = eff_stat_tight
    eff_results['stat_loose'] = eff_stat_loose
    eff_results['syst_bkg']   = eff_syst_bkg
    eff_results['syst_temp']  = eff_syst_temp

    return eff_results

def generate_2d_corr_efficiencies( templates, cuts, lead_eta, subl_eta, var=None ) :
    # integral each region of the 2-d templates

    # lead is on y axis, subl on x
    bin_subl_tight = ( templates['leadPass'].GetXaxis().FindBin(cuts['subl']['tight'][0]), templates['leadPass'].GetXaxis().FindBin(cuts['subl']['tight'][1]) )
    bin_subl_loose = ( templates['leadPass'].GetXaxis().FindBin(cuts['subl']['loose'][0]), templates['leadPass'].GetXaxis().FindBin(cuts['subl']['loose'][1]) )

    err_lead_tight_subl_tight = ROOT.Double()
    int_lead_tight_subl_tight = templates['leadPass'].IntegralAndError(bin_subl_tight[0], bin_subl_tight[1], err_lead_tight_subl_tight)

    err_lead_tight_subl_loose = ROOT.Double()
    int_lead_tight_subl_loose = templates['leadPass'].IntegralAndError(bin_subl_loose[0], bin_subl_loose[1], err_lead_tight_subl_loose)

    err_lead_loose_subl_tight = ROOT.Double()
    int_lead_loose_subl_tight = templates['leadFail'].IntegralAndError(bin_subl_tight[0], bin_subl_tight[1], err_lead_loose_subl_tight)

    err_lead_loose_subl_loose = ROOT.Double()
    int_lead_loose_subl_loose = templates['leadFail'].IntegralAndError(bin_subl_loose[0], bin_subl_loose[1], err_lead_loose_subl_loose)

    int_tt = ufloat( int_lead_tight_subl_tight, err_lead_tight_subl_tight )
    int_tl = ufloat( int_lead_tight_subl_loose, err_lead_tight_subl_loose )
    int_lt = ufloat( int_lead_loose_subl_tight, err_lead_loose_subl_tight )
    int_ll = ufloat( int_lead_loose_subl_loose, err_lead_loose_subl_loose)

    print '2d FF template N TT = ', int_tt
    print '2d FF template N TL = ', int_tl
    print '2d FF template N LT = ', int_lt
    print '2d FF template N LL = ', int_ll

    denominator = int_tt + int_tl + int_lt + int_ll

    if denominator != 0 :
        frac_tt = int_tt/( denominator )
        frac_tl = int_tl/( denominator )
        frac_lt = int_lt/( denominator )
        frac_ll = int_ll/( denominator )
    else :
        frac_tt = ufloat(0,0)
        frac_tl = ufloat(0,0)
        frac_lt = ufloat(0,0)
        frac_ll = ufloat(0,0)

    eff_stat = {}
    eff_syst = {}

    eff_stat['eff_FF_TT'] = frac_tt
    eff_stat['eff_FF_TL'] = frac_tl
    eff_stat['eff_FF_LT'] = frac_lt
    eff_stat['eff_FF_LL'] = frac_ll

    eff_syst['eff_FF_TT'] = frac_tt
    eff_syst['eff_FF_TL'] = frac_tl
    eff_syst['eff_FF_LT'] = frac_lt
    eff_syst['eff_FF_LL'] = frac_ll

    return eff_stat, eff_syst


def get_1d_loose_efficiencies( int_stat, int_syst, lead_eta, subl_eta, var=None) :

    eff_stat = {}
    eff_stat_tight = {}
    eff_stat_loose = {}
    eff_syst_int = {}
    eff_syst_temp = {}

    eff_stat['eff_R_L_lead'] = int_stat['lead']['real']['loose'] / (int_stat['lead']['real']['tight']+int_stat['lead']['real']['loose'])
    eff_stat['eff_F_L_lead'] = int_stat['lead']['fake']['loose'] / (int_stat['lead']['fake']['tight']+int_stat['lead']['fake']['loose'])
    eff_stat['eff_R_L_subl'] = int_stat['subl']['real']['loose'] / (int_stat['subl']['real']['tight']+int_stat['subl']['real']['loose'])
    eff_stat['eff_F_L_subl'] = int_stat['subl']['fake']['loose'] / (int_stat['subl']['fake']['tight']+int_stat['subl']['fake']['loose'])

    int_stat_notunc = {'lead' : {}, 'subl' : {}}
    int_stat_notunc['lead'] = { 'real' : {}, 'fake' : {} }
    int_stat_notunc['subl'] = { 'real' : {}, 'fake' : {} }

    int_stat_notunc['lead']['real']['tight'] = ufloat( int_stat['lead']['real']['tight'].n, 0.0 )
    int_stat_notunc['lead']['fake']['tight'] = ufloat( int_stat['lead']['fake']['tight'].n, 0.0 )
    int_stat_notunc['subl']['real']['tight'] = ufloat( int_stat['subl']['real']['tight'].n, 0.0 )
    int_stat_notunc['subl']['fake']['tight'] = ufloat( int_stat['subl']['fake']['tight'].n, 0.0 )

    int_stat_nolunc = {'lead' : {}, 'subl' : {}}
    int_stat_nolunc['lead'] = { 'real' : {}, 'fake' : {} }
    int_stat_nolunc['subl'] = { 'real' : {}, 'fake' : {} }

    int_stat_nolunc['lead']['real']['loose'] = ufloat( int_stat['lead']['real']['loose'].n, 0.0 )
    int_stat_nolunc['lead']['fake']['loose'] = ufloat( int_stat['lead']['fake']['loose'].n, 0.0 )
    int_stat_nolunc['subl']['real']['loose'] = ufloat( int_stat['subl']['real']['loose'].n, 0.0 )
    int_stat_nolunc['subl']['fake']['loose'] = ufloat( int_stat['subl']['fake']['loose'].n, 0.0 )

    if int_stat['lead']['real']['loose'].n == 0 :
        eff_stat_loose['eff_R_L_lead'] = ufloat( 0.0, 0.0 ) 
    else :
        eff_stat_loose['eff_R_L_lead'] = 1.0 / (1.0 + (int_stat_notunc['lead']['real']['tight']/int_stat['lead']['real']['loose']) )
    if int_stat['lead']['fake']['loose'].n == 0 :
        eff_stat_loose['eff_F_L_lead'] = ufloat( 0.0, 0.0 ) 
    else :
        eff_stat_loose['eff_F_L_lead'] = 1.0 / (1.0 + (int_stat_notunc['lead']['fake']['tight']/int_stat['lead']['fake']['loose']) )
    if int_stat['subl']['real']['loose'].n == 0 :
        eff_stat_loose['eff_R_L_subl'] = ufloat( 0.0, 0.0 ) 
    else :
        eff_stat_loose['eff_R_L_subl'] = 1.0 / (1.0 + (int_stat_notunc['subl']['real']['tight']/int_stat['subl']['real']['loose']) )
    if int_stat['subl']['fake']['loose'].n == 0 :
        eff_stat_loose['eff_F_L_subl'] = ufloat( 0.0, 0.0 ) 
    else :
        eff_stat_loose['eff_F_L_subl'] = 1.0 / (1.0 + (int_stat_notunc['subl']['fake']['tight']/int_stat['subl']['fake']['loose']) )

    if int_stat_nolunc['lead']['real']['loose'].n == 0 :
        eff_stat_tight['eff_R_L_lead'] = ufloat( 0.0, 0.0 ) 
    else :
        eff_stat_tight['eff_R_L_lead'] = 1.0 / (1.0 + (int_stat['lead']['real']['tight']/int_stat_nolunc['lead']['real']['loose']) )
    if int_stat_nolunc['lead']['fake']['loose'].n == 0 :
        eff_stat_tight['eff_F_L_lead'] = ufloat( 0.0, 0.0 ) 
    else :
        eff_stat_tight['eff_F_L_lead'] = 1.0 / (1.0 + (int_stat['lead']['fake']['tight']/int_stat_nolunc['lead']['fake']['loose']) )
    if int_stat_nolunc['subl']['real']['loose'].n == 0 :
        eff_stat_tight['eff_R_L_subl'] = ufloat( 0.0, 0.0 ) 
    else :
        eff_stat_tight['eff_R_L_subl'] = 1.0 / (1.0 + (int_stat['subl']['real']['tight']/int_stat_nolunc['subl']['real']['loose']) )
    if int_stat_nolunc['subl']['fake']['loose'].n == 0 :
        eff_stat_tight['eff_F_L_subl'] = ufloat( 0.0, 0.0 ) 
    else :
        eff_stat_tight['eff_F_L_subl'] = 1.0 / (1.0 + (int_stat['subl']['fake']['tight']/int_stat_nolunc['subl']['fake']['loose']) )

    # make the systematic efficiencies
    # based on the systematics from the
    # input integrals
    if int_syst['lead']['real']['loose'].n == 0 :
        eff_syst_int['eff_R_L_lead'] = ufloat( 0.0, 0.0 ) 
    else :
        eff_syst_int['eff_R_L_lead'] = 1.0 / ( 1.0 + (int_syst['lead']['real']['tight']/int_syst['lead']['real']['loose']) )
    if int_syst['lead']['fake']['loose'].n == 0 :
        eff_syst_int['eff_F_L_lead'] = ufloat( 0.0, 0.0 ) 
    else :
        eff_syst_int['eff_F_L_lead'] = 1.0 / ( 1.0 + (int_syst['lead']['fake']['tight']/int_syst['lead']['fake']['loose']) )
    if int_syst['subl']['real']['loose'].n == 0 :
        eff_syst_int['eff_R_L_subl'] = ufloat( 0.0, 0.0 ) 
    else :
        eff_syst_int['eff_R_L_subl'] = 1.0 / ( 1.0 + (int_syst['subl']['real']['tight']/int_syst['subl']['real']['loose']) )
    if int_syst['subl']['fake']['loose'].n == 0 :
        eff_syst_int['eff_F_L_subl'] = ufloat( 0.0, 0.0 ) 
    else :
        eff_syst_int['eff_F_L_subl'] = 1.0 / ( 1.0 + (int_syst['subl']['fake']['tight']/int_syst['subl']['fake']['loose']) )

    # make the systematic efficiencies
    # based on the systematics from the
    # templates which are set below
    # first get the ratios
    if int_syst['lead']['real']['loose'].n == 0 :
        eff_syst_temp['eff_R_L_lead'] = ufloat( 0.0, 0.0 ) 
    else :
        eff_syst_temp['eff_R_L_lead'] = 1.0 / (1.0 + (int_syst['lead']['real']['tight']/int_syst['lead']['real']['loose']) )
    if int_syst['lead']['fake']['loose'].n == 0 :
        eff_syst_temp['eff_F_L_lead'] = ufloat( 0.0, 0.0 ) 
    else :
        eff_syst_temp['eff_F_L_lead'] = 1.0 / (1.0 + (int_syst['lead']['fake']['tight']/int_syst['lead']['fake']['loose']) )
    if int_syst['subl']['real']['loose'].n == 0 :
        eff_syst_temp['eff_R_L_subl'] = ufloat( 0.0, 0.0 ) 
    else :
        eff_syst_temp['eff_R_L_subl'] = 1.0 / (1.0 + (int_syst['subl']['real']['tight']/int_syst['subl']['real']['loose']) )
    if int_syst['subl']['fake']['loose'].n == 0 :
        eff_syst_temp['eff_F_L_subl'] = ufloat( 0.0, 0.0 ) 
    else :
        eff_syst_temp['eff_F_L_subl'] = 1.0 / (1.0 + (int_syst['subl']['fake']['tight']/int_syst['subl']['fake']['loose']) )

    # get the template uncertainties
    # simply overwrite the current
    # uncertainties so that we have
    # only the template uncertainties
    # eff_syst_temp['eff_R_L_lead'] = ufloat( eff_syst_temp['eff_R_L_lead'].n , math.fabs(eff_syst_temp['eff_R_L_lead'].n)*get_syst_uncertainty( var, 'RealTemplate%s'%systematics, lead_reg, lead_ptrange, 'real', 'loose' ), 'Template_lead_real_loose')
    # eff_syst_temp['eff_F_L_lead'] = ufloat( eff_syst_temp['eff_F_L_lead'].n , math.fabs(eff_syst_temp['eff_F_L_lead'].n)*get_syst_uncertainty( var, 'FakeTemplate%s'%systematics, lead_reg, lead_ptrange, 'fake', 'loose' ), 'Template_lead_fake_loose' )
    # eff_syst_temp['eff_R_L_subl'] = ufloat( eff_syst_temp['eff_R_L_subl'].n , math.fabs(eff_syst_temp['eff_R_L_subl'].n)*get_syst_uncertainty( var, 'RealTemplate%s'%systematics, subl_reg, subl_ptrange, 'real', 'loose' ), 'Template_subl_real_loose' )
    # eff_syst_temp['eff_F_L_subl'] = ufloat( eff_syst_temp['eff_F_L_subl'].n , math.fabs(eff_syst_temp['eff_F_L_subl'].n)*get_syst_uncertainty( var, 'FakeTemplate%s'%systematics, subl_reg, subl_ptrange, 'fake', 'loose' ), 'Template_subl_fake_loose' )

    return eff_stat_tight, eff_stat_loose, eff_syst_int, eff_syst_temp

def get_integral_and_error( hist, bins=None, name='' ) :

    err = ROOT.Double()
    if bins is None :
        val = hist.IntegralAndError( 1, hist.GetNbinsX(), err )
    else :
        if bins[1] is None :
            val = hist.IntegralAndError( bins[0], hist.GetNbinsX(), err )
        else :
            val = hist.IntegralAndError( bins[0], bins[1], err )

    return ufloat( val, err, name )


def get_template_integrals( templates, eff_cuts, lead_eta, subl_eta, var=None) :

    int_stat = {}
    int_stat['lead']={}
    int_stat['subl']={}
    int_stat['lead']['real']={}
    int_stat['subl']['real']={}
    int_stat['lead']['fake']={}
    int_stat['subl']['fake']={}
    
    int_syst = {}
    int_syst['lead']={}
    int_syst['subl']={}
    int_syst['lead']['real']={}
    int_syst['subl']['real']={}
    int_syst['lead']['fake']={}
    int_syst['subl']['fake']={}
    
    bins_lead_real_tight = ( templates['lead']['real'].GetXaxis().FindBin(eff_cuts['lead']['tight'][0]), templates['lead']['real'].GetXaxis().FindBin(eff_cuts['lead']['tight'][1]) )
    bins_lead_real_loose = ( templates['lead']['real'].GetXaxis().FindBin(eff_cuts['lead']['loose'][0]), templates['lead']['real'].GetXaxis().FindBin(eff_cuts['lead']['loose'][1]) )
    bins_lead_fake_tight = ( templates['lead']['fake'].GetXaxis().FindBin(eff_cuts['lead']['tight'][0]), templates['lead']['fake'].GetXaxis().FindBin(eff_cuts['lead']['tight'][1]) )
    bins_lead_fake_loose = ( templates['lead']['fake'].GetXaxis().FindBin(eff_cuts['lead']['loose'][0]), templates['lead']['fake'].GetXaxis().FindBin(eff_cuts['lead']['loose'][1]) )
    bins_subl_real_tight = ( templates['subl']['real'].GetXaxis().FindBin(eff_cuts['subl']['tight'][0]), templates['subl']['real'].GetXaxis().FindBin(eff_cuts['subl']['tight'][1]) )
    bins_subl_real_loose = ( templates['subl']['real'].GetXaxis().FindBin(eff_cuts['subl']['loose'][0]), templates['subl']['real'].GetXaxis().FindBin(eff_cuts['subl']['loose'][1]) )
    bins_subl_fake_tight = ( templates['subl']['fake'].GetXaxis().FindBin(eff_cuts['subl']['tight'][0]), templates['subl']['fake'].GetXaxis().FindBin(eff_cuts['subl']['tight'][1]) )
    bins_subl_fake_loose = ( templates['subl']['fake'].GetXaxis().FindBin(eff_cuts['subl']['loose'][0]), templates['subl']['fake'].GetXaxis().FindBin(eff_cuts['subl']['loose'][1]) )
    
    int_stat['lead']['real']['tight'] = get_integral_and_error( templates['lead']['real'], bins_lead_real_tight, 'Data_lead_real_tight' )
    int_stat['lead']['real']['loose'] = get_integral_and_error( templates['lead']['real'], bins_lead_real_loose, 'Data_lead_real_loose' )
    int_stat['lead']['fake']['tight'] = get_integral_and_error( templates['lead']['fake'], bins_lead_fake_tight, 'Data_lead_fake_tight' )
    int_stat['lead']['fake']['loose'] = get_integral_and_error( templates['lead']['fake'], bins_lead_fake_loose, 'Data_lead_fake_loose' )
    int_stat['subl']['real']['tight'] = get_integral_and_error( templates['subl']['real'], bins_subl_real_tight, 'Data_subl_real_tight' )
    int_stat['subl']['real']['loose'] = get_integral_and_error( templates['subl']['real'], bins_subl_real_loose, 'Data_subl_real_loose' )
    int_stat['subl']['fake']['tight'] = get_integral_and_error( templates['subl']['fake'], bins_subl_fake_tight, 'Data_subl_fake_tight' )
    int_stat['subl']['fake']['loose'] = get_integral_and_error( templates['subl']['fake'], bins_subl_fake_loose, 'Data_subl_fake_loose' )

    # If running with systematics, set the data systs to zero
    # May need to implement non-zero systematics for data in the future
    # The overall template systematics should not be set here
    int_syst['lead']['real']['tight'] = ufloat(int_stat['lead']['real']['tight'].n, 0.0 , 'Data_lead_real_tight' )
    int_syst['lead']['real']['loose'] = ufloat(int_stat['lead']['real']['loose'].n, 0.0 , 'Data_lead_real_loose' )
    int_syst['lead']['fake']['tight'] = ufloat(int_stat['lead']['fake']['tight'].n, 0.0 , 'Data_lead_fake_tight' )
    int_syst['lead']['fake']['loose'] = ufloat(int_stat['lead']['fake']['loose'].n, 0.0 , 'Data_lead_fake_loose' )
    int_syst['subl']['real']['tight'] = ufloat(int_stat['subl']['real']['tight'].n, 0.0 , 'Data_subl_real_tight' )
    int_syst['subl']['real']['loose'] = ufloat(int_stat['subl']['real']['loose'].n, 0.0 , 'Data_subl_real_loose' )
    int_syst['subl']['fake']['tight'] = ufloat(int_stat['subl']['fake']['tight'].n, 0.0 , 'Data_subl_fake_tight' )
    int_syst['subl']['fake']['loose'] = ufloat(int_stat['subl']['fake']['loose'].n, 0.0 , 'Data_subl_fake_loose' )

    # subtract Zg events from fake photon template
    if templates['lead']['bakg'] is not None :
        bkg_int_tight = get_integral_and_error(templates['lead']['bakg'], bins_lead_fake_tight, 'Background_lead_fake_tight')
        bkg_int_loose = get_integral_and_error(templates['lead']['bakg'], bins_lead_fake_loose, 'Background_lead_fake_loose')
        
        syst_bkg_int_tight = ufloat(0.0, math.fabs(bkg_int_tight.n,)*get_syst_uncertainty('bakg'), 'Background_lead_fake_tight')
        syst_bkg_int_loose = ufloat(0.0, math.fabs(bkg_int_loose.n,)*get_syst_uncertainty('bakg'), 'Background_lead_fake_loose')
    
        int_syst['lead']['fake']['tight'] = int_syst['lead']['fake']['tight'] - syst_bkg_int_tight
        int_syst['lead']['fake']['loose'] = int_syst['lead']['fake']['loose'] - syst_bkg_int_loose
        
    if templates['subl']['bakg'] is not None :
        bkg_int_tight = get_integral_and_error(templates['subl']['bakg'], bins_subl_fake_tight, 'Background_subl_fake_tight')
        bkg_int_loose = get_integral_and_error(templates['subl']['bakg'], bins_subl_fake_loose, 'Background_subl_fake_loose')
        
        syst_bkg_int_tight = ufloat(0.0, math.fabs(bkg_int_tight.n,)*get_syst_uncertainty('bakg'), 'Background_subl_fake_tight')
        syst_bkg_int_loose = ufloat(0.0, math.fabs(bkg_int_loose.n,)*get_syst_uncertainty('bakg'), 'Background_subl_fake_loose')
        
        int_syst['subl']['fake']['tight'] = int_syst['subl']['fake']['tight'] - syst_bkg_int_tight
        int_syst['subl']['fake']['loose'] = int_syst['subl']['fake']['loose'] - syst_bkg_int_loose
        
    return int_stat, int_syst

def run_fit( data, efficiencies ) :

    # make the matrix
    matrix = generate_eff_matrix( efficiencies, ndim=len(data) )

    #do the fit!  Invert the matrix and multiply the by counts vectors
    if len( data ) == 3 :
        results = solve_matrix_eq( matrix, [data['TL'], data['LT'], data['LL']] )
    elif len(data) == 4 :
        results = solve_matrix_eq( matrix, [data['TT'],data['TL'], data['LT'], data['LL']] )

    return results


def generate_eff_matrix( eff_dic, ndim=4 ) :
#    print 'ndim = ', ndim

    eff_matrix = [ [ eff_dic['eff_RF_TL'], eff_dic['eff_FR_TL'], eff_dic['eff_FF_TL'] ],
                   [ eff_dic['eff_RF_LT'], eff_dic['eff_FR_LT'], eff_dic['eff_FF_LT'] ], 
                   [ eff_dic['eff_RF_LL'], eff_dic['eff_FR_LL'], eff_dic['eff_FF_LL'] ] ] 
    
    if ndim == 4 :
        eff_matrix = [ [ eff_dic['eff_RR_TT'], eff_dic['eff_RF_TT'], eff_dic['eff_FR_TT'], eff_dic['eff_FF_TT'] ], 
                       [ eff_dic['eff_RR_TL'], eff_matrix[0][0]    , eff_matrix[0][1]    , eff_matrix[0][2]     ],
                       [ eff_dic['eff_RR_LT'], eff_matrix[1][0]    , eff_matrix[1][1]    , eff_matrix[1][2]     ],
                       [ eff_dic['eff_RR_LL'], eff_matrix[2][0]    , eff_matrix[2][1]    , eff_matrix[2][2]     ] ]

    elif ndim != 3 :
        print 'Only Dim 3 and 4 implemented'
        return None

    return eff_matrix

##OKAY
def solve_matrix_eq( matrix_ntries, vector_entries ) :

    ms = []
    mn = []
    for row in matrix_ntries :
        ms_row = []
        mn_row = []
        for col in row :
            ms_row.append( col.s )
            mn_row.append( col.n )
        ms.append( ms_row )
        mn.append( mn_row )

    #Matrices of numbers with uncertainties are best created in one of two ways.
    #The first way is similar to using uarray():
    matrix = unumpy.umatrix( mn, ms )

#    print matrix

    vs = []
    vn = []
    for row in vector_entries :
        vn.append( [ row.n ] )
        vs.append( [ row.s ] )

    vector = unumpy.umatrix( vn, vs )
    
    inv_matrix = None
    try :
        inv_matrix = matrix.getI()
    except :
        print 'Failed to invert matrix, aborting'
        return unumpy.umatrix( [ [1]*len(vs) ], [ [0]*len(vs) ] )

    return inv_matrix*vector

def run_fit_manual(data, eff) :

    alpha_rf = ( (1.0/( (1-eff['eff_R_L_lead'])*eff['eff_F_L_subl']*eff['eff_F_L_lead']*(1-eff['eff_R_L_subl'])*eff['eff_F_L_lead']*eff['eff_F_L_subl'] 
                      + (1-eff['eff_F_L_lead'])*eff['eff_R_L_subl']*eff['eff_F_L_lead']*(1-eff['eff_F_L_subl'])*eff['eff_R_L_lead']*eff['eff_F_L_subl']
                      + (1-eff['eff_F_L_lead'])*eff['eff_F_L_subl']*eff['eff_R_L_lead']*(1-eff['eff_F_L_subl'])*eff['eff_F_L_lead']*eff['eff_R_L_subl']
                      - (1-eff['eff_F_L_lead'])*eff['eff_F_L_subl']*eff['eff_F_L_lead']*(1-eff['eff_R_L_subl'])*eff['eff_R_L_lead']*eff['eff_F_L_subl']
                      - (1-eff['eff_F_L_lead'])*eff['eff_R_L_subl']*eff['eff_R_L_lead']*(1-eff['eff_F_L_subl'])*eff['eff_F_L_lead']*eff['eff_F_L_subl']
                      - (1-eff['eff_R_L_lead'])*eff['eff_F_L_subl']*eff['eff_F_L_lead']*(1-eff['eff_F_L_subl'])*eff['eff_F_L_lead']*eff['eff_R_L_subl'] ) )
              * ( (  eff['eff_F_L_lead']*(1-eff['eff_R_L_subl'])*eff['eff_F_L_lead']*eff['eff_F_L_subl']
                   - eff['eff_F_L_lead']*(1-eff['eff_F_L_subl'])*eff['eff_F_L_lead']*eff['eff_R_L_subl'] )*data['TL']
                + (  (1-eff['eff_F_L_lead'])*eff['eff_F_L_subl']*eff['eff_F_L_lead']*eff['eff_R_L_subl']
                   - (1-eff['eff_F_L_lead'])*eff['eff_R_L_subl']*eff['eff_F_L_lead']*eff['eff_F_L_subl'] )*data['LT']
                + (   (1-eff['eff_F_L_lead'])*eff['eff_R_L_subl']*eff['eff_F_L_lead']*(1-eff['eff_F_L_subl'])
                    - (1-eff['eff_F_L_lead'])*eff['eff_F_L_subl']*eff['eff_F_L_lead']*(1-eff['eff_R_L_subl']) )*data['LL']
              ) )

    alpha_fr = ( (1.0/( (1-eff['eff_R_L_lead'])*eff['eff_F_L_subl']*eff['eff_F_L_lead']*(1-eff['eff_R_L_subl'])*eff['eff_F_L_lead']*eff['eff_F_L_subl']     
                    + (1-eff['eff_F_L_lead'])*eff['eff_R_L_subl']*eff['eff_F_L_lead']*(1-eff['eff_F_L_subl'])*eff['eff_R_L_lead']*eff['eff_F_L_subl']
                    + (1-eff['eff_F_L_lead'])*eff['eff_F_L_subl']*eff['eff_R_L_lead']*(1-eff['eff_F_L_subl'])*eff['eff_F_L_lead']*eff['eff_R_L_subl']
                    - (1-eff['eff_F_L_lead'])*eff['eff_F_L_subl']*eff['eff_F_L_lead']*(1-eff['eff_R_L_subl'])*eff['eff_R_L_lead']*eff['eff_F_L_subl']
                    - (1-eff['eff_F_L_lead'])*eff['eff_R_L_subl']*eff['eff_R_L_lead']*(1-eff['eff_F_L_subl'])*eff['eff_F_L_lead']*eff['eff_F_L_subl']
                    - (1-eff['eff_R_L_lead'])*eff['eff_F_L_subl']*eff['eff_F_L_lead']*(1-eff['eff_F_L_subl'])*eff['eff_F_L_lead']*eff['eff_R_L_subl'] ) )
            * ( (   eff['eff_F_L_lead']*(1-eff['eff_F_L_subl'])*eff['eff_R_L_lead']*eff['eff_F_L_subl']
                  - eff['eff_R_L_lead']*(1-eff['eff_F_L_subl'])*eff['eff_F_L_lead']*eff['eff_F_L_subl'] )*data['TL'] 
              + (   (1-eff['eff_R_L_lead'])*eff['eff_F_L_subl']*eff['eff_F_L_lead']*eff['eff_F_L_subl']
                  - (1-eff['eff_F_L_lead'])*eff['eff_F_L_subl']*eff['eff_R_L_lead']*eff['eff_F_L_subl'] )*data['LT']
              + (   (1-eff['eff_F_L_lead'])*eff['eff_F_L_subl']*eff['eff_R_L_lead']*(1-eff['eff_F_L_subl'])
                  - (1-eff['eff_R_L_lead'])*eff['eff_F_L_subl']*eff['eff_F_L_lead']*(1-eff['eff_F_L_subl']) ) *data['LL']
              ) )
                
    alpha_ff = ( (1.0/( (1-eff['eff_R_L_lead'])*eff['eff_F_L_subl']*eff['eff_F_L_lead']*(1-eff['eff_R_L_subl'])*eff['eff_F_L_lead']*eff['eff_F_L_subl']    
                    + (1-eff['eff_F_L_lead'])*eff['eff_R_L_subl']*eff['eff_F_L_lead']*(1-eff['eff_F_L_subl'])*eff['eff_R_L_lead']*eff['eff_F_L_subl']
                    + (1-eff['eff_F_L_lead'])*eff['eff_F_L_subl']*eff['eff_R_L_lead']*(1-eff['eff_F_L_subl'])*eff['eff_F_L_lead']*eff['eff_R_L_subl']
                    - (1-eff['eff_F_L_lead'])*eff['eff_F_L_subl']*eff['eff_F_L_lead']*(1-eff['eff_R_L_subl'])*eff['eff_R_L_lead']*eff['eff_F_L_subl']
                    - (1-eff['eff_F_L_lead'])*eff['eff_R_L_subl']*eff['eff_R_L_lead']*(1-eff['eff_F_L_subl'])*eff['eff_F_L_lead']*eff['eff_F_L_subl']
                    - (1-eff['eff_R_L_lead'])*eff['eff_F_L_subl']*eff['eff_F_L_lead']*(1-eff['eff_F_L_subl'])*eff['eff_F_L_lead']*eff['eff_R_L_subl'] ) )
            * ( (   eff['eff_R_L_lead']*(1-eff['eff_F_L_subl'])*eff['eff_F_L_lead']*eff['eff_R_L_subl']
                  - eff['eff_F_L_lead']*(1-eff['eff_R_L_subl'])*eff['eff_R_L_lead']*eff['eff_F_L_subl'] )*data['TL']
              + (  (1-eff['eff_F_L_lead'])*eff['eff_R_L_subl']*eff['eff_R_L_lead']*eff['eff_F_L_subl']
                 - (1-eff['eff_R_L_lead'])*eff['eff_F_L_subl']*eff['eff_F_L_lead']*eff['eff_R_L_subl'] ) * data['LT']
              + (   (1-eff['eff_R_L_lead'])*eff['eff_F_L_subl']*eff['eff_F_L_lead']*(1-eff['eff_R_L_subl'])
                  - (1-eff['eff_F_L_lead'])*eff['eff_R_L_subl']*eff['eff_R_L_lead']*(1-eff['eff_F_L_subl']) )*data['LL']
              ) )


    nPred_RF_TT = ( (1-eff['eff_R_L_lead'])*(1-eff['eff_F_L_subl'])* 
                      ( (1.0/( (1-eff['eff_R_L_lead'])*eff['eff_F_L_subl']*eff['eff_F_L_lead']*(1-eff['eff_R_L_subl'])*eff['eff_F_L_lead']*eff['eff_F_L_subl'] 
                      + (1-eff['eff_F_L_lead'])*eff['eff_R_L_subl']*eff['eff_F_L_lead']*(1-eff['eff_F_L_subl'])*eff['eff_R_L_lead']*eff['eff_F_L_subl']
                      + (1-eff['eff_F_L_lead'])*eff['eff_F_L_subl']*eff['eff_R_L_lead']*(1-eff['eff_F_L_subl'])*eff['eff_F_L_lead']*eff['eff_R_L_subl']
                      - (1-eff['eff_F_L_lead'])*eff['eff_F_L_subl']*eff['eff_F_L_lead']*(1-eff['eff_R_L_subl'])*eff['eff_R_L_lead']*eff['eff_F_L_subl']
                      - (1-eff['eff_F_L_lead'])*eff['eff_R_L_subl']*eff['eff_R_L_lead']*(1-eff['eff_F_L_subl'])*eff['eff_F_L_lead']*eff['eff_F_L_subl']
                      - (1-eff['eff_R_L_lead'])*eff['eff_F_L_subl']*eff['eff_F_L_lead']*(1-eff['eff_F_L_subl'])*eff['eff_F_L_lead']*eff['eff_R_L_subl'] ) )
              * ( (  eff['eff_F_L_lead']*(1-eff['eff_R_L_subl'])*eff['eff_F_L_lead']*eff['eff_F_L_subl']
                   - eff['eff_F_L_lead']*(1-eff['eff_F_L_subl'])*eff['eff_F_L_lead']*eff['eff_R_L_subl'] )*data['TL']
                + (  (1-eff['eff_F_L_lead'])*eff['eff_F_L_subl']*eff['eff_F_L_lead']*eff['eff_R_L_subl']
                   - (1-eff['eff_F_L_lead'])*eff['eff_R_L_subl']*eff['eff_F_L_lead']*eff['eff_F_L_subl'] )*data['LT']
                + (   (1-eff['eff_F_L_lead'])*eff['eff_R_L_subl']*eff['eff_F_L_lead']*(1-eff['eff_F_L_subl'])
                    - (1-eff['eff_F_L_lead'])*eff['eff_F_L_subl']*eff['eff_F_L_lead']*(1-eff['eff_R_L_subl']) )*data['LL']
              ) ) )

    nPred_FR_TT = ( (1-eff['eff_F_L_lead'])*(1-eff['eff_R_L_subl'])* 
                   ( (1.0/( (1-eff['eff_R_L_lead'])*eff['eff_F_L_subl']*eff['eff_F_L_lead']*(1-eff['eff_R_L_subl'])*eff['eff_F_L_lead']*eff['eff_F_L_subl']     
                    + (1-eff['eff_F_L_lead'])*eff['eff_R_L_subl']*eff['eff_F_L_lead']*(1-eff['eff_F_L_subl'])*eff['eff_R_L_lead']*eff['eff_F_L_subl']
                    + (1-eff['eff_F_L_lead'])*eff['eff_F_L_subl']*eff['eff_R_L_lead']*(1-eff['eff_F_L_subl'])*eff['eff_F_L_lead']*eff['eff_R_L_subl']
                    - (1-eff['eff_F_L_lead'])*eff['eff_F_L_subl']*eff['eff_F_L_lead']*(1-eff['eff_R_L_subl'])*eff['eff_R_L_lead']*eff['eff_F_L_subl']
                    - (1-eff['eff_F_L_lead'])*eff['eff_R_L_subl']*eff['eff_R_L_lead']*(1-eff['eff_F_L_subl'])*eff['eff_F_L_lead']*eff['eff_F_L_subl']
                    - (1-eff['eff_R_L_lead'])*eff['eff_F_L_subl']*eff['eff_F_L_lead']*(1-eff['eff_F_L_subl'])*eff['eff_F_L_lead']*eff['eff_R_L_subl'] ) )
            * ( (   eff['eff_F_L_lead']*(1-eff['eff_F_L_subl'])*eff['eff_R_L_lead']*eff['eff_F_L_subl']
                  - eff['eff_R_L_lead']*(1-eff['eff_F_L_subl'])*eff['eff_F_L_lead']*eff['eff_F_L_subl'] )*data['TL'] 
              + (   (1-eff['eff_R_L_lead'])*eff['eff_F_L_subl']*eff['eff_F_L_lead']*eff['eff_F_L_subl']
                  - (1-eff['eff_F_L_lead'])*eff['eff_F_L_subl']*eff['eff_R_L_lead']*eff['eff_F_L_subl'] )*data['LT']
              + (   (1-eff['eff_F_L_lead'])*eff['eff_F_L_subl']*eff['eff_R_L_lead']*(1-eff['eff_F_L_subl'])
                  - (1-eff['eff_R_L_lead'])*eff['eff_F_L_subl']*eff['eff_F_L_lead']*(1-eff['eff_F_L_subl']) ) *data['LL']
              ) ) )

    nPred_FF_TT = ( (1-eff['eff_F_L_lead'])*(1-eff['eff_F_L_subl'])* 
                  ( (1.0/( (1-eff['eff_R_L_lead'])*eff['eff_F_L_subl']*eff['eff_F_L_lead']*(1-eff['eff_R_L_subl'])*eff['eff_F_L_lead']*eff['eff_F_L_subl']    
                    + (1-eff['eff_F_L_lead'])*eff['eff_R_L_subl']*eff['eff_F_L_lead']*(1-eff['eff_F_L_subl'])*eff['eff_R_L_lead']*eff['eff_F_L_subl']
                    + (1-eff['eff_F_L_lead'])*eff['eff_F_L_subl']*eff['eff_R_L_lead']*(1-eff['eff_F_L_subl'])*eff['eff_F_L_lead']*eff['eff_R_L_subl']
                    - (1-eff['eff_F_L_lead'])*eff['eff_F_L_subl']*eff['eff_F_L_lead']*(1-eff['eff_R_L_subl'])*eff['eff_R_L_lead']*eff['eff_F_L_subl']
                    - (1-eff['eff_F_L_lead'])*eff['eff_R_L_subl']*eff['eff_R_L_lead']*(1-eff['eff_F_L_subl'])*eff['eff_F_L_lead']*eff['eff_F_L_subl']
                    - (1-eff['eff_R_L_lead'])*eff['eff_F_L_subl']*eff['eff_F_L_lead']*(1-eff['eff_F_L_subl'])*eff['eff_F_L_lead']*eff['eff_R_L_subl'] ) )
            * ( (   eff['eff_R_L_lead']*(1-eff['eff_F_L_subl'])*eff['eff_F_L_lead']*eff['eff_R_L_subl']
                  - eff['eff_F_L_lead']*(1-eff['eff_R_L_subl'])*eff['eff_R_L_lead']*eff['eff_F_L_subl'] )*data['TL']
              + (  (1-eff['eff_F_L_lead'])*eff['eff_R_L_subl']*eff['eff_R_L_lead']*eff['eff_F_L_subl']
                 - (1-eff['eff_R_L_lead'])*eff['eff_F_L_subl']*eff['eff_F_L_lead']*eff['eff_R_L_subl'] ) * data['LT']
              + (   (1-eff['eff_R_L_lead'])*eff['eff_F_L_subl']*eff['eff_F_L_lead']*(1-eff['eff_R_L_subl'])
                  - (1-eff['eff_F_L_lead'])*eff['eff_R_L_subl']*eff['eff_R_L_lead']*(1-eff['eff_F_L_subl']) )*data['LL']
              ) ) )
       
    return {'alpha_RF' : alpha_rf, 'alpha_FR' : alpha_fr, 'alpha_FF' : alpha_ff, 'pred_RF_TT' : nPred_RF_TT, 'pred_FR_TT' : nPred_FR_TT, 'pred_FF_TT' : nPred_FF_TT }


def save_results( results, outputDir, namePostfix='' ) :

    if outputDir is None :
        return

    fname = outputDir + '/results%s.pickle' %namePostfix

    if not os.path.isdir( os.path.dirname( fname ) ) :
        os.makedirs( os.path.dirname( fname ) )
        
    file = open( fname, 'w' )
    pickle.dump( results, file )
    file.close()

######## Create 2x2 histogram ######
def create_2dEff_hist(eff_2d_stat, lead_eta, subl_eta, outputDir, ndim) :

    eff_hist = {}
    eff_can = {}
    eff_latex = {}

    if ndim == 4 :
        photons = ['RR', 'RF', 'FR', 'FF']
    else :
        photons = ['RF', 'FR', 'FF']

    for pho in photons :
        eff_hist[pho] = ROOT.TH2F(pho, '', 2, 0, 2, 2, 0, 2)
        eff_can[pho] = ROOT.TCanvas(pho+'_can', '', 800, 800)
        rname = ''

        if pho[0] == 'R' :
            rname += 'real+'
        if pho[0] == 'F' :
            rname += 'fake+'
        if pho[1] == 'R' :
            rname += 'real'
        if pho[1] == 'F' :
            rname += 'fake'

        eff_latex[pho] = ROOT.TLatex(-0.5, 0.9, rname + ' template')
        eff_latex[pho].SetNDC()
        eff_latex[pho].SetX(0.05)
        eff_latex[pho].SetY(0.95)
        eff_latex[pho].SetTextSize(0.06)

    for val in eff_can.values() :
        val.SetLeftMargin(0.15)
        val.SetTopMargin(0.05)

    for val in eff_hist.values() :
        val.GetXaxis().SetBinLabel( 2, 'P' )
        val.GetXaxis().SetBinLabel( 1, 'F' )
        val.GetYaxis().SetBinLabel( 2, 'P' )
        val.GetYaxis().SetBinLabel( 1, 'F' )
        val.GetXaxis().SetTitle( 'Leading photon' )
        val.GetYaxis().SetTitle( 'Trailing photon' )
        val.GetZaxis().SetTitle( 'Probability' )
        val.GetXaxis().SetLabelSize( 0.07 )
        val.GetYaxis().SetLabelSize( 0.07 )
        val.GetZaxis().SetLabelSize( 0.054 )
        val.GetXaxis().SetTitleSize( 0.07 )
        val.GetYaxis().SetTitleSize( 0.07 )
        val.GetZaxis().SetTitleSize( 0.060 )
        val.GetXaxis().CenterTitle()
        val.GetYaxis().CenterTitle()
        val.GetXaxis().SetTitleOffset( 1.0 )
        val.GetYaxis().SetTitleOffset( 1.0 )
        val.GetZaxis().SetTitleOffset( 1.1 )
        val.SetLineColor(6)
        val.SetLineWidth(2)
        val.SetMinimum(0)
        val.SetMaximum(0.8)
        val.SetStats(0)

    for pho in photons :
        
        eff_can[pho].cd()

        eff_hist[pho].SetBinContent(1, 1, eff_2d_stat['eff_%s_LL'%pho].n)
        eff_hist[pho].SetBinContent(1, 2, eff_2d_stat['eff_%s_LT'%pho].n)
        eff_hist[pho].SetBinContent(2, 1, eff_2d_stat['eff_%s_TL'%pho].n)
        eff_hist[pho].SetBinContent(2, 2, eff_2d_stat['eff_%s_TT'%pho].n)
        eff_hist[pho].Draw('lego4FB')

        eff_can[pho].SaveAs(outputDir+'/%s_%s_%s.pdf' %(lead_eta, subl_eta, pho) )
####################################


####################################
_fitVar = ['Sieie', 'ChIso']
fitVar = _fitVar[1]
iso = 'nomi'

## do_nominal_fit ##
path = '/uscmst1b_scratch/lpc1/old_scratch/lpceg/yurii/EnSc/sachiko'
dir_tmp = '/Zgg/NEW_WORK_09010/Zgg/Muon/TempHist/TEMP_ROOT_0111/'

# Input_Zgg       = path + dir_tmp + 'ROOT_Zgg/muon_Zgg_'+fitVar+'_533.root'
# Input_Zg        = path + dir_tmp + 'ROOT_Sig/muon_Sig_'+fitVar+'_'+iso+'.root'
# Input_data_Zjet = path + dir_tmp + 'ROOT_Bkg/muon_subtrBkg_'+fitVar+'_'+iso+'.root'
# Input_FF_corr   = path + dir_tmp + 'ROOT_CorrFF/Temp_Corr_FF_nomi.root'
# Input_mc_Zjet   = path + dir_tmp + 'ROOT_BkgMC/muon_Bkg_'+fitVar+'_nomi.root'

Input_Zgg       = path + dir_tmp + 'ROOT_allPt/muon_allpT_Zgg_'+fitVar+'_nomi.root'
Input_Zg        = path + dir_tmp + 'ROOT_allPt/muon_allpT_Sig_'+fitVar+'_nomi.root'
Input_data_Zjet = path + dir_tmp + 'ROOT_allPt/muon_allpT_subtrBkg_'+fitVar+'_nomi.root'
Input_mc_Zjet   = path + dir_tmp + 'ROOT_allPt/muon_allpT_Bkg_'+fitVar+'_nomi.root'
Input_FF_corr   = path + dir_tmp + 'ROOT_allPt/Temp_allPt_Corr_FF_nomi.root'

#pt_bins = ['1', '2', '3']
pt_bins = ['1']
for bin in pt_bins :

    templates_reg = {}
    templates_reg['EB'] = {}
    templates_reg['EE'] = {}
    templates_reg['EB']['real'] = getWeightedHist(Input_Zg, 'h_'+fitVar+'_EB_'+bin)
    templates_reg['EE']['real'] = getWeightedHist(Input_Zg, 'h_'+fitVar+'_EE_'+bin)
    templates_reg['EB']['fake'] = getWeightedHist(Input_data_Zjet, 'h_'+fitVar+'_EB_'+bin)
    templates_reg['EE']['fake'] = getWeightedHist(Input_data_Zjet, 'h_'+fitVar+'_EE_'+bin)
    templates_reg['EB']['bakg'] = getWeightedHist(Input_mc_Zjet, 'h_'+fitVar+'_EB_'+bin)
    templates_reg['EE']['bakg'] = getWeightedHist(Input_mc_Zjet, 'h_'+fitVar+'_EE_'+bin)

    templatrs_corr = None
    ffcorr = get_default_cuts(fitVar)

    templates_corr = {'leadPass' : {}, 'leadFail' : {}}
    templates_corr['leadFail'][('EB','EB')] = getWeightedHist(Input_FF_corr, 'h_'+fitVar+'_LeadFail_EB_EB_'+bin)
    templates_corr['leadPass'][('EB','EB')] = getWeightedHist(Input_FF_corr, 'h_'+fitVar+'_LeadPass_EB_EB_'+bin)
    templates_corr['leadFail'][('EB','EE')] = getWeightedHist(Input_FF_corr, 'h_'+fitVar+'_LeadFail_EB_EE_'+bin)
    templates_corr['leadPass'][('EB','EE')] = getWeightedHist(Input_FF_corr, 'h_'+fitVar+'_LeadPass_EB_EE_'+bin)
    templates_corr['leadFail'][('EE','EB')] = getWeightedHist(Input_FF_corr, 'h_'+fitVar+'_LeadFail_EE_EB_'+bin)
    templates_corr['leadPass'][('EE','EB')] = getWeightedHist(Input_FF_corr, 'h_'+fitVar+'_LeadPass_EE_EB_'+bin)

    # templates_corr['leadFail'][('EB','EB')] = getWeightedHist(Input_FF_corr, 'h_'+fitVar+'_EB_EB_'+bin,ffcorr['EB']['loose'][0],ffcorr['EB']['loose'][1])
    # templates_corr['leadPass'][('EB','EB')] = getWeightedHist(Input_FF_corr, 'h_'+fitVar+'_EB_EB_'+bin,ffcorr['EB']['tight'][0],ffcorr['EB']['tight'][1])
    # templates_corr['leadFail'][('EB','EE')] = getWeightedHist(Input_FF_corr, 'h_'+fitVar+'_EB_EE_'+bin,ffcorr['EB']['loose'][0],ffcorr['EB']['loose'][1])
    # templates_corr['leadPass'][('EB','EE')] = getWeightedHist(Input_FF_corr, 'h_'+fitVar+'_EB_EE_'+bin,ffcorr['EB']['tight'][0],ffcorr['EB']['tight'][1])
    # templates_corr['leadFail'][('EE','EB')] = getWeightedHist(Input_FF_corr, 'h_'+fitVar+'_EE_EB_'+bin,ffcorr['EE']['loose'][0],ffcorr['EE']['loose'][1])
    # templates_corr['leadPass'][('EE','EB')] = getWeightedHist(Input_FF_corr, 'h_'+fitVar+'_EE_EB_'+bin,ffcorr['EE']['tight'][0],ffcorr['EE']['tight'][1])

    gg_hist_bothiso_name = {'leadPass' : {}, 'leadFail' : {}}
    gg_hist_bothiso_name['leadPass'][('EB','EB')] = getWeightedHist(Input_Zgg, 'h_'+fitVar+'_LeadPass_EB_EB_'+bin)
    gg_hist_bothiso_name['leadFail'][('EB','EB')] = getWeightedHist(Input_Zgg, 'h_'+fitVar+'_LeadFail_EB_EB_'+bin)
    gg_hist_bothiso_name['leadPass'][('EB','EE')] = getWeightedHist(Input_Zgg, 'h_'+fitVar+'_LeadPass_EB_EE_'+bin)
    gg_hist_bothiso_name['leadFail'][('EB','EE')] = getWeightedHist(Input_Zgg, 'h_'+fitVar+'_LeadFail_EB_EE_'+bin)
    gg_hist_bothiso_name['leadPass'][('EE','EB')] = getWeightedHist(Input_Zgg, 'h_'+fitVar+'_LeadPass_EE_EB_'+bin)
    gg_hist_bothiso_name['leadFail'][('EE','EB')] = getWeightedHist(Input_Zgg, 'h_'+fitVar+'_LeadFail_EE_EB_'+bin)

    regions = [ ('EB', 'EB'), ('EB', 'EE'), ('EE', 'EB')]
    for reg in regions :

        print '######################################################################'
        print 'pt bin = %s,   lead eta = %s,   subl eta = %s' % (bin, reg[0], reg[1])
        print '######################################################################'

        templates = {}
        templates['lead'] = {}
        templates['subl'] = {}
        templates['lead']['real'] = templates_reg[reg[0]]['real']
        templates['subl']['real'] = templates_reg[reg[1]]['real']
        templates['lead']['fake'] = templates_reg[reg[0]]['fake']
        templates['subl']['fake'] = templates_reg[reg[1]]['fake']
        templates['lead']['bakg'] = templates_reg[reg[0]]['bakg']
        templates['subl']['bakg'] = templates_reg[reg[1]]['bakg']

        # templates for colleration of fake-fake
        templates_corr_inc = None
        if templates_corr is not None :
            templates_corr_inc = { }
            templates_corr_inc['leadPass'] = templates_corr['leadPass'][(reg[0], reg[1])]
            templates_corr_inc['leadFail'] = templates_corr['leadFail'][(reg[0], reg[1])]
        
        # gg_histogram (2D histogram between lead and sublead photon with no iso cut but eta range for two photons)
        if gg_hist_bothiso_name is not None :
            gg_hist_bothiso = { }
            gg_hist_bothiso['leadPass'] = gg_hist_bothiso_name['leadPass'][(reg[0], reg[1])]
            gg_hist_bothiso['leadFail'] = gg_hist_bothiso_name['leadFail'][(reg[0], reg[1])]

        # gg_hist_name = 'h_'+fitVar+'_'+reg[0]+'_'+reg[1]+'_'+bin
        # gg_hist = get2PhoHist(Input_Zgg, gg_hist_name)

        namePostfix = '_%s_%s' %( reg[0], reg[1] )

        outputDir = '/uscmst1b_scratch/lpc1/old_scratch/lpceg/yurii/EnSc/sachiko/Zgg/NEW_WORK_09010/Zgg/Muon/Fitting/FIT_1023/' + fitVar + '_' + reg[0] + reg[1] + '_' + bin
        
#        (results_stat_data, results_stat_temp, results_syst_bkg, results_syst_temp) = run_diphoton_fit(templates, gg_hist, reg[0], reg[1], templates_corr_inc, 4, fitVar, outputDir)
        (text_results_stat_dataSR,text_results_stat_dataSB, text_results_stat_temp_tight, text_results_stat_temp_loose, text_results_stat_ff,text_results_syst_bkg, text_results_syst_temp) = run_diphoton_fit(templates, gg_hist_bothiso, reg[0], reg[1], templates_corr_inc, 4, fitVar, outputDir)
