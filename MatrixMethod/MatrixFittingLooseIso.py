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

_real_syst = { 'EB' : (0.082641), 'EE' : (0.247633) } 
_fake_syst = { 'EB' : (0.025057), 'EE' : (0.078972) }

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


def run_corrected_diphoton_fit(templates_leadiso, templates_subliso, gg_hist_leadiso, gg_hist_subliso, gg_hist_bothiso, lead_eta, subl_eta, templates_corr, fitvar, ndim, outputDir=None) :
    
    accept_reg = ['EB', 'EE']
    if lead_eta not in accept_reg :
        print 'Lead region does not make sense'
        return
    if subl_eta not in accept_reg :
        print 'Subl region does not make sense'
        return

    # get the defaults
    cuts = get_default_cuts(var=fitvar)
    
    bins_lead_loose = ( None, None )
    bins_lead_tight = ( None, None )
    if isinstance( gg_hist_leadiso, dict ) :

        bins_subl_tight = ( gg_hist_leadiso['leadPass'].GetXaxis().FindBin( cuts[subl_eta]['tight'][0] ), gg_hist_leadiso['leadPass'].GetXaxis().FindBin( cuts[subl_eta]['tight'][1] ) )
        bins_subl_loose = ( gg_hist_leadiso['leadPass'].GetXaxis().FindBin( cuts[subl_eta]['loose'][0] ), gg_hist_leadiso['leadPass'].GetXaxis().FindBin( cuts[subl_eta]['loose'][1] ) )

        # Integrate to the the data in the four regions
        Ndata_TT_leadiso = gg_hist_leadiso['leadPass'].Integral( bins_subl_tight[0], bins_subl_tight[1] )
        Ndata_TL_leadiso = gg_hist_leadiso['leadPass'].Integral( bins_subl_loose[0], bins_subl_loose[1] )
        Ndata_LT_leadiso = gg_hist_leadiso['leadFail'].Integral( bins_subl_tight[0], bins_subl_tight[1] )
        Ndata_LL_leadiso = gg_hist_leadiso['leadFail'].Integral( bins_subl_loose[0], bins_subl_loose[1] )

        Ndata_TT_subliso = gg_hist_subliso['leadPass'].Integral( bins_subl_tight[0], bins_subl_tight[1] )
        Ndata_TL_subliso = gg_hist_subliso['leadPass'].Integral( bins_subl_loose[0], bins_subl_loose[1] )
        Ndata_LT_subliso = gg_hist_subliso['leadFail'].Integral( bins_subl_tight[0], bins_subl_tight[1] )
        Ndata_LL_subliso = gg_hist_subliso['leadFail'].Integral( bins_subl_loose[0], bins_subl_loose[1] )

        Ndata_TT_bothiso = gg_hist_bothiso['leadPass'].Integral( bins_subl_tight[0], bins_subl_tight[1] )
        
    else:
        # Find the bins corresponding to the cuts
        # lead photon on X axis, subl on Y axis
        bins_lead_tight = ( gg_hist_leadiso.GetXaxis().FindBin( cuts[lead_eta]['tight'][0] ), gg_hist_leadiso.GetXaxis().FindBin( cuts[lead_eta]['tight'][1] ) )
        bins_lead_loose = ( gg_hist_leadiso.GetXaxis().FindBin( cuts[lead_eta]['loose'][0] ), gg_hist_leadiso.GetXaxis().FindBin( cuts[lead_eta]['loose'][1] ) )
        bins_subl_tight = ( gg_hist_subliso.GetYaxis().FindBin( cuts[subl_eta]['tight'][0] ), gg_hist_subliso.GetYaxis().FindBin( cuts[subl_eta]['tight'][1] ) )
        bins_subl_loose = ( gg_hist_subliso.GetYaxis().FindBin( cuts[subl_eta]['loose'][0] ), gg_hist_subliso.GetYaxis().FindBin( cuts[subl_eta]['loose'][1] ) )
        
        # Integrate to the the data in the four regions
        Ndata_TT_leadiso = gg_hist_leadiso.Integral( bins_lead_tight[0], bins_lead_tight[1], bins_subl_tight[0], bins_subl_tight[1] )
        Ndata_TL_leadiso = gg_hist_leadiso.Integral( bins_lead_tight[0], bins_lead_tight[1], bins_subl_loose[0], bins_subl_loose[1] )
        Ndata_LT_leadiso = gg_hist_leadiso.Integral( bins_lead_loose[0], bins_lead_loose[1], bins_subl_tight[0], bins_subl_tight[1] )
        Ndata_LL_leadiso = gg_hist_leadiso.Integral( bins_lead_loose[0], bins_lead_loose[1], bins_subl_loose[0], bins_subl_loose[1] )
        
        Ndata_TT_subliso = gg_hist_subliso.Integral( bins_lead_tight[0], bins_lead_tight[1], bins_subl_tight[0], bins_subl_tight[1] )
        Ndata_TL_subliso = gg_hist_subliso.Integral( bins_lead_tight[0], bins_lead_tight[1], bins_subl_loose[0], bins_subl_loose[1] )
        Ndata_LT_subliso = gg_hist_subliso.Integral( bins_lead_loose[0], bins_lead_loose[1], bins_subl_tight[0], bins_subl_tight[1] )
        Ndata_LL_subliso = gg_hist_subliso.Integral( bins_lead_loose[0], bins_lead_loose[1], bins_subl_loose[0], bins_subl_loose[1] )
        
        Ndata_TT_bothiso = gg_hist_bothiso.Integral( bins_lead_tight[0], bins_lead_tight[1], bins_subl_tight[0], bins_subl_tight[1] )

    # arragnge the cuts by region
    eff_cuts = {}
    eff_cuts['lead'] = {}
    eff_cuts['subl'] = {}
    eff_cuts['lead']['tight'] = cuts[lead_eta]['tight']
    eff_cuts['lead']['loose'] = cuts[lead_eta]['loose']
    eff_cuts['subl']['tight'] = cuts[subl_eta]['tight']
    eff_cuts['subl']['loose'] = cuts[subl_eta]['loose']

    # get template integrals
    #int_leadiso = get_template_integrals( templates_leadiso, eff_cuts )
    #int_subliso = get_template_integrals( templates_subliso, eff_cuts )
    (stat_int_leadiso, syst_int_leadiso) = get_template_integrals( templates_leadiso, eff_cuts, lead_eta, subl_eta, var=fitvar)
    (stat_int_subliso, syst_int_subliso) = get_template_integrals( templates_subliso, eff_cuts, lead_eta, subl_eta, var=fitvar)

    if stat_int_leadiso['subl']['fake']['tight'] != 0 :
        iso_eff_subl_tight = stat_int_subliso['subl']['fake']['tight'] / stat_int_leadiso['subl']['fake']['tight']
    else :
        iso_eff_subl_tight = ufloat( 0, 0)
    if stat_int_leadiso['subl']['fake']['loose'] != 0 :
        iso_eff_subl_loose = stat_int_subliso['subl']['fake']['loose'] / stat_int_leadiso['subl']['fake']['loose']
    else :
        iso_eff_subl_loose = ufloat( 0, 0)
    if stat_int_subliso['lead']['fake']['tight'] != 0 :
        iso_eff_lead_tight = stat_int_leadiso['lead']['fake']['tight'] / stat_int_subliso['lead']['fake']['tight']
    else :
        iso_eff_lead_tight = ufloat( 0, 0)
    if stat_int_subliso['lead']['fake']['loose'] != 0 :
        iso_eff_lead_loose = stat_int_leadiso['lead']['fake']['loose'] / stat_int_subliso['lead']['fake']['loose']
    else :
        iso_eff_lead_loose = ufloat( 0, 0)

    #-----------------------------------------
    # Use data with loosened iso on the Loose photon
    # Multiply by the efficiency of the loosened iso
    #-----------------------------------------
     
    # lead has loosened iso
    # Correct data in LT region by loosening isolation on the lead photon 
    # and correct by the efficiency of the loosened selection
    Ncorr_LT         = Ndata_LT_subliso * iso_eff_lead_loose
    Ncorr_LL_subliso = Ndata_LL_subliso * iso_eff_lead_loose
    Ncorr_TT_subliso = Ndata_TT_subliso * iso_eff_lead_loose

    # subl has loosened iso
    # correct data in TL region by loosening isolation on the subl photon
    # and correct by the efficiency of the loosened selection
    Ncorr_TL         = Ndata_TL_leadiso * iso_eff_subl_loose
    Ncorr_LL_leadiso = Ndata_LL_leadiso * iso_eff_subl_loose
    Ncorr_TT_leadiso = Ndata_TT_leadiso * iso_eff_subl_loose

    # use the average of the two
    Ncorr_LL = ( Ncorr_LL_leadiso + Ncorr_LL_subliso )/2.
    Ncorr_TT = ( Ncorr_TT_leadiso + Ncorr_TT_subliso )/2.

    print 'NData both iso TT = %d' %Ndata_TT_bothiso
    print 'NData orig leadiso , TT = %d, TL = %d, LT = %d, LL = %d' %( Ndata_TT_leadiso, Ndata_TL_leadiso, Ndata_LT_leadiso, Ndata_LL_leadiso )
    print 'NData orig subliso , TT = %d, TL = %d, LT = %d, LL = %d' %( Ndata_TT_subliso, Ndata_TL_subliso, Ndata_LT_subliso, Ndata_LL_subliso )
    print 'iso_eff_subl_tight = %s, iso_eff_subl_loose = %s, iso_eff_lead_tight = %s, iso_eff_lead_loose= %s' %( iso_eff_subl_tight, iso_eff_subl_loose, iso_eff_lead_tight, iso_eff_subl_loose )
    print 'NData corr, TL = %s, LT = %s, LLlead = %s, LLsubl = %s, LL = %s, TTlead = %s, TTsubl = %s, TT = %s ' %( Ncorr_TL, Ncorr_LT, Ncorr_LL_leadiso, Ncorr_LL_subliso, Ncorr_LL, Ncorr_TT_leadiso, Ncorr_TT_subliso, Ncorr_TT )

    templates_comb = {}
    templates_comb['lead'] = {}
    templates_comb['subl'] = {}
    templates_comb['lead']['real'] = templates_leadiso['lead']['real']
    templates_comb['lead']['fake'] = templates_leadiso['lead']['fake']
    templates_comb['lead']['bakg'] = templates_leadiso['lead']['bakg']
    templates_comb['subl']['real'] = templates_subliso['subl']['real']
    templates_comb['subl']['fake'] = templates_subliso['subl']['fake']
    templates_comb['subl']['bakg'] = templates_leadiso['subl']['bakg']

    # get 2-d efficiencies from 1-d inputs
    eff_results = generate_2d_efficiencies(templates_comb, eff_cuts, lead_eta, subl_eta, var=fitvar)

    if templates_corr is not None :

        eff_ff_2d_stat, eff_ff_2d_syst = generate_2d_corr_efficiencies( templates_corr, eff_cuts, lead_eta, subl_eta, var=fitvar )

        print 'REPLACE FF with correlated templates'
        print 'eff_FF_TT before = %s, after = %s (tight unc)' %( eff_results['stat_tight']['eff_FF_TT'], eff_ff_2d_stat['eff_FF_TT'] )
        print 'eff_FF_TL before = %s, after = %s (tight unc)' %( eff_results['stat_tight']['eff_FF_TL'], eff_ff_2d_stat['eff_FF_TL'] )
        print 'eff_FF_LT before = %s, after = %s (tight unc)' %( eff_results['stat_tight']['eff_FF_LT'], eff_ff_2d_stat['eff_FF_LT'] )
        print 'eff_FF_LL before = %s, after = %s (tight unc)' %( eff_results['stat_tight']['eff_FF_LL'], eff_ff_2d_stat['eff_FF_LL'] )

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

    else:
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
            eff_2d_nouncert[key] = ufloat( val.n, 0 )
        else :
            eff_2d_nouncert[key] =  val

    #Broken
    if ndim == 3 :
        datacorr = {}
        datacorr['TL'] = Ncorr_TL
        datacorr['LT'] = Ncorr_LT
        datacorr['LL'] = Ncorr_LL

        datacorr_nostat = {}
        datacorr_nostat['TL'] = ufloat( Ncorr_TL.n, 0.0 )
        datacorr_nostat['LT'] = ufloat( Ncorr_LT.n, 0.0 )
        datacorr_nostat['LL'] = ufloat( Ncorr_LL.n, 0.0 )

        results_stat_data       = run_fit( datacorr, eff_2d_nouncert )
        results_stat_temp_tight = run_fit( datacorr_nostat, eff_results['stat_tight'] )
        results_stat_temp_loose = run_fit( datacorr_nostat, eff_results['stat_loose'] )
        results_stat_ff         = run_fit( datacorr_nostat, eff_results['stat_ff'] )
        results_syst_bkg        = run_fit( datacorr_nostat, eff_results['syst_bkg'] )
        results_syst_temp       = run_fit( datacorr_nostat, eff_results['syst_temp'] )

        datacorr['TT'] = ufloat(0, 0)

        text_results_stat_data       = collect_results( results_stat_data      , datacorr, eff_2d_nouncert , templates_comb, eff_cuts, ndim)
        text_results_stat_temp_tight = collect_results( results_stat_temp_tight, datacorr, eff_results['stat_tight']     , templates_comb, eff_cuts, ndim)
        text_results_stat_temp_loose = collect_results( results_stat_temp_loose, datacorr, eff_results['stat_loose']     , templates_comb, eff_cuts, ndim)
        text_results_stat_ff         = collect_results( results_stat_ff        , datacorr, eff_results['stat_ff'   ]     , templates_comb, eff_cuts, ndim)
        text_results_syst_bkg        = collect_results( results_syst_bkg       , datacorr, eff_results['syst_bkg'] , templates_comb, eff_cuts, ndim)
        text_results_syst_temp       = collect_results( results_syst_temp      , datacorr, eff_results['syst_temp'], templates_comb, eff_cuts, ndim)

        #Broken
        return text_results_stat_data, text_results_stat_temp_tight, text_results_stat_temp_loose, text_results_stat_ff, text_results_syst_bkg, text_results_syst_temp

    else :

        # keep all data uncertainties for storing
        datacorr = {}
        datacorr['TT'] = ufloat( Ndata_TT_bothiso, math.sqrt(Ndata_TT_bothiso))
        datacorr['TL'] = Ncorr_TL
        datacorr['LT'] = Ncorr_LT
        datacorr['LL'] = Ncorr_LL

        # only consider uncertainty from the SR
        datacorrSR = {}
        datacorrSR['TT'] = ufloat( Ndata_TT_bothiso, math.sqrt(Ndata_TT_bothiso))
        datacorrSR['TL'] = ufloat(Ncorr_TL.n, 0.0 )
        datacorrSR['LT'] = ufloat(Ncorr_LT.n, 0.0 )
        datacorrSR['LL'] = ufloat(Ncorr_LL.n, 0.0 )

        # only consdier uncertainty from SB
        datacorrSB = {}
        datacorrSB['TT'] = ufloat( Ndata_TT_bothiso, 0.0)
        datacorrSB['TL'] = Ncorr_TL
        datacorrSB['LT'] = Ncorr_LT
        datacorrSB['LL'] = Ncorr_LL

        datacorr_nostat = {}
        datacorr_nostat['TT'] = ufloat( Ndata_TT_bothiso, 0.0 )
        datacorr_nostat['TL'] = ufloat( Ncorr_TL.n, 0.0 )
        datacorr_nostat['LT'] = ufloat( Ncorr_LT.n, 0.0 )
        datacorr_nostat['LL'] = ufloat( Ncorr_LL.n, 0.0 )

        results_stat_dataSR     = run_fit( datacorrSR       , eff_2d_nouncert)
        results_stat_dataSB     = run_fit( datacorrSB       , eff_2d_nouncert)
        results_stat_temp_tight = run_fit( datacorr_nostat, eff_results['stat_tight'] )
        results_stat_temp_loose = run_fit( datacorr_nostat, eff_results['stat_loose'] )
        results_stat_ff         = run_fit( datacorr_nostat, eff_results['stat_ff'] )
        results_syst_bkg        = run_fit( datacorr_nostat, eff_results['syst_bkg'] )
        results_syst_temp       = run_fit( datacorr_nostat, eff_results['syst_temp'] )

        text_results_stat_dataSR     = collect_results( results_stat_dataSR      , datacorr, eff_2d_nouncert  , templates_comb, eff_cuts, ndim)
        text_results_stat_dataSB     = collect_results( results_stat_dataSB      , datacorr, eff_2d_nouncert  , templates_comb, eff_cuts, ndim)
        text_results_stat_temp_tight = collect_results( results_stat_temp_tight, datacorr, eff_results['stat_tight'], templates_comb, eff_cuts, ndim)
        text_results_stat_temp_loose = collect_results( results_stat_temp_loose, datacorr, eff_results['stat_loose'], templates_comb, eff_cuts, ndim)
        text_results_stat_ff         = collect_results( results_stat_ff        , datacorr, eff_results['stat_ff'   ]     , templates_comb, eff_cuts, ndim)
        text_results_syst_bkg        = collect_results( results_syst_bkg       , datacorr, eff_results['syst_bkg']  , templates_comb, eff_cuts, ndim)
        text_results_syst_temp       = collect_results( results_syst_temp      , datacorr, eff_results['syst_temp'] , templates_comb, eff_cuts, ndim)

        # print 'text_results_syst_bkg'
        # print text_results_syst_bkg

        print '*******For SignalRegion Data Results***********************'
        print 'Npred_RF_TT = ', text_results_stat_dataSR['Npred_RF_TT']
        print 'Npred_FR_TT = ', text_results_stat_dataSR['Npred_FR_TT']
        print 'Npred_FF_TT = ', text_results_stat_dataSR['Npred_FF_TT']
        print 'Sum = ', (text_results_stat_dataSR['Npred_RF_TT']+text_results_stat_dataSR['Npred_FR_TT']+text_results_stat_dataSR['Npred_FF_TT'])

        print '*******For Sideband Data Results***********************'
        print 'Npred_RF_TT = ', text_results_stat_dataSB['Npred_RF_TT']
        print 'Npred_FR_TT = ', text_results_stat_dataSB['Npred_FR_TT']
        print 'Npred_FF_TT = ', text_results_stat_dataSB['Npred_FF_TT']
        print 'Sum = ', (text_results_stat_dataSB['Npred_RF_TT']+text_results_stat_dataSB['Npred_FR_TT']+text_results_stat_dataSB['Npred_FF_TT'])

        print '*******For TightTemp Data Results***********************'
        print 'Npred_RF_TT = ', text_results_stat_temp_tight['Npred_RF_TT']
        print 'Npred_FR_TT = ', text_results_stat_temp_tight['Npred_FR_TT']
        print 'Npred_FF_TT = ', text_results_stat_temp_tight['Npred_FF_TT']
        print 'Sum = ', (text_results_stat_temp_tight['Npred_RF_TT']+text_results_stat_temp_tight['Npred_FR_TT']+text_results_stat_temp_tight['Npred_FF_TT'])

        print '*******For LooseTemp Data Results***********************'
        print 'Npred_RF_TT = ', text_results_stat_temp_loose['Npred_RF_TT']
        print 'Npred_FR_TT = ', text_results_stat_temp_loose['Npred_FR_TT']
        print 'Npred_FF_TT = ', text_results_stat_temp_loose['Npred_FF_TT']
        print 'Sum = ', (text_results_stat_temp_loose['Npred_RF_TT']+text_results_stat_temp_loose['Npred_FR_TT']+text_results_stat_temp_loose['Npred_FF_TT'])

        print '*******For FF Temp Data Results*************************'
        print 'Npred_RF_TT = ', text_results_stat_ff['Npred_RF_TT']
        print 'Npred_FR_TT = ', text_results_stat_ff['Npred_FR_TT']
        print 'Npred_FF_TT = ', text_results_stat_ff['Npred_FF_TT']
        print 'Sum = ', (text_results_stat_ff['Npred_RF_TT']+text_results_stat_ff['Npred_FR_TT']+text_results_stat_ff['Npred_FF_TT'])

        print '*******For Bkg Temp Data Results*************************'
        print 'Npred_RF_TT = ', text_results_syst_bkg['Npred_RF_TT']
        print 'Npred_FR_TT = ', text_results_syst_bkg['Npred_FR_TT']
        print 'Npred_FF_TT = ', text_results_syst_bkg['Npred_FF_TT']
        print 'Sum = ', (text_results_syst_bkg['Npred_RF_TT']+text_results_syst_bkg['Npred_FR_TT']+text_results_syst_bkg['Npred_FF_TT'])

        print '*******For SysTemp Data Results*************************'
        print 'Npred_RF_TT = ', text_results_syst_temp['Npred_RF_TT']
        print 'Npred_FR_TT = ', text_results_syst_temp['Npred_FR_TT']
        print 'Npred_FF_TT = ', text_results_syst_temp['Npred_FF_TT']
        print 'Sum = ', (text_results_syst_temp['Npred_RF_TT']+text_results_syst_temp['Npred_FR_TT']+text_results_syst_temp['Npred_FF_TT'])
        print '********************************************************'

        return text_results_stat_dataSR,text_results_stat_dataSB, text_results_stat_temp_tight, text_results_stat_temp_loose, text_results_stat_ff, text_results_syst_bkg, text_results_syst_temp

def generate_1d_efficiencies( templates, eff_cuts, lead_eta, subl_eta, var=None ) :

    (int_stat, int_syst) = get_template_integrals( templates, eff_cuts, lead_eta, subl_eta, var=var)

    (eff_1d_stat_tight, eff_1d_stat_loose, eff_1d_syst_bkg, eff_1d_syst_temp) = get_1d_loose_efficiencies( int_stat, int_syst, lead_eta, subl_eta, var=var)

    #    return eff_1d_stat, eff_1d_syst
    return eff_1d_stat_tight, eff_1d_stat_loose, eff_1d_syst_bkg, eff_1d_syst_temp

def generate_2d_efficiencies( templates, eff_cuts, lead_eta, subl_eta, var=None) :

    (int_stat, int_syst) = get_template_integrals( templates, eff_cuts, lead_eta, subl_eta, var=var)

    # get efficiencies with three sources of uncertainty
    # 1. statistical uncertainty on templates
    # 2. systematics uncertainty on background subtraction
    # 3. systematic uncertainty on template shapes
    # The integrals from get_template_integrals
    # give 1 and 2

    eff_results = {}

    (eff_1d_stat_tight, eff_1d_stat_loose, eff_1d_syst_bkg, eff_1d_syst_temp) = get_1d_loose_efficiencies( int_stat, int_syst, lead_eta, subl_eta, var=var)

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

    # eff_stat
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
    
    # eff_syst_bkg
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

    # eff_syst_temp
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
    eff_syst_temp['eff_R_L_lead'] = ufloat( eff_syst_temp['eff_R_L_lead'].n , math.fabs(eff_syst_temp['eff_R_L_lead'].n)*_real_syst[lead_eta], 'Template_lead_real_loose')
    eff_syst_temp['eff_F_L_lead'] = ufloat( eff_syst_temp['eff_F_L_lead'].n , math.fabs(eff_syst_temp['eff_F_L_lead'].n)*_fake_syst[lead_eta], 'Template_lead_fake_loose' )
    eff_syst_temp['eff_R_L_subl'] = ufloat( eff_syst_temp['eff_R_L_subl'].n , math.fabs(eff_syst_temp['eff_R_L_subl'].n)*_real_syst[subl_eta], 'Template_subl_real_loose' )
    eff_syst_temp['eff_F_L_subl'] = ufloat( eff_syst_temp['eff_F_L_subl'].n , math.fabs(eff_syst_temp['eff_F_L_subl'].n)*_fake_syst[subl_eta], 'Template_subl_fake_loose' )

    return eff_stat_tight, eff_stat_loose, eff_syst_int, eff_syst_temp

def run_fit( data, efficiencies ) :

    # make the matrix
    matrix = generate_eff_matrix( efficiencies, ndim=len(data) )

    #do the fit!  Invert the matrix and multiply the by counts vectors
    if len( data ) == 3 :
        results = solve_matrix_eq( matrix, [data['TL'], data['LT'], data['LL']] )
    elif len(data) == 4 :
        results = solve_matrix_eq( matrix, [data['TT'],data['TL'], data['LT'], data['LL']] )

    return results

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


    if templates['lead']['bakg'] is not None :
        bkg_int_tight = get_integral_and_error( templates['lead']['bakg'], bins_lead_fake_tight, 'Background_lead_fake_tight' ) 
        bkg_int_loose = get_integral_and_error( templates['lead']['bakg'], bins_lead_fake_loose, 'Background_lead_fake_loose' )
        
        # print 'bkg_int_tight = ', bkg_int_tight
        # print 'bkg_int_loose = ', bkg_int_loose
        
        syst_bkg_int_tight = ufloat(0.0, math.fabs(bkg_int_tight.n)*0.20, 'Background_lead_fake_tight')
        syst_bkg_int_loose = ufloat(0.0, math.fabs(bkg_int_loose.n)*0.20, 'Background_lead_fake_loose')

        int_syst['lead']['fake']['tight'] = int_syst['lead']['fake']['tight'] + syst_bkg_int_tight
        int_syst['lead']['fake']['loose'] = int_syst['lead']['fake']['loose'] + syst_bkg_int_loose
    
    if templates['subl']['bakg'] is not None :
        bkg_int_tight = get_integral_and_error( templates['subl']['bakg'], bins_subl_fake_tight, 'Background_subl_fake_tight' ) 
        bkg_int_loose = get_integral_and_error( templates['subl']['bakg'], bins_subl_fake_loose, 'Background_subl_fake_loose' )
        
        # print 'bkg_int_tight = ', bkg_int_tight
        # print 'bkg_int_loose = ', bkg_int_loose
        
        syst_bkg_int_tight = ufloat(0.0, math.fabs(bkg_int_tight.n)*0.20, 'Background_subl_fake_tight')
        syst_bkg_int_loose = ufloat(0.0, math.fabs(bkg_int_loose.n)*0.20, 'Background_subl_fake_loose')
        
        int_syst['subl']['fake']['tight'] = int_syst['subl']['fake']['tight'] + syst_bkg_int_tight
        int_syst['subl']['fake']['loose'] = int_syst['subl']['fake']['loose'] + syst_bkg_int_loose

    return int_stat, int_syst

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


def generate_2d_corr_efficiencies( templates, cuts, lead_eta, subl_eta, var=None ) :
    # integral each region of the 2-d templates

    # lead is on y axis, subl on x

    new_cuts = get_default_cuts('ChIso')
#    print 'new_cuts = ', new_cuts
    new_eff_cuts = {}
    new_eff_cuts['lead'] = {}
    new_eff_cuts['subl'] = {}
    new_eff_cuts['lead']['tight'] = new_cuts[lead_eta]['tight']
    new_eff_cuts['lead']['loose'] = new_cuts[lead_eta]['loose']
    new_eff_cuts['subl']['tight'] = new_cuts[subl_eta]['tight']
    new_eff_cuts['subl']['loose'] = new_cuts[subl_eta]['loose']

    bin_subl_tight = (templates['leadPass'].GetXaxis().FindBin(new_eff_cuts['subl']['tight'][0]), templates['leadPass'].GetXaxis().FindBin(new_eff_cuts['subl']['tight'][1]) )
    bin_subl_loose = (templates['leadPass'].GetXaxis().FindBin(new_eff_cuts['subl']['loose'][0]), templates['leadPass'].GetXaxis().FindBin(new_eff_cuts['subl']['loose'][1]) )
    # bin_subl_tight = ( templates['leadPass'].GetXaxis().FindBin(cuts['subl']['tight'][0]), templates['leadPass'].GetXaxis().FindBin(cuts['subl']['tight'][1]) )
    # bin_subl_loose = ( templates['leadPass'].GetXaxis().FindBin(cuts['subl']['loose'][0]), templates['leadPass'].GetXaxis().FindBin(cuts['subl']['loose'][1]) )

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
    int_ll = ufloat( int_lead_loose_subl_loose, err_lead_loose_subl_loose )

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

    # Matrices of numbers with uncertainties are best created in one of two ways. The first way is similar to using uarray():
    matrix = unumpy.umatrix( mn, ms )

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


####################################
_fitVar = ['Sieie', 'ChIso']
_iso = ['533', '855', '1077', '1299', '151111', '201616']
iso = _iso[0]
fitVar = _fitVar[0]

## do_nominal_fit ##
path          = '/uscmst1b_scratch/lpc1/old_scratch/lpceg/yurii/EnSc/sachiko'
dir_tmp       = '/Zgg/NEW_WORK_09010/Zgg/Muon/TempHist/TEMP_ROOT_0111/'
dir_tmp_loose       = '/Zgg/NEW_WORK_09010/Zgg/Muon/TempHist/TEMP_ROOT_1023/'

# Input_Zg        = path + dir_tmp + 'ROOT_Sig/muon_Sig_'+fitVar+'_nomi.root'
# Input_data_Zjet = path + dir_tmp + 'ROOT_Bkg/muon_subtrBkg_'+fitVar+'_nomi.root'
# Input_mc_Zjet   = path + dir_tmp + 'ROOT_BkgMC/muon_Bkg_'+fitVar+'_nomi.root'
# Input_FF_corr   = path + dir_tmp + 'ROOT_CorrFF/Temp_Corr_FF_nomi.root'
# #Input_Zgg       = path + dir_tmp + 'ROOT_Zgg/muon_Zgg_'+fitVar+'_'+iso+'.root'
# Input_Zgg       = path + dir_tmp + 'ROOT_Zgg/electron_Zgg_'+fitVar+'_'+iso+'.root'
# Input_Zg_loose        = path + dir_tmp + 'ROOT_Sig/muon_Sig_'+fitVar+'_'+iso+'.root'
# Input_data_Zjet_loose = path + dir_tmp + 'ROOT_Bkg/muon_subtrBkg_'+fitVar+'_'+iso+'.root'
# Input_mc_Zjet_loose   = path + dir_tmp + 'ROOT_BkgMC/muon_Bkg_'+fitVar+'_'+iso+'.root'
# pt_reg = ['1', '2', '3']

Input_Zg              = path + dir_tmp + 'ROOT_allPt/muon_allpT_Sig_'+fitVar+'_nomi.root'
Input_data_Zjet       = path + dir_tmp + 'ROOT_allPt/muon_allpT_subtrBkg_'+fitVar+'_nomi.root'
Input_mc_Zjet         = path + dir_tmp + 'ROOT_allPt/muon_allpT_Bkg_'+fitVar+'_nomi.root'
Input_FF_corr         = path + dir_tmp + 'ROOT_allPt/Temp_allPt_Corr_FF_nomi.root'
#Input_Zgg             = path + dir_tmp + 'ROOT_allPt/muon_allpT_Zgg_'+fitVar+'_'+iso+'.root'
Input_Zgg             = path + dir_tmp + 'ROOT_allPt/electron_allpT_Zgg_'+fitVar+'_'+iso+'.root'
Input_Zg_loose        = path + dir_tmp + 'ROOT_allPt/muon_allpT_Sig_'+fitVar+'_'+iso+'.root'
Input_data_Zjet_loose = path + dir_tmp + 'ROOT_allPt/muon_allpT_subtrBkg_'+fitVar+'_'+iso+'.root'
Input_mc_Zjet_loose   = path + dir_tmp + 'ROOT_allPt/muon_allpT_Bkg_'+fitVar+'_'+iso+'.root'
pt_reg = ['1']

for bin in pt_reg :

    # templates_reg = {}
    # templates_reg['EB'] = {}
    # templates_reg['EE'] = {}
    # templates_reg['EB']['real'] = getWeightedHist(Input_Zg, 'h_'+fitVar+'_EB_'+bin)
    # templates_reg['EE']['real'] = getWeightedHist(Input_Zg, 'h_'+fitVar+'_EE_'+bin)
    # templates_reg['EB']['fake'] = getWeightedHist(Input_data_Zjet, 'h_'+fitVar+'_EB_'+bin)
    # templates_reg['EE']['fake'] = getWeightedHist(Input_data_Zjet, 'h_'+fitVar+'_EE_'+bin)
    # templates_reg['EB']['bakg'] = getWeightedHist(Input_mc_Zjet, 'h_'+fitVar+'_EB_'+bin)
    # templates_reg['EE']['bakg'] = getWeightedHist(Input_mc_Zjet, 'h_'+fitVar+'_EE_'+bin)
 
    templates_iso_reg = {}
    templates_iso_reg['EB'] = {}
    templates_iso_reg['EE'] = {}
    templates_iso_reg['EB']['real'] = getWeightedHist(Input_Zg, 'h_'+fitVar+'_EB_'+bin)
    templates_iso_reg['EE']['real'] = getWeightedHist(Input_Zg, 'h_'+fitVar+'_EE_'+bin)
    templates_iso_reg['EB']['fake'] = getWeightedHist(Input_data_Zjet, 'h_'+fitVar+'_EB_'+bin)
    templates_iso_reg['EE']['fake'] = getWeightedHist(Input_data_Zjet, 'h_'+fitVar+'_EE_'+bin)
    templates_iso_reg['EB']['bakg'] = getWeightedHist(Input_mc_Zjet, 'h_'+fitVar+'_EB_'+bin)
    templates_iso_reg['EE']['bakg'] = getWeightedHist(Input_mc_Zjet, 'h_'+fitVar+'_EE_'+bin)

    templates_noiso_reg = {}
    templates_noiso_reg['EB'] = {}
    templates_noiso_reg['EE'] = {}
    templates_noiso_reg['EB']['real'] = getWeightedHist(Input_Zg_loose, 'h_'+fitVar+'_EB_'+bin)
    templates_noiso_reg['EE']['real'] = getWeightedHist(Input_Zg_loose, 'h_'+fitVar+'_EE_'+bin)
    templates_noiso_reg['EB']['fake'] = getWeightedHist(Input_data_Zjet_loose, 'h_'+fitVar+'_EB_'+bin)
    templates_noiso_reg['EE']['fake'] = getWeightedHist(Input_data_Zjet_loose, 'h_'+fitVar+'_EE_'+bin)
    templates_noiso_reg['EB']['bakg'] = getWeightedHist(Input_mc_Zjet_loose, 'h_'+fitVar+'_EB_'+bin)
    templates_noiso_reg['EE']['bakg'] = getWeightedHist(Input_mc_Zjet_loose, 'h_'+fitVar+'_EE_'+bin)

    templatrs_corr = None
    ffcorr = get_default_cuts(fitVar)

    templates_corr = {'leadPass' : {}, 'leadFail' : {}}

    templates_corr['leadFail'][('EB','EB')] = getWeightedHist(Input_FF_corr, 'h_ChIso_LeadFail_EB_EB_'+bin)
    templates_corr['leadPass'][('EB','EB')] = getWeightedHist(Input_FF_corr, 'h_ChIso_LeadPass_EB_EB_'+bin)
    templates_corr['leadFail'][('EB','EE')] = getWeightedHist(Input_FF_corr, 'h_ChIso_LeadFail_EB_EE_'+bin)
    templates_corr['leadPass'][('EB','EE')] = getWeightedHist(Input_FF_corr, 'h_ChIso_LeadPass_EB_EE_'+bin)
    templates_corr['leadFail'][('EE','EB')] = getWeightedHist(Input_FF_corr, 'h_ChIso_LeadFail_EE_EB_'+bin)
    templates_corr['leadPass'][('EE','EB')] = getWeightedHist(Input_FF_corr, 'h_ChIso_LeadPass_EE_EB_'+bin)

    # templates_corr['leadFail'][('EB','EB')] = getWeightedHist(Input_FF_corr, 'h_'+fitVar+'_LeadFail_EB_EB_'+bin)
    # templates_corr['leadPass'][('EB','EB')] = getWeightedHist(Input_FF_corr, 'h_'+fitVar+'_LeadPass_EB_EB_'+bin)
    # templates_corr['leadFail'][('EB','EE')] = getWeightedHist(Input_FF_corr, 'h_'+fitVar+'_LeadFail_EB_EE_'+bin)
    # templates_corr['leadPass'][('EB','EE')] = getWeightedHist(Input_FF_corr, 'h_'+fitVar+'_LeadPass_EB_EE_'+bin)
    # templates_corr['leadFail'][('EE','EB')] = getWeightedHist(Input_FF_corr, 'h_'+fitVar+'_LeadFail_EE_EB_'+bin)
    # templates_corr['leadPass'][('EE','EB')] = getWeightedHist(Input_FF_corr, 'h_'+fitVar+'_LeadPass_EE_EB_'+bin)

    gg_hist_leadiso_name = {'leadPass' : {}, 'leadFail' : {}}
    gg_hist_leadiso_name['leadPass'][('EB','EB')] = getWeightedHist(Input_Zgg, 'h_'+fitVar+'_LeadIso_LeadPass_EB_EB_'+bin)
    gg_hist_leadiso_name['leadFail'][('EB','EB')] = getWeightedHist(Input_Zgg, 'h_'+fitVar+'_LeadIso_LeadFail_EB_EB_'+bin)
    gg_hist_leadiso_name['leadPass'][('EB','EE')] = getWeightedHist(Input_Zgg, 'h_'+fitVar+'_LeadIso_LeadPass_EB_EE_'+bin)
    gg_hist_leadiso_name['leadFail'][('EB','EE')] = getWeightedHist(Input_Zgg, 'h_'+fitVar+'_LeadIso_LeadFail_EB_EE_'+bin)
    gg_hist_leadiso_name['leadPass'][('EE','EB')] = getWeightedHist(Input_Zgg, 'h_'+fitVar+'_LeadIso_LeadPass_EE_EB_'+bin)
    gg_hist_leadiso_name['leadFail'][('EE','EB')] = getWeightedHist(Input_Zgg, 'h_'+fitVar+'_LeadIso_LeadFail_EE_EB_'+bin)

    gg_hist_subliso_name = {'leadPass' : {}, 'leadFail' : {}}
    gg_hist_subliso_name['leadPass'][('EB','EB')] = getWeightedHist(Input_Zgg, 'h_'+fitVar+'_SublIso_LeadPass_EB_EB_'+bin)
    gg_hist_subliso_name['leadFail'][('EB','EB')] = getWeightedHist(Input_Zgg, 'h_'+fitVar+'_SublIso_LeadFail_EB_EB_'+bin)
    gg_hist_subliso_name['leadPass'][('EB','EE')] = getWeightedHist(Input_Zgg, 'h_'+fitVar+'_SublIso_LeadPass_EB_EE_'+bin)
    gg_hist_subliso_name['leadFail'][('EB','EE')] = getWeightedHist(Input_Zgg, 'h_'+fitVar+'_SublIso_LeadFail_EB_EE_'+bin)
    gg_hist_subliso_name['leadPass'][('EE','EB')] = getWeightedHist(Input_Zgg, 'h_'+fitVar+'_SublIso_LeadPass_EE_EB_'+bin)
    gg_hist_subliso_name['leadFail'][('EE','EB')] = getWeightedHist(Input_Zgg, 'h_'+fitVar+'_SublIso_LeadFail_EE_EB_'+bin)

    gg_hist_bothiso_name = {'leadPass' : {}, 'leadFail' : {}}
    gg_hist_bothiso_name['leadPass'][('EB','EB')] = getWeightedHist(Input_Zgg, 'h_'+fitVar+'_BothIso_LeadPass_EB_EB_'+bin)
    gg_hist_bothiso_name['leadFail'][('EB','EB')] = getWeightedHist(Input_Zgg, 'h_'+fitVar+'_BothIso_LeadFail_EB_EB_'+bin)
    gg_hist_bothiso_name['leadPass'][('EB','EE')] = getWeightedHist(Input_Zgg, 'h_'+fitVar+'_BothIso_LeadPass_EB_EE_'+bin)
    gg_hist_bothiso_name['leadFail'][('EB','EE')] = getWeightedHist(Input_Zgg, 'h_'+fitVar+'_BothIso_LeadFail_EB_EE_'+bin)
    gg_hist_bothiso_name['leadPass'][('EE','EB')] = getWeightedHist(Input_Zgg, 'h_'+fitVar+'_BothIso_LeadPass_EE_EB_'+bin)
    gg_hist_bothiso_name['leadFail'][('EE','EB')] = getWeightedHist(Input_Zgg, 'h_'+fitVar+'_BothIso_LeadFail_EE_EB_'+bin)

    regions = [ ('EB', 'EB'), ('EB', 'EE'), ('EE', 'EB')]
    for reg in regions :
    
        print '######################################################################'
        print 'pt bin = %s,   lead eta = %s,   subl eta = %s' % (bin, reg[0], reg[1])
        print '######################################################################'

        # templates_leadiso = {'lead' : { 'real' : {}, 'fake' : {} }, 'subl' : {'real' : {}, 'fake' : {} } }
        # templates_subliso = {'lead' : { 'real' : {}, 'fake' : {} }, 'subl' : {'real' : {}, 'fake' : {} } }
        
        # templates_leadiso['lead']['real'] = templates_reg[reg[0]]['real']
        # templates_leadiso['subl']['real'] = templates_reg[reg[1]]['real']
        # templates_leadiso['lead']['fake'] = templates_reg[reg[0]]['fake']
        # templates_leadiso['subl']['fake'] = templates_reg[reg[1]]['fake']
        # templates_leadiso['lead']['bakg'] = templates_reg[reg[0]]['bakg']
        # templates_leadiso['subl']['bakg'] = templates_reg[reg[1]]['bakg']

        # templates_subliso['lead']['real'] = templates_reg[reg[0]]['real']
        # templates_subliso['subl']['real'] = templates_reg[reg[1]]['real']
        # templates_subliso['lead']['fake'] = templates_reg[reg[0]]['fake']
        # templates_subliso['subl']['fake'] = templates_reg[reg[1]]['fake']
        # templates_subliso['lead']['bakg'] = templates_reg[reg[0]]['bakg']
        # templates_subliso['subl']['bakg'] = templates_reg[reg[1]]['bakg']

        templates_leadiso = {}
        templates_leadiso['lead'] = {}
        templates_leadiso['subl'] = {}
        templates_leadiso['lead']['real'] = templates_iso_reg[reg[0]]['real']
        templates_leadiso['subl']['real'] = templates_noiso_reg[reg[1]]['real']
        templates_leadiso['lead']['fake'] = templates_iso_reg[reg[0]]['fake']
        templates_leadiso['subl']['fake'] = templates_noiso_reg[reg[1]]['fake']
        templates_leadiso['lead']['bakg'] = templates_iso_reg[reg[0]]['bakg']
        templates_leadiso['subl']['bakg'] = templates_noiso_reg[reg[1]]['bakg']

        templates_subliso = {}
        templates_subliso['lead'] = {}
        templates_subliso['subl'] = {}
        templates_subliso['lead']['real'] = templates_noiso_reg[reg[0]]['real']
        templates_subliso['subl']['real'] = templates_iso_reg[reg[1]]['real']
        templates_subliso['lead']['fake'] = templates_noiso_reg[reg[0]]['fake']
        templates_subliso['subl']['fake'] = templates_iso_reg[reg[1]]['fake']
        templates_subliso['lead']['bakg'] = templates_noiso_reg[reg[0]]['bakg']
        templates_subliso['subl']['bakg'] = templates_iso_reg[reg[1]]['bakg']

        templates_nom = {}
        templates_nom['lead'] = {}
        templates_nom['subl'] = {}
        templates_nom['subl']['real'] = templates_subliso['subl']['real']
        templates_nom['subl']['fake'] = templates_subliso['subl']['fake']
        templates_nom['subl']['bakg'] = templates_subliso['subl']['bakg']
        templates_nom['lead']['real'] = templates_leadiso['lead']['real']
        templates_nom['lead']['fake'] = templates_leadiso['lead']['fake']
        templates_nom['lead']['bakg'] = templates_leadiso['lead']['bakg']

        # templates for colleration of fake-fake
        templates_corr_inc = None
        if templates_corr is not None :
            templates_corr_inc = { }
            templates_corr_inc['leadPass'] = templates_corr['leadPass'][(reg[0], reg[1])]
            templates_corr_inc['leadFail'] = templates_corr['leadFail'][(reg[0], reg[1])]
        
        gg_hist_leadiso = None
        gg_hist_subliso = None
        gg_hist_bothiso = None
        if gg_hist_leadiso_name is not None :
            gg_hist_leadiso = { }
            gg_hist_leadiso['leadPass'] = gg_hist_leadiso_name['leadPass'][(reg[0], reg[1])]
            gg_hist_leadiso['leadFail'] = gg_hist_leadiso_name['leadFail'][(reg[0], reg[1])]

        if gg_hist_subliso_name is not None :
            gg_hist_subliso = { }
            gg_hist_subliso['leadPass'] = gg_hist_subliso_name['leadPass'][(reg[0], reg[1])]
            gg_hist_subliso['leadFail'] = gg_hist_subliso_name['leadFail'][(reg[0], reg[1])]

        if gg_hist_bothiso_name is not None :
            gg_hist_bothiso = { }
            gg_hist_bothiso['leadPass'] = gg_hist_bothiso_name['leadPass'][(reg[0], reg[1])]
            gg_hist_bothiso['leadFail'] = gg_hist_bothiso_name['leadFail'][(reg[0], reg[1])]
        
        outputDir = '/uscmst1b_scratch/lpc1/old_scratch/lpceg/yurii/EnSc/sachiko/Zgg/NEW_WORK_09010/Zgg/Muon/Fitting/FIT_1023/Loose0114/' + fitVar + '_' + reg[0] + reg[1] + '_' + bin


        (text_results_stat_dataSR,text_results_stat_dataSB, text_results_stat_temp_tight, text_results_stat_temp_loose, text_results_stat_ff, text_results_syst_bkg, text_results_syst_temp) = run_corrected_diphoton_fit(templates_leadiso, templates_subliso, gg_hist_leadiso, gg_hist_subliso, gg_hist_bothiso, reg[0], reg[1], templates_corr_inc, fitVar, 4, outputDir)
