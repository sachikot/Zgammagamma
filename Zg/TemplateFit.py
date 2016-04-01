from ROOT import TH1, TH1F, TH2, TH2D, TFile, TCanvas, THStack, TLegend, TMath
from ROOT import gROOT, gStyle
from ROOT import *
import array, glob, ROOT

ROOT.gROOT.SetBatch()
ROOT.gStyle.SetOptStat(0)

def getTemp1DHist(filename, histname):
    file = TFile.Open(filename)    
    temp1DHist = file.Get(histname)
    temp1DHist.SetDirectory(0)
    temp1DHist.SetFillColor(0)

    return temp1DHist

def getTemp2DHist(filename, histname, xMin, xMax, projAxis):
    file = TFile.Open(filename)
    #    print "histoname=%s" % histname
    #    file.Print()
    Hist2D = file.Get(histname)
    #    print type (Hist2D)
    Hist2D.SetDirectory(0)

    if projAxis == 'Y':
        firstBin = Hist2D.GetYaxis().FindBin(xMin)
        endBin = Hist2D.GetYaxis().FindBin(xMax)
    elif projAxis == 'X':
        firstBin = Hist2D.GetXaxis().FindBin(xMin)
        endBin = Hist2D.GetXaxis().FindBin(xMax)
    else:
        print 'Error!!! I need to put X or Y on projAxis!!'
        return

    # project a 2D histogram into a 1D histogram
    if projAxis == 'Y':
        temp2DHist = Hist2D.ProjectionX('projX_'+histname+'_'+str(xMin)+'_'+str(xMax), firstBin, endBin) # str(number): transfer the number to string
    else:
        temp2DHist = Hist2D.ProjectionY('projY_'+histname+'_'+str(xMin)+'_'+str(xMax), firstBin, endBin)

    temp2DHist.SetDirectory(0)
    temp2DHist.SetFillColor(0)
    
    return temp2DHist

def getWeightedHist(filename, histName, leftCut=None, rightCut=None):

    if leftCut==None and rightCut==None:
        # template is 1D histogram
        sumHist = getTemp1DHist(filename, histName)
        return sumHist
    else:
        # template is projected 2D histogram
        sumHistp = getTemp2DHist(filename, histName, leftCut, rightCut, 'X')
        return sumHistp
    
def makeFit(varname, varMin, varMax, signalHist, bkgHist, dataHist, plotName):
    # PhoSCRChHadIso will be input into varname 
    #    print 'variname = ', varname
    # RooFit variables
    argVar = RooRealVar(varname, varname, varMin, varMax)
    argList = RooArgList()
    argList.add(argVar)
    argSet = RooArgSet()
    argSet.add(argVar)

    ############ create PDF files ################
    # signal
    signalDataHist = RooDataHist('signalDataHist', 'signal RooDataHist', argList, signalHist)
    signalPdf = RooHistPdf('signalPdf', varname+' of Signal', argSet, signalDataHist)
    print 'signalPdf = ', signalPdf
    # background
    bkgDataHist = RooDataHist('BackgroundDataHist', 'BackgroundDataHist', argList, bkgHist)
    backgroundPdf = RooHistPdf('backgroundPdf', varname+' of background', argSet, bkgDataHist)
    # data
    dataDataHist = RooDataHist('data '+varname, varname+' in Data', argList, dataHist)

    # signal fraction parameter
    fractionSignal = RooRealVar('signal fraction', 'signal fraction', 0.5, 0.0, 1.0)

    # total pdf of signal and background
    # model(x) = fsig*sig(x) + (1-fsig)*bkg(x)
    sumPDF = RooAddPdf('sumPdf', 'sig and bkg pdf', signalPdf, backgroundPdf, fractionSignal)

    ########## Now do the fit and plot the result #########
    # fit
    sumPDF.fitTo(dataDataHist, RooFit.SumW2Error(kFALSE), RooFit.PrintLevel(-1))

    # plot PDF and toy data overlaid
    if plotName!='':
        c1 = TCanvas('c1', 'c1', 700, 600)
        plotter = RooPlot(argVar, varMin, varMax, 20) # nBins = 20 is dummy
        dataDataHist.plotOn(plotter, RooFit.Name('Data'))
        sumPDF.plotOn(plotter, RooFit.Name('sum'), RooFit.LineColor(kOrange+10))
        sumPDF.plotOn(plotter, RooFit.Components('signalPdf'), RooFit.Name('signal'), RooFit.LineColor(kGreen-4))
        sumPDF.plotOn(plotter, RooFit.Components('backgroundPdf'), RooFit.Name('background'), RooFit.LineColor(kAzure))
        sumPDF.paramOn(plotter)

        plotter.Draw()
        c1.SaveAs(plotName)

    print 'fit returned value ',fractionSignal.getVal(), ' +- ',fractionSignal.getError()
    return (fractionSignal.getVal(),fractionSignal.getError())

def optimizeBinBoundaries(templ, minAccumulate, firstBinValue = -9999, endBinValue = 9999):
    nbins = templ.GetNbinsX()
    accumulate = 0
    binlist = []

    # ignore empty bins in the beginning
    for firstIndex in xrange(1, nbins+1):
        if templ.GetBinLowEdge(firstIndex)>=firstBinValue:
            binlist.append(templ.GetBinLowEdge(firstIndex))
            break

    for Index in xrange(firstIndex, nbins+1):
        if templ.GetBinLowEdge(Index)>=endBinValue:
            binlist.append(endBinValue)
            break
        if accumulate > minAccumulate:
            #            print "accumulate=",accumulate,"minAccumulate",minAccumulate
            accumulate=0
            binlist.append(templ.GetBinLowEdge(Index))
            #        print Index,"th bincontent=",templ.GetBinContent(Index)
        accumulate+=templ.GetBinContent(Index)

    binlist.append(templ.GetBinLowEdge(nbins+1)) # bin end

    newarray = array.array('d') 
    newarray.fromlist(binlist)
    return newarray

def makeUniformHist(hist):
    nbins = hist.GetNbinsX()
    print 'nbins = ', nbins

    #    print '####################################'
    #    print "histoname for makeUniformHist=%s" % hist
    #    print '####################################'    
    
    outHist = TH1F(hist.GetName()+'_uni',hist.GetName()+'_uni',nbins,1,nbins+1)
    for index in xrange(1, nbins+1):
        outHist.SetBinContent(index, hist.GetBinContent(index))
        outHist.SetBinError(index, hist.GetBinError(index))
    return outHist

#################### MAIN CODE: Do The Fitting #######################
InputFile = 'Histogram_SCRIso_fit.root'
Sig_temp = "pho1_RandomCone"
histo2d_temp = "pho1_Sihih_ChIso"
FitVariable = 'SCRChIsolation'

# sideband region for bkg template
sb_left = 0.011
sb_right = 0.014

# selection cut of sigmaIetaIeta
SihihSelCut = 0.011

# Fit range
LowFitRange = -1.0
HighFitRange = 20

# signal template
sig_templ = getWeightedHist(InputFile, Sig_temp)
# background template
bkg_templ = getWeightedHist(InputFile, histo2d_temp, sb_left, sb_right)
# pseudo data
pseudo_data_templ = getWeightedHist(InputFile, histo2d_temp, 0.0, SihihSelCut)

print 'Signal Integral = ',sig_templ.Integral()
print 'Background Integral = ',bkg_templ.Integral()
print 'Pseudo Data Integral = ',pseudo_data_templ.Integral()

# Rebin the histograms to make equal bins
bkg_templ.Rebin(2)
pseudo_data_templ.Rebin(2)

boundaries = optimizeBinBoundaries(bkg_templ, 500, LowFitRange, HighFitRange)

# binlist = [-1.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 10.0, 15.0, 20.0]
# print 'using custom bins for Data:',binlist
# boundaries = array.array('d')
# boundaries.fromlist(binlist)

# I can set the bin boundaries manually for fitting data.
print 'new number of boundaries = ',len(boundaries)
pseudo_data_templRbin = pseudo_data_templ.Rebin(len(boundaries)-1, pseudo_data_templ.GetName()+'rebin', boundaries)
bkg_templRbin = bkg_templ.Rebin(len(boundaries)-1, bkg_templ.GetName()+'rebin', boundaries)
sig_templRbin = sig_templ.Rebin(len(boundaries)-1, sig_templ.GetName()+'rebin', boundaries)

# to avoid signal going below fit range in the plot
print 'sig_templRbin underflow = ', sig_templRbin.GetBinContent(0)

sig_templRbin.SetBinContent(0, 0)

pseudo_data_templUni = makeUniformHist(pseudo_data_templRbin)
bkg_templUni = makeUniformHist(bkg_templRbin)
sig_templUni = makeUniformHist(sig_templRbin)

# using LowFitRange and HighFitRange to fo the fit doesn't work well
# to get the good fit result, do fit for fixed size histograms
LowUniFitRange = pseudo_data_templRbin.FindBin(LowFitRange)
if LowUniFitRange == 0:
#    print 'lower bin is underflow, setting to 1'
    LowUniFitRange = 1

HighUniFitRange = pseudo_data_templRbin.FindBin(HighFitRange) + 1
if pseudo_data_templRbin.GetBinLowEdge(HighUniFitRange-1) == HighFitRange:
#    print 'upper bin in on the border of fit range, reducing'
    HighUniFitRange -= 1

print 'fitting in the bin range = ',LowUniFitRange,': ', HighUniFitRange

# Fit results: overlap signal, bkg, sum of sig and bkg, and pseudo data
(fitSigFrac,fitSigFracErr) = makeFit(FitVariable + ' bin number', LowUniFitRange, HighUniFitRange, sig_templUni, bkg_templUni, pseudo_data_templUni, 'fit_'+FitVariable+'_uni.png')

#############################
leftbin = LowUniFitRange
rightbin = HighUniFitRange

c2 = TCanvas('c2', 'c2', 700, 600)

leg = TLegend(0.55, 0.65, 0.85, 0.85)
leg.SetFillColor(0)
leg.SetBorderSize(0)

print 'bins for normalized = ', leftbin, ' : ', rightbin
sig_templRbin.Scale(fitSigFrac*pseudo_data_templRbin.Integral(leftbin,rightbin)/sig_templRbin.Integral(leftbin,rightbin))
bkg_templRbin.Scale((1.0-fitSigFrac)*pseudo_data_templRbin.Integral(leftbin,rightbin)/bkg_templRbin.Integral(leftbin,rightbin))

sig_templRbin.SetLineColor(kGreen-4)
sig_templRbin.SetFillColor(kGreen-4)
sig_templRbin.GetXaxis().SetRangeUser(LowFitRange, HighFitRange)

bkg_templRbin.SetLineColor(kAzure)
bkg_templRbin.SetFillColor(kAzure)
bkg_templRbin.GetXaxis().SetRangeUser(LowFitRange, HighFitRange)

stack = THStack(FitVariable+'_stack', FitVariable+'_stack')
stack.Add(bkg_templRbin)
stack.Add(sig_templRbin)
stack.Draw('hist')
#stack.SetMaximum(120)
stack.GetXaxis().SetTitle('photon '+FitVariable + ' [GeV]')

pseudo_data_templRbin.SetMarkerColor(1)
pseudo_data_templRbin.SetMarkerStyle(8)
pseudo_data_templRbin.Draw('esame')

leg.AddEntry(pseudo_data_templRbin, 'Pseudo Data', 'lp')
leg.AddEntry(sig_templRbin, 'Signal', 'f')
leg.AddEntry(bkg_templRbin, 'Background', 'f')
leg.Draw('same')

c2.SaveAs(FitVariable+'_overlap.png')

leg.Clear()
