#! /bin/env python
'''

Subtract expected WZ and ZZ contamination from FR numerator and denominators.

Author: Evan K. Frii

'''

import logging
import sys
logging.basicConfig(stream=sys.stderr, level=logging.INFO)
from RecoLuminosity.LumiDB import argparse
import fnmatch
from FinalStateAnalysis.PlotTools.RebinView import RebinView
from FinalStateAnalysis.PlotTools.SubtractionView import SubtractionView
import glob
import os

log = logging.getLogger("CorrectFakeRateData")
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--files', nargs='+')
    parser.add_argument('--lumifiles', nargs='+')
    parser.add_argument('--outputfile', required=True)
    parser.add_argument('--denom', required=True, help='Path to denom')
    parser.add_argument('--numerator', required=True, help='Path to numerator')
    parser.add_argument('--rebin', type=str, default="1")

    args = parser.parse_args()

    from rootpy import io
    import ROOT
    import rootpy.plotting.views as views
    import rootpy.plotting as plotting
    from FinalStateAnalysis.MetaData.data_views import data_views

    files = []
    for pattern in args.files:
        files.extend(glob.glob(pattern))

    log.info("Loading data from %i files", len(files))

    lumifiles = []
    for pattern in args.lumifiles:
        lumifiles.extend(glob.glob(pattern))

    the_views = data_views(files, lumifiles)

    outputdir = os.path.dirname(args.outputfile)
    if outputdir and not os.path.exists(outputdir):
        os.makedirs(outputdir)

    log.info("Rebinning with factor %s", args.rebin)
    def rebin_view(x):
        ''' Make a view which rebins histograms '''
        binning = None
        binning = eval(args.rebin)
        return RebinView(x, binning)

    def round_to_ints(x):
        new = x.Clone()
        new.Reset()
        for bin in range(x.GetNbinsX()+1):
            binsy = range(x.GetNbinsY()+1) if isinstance(x, ROOT.TH2) else [-1]
            for biny in binsy:
                nentries = ROOT.TMath.Nint(x.GetBinContent(bin, biny)) \
                    if isinstance(x, ROOT.TH2) else \
                    ROOT.TMath.Nint(x.GetBinContent(bin))
                centerx = x.GetXaxis().GetBinCenter(bin)
                centery = x.GetYaxis().GetBinCenter(biny) \
                    if isinstance(x, ROOT.TH2) else \
                    0.
                for _ in range(nentries):
                    if isinstance(x, ROOT.TH2):
                        new.Fill(centerx, centery)
                    else:
                        new.Fill(centerx)
        return new

    def int_view(x):
        return views.FunctorView(x, round_to_ints)

    def get_view(sample_pattern):
        for sample, sample_info in the_views.iteritems():
            if fnmatch.fnmatch(sample, sample_pattern):
                return rebin_view(sample_info['view'])
        raise KeyError("I can't find a view that matches %s, I have: %s" % (
            sample_pattern, " ".join(the_views.keys())))

            #from pdb import set_trace; set_trace()
    wz_view = get_view('WZ*')
    zz_view = get_view('ZZ*')
    data = rebin_view(the_views['data']['view'])
    
    corrected_view = int_view(
        SubtractionView(data, wz_view, zz_view, restrict_positive=True))

    log.debug('creating output file')
    output = io.open(args.outputfile, 'RECREATE')
    output.cd()

    log.debug('getting from corrected view')
    corr_numerator = corrected_view.Get(args.numerator)
    corr_denominator = corrected_view.Get(args.denom)

    log.info("Corrected:   %0.2f/%0.2f = %0.1f%%",
             corr_numerator.Integral(),
             corr_denominator.Integral(),
             100*corr_numerator.Integral()/corr_denominator.Integral()
             if corr_denominator.Integral() else 0
            )

    uncorr_numerator = data.Get(args.numerator)
    uncorr_denominator = data.Get(args.denom)

    wz_integral = wz_view.Get(args.numerator).Integral()
    zz_integral = zz_view.Get(args.numerator).Integral()

    log.info("Numerator integrals data: %.2f WZ: %.2f, ZZ: %.2f. Corrected numerator: %.2f",
             uncorr_numerator.Integral(),
             wz_integral,
             zz_integral,
             corr_numerator.Integral()
            )

    log.info("Uncorrected: %0.2f/%0.2f = %0.1f%%",
             uncorr_numerator.Integral(),
             uncorr_denominator.Integral(),
             100*uncorr_numerator.Integral()/uncorr_denominator.Integral()
             if uncorr_denominator.Integral() else 0
            )


    corr_numerator.SetName('numerator')
    corr_denominator.SetName('denominator')

    uncorr_numerator.SetName('numerator_uncorr')
    uncorr_denominator.SetName('denominator_uncorr')

    corr_numerator.Write()
    corr_denominator.Write()
    uncorr_numerator.Write()
    uncorr_denominator.Write()

    #make the unbinned plots
    args.rebin = '1'
    wz_view    = get_view('WZ*')
    zz_view    = get_view('ZZ*')
    data       = the_views['data']['view']
    
    corrected_view = int_view(
        SubtractionView(data, wz_view, zz_view, restrict_positive=True))

    corr_numerator_unrebinned     = corrected_view.Get(args.numerator)
    corr_denominator_unrebinned   = corrected_view.Get(args.denom)
    uncorr_numerator_unrebinned   = data.Get(args.numerator)
    uncorr_denominator_unrebinned = data.Get(args.denom)

    corr_numerator_unrebinned.SetName('numerator_unrebinned')
    corr_denominator_unrebinned.SetName('denominator_unrebinned')
    uncorr_numerator_unrebinned.SetName('numerator_uncorr_unrebinned')
    uncorr_denominator_unrebinned.SetName('denominator_uncorr_unrebinned')

    corr_numerator_unrebinned.Write()
    corr_denominator_unrebinned.Write()
    uncorr_numerator_unrebinned.Write()
    uncorr_denominator_unrebinned.Write()

