#!/usr/bin/env python

'''

Fit the probability of mismeasuring the charge of an electron. Each histogram contains the dielectron mass spectrum inside a bin

Usage:
    fitElectronChargeFlipProbability.py input.root output.root path/to/os path/to/ss "Z peak functional shape" "background functional shape"

Where [efficiency] is a RooFit factory command.

'''

import array
from RecoLuminosity.LumiDB import argparse
from FinalStateAnalysis.PlotTools.RebinView import RebinView
from FinalStateAnalysis.PlotTools.THBin import zipBins
import logging
import sys
args = sys.argv[:]
sys.argv = [sys.argv[0]]

log = logging.getLogger("make_charge_flip_probability_map")

def get_x_binning(histo):
    bin_edges = []
    for i in range(histo.GetNbinsX()+1):
        bin_edges.append(histo.GetXaxis().GetBinLowEdge(i+1))
    return array.array('d', bin_edges)

def get_y_binning(histo):
    bin_edges = []
    for i in range(histo.GetNbinsY()+1):
        bin_edges.append(histo.GetYaxis().GetBinLowEdge(i+1))
    return array.array('d', bin_edges)

## def rebin2d(histo, rebinx, rebiny):
##     # Just merging bins
##     if isinstance(binning, int):
##         return histo.Rebin2D(rebinx, rebiny, histo.GetName()+'_rebin')
##     bin_arrayx = array.array('d', rebinx)
##     bin_arrayy = array.array('d', rebiny)
##     new_histo  = ROOT.TH2(histo.GetName() + '_rebin', histo.GetTitle(), len(binning)-1) histogram.Rebin(len(binning)-1, histogram.GetName() + 'rebin', bin_array)
    


if __name__ == "__main__":
    #################################
    ###Argparse options  ############
    #################################
    parser = argparse.ArgumentParser()
    parser.add_argument('output', metavar='output.root', help='Output root file')
    parser.add_argument('num', metavar='/path/to/numerator',   help='Path to numerator histogram')
    parser.add_argument('den', metavar='/path/to/denominator', help='Path to denominator histogram')
    parser.add_argument('input', nargs='+', metavar='input.root', help='Input root files - will be summed')
    parser.add_argument('--verbose', action='store_true', help='More log output')
    parser.add_argument('--rebinX', metavar='N', type=int, required=False, help='Rebin histograms before fitting. '
                        #'Can be either an integer or a list of bin edges.'
                        #' If variable binning is used, the numbers should '
                        #' be specified as a comma separated list w/o spaces.'
        )
    parser.add_argument('--rebinY', metavar='N', type=int, required=False, help='Rebin histograms before fitting. '
                        #'Can be either an integer or a list of bin edges.'
                        #' If variable binning is used, the numbers should '
                        #' be specified as a comma separated list w/o spaces.'
        )


    args = parser.parse_args(args[1:])

    if args.verbose:
        logging.basicConfig(level=logging.DEBUG, stream=sys.stderr)
    else:
        logging.basicConfig(level=logging.INFO, stream=sys.stderr)

    #################################
    ###main part         ############
    #################################
    from rootpy.plotting import views
    import rootpy.io as io
    import ROOT
    ROOT.gROOT.SetBatch(True)

    # Build view of input histograms
    log.info("Merging input files")
    input_files = [io.open(x) for x in args.input]
    input_view  = views.SumView(*input_files)
    numerator   = input_view.Get(args.num)
    denominator = input_view.Get(args.den)

    #Checks if rebinning is needed
    if (args.rebinX and args.rebinX > 1) or (args.rebinY and args.rebinY > 1):
        binx = args.rebinX if args.rebinX else 1
        biny = args.rebinY if args.rebinY else 1
        numerator.Rebin2D(binx, biny)
        denominator.Rebin2D(binx, biny) 
        
    eff_map          = numerator.Clone("efficiency_map")
    eff_map_statUp   = numerator.Clone("efficiency_map_statUp")
    eff_map_statDown = numerator.Clone("efficiency_map_statDown")
    worse_rel_err    = numerator.Clone("worse_rel_err")

    eff_map.SetTitle( eff_map.GetName() )
    eff_map_statUp.SetTitle( eff_map_statUp.GetName() )
    eff_map_statDown.SetTitle( eff_map_statDown.GetName() )
    worse_rel_err.SetTitle( worse_rel_err.GetName() )
    
    log.info("Making efficiency map")
    for ibinx in range(1,eff_map.GetNbinsX()+1):
        for ibiny in range(eff_map.GetNbinsY()+1):
            inum      = numerator.GetBinContent(ibinx, ibiny)
            iden      = denominator.GetBinContent(ibinx, ibiny)
            ieff       = inum / iden if iden <> 0 else 0
            statUp    = ROOT.TEfficiency.ClopperPearson( int(iden), int(inum), 0.682689492137, True)
            statDown  = ROOT.TEfficiency.ClopperPearson( int(iden), int(inum), 0.682689492137, False)
            worse_err = (abs(statUp - ieff)/ieff if abs(statUp - ieff) > abs(statDown - ieff) else abs(statDown - ieff)/ieff) if ieff else 0.
            eff_map.SetBinContent(ibinx, ibiny, ieff)
            eff_map.SetBinError(ibinx, ibiny, worse_err*ieff)
            eff_map_statUp.SetBinContent(ibinx, ibiny, statUp)
            eff_map_statDown.SetBinContent(ibinx, ibiny, statDown)
            worse_rel_err.SetBinContent(ibinx, ibiny, worse_err)
            
    outFile = ROOT.TFile(args.output,'recreate') #FIXME move to rootpy io
    outFile.cd()
    eff_map.Write()
    eff_map_statUp.Write()
    eff_map_statDown.Write()
    worse_rel_err.Write()
    outFile.Close()
