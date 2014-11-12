#! /bin/env python

import ROOT
import fnmatch
from optparse import OptionParser
import os
import re
import logging
import sys

__author__  = "Mauro Verzetti (mauro.verzetti@cern.ch)"
__doc__ = """checks that the shapes for the limit are looking good"""

parser = OptionParser(description=__doc__)
parser.add_option('--fakes', '-f', type=str, default = '*fakes*',
                  help='pattern to match to fakes',dest='fakes')
parser.add_option('--exclude', '-e', type=str, default = '*_?',
                  help='pattern of dirs to exclude',dest='exclude')
parser.add_option('--signal', '-s', type=str, default = 'WH*',
                  help='pattern of plots to be recognised as signals',dest='signal')
parser.add_option('--shape-unc-matcher', default='*CMS_*',
                    dest='shapematch',
                    help='Shell glob-style matcher for shape errors.  These shapes arent shown.')
parser.add_option('--level', '-l', type=str, default = 'INFO',
                  help='pattern of dirs to excludelogging level', dest='level')

args, shapes_filenames = parser.parse_args()

logging.basicConfig(stream=sys.stderr, level=getattr(logging, args.level))

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat('111111111')

def walk(inputdir):
    ''' Recursive function which generates (path, subdirs, histos) '''
    directories = []
    histos = []
    for key in inputdir.GetListOfKeys():
        # Keep track of stuff we find in this directory
        name = key.GetName()
        classname = key.GetClassName()
        if classname.startswith('TDirectory'):
            directories.append(name)
        elif isinstance(inputdir.Get(name), ROOT.TH1):
            histos.append(name)
        path = inputdir.GetPath().split(':')[1]
        # Yield the stuff in this directory
    yield (path, tuple(directories), tuple(histos))
    # Now get the stuff in sub directories
    for directory in directories:
        for subresult in walk(inputdir.Get(directory)):
            yield subresult

def seperate_histos(histo_list, signal_pattern, err_pattern, fakes_pattern):
    ''' Separate histogram lists into backgrounds, shape_uncs, and signals '''
    signals = []
    shapes = []
    bkgs = []
    fakes = []
    for x in histo_list:
        if 'obs' in x:
            continue
        if fnmatch.fnmatch(x, signal_pattern):
            signals.append(x)
        elif fnmatch.fnmatch(x, err_pattern):
            shapes.append(x)
        elif fnmatch.fnmatch(x, fakes_pattern):
            fakes.append(x)
        else:
            bkgs.append(x)
    return (tuple(signals), tuple(shapes), tuple(fakes), tuple(bkgs))

def count_holes(histo):
    nbins_x = histo.GetNbinsX()
    last_filled = 0
    first_filled = nbins_x
    for i in xrange(nbins_x, 1, -1):
        if histo.GetBinContent(i):
            last_filled = i
            break
    for i in xrange(1, nbins_x+1):
        if histo.GetBinContent(i):
            first_filled = i
            break
    
    if last_filled < first_filled:
        return -1
        #raise Exception("Error in count_holes")
    
    holes = sum(
        1 for i in xrange(first_filled, last_filled) 
        if not histo.GetBinContent(i)
    )
    return holes

def count_negatives(histo):
    nbins_x = histo.GetNbinsX()
    negatives = sum(
        1 for i in xrange(1, nbins_x+1) 
        if (histo.GetBinContent(i) + histo.GetBinError(i)) < 0
    )
    return negatives

def check_empty(histo):
    return histo.Integral() == 0

def check_sig_bkg(signals, backgrounds):
    nbins_x = signals[0].GetNbinsX() #same as 
    bins    = 0
    for ibin in xrange(1, nbins_x+1):
        bkg_sum = sum(
            h.GetBinContent(ibin) for h in backgrounds
            if h.GetBinContent(ibin) >= 0
            )
        if not bkg_sum: #in background is empty
            #check signals
            bins += int(
                any(
                    h.GetBinContent(ibin) > 0 
                    for h in signals
                    )
                )
    return bins

def pair_signals(signals):
    regex = re.compile("[a-zA-Z_]+(?P<mass>\d+)")
    pairs = {}
    for name in signals:
        logging.debug(name)
        mass = regex.match(name).group('mass')
        pairs[mass] = pairs.get(mass,[]) + [name]
    return pairs


good_shapes = []
bad_shapes  = []
ugly_shapes = [] #not used

for shape_filename in shapes_filenames:
    logging.info( "Inspecting %s" % shape_filename)
    shape_file = ROOT.TFile.Open(shape_filename)
    
    info = {}
    for path, subdirs, histos in walk(shape_file):
        if not histos or fnmatch.fnmatch(path, args.exclude):
            continue

        h_dict = {}
        for h in histos:
            h_dict[h] = shape_file.Get(os.path.join(path, h))

        logging.debug('%s, %s, %s' % (path.__repr__(), subdirs.__repr__(), histos.__repr__()))
        signals, unc, fakes, bkgs = seperate_histos(
            histos, args.signal, args.shapematch, args.fakes)
        
        paired_signals = pair_signals(signals)
        logging.debug('paired signals: %s' % paired_signals)
        signals = []
        for mass, names in paired_signals.iteritems():
            if len(names) == 1:
                signals += names
            if len(names) == 2:
                h_dict[names[0]].Add(h_dict[names[1]])
                signals.append(names[0])
        
        logging.debug('new signals: %s' % signals)
        signals = tuple(signals)
        #count holes
        for h in fakes+bkgs:
            hpath = os.path.join(path, h)
            info[hpath] = {}
            info[hpath]['holes'] = count_holes(h_dict[h])

        for h in signals:
            info[hpath]['empty'] = check_empty(h_dict[h])

        #count negatives
        for h in fakes:
            hpath = os.path.join(path, h)
            info[hpath]['negatives'] = count_negatives(h_dict[h])
        
        h_sigs = [ h_dict[h] for h in signals ]
        h_bkgs = [ h_dict[h] for h in fakes+bkgs ]
        info[path] = {
            'signal_only' : check_sig_bkg(h_sigs, h_bkgs)
            }

    has_errors = False
    for path, values in info.iteritems():
        if values.get('holes') or values.get('negatives'):
            has_errors = True
            logging.info( '%s found with %i holes and %i negative bins!' % (path, values['holes'], values.get('negatives', 0)) )
        if values.get('empty'):
            has_errors = True
            logging.info( '%s found to be empty!' % path )            
        if values.get('signal_only'):
            has_errors = True
            logging.info('category %s has %i bins with signal-only contribution!' % (path, values['signal_only']))
    if has_errors:
        bad_shapes.append(shape_filename)
        logging.warning('%s has errors!' % shape_filename)
    else: 
        good_shapes.append(shape_filename)
    logging.info( "\n\n\n")
    shape_file.Close()



print '\n\nBAD SHAPE FILES:\n\n'
print '\n'.join( bad_shapes )


print '\n\nGOOD SHAPE FILES:\n\n'
print '\n'.join( good_shapes )
