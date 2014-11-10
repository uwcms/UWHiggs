import os
import re
import glob
import ROOT
from math import sqrt
from FinalStateAnalysis.StatTools.RooFunctorFromWS import build_roofunctor, make_corrector_from_histo

frfit_dir = os.path.join('results', os.environ['jobid'], 'fakerate_fits')+'/'
        
tau_fake_rate = build_roofunctor(
    frfit_dir+'t_os_tLoose_tTigh_tAbsEta.root', 
    'fit_efficiency', 
    'efficiency'
)

def inflate_systematics(functor, inflation):
    '''Christian recipe: INFLATE ALL THE SYSTEMATICS!'''
    def fcn(*args):
        return functor(*args) * inflation
    return fcn

    #tau_fake_rate_up = inflate_systematics(tau_fake_rate, 1.3)
    #tau_fake_rate_dw = inflate_systematics(tau_fake_rate, 0.7)

###Tau fakes as it should be
tau_fake_rate_up = make_corrector_from_histo(
    frfit_dir+'t_os_tLoose_tTigh_tAbsEta_2.root', 
    'efficiency_up', 
    '1D'
)

tau_fake_rate_dw = make_corrector_from_histo(
    frfit_dir+'t_os_tLoose_tTigh_tAbsEta_2.root', 
    'efficiency_dw', 
    '1D'
)


