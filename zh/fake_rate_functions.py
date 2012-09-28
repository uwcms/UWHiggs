import os
from FinalStateAnalysis.StatTools.RooFunctorFromWS import build_roofunctor

################################################################################
#### Fitted fake rate functions ################################################
################################################################################

# Get fitted fake rate functions
frfit_dir = os.path.join('results', os.environ['jobid'], 'fakerate_fits')
mu_fr = lambda x: 0.01
## build_roofunctor(
##     frfit_dir + '/m_wjets_pt20_pfidiso02_muonJetPt.root',
##     'fit_efficiency', # workspace name
##     'efficiency'
## )
e_fr = lambda x: 0.01
## build_roofunctor(
##     frfit_dir + '/m_wjets_pt10_pfidiso02_muonJetPt.root',
##     'fit_efficiency', # workspace name
##     'efficiency'
## )


tau_fr = lambda x: 0.01
## build_roofunctor(
##     frfit_dir + '/t_ztt_pt20_mvaloose_tauPt.root',
##     'fit_efficiency', # workspace name
##     'efficiency'
## )
