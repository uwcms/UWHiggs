import os
from FinalStateAnalysis.StatTools.RooFunctorFromWS import build_roofunctor

################################################################################
#### Fitted fake rate functions ################################################
################################################################################

# Get fitted fake rate functions
frfit_dir = os.path.join('results', os.environ['jobid'], 'fakerate_fits')
mu_fr = build_roofunctor( #lambda x: 0.1
    frfit_dir + '/m_zlt_pt10_idPassed_muonPt.root',
    'fit_efficiency', # workspace name
    'efficiency'
)

e_fr =build_roofunctor( #lambda x: 0.2
    frfit_dir + '/e_zlt_pt10_idPassed_electronPt.root',
    'fit_efficiency', # workspace name
    'efficiency'
)


tau_fr = build_roofunctor(
    frfit_dir + '/t_ztt_pt10_MediumIso_tauPt.root',
    'fit_efficiency', # workspace name
    'efficiency'
)
