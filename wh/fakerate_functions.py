import os
from FinalStateAnalysis.StatTools.RooFunctorFromWS import build_roofunctor, make_corrector_from_th2
from FinalStateAnalysis.StatTools.VariableScaler import make_scaler
from optimizer import leading_lepton_iso_tag, subleading_lepton_iso_tag

################################################################################
#### Fitted fake rate functions ################################################
################################################################################

# Get fitted fake rate functions
frfit_dir = os.path.join('results', os.environ['jobid'], 'fakerate_fits')


##################
## 1D Muon Func ##
##################

highpt_mu_fr = build_roofunctor(
    frfit_dir + '/m_wjets_pt20_%s_muonJetPt.root' % leading_lepton_iso_tag,
    'fit_efficiency', # workspace name
    'efficiency'
)
lowpt_mu_fr = build_roofunctor(
    frfit_dir + '/m_wjets_pt10_%s_muonJetPt.root' % subleading_lepton_iso_tag,
    'fit_efficiency', # workspace name
    'efficiency'
)

highpt_mu_qcd_fr = build_roofunctor(
    frfit_dir + '/m_qcd_pt20_%s_muonJetPt.root' % leading_lepton_iso_tag,
    'fit_efficiency', # workspace name
    'efficiency'
)

lowpt_mu_qcd_fr = build_roofunctor(
    frfit_dir + '/m_qcd_pt10_%s_muonJetPt.root' % subleading_lepton_iso_tag,
    'fit_efficiency', # workspace name
    'efficiency'
)

#######################
## 1D Electrons Func ##
#######################

lowpt_e_qcd_fr = build_roofunctor(
    frfit_dir + '/e_qcd_pt10_%s_eJetPt.root' % subleading_lepton_iso_tag,
    'fit_efficiency', # workspace name
    'efficiency'
)

lowpt_e_fr = build_roofunctor(
    frfit_dir + '/e_wjets_pt10_%s_eJetPt.root' % subleading_lepton_iso_tag,
    'fit_efficiency', # workspace name
    'efficiency'
)

highpt_e_qcd_fr = build_roofunctor(
    frfit_dir + '/e_qcd_pt20_%s_eJetPt.root' % leading_lepton_iso_tag,
    'fit_efficiency', # workspace name
    'efficiency'
)

highpt_e_fr = build_roofunctor(
    frfit_dir + '/e_wjets_pt20_%s_eJetPt.root' % leading_lepton_iso_tag,
    'fit_efficiency', # workspace name
    'efficiency'
)

highpt_ee_fr = build_roofunctor(
    frfit_dir + '/ee_wjetsNoZmass_pt20_%s_electronJetPt.root' % leading_lepton_iso_tag,
    'fit_efficiency', # workspace name
    'efficiency'
)

lowpt_ee_fr = build_roofunctor(
    frfit_dir + '/ee_wjetsNoZmass_pt10_%s_electronJetPt.root' % subleading_lepton_iso_tag,
    'fit_efficiency', # workspace name
    'efficiency'
)

highpt_ee_qcd_fr = build_roofunctor(
    frfit_dir + '/ee_qcd_pt20_%s_electronJetPt.root' % leading_lepton_iso_tag,
    'fit_efficiency', # workspace name
    'efficiency'
)

lowpt_ee_qcd_fr = build_roofunctor(
    frfit_dir + '/ee_qcd_pt10_%s_electronJetPt.root' % subleading_lepton_iso_tag,
    'fit_efficiency', # workspace name
    'efficiency'
)



##################
## 1D Taus Func ##
##################

tau_fr = build_roofunctor(
    frfit_dir + '/t_ztt_pt20_mvaloose_tauPt.root',
    'fit_efficiency', # workspace name
    'efficiency'
)

tau_qcd_fr = tau_fr 

e_charge_flip      = make_corrector_from_th2(frfit_dir+"/charge_flip_prob_map.root", "efficiency_map")         
e_charge_flip_up   = make_corrector_from_th2(frfit_dir+"/charge_flip_prob_map.root", "efficiency_map_statUp")  
e_charge_flip_down = make_corrector_from_th2(frfit_dir+"/charge_flip_prob_map.root", "efficiency_map_statDown")
mass_scaler        = make_scaler(frfit_dir+"/charge_flip_prob_map.root", 'mass_scale')


