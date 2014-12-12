from ETauTree import ETauTree
import sys
import logging
logging.basicConfig(stream=sys.stderr, level=logging.WARNING)
import os
from pdb import set_trace
import ROOT
import math
import glob
import array
import mcCorrections
import baseSelections as selections
import FinalStateAnalysis.PlotTools.pytree as pytree
from FinalStateAnalysis.PlotTools.decorators import  memo_last
from FinalStateAnalysis.PlotTools.MegaBase import MegaBase
from math import sqrt, pi, cos
from fakerate_functions import tau_fake_rate, tau_fake_rate_up, tau_fake_rate_dw, e_fake_rate, e_fake_rate_up, e_fake_rate_dw
import itertools
import traceback
from FinalStateAnalysis.PlotTools.decorators import memo
import FinalStateAnalysis.PlotTools.pytree as pytree
from FinalStateAnalysis.Utilities.struct import struct
import optimizer

@memo
def getVar(name, var):
    return name+var

met_et  = 'pfMet_Et%s'
met_phi = 'pfMet_Phi%s'
t_pt  = 'tPt%s'
etMass = 'e_t_Mass%s'
@memo
def met(shift=''):
    return met_et % shift

@memo
def metphi(shift=''):
    return met_phi % shift
@memo
def tpt(shift=''):
    return t_pt % shift

@memo
def vismass(shift=''):
    return etMass % shift

@memo
def split(string, separator='#'):
    return tuple(attr.split(separator))

def create_mapper(mapping):
    def _f(path):
        for key, out in mapping.iteritems():
            if key == path:
                path = path.replace(key,out)
                print 'path', path
        return path
    return _f

def attr_getter(attribute):
    '''return a function that gets an attribute'''
    def f(row, weight):
        return (getattr(row,attribute), weight)
    return f

def merge_functions(fcn_1, fcn_2):
    '''merges two functions to become a TH2'''
    def f(row, weight):
        r1, w1 = fcn_1(row, weight)
        r2, w2 = fcn_2(row, weight)
        w = w1 if w1 and w2 else None
        return ((r1, r2), w)
    return f

def collmass(row, met, metPhi):
    ptnu =abs(met*cos(deltaPhi(metPhi, row.tPhi)))
    visfrac = row.tPt/(row.tPt+ptnu)
    #print met, cos(deltaPhi(metPhi, row.tPhi)), ptnu, visfrac
    return (row.e_t_Mass / sqrt(visfrac))

def deltaPhi(phi1, phi2):
    PHI = abs(phi1-phi2)
    if PHI<=pi:
        return PHI
    else:
        return 2*pi-PHI

def deltaR(phi1, ph2, eta1, eta2):
    deta = eta1 - eta2
    dphi = abs(phi1-phi2)
    if (dphi>pi) : dphi = 2*pi-dphi
    return sqrt(deta*deta + dphi*dphi);

def make_collmass_systematics(shift):
    if shift.startswith('tes'):
        ptnu =abs(met*cos(deltaPhi(metPhi, row.tPhi)))
        visfrac = tpt/(tpt+ptnu)
        vis_mass = vismass(shift)
        return (vis_mass / sqrt(visfrac))
    else:
        met_name = met(shift)
        phi_name = metphi(shift)
        def collmass_shifted(row, weight):
            met = getattr(row, met_name)
            phi = getattr(row, phi_name)
            return collmass(row, met, phi), weight
        return collmass_shifted

class LFVHETauAnalyzerMVA(MegaBase):
    tree = 'et/final/Ntuple'
    def __init__(self, tree, outfile, **kwargs):
        logging.debug('LFVHETauAnalyzerMVA constructor')
        self.channel='ET'
        super(LFVHETauAnalyzerMVA, self).__init__(tree, outfile, **kwargs)
        self.tree = ETauTree(tree)
        self.out=outfile
        self.histograms = {}

        #understand what we are running
        target = os.path.basename(os.environ['megatarget'])
        self.is_data = target.startswith('data_')
        self.is_embedded = ('Embedded' in target)
        self.is_mc = not (self.is_data or self.is_embedded)
        self.efake   = e_fake_rate(0.2)
        self.efakeup = e_fake_rate_up(0.2)
        self.efakedw = e_fake_rate_dw(0.2)

        #systematics used
        self.systematics = {
            'trig' : (['', 'trp1s', 'trm1s'] if not self.is_data else []),
            'pu'   : (['', 'p1s', 'm1s'] if self.is_mc else []),
            'eid'  : (['', 'eidp1s','eidm1s']   if False else []),
            'eiso' : (['', 'eisop1s','eisom1s'] if False else []),
            'jes'  : (['', '_jes_plus','_jes_minus'] if self.is_mc else ['']),
            'mvetos': (['', 'mVetoUp', 'mVetoDown'] if False else ['']),
            'tvetos': (['', 'tVetoUp', 'tVetoDown'] if False else ['']),
            'evetos': (['', 'eVetoUp', 'eVetoDown'] if False else ['']),
            'met'  : (["_mes_plus", "_ues_plus", "_mes_minus", "_ues_minus"] if self.is_mc else []),
            'tes'  : (["", "_tes_plus", "_tes_minus"] if not self.is_data else ['']),
            'ees'  : (["", "_ees_plus", '_ees_minus'] if not self.is_data else [''])
        }

        #self filling histograms
        coll_mass = make_collmass_systematics('') #no sys shift
        self.histo_locations = {} #just a mapping of the histograms we have to avoid changing self.histograms indexing an screw other files
        self.hfunc   = { #maps the name of non-trivial histograms to a function to get the proper value, the function MUST have two args (evt and weight). Used in fill_histos later
            'nTruePU' : lambda row, weight: (row.nTruePU,None),
            'weight'  : lambda row, weight: (weight,None) if weight is not None else (1.,None),
            'Event_ID': lambda row, weight: (array.array("f", [row.run,row.lumi,int(row.evt)/10**5,int(row.evt)%10**5] ), None),
            'h_collmass_pfmet' : coll_mass,
            'h_collmass_vs_dPhi_pfmet' : merge_functions(
                attr_getter('tToMETDPhi'),
                coll_mass
            ),
            'MetEt_vs_dPhi' : merge_functions(
                lambda row, weight: (deltaPhi(row.tPhi, getattr(row, metphi())), weight),
                attr_getter('type1_pfMet_Et')
            ),
            'ePFMET_DeltaPhi' : lambda row, weight: (deltaPhi(row.ePhi, getattr(row, metphi())), weight),
            'tPFMET_DeltaPhi' : lambda row, weight: (deltaPhi(row.tPhi, getattr(row, metphi())), weight),
            'evtInfo' : lambda row, weight: (struct(run=row.run,lumi=row.lumi,evt=row.evt,weight=weight), None)
            }
        for shift in self.systematics['met']:
            #patch name
            postfix = shift
            self.hfunc['h_collmass_pfmet%s' % postfix] = make_collmass_systematics(shift)
        for shift in self.systematics['tes']:
            #patch name
            postfix = shift
            self.hfunc['h_collmass_pfmet%s' % postfix] = make_collmass_systematics(shift)
        for shift in self.systematics['jes']:
            #patch name
            postfix = shift
            self.hfunc['h_collmass_pfmet%s' % postfix] = make_collmass_systematics(shift)
        for shift in self.systematics['ees']:
            #patch name
            postfix = shift
            self.hfunc['h_collmass_pfmet%s' % postfix] = make_collmass_systematics(shift)
        

        #PU correctors
        self.pucorrector = mcCorrections.make_shifted_weights(
            mcCorrections.make_puCorrector('singlee'),
            ['p1s', 'm1s'],
            [mcCorrections.make_puCorrectorUp('singlee'), mcCorrections.make_puCorrectorDown('singlee')]
        )     
        self.trig_weight = mcCorrections.trig_efficiency if self.is_embedded else mcCorrections.trig_correction

    @staticmethod 
    def tau_veto(row):
        if not row.tAntiMuonLoose2 or not row.tAntiElectronMVA5Tight or not row.tDecayFinding :
            return False

    @staticmethod
    def obj1_matches_gen(row):
        return row.eGenPdgId == -1*row.eCharge*11

    @staticmethod 
    def obj3_matches_gen(row):
        return t.genDecayMode != -2 

    def event_weight(self, row, sys_shifts):
        if self.is_data:
            return {'' : 1.}

        weights = {}
        for shift in sys_shifts:
            embedded_weight = row.EmbPtWeight*mcCorrections.eEmb_correction( row, 'e', shift=shift) if self.is_embedded else 1.

            weights[shift] = embedded_weight *\
                             mcCorrections.eid_correction( row, 'e', shift=shift) * \
                             mcCorrections.eiso_correction(row, 'e', shift=shift) * \
                             self.trig_weight(row, 'e', shift=shift) * \
                             self.pucorrector(row.nTruePU, shift=shift)
                       
        return weights
## 
    def begin(self):
        logging.debug('Booking histograms directory tree')
        sys_shifts = self.systematics['trig'] + \
                     self.systematics['pu'] + \
                     self.systematics['eid'] + \
                     self.systematics['eiso'] + \
                     self.systematics['mvetos'] + \
                     self.systematics['tvetos'] + \
                     self.systematics['evetos'] + \
                     [i.strip('_') for i in self.systematics['jes']] + \
                     [i.strip('_') for i in self.systematics['tes']] + \
                     [i.strip('_') for i in self.systematics['ees']] + \
                     ['tLoose', 'tLoose/Up', 'tLoose/Down', 'tLooseUnweight'] + \
                     ['eLoose', 'eLoose/Up', 'eLoose/Down'] +\
                     ['etLoose', 'etLoose/Up', 'etLoose/Down'] 
        sys_shifts = list( set( sys_shifts ) ) #remove double dirs
        processtype=['gg']
        threshold=['ept30']
        signs =['os', 'ss']
        jetN = ['0', '1', '2', '3'] #[''.join(i) for i in itertools.product(['0', '1', '2', '3'], self.systematics['jes'])]
        full_met_systematics = self.systematics['met']+self.systematics['jes']+ \
            self.systematics['tes']+self.systematics['ees']
        #remove empty string
        full_met_systematics = [i for i in full_met_systematics if i]
        folder=[]

        for tuple_path in itertools.product(sys_shifts, signs, processtype, threshold, jetN):
             folder.append(os.path.join(*tuple_path))
             path = list(tuple_path)
             path.append('selected')
             folder.append(os.path.join(*path))
             prefix_path = os.path.join(*tuple_path)
             for region in optimizer.regions[tuple_path[-1]]:
                 folder.append(
                     os.path.join(prefix_path, region)
                     )

        def book_with_sys(location, name, *args, **kwargs):
            postfixes = kwargs['postfixes']
            del kwargs['postfixes']
            self.book(location, name, *args, **kwargs)
            for postfix in postfixes:
                #patch name to be removed
                fix = postfix
                self.book(location, name+fix, *args, **kwargs)

        self.book('os/gg/ept30/', "h_collmass_pfmet" , "h_collmass_pfmet",  32, 0, 320)
        self.book('os/gg/ept30/', "e_t_Mass",  "h_vismass",  32, 0, 320)

        for f in folder: 
            #self.book(
            #    f,
            #    'evtInfo', 'evtInfo',
            #    'run/l:lumi/l:evt/l:weight/D',
            #    type=pytree.PyTree
            #)
            #@# self.book(f,"weight", "weight", 100, 0, 10)
            #@# self.book(f,"tPt", "tau p_{T}", 40, 0, 200)
            #@# self.book(f,"tPt_tes_plus", "tau p_{T} (tes+)", 40, 0, 200)
            #@# self.book(f,"tPt_tes_minus", "tau p_{T} (tes-)",40, 0, 200)
            #@# 
            #@# self.book(f,"tPhi", "tau phi", 26, -3.25, 3.25)
            #@# self.book(f,"tEta", "tau eta",  10, -2.5, 2.5)
            #@# 
            #@# self.book(f,"ePt", "e p_{T}", 40, 0, 200)
            #@# self.book(f,"ePt_ees_plus", "e p_{T} (ees+)", 40, 0, 200)
            #@# self.book(f,"ePt_ees_minus", "e p_{T} (ees-)",40, 0, 200)
            #@# 
            #@# self.book(f,"ePhi", "e phi",  26, -3.2, 3.2)
            #@# self.book(f,"eEta", "e eta", 10, -2.5, 2.5)
            #@# 
            #@# self.book(f, "e_t_DPhi", "e-tau DeltaPhi" , 20, 0, 3.2)
            #@# self.book(f, "e_t_DR", "e-tau DeltaR" , 20, 0, 3.2)
            
            #self.book(f, "h_collmass_pfmet",  "h_collmass_pfmet",  32, 0, 320)
            book_with_sys(f, "h_collmass_pfmet",  "h_collmass_pfmet",  40, 0, 400, 
                          postfixes=full_met_systematics)

            #@# self.book(f, "h_collmass_vs_dPhi_pfmet",  "h_collmass_vs_dPhi_pfmet", 20, 0, 3.2, 40, 0, 400, type=ROOT.TH2F)
            #@# 
            #@# self.book(f, "e_t_Mass",  "h_vismass",  40, 0, 400)
            #@# self.book(f, "e_t_Mass_tes_plus" ,  "h_vismass_tes_plus", 40 , 0, 400)
            #@# self.book(f, "e_t_Mass_tes_minus",  "h_vismass_tes_minus",40 , 0, 400)
            #@# self.book(f, "e_t_Mass_ees_plus" ,  "h_vismass_ees_plus", 40 , 0, 400)
            #@# self.book(f, "e_t_Mass_ees_minus",  "h_vismass_ees_minus",40 , 0, 400)
            #@# 
            #@# self.book(f, "MetEt_vs_dPhi", "PFMet vs #Delta#phi(#tau,PFMet)", 20, 0, 3.2, 40, 0, 400, type=ROOT.TH2F)
            #@# 
            #@# self.book(f, "tPFMET_DeltaPhi", "tau-type1PFMET DeltaPhi" , 20, 0, 3.2)
            #@# 
            #@# self.book(f, "ePFMET_DeltaPhi", "e-PFMET DeltaPhi" , 20, 0, 3.2)
            #@# 
            #@# #self.book(f,"tMtToPFMET", "tau-PFMET M_{T}" , 200, 0, 200)
            #@# book_with_sys(f, "tMtToPfMet", "tau-PFMET M_{T}" , 40, 0, 200,
            #@#               postfixes=full_met_systematics)
            #@# #self.book(f,"eMtToPFMET", "e-PFMET M_{T}" , 200, 0, 200)
            #@# book_with_sys(f, "eMtToPfMet", "e-PFMET M_{T}" , 40, 0, 200,
            #@#               postfixes=full_met_systematics)
            #@# 
            #@# #self.book(f, "pfMetEt",  "pfMetEt",  200, 0, 200)
            #@# book_with_sys(f, "pfMet_Et",  "pfMet_Et",  40, 0, 200, postfixes=full_met_systematics)
            #@# 
            #@# #self.book(f, "pfMetPhi",  "pfMetPhi", 100, -3.2, 3.2)
            #@# book_with_sys(f, "pfMet_Phi",  "pfMet_Phi", 26, -3.2, 3.2, postfixes=full_met_systematics)
            #@#  
            #@# self.book(f, "jetVeto20", "Number of jets, p_{T}>20", 5, -0.5, 4.5) 
            #@# self.book(f, "jetVeto30", "Number of jets, p_{T}>30", 5, -0.5, 4.5) 
        
        #index dirs and histograms
        for key, value in self.histograms.iteritems():
            location = os.path.dirname(key)
            name     = os.path.basename(key)
            if location not in self.histo_locations:
                self.histo_locations[location] = {name : value}
            else:
                self.histo_locations[location][name] = value

    def fakerate_weights(self, tEta): 
        tLoose    = tau_fake_rate(tEta)
        tLooseUp  = tau_fake_rate_up(tEta) 
        tLooseDown= tau_fake_rate_dw(tEta) 
        
        tLoose    = tLoose     / (1. - tLoose    ) 
        tLooseUp  = tLooseUp   / (1. - tLooseUp  ) 
        tLooseDown= tLooseDown / (1. - tLooseDown) 
        
 
        eLoose    = self.efake   
        eLooseUp  = self.efakeup
        eLooseDown= self.efakedw
        
        eLoose    = eLoose     / (1. - eLoose    ) 
        eLooseUp  = eLooseUp   / (1. - eLooseUp  ) 
        eLooseDown= eLooseDown / (1. - eLooseDown) 


        etLoose    = tLoose * eLoose
        etLooseUp  = tLooseUp * eLooseUp
        etLooseDown= tLooseDown *eLooseDown

        frweight = {
            'tLoose'     : tLoose    ,
            'tLoose/Up'   : tLooseUp  ,
            'tLoose/Down' : tLooseDown,
            'eLoose'      : eLoose    ,
            'eLoose/Up'   : eLooseUp  ,
            'eLoose/Down' : eLooseDown,
            'etLoose'      : etLoose    ,
            'etLoose/Up'   : etLooseUp  ,
            'etLoose/Down' : etLooseDown,
            'tLooseUnweight' : 1.,
        }

        return  frweight;

    def fill_histos(self, folder_str, row, weight, filter_label = ''):
        '''fills histograms'''
        #find all keys matching
        for attr, value in self.histo_locations[folder_str].iteritems():
            name = attr
            #if attr=='DEBUG':
            #    set_trace()
            if filter_label:
                if not attr.startswith(filter_label+'$'):
                    continue
                attr = attr.replace(filter_label+'$', '')
            if value.InheritsFrom('TH2'):
                if attr in self.hfunc:
                    try:
                        result, out_weight = self.hfunc[attr](row, weight)
                    except Exception as e:
                        raise RuntimeError("Error running function %s. Error: \n\n %s" % (attr, str(e)))
                    r1, r2 = result
                    if out_weight is None:
                        value.Fill( r1, r2 ) #saves you when filling NTuples!
                    else:
                        value.Fill( r1, r2, out_weight )
                else:
                    attr1, attr2 = split(attr)
                    v1 = getattr(row,attr1)
                    v2 = getattr(row,attr2)
                    value.Fill( v1, v2, weight ) if weight is not None else value.Fill( v1, v2 )
            else:
                if attr in self.hfunc:
                    try:
                        result, out_weight = self.hfunc[attr](row, weight)
                    except Exception as e:
                        raise RuntimeError("Error running function %s. Error: \n\n %s" % (attr, str(e)))
                    if out_weight is None:
                        value.Fill( result ) #saves you when filling NTuples!
                    else:
                        value.Fill( result, out_weight )
                else:
                    value.Fill( getattr(row,attr), weight ) if weight is not None else value.Fill( getattr(row,attr) )
        return None
    
    def process(self):
        logging.debug('Starting processing')
        systematics = self.systematics
        
        frw = []
        lock =()
        ievt = 0
        logging.debug('Starting evt loop')

        #pre-compute static things
        sys_shifts = systematics['trig'] + \
            systematics['pu'] + \
            systematics['eid'] + \
            systematics['eiso'] #+ \ 

        jes_dirs = [i.strip('_') for i in systematics['jes']]
        #anyway the first is always ""
        tes_dirs = [i.strip('_') for i in systematics['tes']][1:]
        ees_dirs = [i.strip('_') for i in systematics['ees']][1:]

        for row in self.tree:
            if (ievt % 100) == 0:
                logging.debug('New event')
            ievt += 1
            #avoid double counting events!
            evt_id = (row.run, row.lumi, row.evt)
            if evt_id == lock: continue
            if lock != () and evt_id == lock:
                logging.info('Removing duplicate of event: %d %d %d' % evt_id)

            #
            #preselection, common to everything and everyone
            #
            #trigger
            if self.is_embedded :
                if not bool(row.doubleMuPass) : continue
            else: 
                if not bool(row.singleE27WP80Pass) : continue
                if not bool(row.eMatchesSingleE27WP80): continue

            #objects
            if not selections.eSelection(row, 'e'): continue
            #if row.ePt < 30 : continue
            if not selections.tauSelection(row, 't'): continue
            if not row.tAntiElectronMVA5Tight : continue
            if not row.tAntiMuon2Loose : continue
            if not row.tLooseIso3Hits : continue
            logging.debug('object selection passed')
            #e ID/ISO
            if not selections.lepton_id_iso(row, 'e', 'eid13Tight_idiso02'): continue
            logging.debug('Passed preselection')

            #
            # Compute event weight
            #
            #event weight
            #sys_shifts is precomputed

            #set_trace()
            weight_map = self.event_weight(row, sys_shifts)

            #Fill embedded sample normalization BEFORE the vetoes
            if not row.e_t_SS:
                self.fill_histos('os/gg/ept30', row, weight_map[''])

            # it is better vetoing on b-jets  after the histo for the DY embedded
            #bjet veto
            if row.bjetCSVVeto30!=0 : continue

            #tau ID, id Tau is tight then go in full selection, otherwise use for fakes
            isTauTight = bool(row.tTightIso3Hits)
            isETight = bool(selections.lepton_id_iso(row, 'e', 'eid13Tight_etauiso01'))
            etau_category = ['']
            if (not isETight) and (not isTauTight):
                etau_category = ['etLoose', 'etLoose/Up', 'etLoose/Down']
            elif (not isTauTight):
                etau_category = ['tLoose', 'tLoose/Up', 'tLoose/Down', 'tLooseUnweight']
            elif (not isETight):
                etau_category = ['eLoose', 'eLoose/Up', 'eLoose/Down']

            #jet category
            central = struct(
                njets = min(row.jetVeto30, 3),
                tPt = row.tPt,
                ePt = row.ePt
                )
            jets = [min(row.jetVeto30, 3), min(row.jetVeto30jes_plus, 3),  min(row.jetVeto30jes_minus, 3)]
            tpts = [row.tPt_tes_plus, row.tPt_tes_minus]
            epts = [row.ePt_ees_plus, row.ePt_ees_minus]
            sys_effects = [(name, central.clone(njets = jnum)) for name, jnum in zip(jes_dirs, jets)]
            sys_effects.extend(
                [(name, central.clone(tPt = pt)) for name, pt in zip(tes_dirs, tpts)]
                )
            sys_effects.extend(
                [(name, central.clone(ePt = pt)) for name, pt in zip(ees_dirs, epts)]
                )

            #
            # Full tight selection
            #
            passes_full_selection = False
            selection_categories = []
            for name, shifted in sys_effects:
                #preselection. flat pt values
                if shifted.ePt > 30 and shifted.tPt > 30:
                    selection_categories.append((name, '%i' % shifted.njets, ''))
                else:
                    continue

                if shifted.njets == 0 :
                    selection_categories.extend([
                        (name, '0', i) for i in optimizer.compute_regions_0jet(
                            shifted.tPt, shifted.ePt, deltaPhi(row.ePhi, row.tPhi),
                            row.tMtToPfMet)
                        ])
                    
                    if shifted.tPt < 30: continue #was 35
                    if shifted.ePt < 50: continue #was 40
                    if deltaPhi(row.ePhi, row.tPhi) < 2.5 : continue #was 2.7
                    if row.tMtToPfMet > 70 : continue #was 50
                    selection_categories.append((name, '0', 'selected'))
                    passes_full_selection = True 
                elif shifted.njets == 1 :
                    selection_categories.extend([
                        (name, '1', i) for i in optimizer.compute_regions_1jet(
                            shifted.tPt, shifted.ePt, row.tMtToPfMet)
                        ])

                    if shifted.tPt < 30: continue #was 40
                    if shifted.ePt < 40: continue #was 35
                    if row.tMtToPfMet > 45 : continue #was 35
                    selection_categories.append((name, '1', 'selected'))
                    passes_full_selection = True 
                elif shifted.njets == 2 :
                    selection_categories.extend([
                        (name, '2', i) for i in optimizer.compute_regions_2jet(
                            shifted.tPt, shifted.ePt, row.tMtToPfMet, row.vbfMass,
                            row.vbfDeta)
                        ])

                    if shifted.tPt < 30: continue #was 40
                    if shifted.ePt < 40: continue #was 30
                    if row.tMtToPfMet > 55 : continue #was 35
                    if row.vbfMass < 500 : continue #was 550
                    if row.vbfDeta < 2.5 : continue #was 3.5
                    selection_categories.append((name, '2', 'selected'))
                    passes_full_selection = True 

            if passes_full_selection:
                logging.debug('Passed full selection')
                
            #
            #different selections
            #
            sign = 'ss' if row.e_t_SS else 'os'
            processtype ='gg'
            ptthreshold = ['ept30']

            #
            # Lepton vetoes
            #
            tvetoes = [row.tauVetoPt20EleTight3MuLoose, row.tauVetoPt20EleTight3MuLoose_tes_plus, row.tauVetoPt20EleTight3MuLoose_tes_minus]
            mvetoes = [row.muVetoPt5IsoIdVtx          , row.muVetoPt5IsoIdVtx_mes_plus          , row.muVetoPt5IsoIdVtx_mes_minus          ]
            evetoes = [row.eVetoCicLooseIso           , row.eVetoCicLooseIso_ees_plus           , row.eVetoCicLooseIso_ees_minus           ]
            
            tdirs = [ i for i, j in zip( systematics['tvetos'], tvetoes) if not j]
            mdirs = [ i for i, j in zip( systematics['mvetos'], mvetoes) if not j]
            edirs = [ i for i, j in zip( systematics['evetos'], evetoes) if not j]

            #if any of the lists is empty
            #set_trace()
            if not tdirs or not mdirs or not edirs:
                continue
            logging.debug('Passed Vetoes')

            #make all possible veto combos...
            all_dirs = [''.join(i) for i in itertools.product(tdirs, mdirs, edirs)]
            #...and choose only the meaningful ones
            veto_sys = set(systematics['tvetos']+systematics['mvetos']+systematics['evetos'])
            all_dirs = [i for i in all_dirs if i in veto_sys]
                        
            sys_directories = all_dirs + sys_shifts 
            #remove duplicates
            sys_directories = list(set(sys_directories))
            
            #at least one loose object
            if (not isETight) or (not isTauTight):
                #if is a loose tau just compute the fakes!                            
                sys_directories = etau_category
                
                #gather the one and only weight we do care about
                mc_weight = weight_map['']

                #weights are the fr ones...
                weight_map = self.fakerate_weights(row.tEta)
                for i in weight_map:
                    #...times the mc weight (if any)
                    weight_map[i] *= mc_weight


            #Fill histograms in appropriate direcotries
            #if passes_full_selection:
            #dirs = [os.path.join(sys, sign, processtype, e_thr, jet_dir) for sys, e_thr, jet_dir in itertools.product(sys_directories, ptthreshold, jet_directories)]
            #if len(dirs) <> len(set(dirs)):
            #    set_trace()
            for sys, e_thr, selection in itertools.product(sys_directories, ptthreshold, selection_categories):
                selection_sys, jet_dir, selection_step = selection
                #if we have multiple systematic shifts applied
                #reject the combination
                if selection_sys and sys:
                    continue

                #if we fill a histogram, lock the event
                lock = evt_id
                dir_name = os.path.join(sys, selection_sys, sign, processtype, 
                                        e_thr, jet_dir, selection_step)
                if dir_name[-1] == '/':
                    dir_name = dir_name[:-1]
                if passes_full_selection:
                    logging.debug('Filling %s' % dir_name)
                #fill them!
                weight_to_use = weight_map[sys] if sys in weight_map else weight_map['']
                self.fill_histos(dir_name, row, weight_to_use)
             
            
    def finish(self):
        self.write_histos()


