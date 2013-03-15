'''

Make control plots of Z->mumu events.

Author: Evan K. Friis, UW

'''

import EETree
from FinalStateAnalysis.PlotTools.MegaBase import MegaBase
import glob
import os
import mcCorrectors
#from FinalStateAnalysis.PlotTools.decorators import decorator
import fakerate_functions as frfits
import baseSelections as selections
import ROOT
from functools import wraps

math = ROOT.TMath

def sc_inv_mass(row):
    pt1 = row.e1SCEnergy / math.CosH(row.e1SCEta)
    pt2 = row.e2SCEnergy / math.CosH(row.e2SCEta)
    return math.Sqrt(2*pt1*pt2*( math.CosH(row.e1SCEta - row.e2SCEta) - math.Cos(row.e1SCPhi - row.e2SCPhi) ) )

#@decorator
def make_weight(fcn):
    #@wraps(fcn)
    def f_(*args):
        prob = fcn(*args)
        if prob:
            return prob / (1-prob)
        else:
            raise Exception("error in calulating weight")
    return f_

#SUMMED PROBABILITIES
charge_flip       = make_weight(frfits.e_charge_flip      )
charge_flip_up    = make_weight(frfits.e_charge_flip_up   )
charge_flip_down  = make_weight(frfits.e_charge_flip_down )

highpt_e_charge_flip      = make_weight( frfits.highpt_e_charge_flip      )
highpt_e_charge_flip_up   = make_weight( frfits.highpt_e_charge_flip_up   )
highpt_e_charge_flip_down = make_weight( frfits.highpt_e_charge_flip_down )

lowpt_e_charge_flip       = make_weight( frfits.lowpt_e_charge_flip       )
lowpt_e_charge_flip_up    = make_weight( frfits.lowpt_e_charge_flip_up    )
lowpt_e_charge_flip_down  = make_weight( frfits.lowpt_e_charge_flip_down  )

leading_e_fr_wjets    = make_weight(frfits.highpt_ee_fr )
leading_e_fr_qcd      = make_weight(frfits.highpt_ee_qcd_fr )
subleading_e_fr_wjets = make_weight(frfits.lowpt_ee_fr )
subleading_e_fr_qcd   = make_weight(frfits.lowpt_ee_qcd_fr )

def assign_charge_weight_two_maps(dir_, row):
    if "charge_weightSysUp"  in dir_: return highpt_e_charge_flip_up(row.e1AbsEta,row.e1Pt) + lowpt_e_charge_flip_up(row.e2AbsEta,row.e2Pt)
    if "charge_weightSysDwn" in dir_: return highpt_e_charge_flip_down(row.e1AbsEta,row.e1Pt) + lowpt_e_charge_flip_down(row.e2AbsEta,row.e2Pt)
    if "charge_weight"       in dir_:
        return highpt_e_charge_flip(row.e1AbsEta,row.e1Pt) + lowpt_e_charge_flip(row.e2AbsEta,row.e2Pt)
    return 1 #No Charge w to be applied!

def assign_charge_weight_one_maps(dir_, row):
    if "charge_weightSysUp"  in dir_: return charge_flip_up(row.e1AbsEta,row.e1Pt)   + charge_flip_up(row.e2AbsEta,row.e2Pt)
    if "charge_weightSysDwn" in dir_: return charge_flip_down(row.e1AbsEta,row.e1Pt) + charge_flip_down(row.e2AbsEta,row.e2Pt)
    if "charge_weight"       in dir_:
        return charge_flip(row.e1AbsEta,row.e1Pt) + charge_flip(row.e2AbsEta,row.e2Pt)
    return 1 #No Charge w to be applied!

assign_charge_weight = assign_charge_weight_two_maps

def assign_id_weight(dir_, row):
    if 'wjet_w' not in dir_ and 'qcd_w' not in dir_: #NO Weights to be applied!
        return 1
    leading_fr = leading_e_fr_wjets if 'wjet_w' in dir_ else leading_e_fr_qcd
    sublead_fr = subleading_e_fr_wjets if 'wjet_w' in dir_ else subleading_e_fr_qcd
    if "f1f2"  in dir_: return leading_fr( max(row.e1JetPt, row.e1Pt) ) * sublead_fr( max(row.e2JetPt, row.e2Pt) )
    if "f2"    in dir_: return sublead_fr( max(row.e2JetPt, row.e2Pt) )
    if "f1"    in dir_: return leading_fr( max(row.e1JetPt, row.e1Pt) )

class ControlZEE(MegaBase):
    tree = 'ee/final/Ntuple'
    def __init__(self, tree, outfile, **kwargs):
        super(ControlZEE, self).__init__(tree, outfile, **kwargs)
        # Use the cython wrapper
        self.tree = EETree.EETree(tree)
        self.out = outfile
        # Histograms for each category
        self.dir_based_histograms = {}
        self.is7TeV = '7TeV' in os.environ['jobid']
        self.pucorrector = mcCorrectors.make_puCorrector('doublee')

    def begin(self):
        def book_(dirname):
            self.book(dirname, 'ePt'     , 'e Pt', 100, 0, 100)
            self.book(dirname, 'eAbsEta' , 'e Abs Eta', 80, 0, 2.5)
            self.book(dirname, 'SCEnergy', 'electron Super Cluster energy', 500, 0, 1000)
            self.book(dirname, 'SCDPhi'  , 'electron Super Cluster DeltaPhi', 180, 0, math.Pi())
            self.book(dirname, 'TrkMass' , 'Dielectrons invariant mass; M_{ee} [GeV];counts', 110, 40, 150)
            self.book(dirname, 'TrkMass_NOSCALE' , 'Dielectrons invariant mass; M_{ee} [GeV];counts', 110, 40, 150)
            self.book(dirname, 'SCMass'  , 'Dielectrons Super Cluster invariant mass; M_{ee} [GeV];counts', 110, 40, 150)
            self.book(dirname, "e1Pt"    , "electron 1 Pt", 100, 0, 100)
            self.book(dirname, "e2Pt"    , "electron 2 Pt", 100, 0, 100)
            self.book(dirname, "e1AbsEta", "electron 1 abseta", 100, 0., 2.5)
            self.book(dirname, "e2AbsEta", "electron 2 abseta", 100, 0., 2.5)

        self.dirs = ['/'.join([sign,id_,weight,ch_weight]) for sign in ['os','ss']
                for id_ in [h+k for h in ['p1','f1'] for k in ['p2','f2'] ]
                for weight in ['','wjet_w','qcd_w']
                for ch_weight in ["","charge_weight","charge_weightSysUp","charge_weightSysDwn"]
                if 'f' in id_ or weight == ''
                if ch_weight == '' or sign == 'os'
            ] #builds all possible dirs o(s)s/p(f)1p(f)2(/wjet(qcd)_w)(/charge_weight)
        self.dirs = map(lambda x: x.replace('//','/') , self.dirs)
        self.dirs = map(lambda x: x[:-1] if x[-1] == '/' else x, self.dirs) #removes trailing / when needed
                                                                  # There are actually more dirs than what I need, but who cares?
        for dirname in self.dirs:
            book_(dirname)

        for key, hist in self.histograms.iteritems():
            charpos  = key.rfind('/')
            location = key[ : charpos]
            name     = key[ charpos + 1 :]
            if location not in self.dir_based_histograms:
                self.dir_based_histograms[location] = {}
            self.dir_based_histograms[location][name] = hist

    def evt_weight(self, row):
        if row.run > 2:
            return 1.
        else:
            return self.pucorrector(row.nTruePU) * \
                mcCorrectors.get_electron_corrections(row,'e1','e2')
        

    def preselection(self, row):
        ''' Preselection applied to events.

        Excludes FR object IDs and sign cut.
        '''
        if not row.doubleEPass: return False
        if row.e1Pt < 20: return False
        if not selections.eSelection(row, 'e1'): return False
        if not selections.eSelection(row, 'e2'): return False
        if not selections.vetos(row): return False
            #region = selections.control_region_ee(row)
        return True

    def obj1_id(self, row):
        return selections.h2tau_eid(row, 'e1')

    def obj2_id(self, row):
        return selections.h2tau_eid(row, 'e2')

    def process(self):

        histos = self.dir_based_histograms
        def fill_histos(dirname, row, weight):
            mass = frfits.mass_scaler(row.e1_e2_Mass) if "charge_weight" in dirname else row.e1_e2_Mass
            histos[dirname]['ePt'     ].Fill(row.e1Pt,weight)
            histos[dirname]['eAbsEta' ].Fill(row.e1AbsEta,weight)
            histos[dirname]['SCEnergy'].Fill(row.e1SCEnergy,weight)
            histos[dirname]['ePt'     ].Fill(row.e2Pt,weight)
            histos[dirname]['eAbsEta' ].Fill(row.e2AbsEta,weight)
            histos[dirname]['SCEnergy'].Fill(row.e2SCEnergy,weight)
            histos[dirname]['TrkMass' ].Fill(mass, weight)
            #histos[dirname]['SCMass'  ].Fill(sc_inv_mass(row),weight)
            histos[dirname]["e1Pt"    ].Fill(row.e1Pt,weight)
            histos[dirname]["e2Pt"    ].Fill(row.e2Pt,weight)
            histos[dirname]["e1AbsEta"].Fill(row.e1AbsEta,weight)
            histos[dirname]["e2AbsEta"].Fill(row.e2AbsEta,weight)
            histos[dirname]['TrkMass_NOSCALE'].Fill(row.e1_e2_Mass, weight)

        for row in self.tree:
            if not self.preselection(row):
                continue
            pair_sign = 'ss' if bool(row.e1_e2_SS) else 'os'
            id1       = 'p1' if self.obj1_id(row)  else 'f1'
            id2       = 'p2' if self.obj2_id(row)  else 'f2'
            to_match  = pair_sign+'/'+id1+id2
            matchingd = (i for i in self.dirs if i.startswith(to_match)) # filter(lambda x: x.startswith(pair_sign+'/'+id1+id2), self.dirs)
            for dir_ in matchingd:
                weight = 1
                weight *= assign_id_weight(dir_, row)
                weight *= assign_charge_weight(dir_, row)
                fill_histos(dir_, row, weight)

    def finish(self):
        self.write_histos()
