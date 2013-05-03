'''

Measures probability for a gen electron to be reconstructed with different charge

Author: M. Verzetti, UZH

'''

import EETree
#from EETauTree import EETauTree
from FinalStateAnalysis.PlotTools.MegaBase import MegaBase
import optimizer as selections
import math
import os
import ROOT
import array
from FinalStateAnalysis.PlotTools.decorators import memo

@memo
def name_mapping(name, var):
    return name+var

barrelThr    = 1.48
accThr       = 2.5
binsEndcap   = 5
binsBarrel   = 4
ptbins       = [10.,30.,50.,130.]

class ChargeFlipProbabilityEE(MegaBase):
    tree = 'ee/final/Ntuple'
    #tree = 'eet/final/Ntuple'
    def __init__(self, tree, outfile, **kwargs):
        super(ChargeFlipProbabilityEE, self).__init__(tree, outfile, **kwargs)
        # Use the cython wrapper
        self.tree = EETree.EETree(tree)
        self.out = outfile
        # Histograms for each category
        self.histograms = {}
        self.is7TeV = '7TeV' in os.environ['jobid']

    def begin(self):
        # Charge mis-ID measurements
        etabins    = [(i*barrelThr/binsBarrel) for i in range(binsBarrel+1)] + \
            [ (barrelThr+i*(accThr-barrelThr)/binsEndcap) for i in range(1,binsEndcap+1)] 
        self.book('charge', 'flipped_electrons', 'Flipped electrons distribution; |#eta|; p_{T} [GeV]',
                  len(etabins)-1, array.array('f',etabins), len(ptbins) -1, array.array('f',ptbins), type=ROOT.TH2F)
        self.book('charge', 'matched_electrons', 'Gen-matched electrons distribution; |#eta|; p_{T} [GeV]',
                  len(etabins)-1, array.array('f',etabins), len(ptbins) -1, array.array('f',ptbins), type=ROOT.TH2F)
        self.book('charge', 'flipped_e1', 'Flipped electrons distribution; |#eta|; p_{T} [GeV]',
                  len(etabins)-1, array.array('f',etabins), len(ptbins) -1, array.array('f',ptbins), type=ROOT.TH2F)
        self.book('charge', 'matched_e1', 'Gen-matched electrons distribution; |#eta|; p_{T} [GeV]',
                  len(etabins)-1, array.array('f',etabins), len(ptbins) -1, array.array('f',ptbins), type=ROOT.TH2F)
        self.book('charge', 'flipped_e2', 'Flipped electrons distribution; |#eta|; p_{T} [GeV]',
                  len(etabins)-1, array.array('f',etabins), len(ptbins) -1, array.array('f',ptbins), type=ROOT.TH2F)
        self.book('charge', 'matched_e2', 'Gen-matched electrons distribution; |#eta|; p_{T} [GeV]',
                  len(etabins)-1, array.array('f',etabins), len(ptbins) -1, array.array('f',ptbins), type=ROOT.TH2F)
        self.book('charge', 'pt_in_barrel', 'Gen-matched pt; p_{T} [GeV]', 130,0.,130.)
        self.book('charge', 'pt_in_endcap', 'Gen-matched pt; p_{T} [GeV]', 130,0.,130.)
        self.book('charge', 'os_trkMass', '', 130,0.,130.)
        self.book('charge', 'ss_trkMass', '', 130,0.,130.)
        ## self.book('charge', 'trkMass_endcap', '', 130,0.,130.)
        ## self.book('charge', 'trkMass_barrel', '', 130,0.,130.)
            
        #print self.histograms.keys()

    def process(self):

        def preselection(row):
            if not row.doubleEPass:  return False
            #if (row.e1GenPdgId == -999) or (row.e2GenPdgId == -999): return False
            if row.e1Pt < 20:        return False
            if not selections.eSelection(row, 'e1'): return False
            if not selections.eSelection(row, 'e2'): return False
            if any([ row.muVetoPt5,
                      row.tauVetoPt20,
                      row.eVetoCicTightIso]):        return False
            if not selections.leading_lepton_id_iso(row, 'e1'):    return False
            if not selections.subleading_lepton_id_iso(row, 'e2'): return False
            if not (row.jetVeto40 >= 1):              return False
            return True
            ## return bool(selections.control_region_ee(row) == 'zee')

        #def fill(the_histos, row):

        histos = self.histograms
        for row in self.tree:
            if not preselection(row):
                continue

            if row.e1GenPdgId != -999 and row.e2GenPdgId != -999:
                if row.e1_e2_SS:
                    histos['charge/ss_trkMass'].Fill(row.e1_e2_Mass)
                else:
                    histos['charge/os_trkMass'].Fill(row.e1_e2_Mass)
                    
            for el in ['e1','e2']:
                pdgid = getattr(row, name_mapping(el,'GenPdgId'))
                if pdgid != -999:
                    eta    = getattr(row, name_mapping(el,'AbsEta') )
                    pt     = getattr(row, name_mapping(el,'Pt') )
                    charge = getattr(row, name_mapping(el,'Charge') )
                    if eta < barrelThr:
                        histos['charge/pt_in_barrel'].Fill(pt)
                    elif eta < accThr:
                        histos['charge/pt_in_endcap'].Fill(pt)
                    histos['charge/matched_electrons'].Fill( eta, pt )
                    single_hist = 'charge/matched_%s' % el
                    histos[single_hist].Fill( eta, pt )
                    if charge*(11) == pdgid: ##e- --> -11; e+ --> 11 MISMEASURED ELECTRON
                        single_hist = 'charge/flipped_%s' % el
                        histos[single_hist].Fill( eta, pt )
                        histos['charge/flipped_electrons'].Fill( eta, pt )

    def finish(self):
        self.write_histos()
