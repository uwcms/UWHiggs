'''

Measure fake rates in dielectron events.

NB that the triggers we are using always have the same iso on the leading
and subleading legs.

We measure W+jet (iso electron) control regions.

The layout of output is:

    region/denom_tag/var1
    region/denom_tag/var2
    region/denom_tag/num_tag/var1
    region/denom_tag/num_tag/var2

Author: Evan K. Friis, UW

'''

import EETree
from FinalStateAnalysis.PlotTools.MegaBase import MegaBase
import baseSelections as selections
import math
import os
import ROOT
import array

barrelThr = 1.48
accThr    = 2.5
binsEndcap = 5
binsBarrel = 4
ptbins     = [10.,30.,50.,130.]

class ChargeFlipProbabilityEE(MegaBase):
    tree = 'ee/final/Ntuple'
    def __init__(self, tree, outfile, **kwargs):
        super(ChargeFlipProbabilityEE, self).__init__(tree, outfile, **kwargs)
        # Use the cython wrapper
        self.tree = EETree.EETree(tree)
        self.out = outfile
        # Histograms for each category
        self.histograms = {}
        self.is7TeV = '7TeV' in os.environ['jobid']
        self.charge_ids = ['ChargeIdTight']

    def begin(self):
        # Charge mis-ID measurements
        for charge_id in self.charge_ids:
            self.book( charge_id, 'os_trkMass', '', 130,0.,130.)
            for electron in ['leading_electron','subleading_electron']:
                etabins    = [(i*barrelThr/binsBarrel) for i in range(binsBarrel+1)]+[ (barrelThr+i*(accThr-barrelThr)/binsEndcap) for i in range(1,binsEndcap+1)]
                folder     = '/'.join([charge_id,electron])
                self.book( folder, 'flipped_electrons', 'Flipped electrons distribution; |#eta|; p_{T} [GeV]',     len(etabins)-1, array.array('f',etabins), len(ptbins) -1, array.array('f',ptbins), type=ROOT.TH2F)
                self.book( folder, 'matched_electrons', 'Gen-matched electrons distribution; |#eta|; p_{T} [GeV]', len(etabins)-1, array.array('f',etabins), len(ptbins) -1, array.array('f',ptbins), type=ROOT.TH2F)
                self.book( folder, 'pt_in_barrel', 'Gen-matched pt; p_{T} [GeV]', 130,0.,130.)
                self.book( folder, 'pt_in_endcap', 'Gen-matched pt; p_{T} [GeV]', 130,0.,130.)
                self.book( folder, 'ss_trkMass', '', 130,0.,130.)
        ## self.book('charge', 'trkMass_endcap', '', 130,0.,130.)
        ## self.book('charge', 'trkMass_barrel', '', 130,0.,130.)
            
        #print self.histograms.keys()

    def process(self):

        def preselection(row):
            if not row.doubleEPass: return False
            if not row.e1Pt > 20: return False
            if not selections.eSelection(row, 'e1'): return False #already applies charge id tight, FIXME if needed
            if not selections.eSelection(row, 'e2'): return False #already applies charge id tight, FIXME if needed
            if not selections.vetos(row): return False
            return bool(selections.control_region_ee(row) == 'zee')

        #def fill(the_histos, row):

        histos = self.histograms
        for row in self.tree:
            if not preselection(row):
                continue

            for charge_id in self.charge_ids:
                #No need to check if the electrons pass the id, it is already applied in preselection (we currently have only one!)
                if row.e1GenPdgId != -999 and row.e2GenPdgId != -999 and not row.e1_e2_SS:
                    histos[charge_id+'/os_trkMass'].Fill(row.e1_e2_Mass)

                for el, electron in zip(['e1','e2'],['leading_electron','subleading_electron']):
                    folder     = '/'.join([charge_id,electron,''])
                    pdgid = getattr(row,el+'GenPdgId')
                    if pdgid != -999:
                        eta = getattr(row,el+'AbsEta')
                        pt  = getattr(row,el+'Pt')
                        if eta < barrelThr:
                            histos[folder+'pt_in_barrel'].Fill(pt)
                        elif eta < accThr:
                            histos[folder+'pt_in_endcap'].Fill(pt)
                        histos[folder+'matched_electrons'].Fill( eta, pt )
                        if getattr(row,el+'Charge')*(11) == pdgid: ##e- --> -11; e+ --> 11 MISMEASURED ELECTRON 
                            histos[folder+'flipped_electrons'].Fill( eta, pt )
                            if row.e1GenPdgId != -999 and row.e2GenPdgId != -999 and row.e1_e2_SS: #one of the two is certainly tru, but is faster this way
                                histos[folder+'ss_trkMass'].Fill(row.e1_e2_Mass)

    def finish(self):
        self.write_histos()


 
 
