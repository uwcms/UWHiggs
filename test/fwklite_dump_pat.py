#! /bin/env python

import ROOT
import pprint
from DataFormats.FWLite import Events, Handle
ROOT.gROOT.SetBatch()


events = Events('/hdfs/store/user/tapas/2012-07-27-7TeV-PatTuple/data_DoubleMu_Run2011B_PromptReco_v1/1/patTuple_cfg-B21BFB50-B0F3-E011-917B-BCAEC5329716.root')
handle = Handle('pat::TriggerEvent')
label  = "patTriggerEvent"

evt = events.__iter__().next()
evt.getByLabel(label, handle)
obj = handle.product()
algo_names = [i.name() for i in obj.paths()]
print filter(lambda x: 'Mu8_Ele17' in x, algo_names)
        
