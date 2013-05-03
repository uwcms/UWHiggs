'small interface module to deal with optimizization'

import os
from baseSelections import *

leading_lepton_iso_tag    = os.environ['leadLeptonIsoTag'] if 'leadLeptonIsoTag' in os.environ else 'h2taucuts'
subleading_lepton_iso_tag = os.environ['subleadLeptonIsoTag'] if 'subleadLeptonIsoTag' in os.environ else 'h2taucuts'
lt_lower_threshold        = float(os.environ['ltThreshold']) if 'ltThreshold' in os.environ else 80.

leading_lepton_id_iso    = lepton_ids[leading_lepton_iso_tag]
subleading_lepton_id_iso = lepton_ids[leading_lepton_iso_tag]
