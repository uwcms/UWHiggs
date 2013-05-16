'''

Generic base class for WH leptonic analysis.

WH has three objects:
    obj1
    obj2
    obj3

Objects 1 & 2 must be SS in the signal region.   The OS region is kept as a
control.  The subclasses must define the following functions:

    self.preselection - includes loose preselection on objects
    slef.sign_cut - returns true for SS (signal), false for OS
    self.obj1_id(row)
    self.obj2_id(row)
    self.obj3_id(row)

    self.event_weight(row) - general corrections
    self.obj1_weight(row) - returns double
    self.obj2_weight(row)
    self.obj3_weight(row)

Additionally, the fake rate weight in QCD events must be defined.

    self.obj1_qcd_weight(row) - returns double
    self.obj2_qcd_weight(row)
    self.obj3_qcd_weight(row)

    self.book_histos(folder) # books histograms in a given folder (region)
    self.fill_histos(histos, row, weight) # fill histograms in a given region

The output histogram has the following structure:

    ss/
        p1p2p3/   - signal region
        == Type1 FRs ==
        p1p2f3/   - object 3 fails
            w3/   - with weight 3 applied
            q3/   - with weight 3 (from QCD events) applied
        ... same for 1 and 2
        == Type2 FRs ==
        f1f2p3/   - objects 1 & 2 fail
            w1w2/
            q1q2/
        == Type3 FRs ==
        f1f2f3/   - everything fails
            w3/   - extrapolate triple fakes to type2 region
            q3/

    os/           - OS control region

'''

from FinalStateAnalysis.PlotTools.MegaBase import MegaBase
import optimizer
import os
import ROOT
import math
#from pudb import set_trace; set_trace()
def quad(*xs):
    return math.sqrt(sum(x * x for x in xs))

def inv_mass(*args):
    lorentz_vecs = [ROOT.TLorentzVector() for i in xrange(len(args))]
    for v,i in zip(lorentz_vecs,args):
        v.SetPtEtaPhiM(*i)
    return sum(lorentz_vecs,ROOT.TLorentzVector()).M() #We need to give the staring point otherwise it starts from an int and it does not work

class cut_flow_tracker(object):
    def __init__(self, hist):
        self.labels   = [hist.GetXaxis().GetBinLabel(i+1) for i in range(hist.GetNbinsX())]
        self.cut_flow = dict([ (i, False) for i in self.labels])
        self.hist     = hist
        self.evt_info = [-1, -1, -1]
        self.disabled = 'CutFlow' not in os.environ
    
    def fill(self, label):
        self.cut_flow[label] = True

    def Fill(self, *labels):
        if self.disabled:
            return
        for label in labels:
            self.fill(label)

    def flush(self):
        if self.disabled:
            return
        for i, label in enumerate(self.labels):
            val = self.cut_flow[label]
            if val:
                self.hist.Fill(i+0.5)

    def new_row(self, *args):
        if self.disabled:
            return
        if self.evt_info != list(args):
            self.flush()
            self.evt_info = list(args)
            self.cut_flow = dict([ (i, False) for i in self.labels])
                        
class WHAnalyzerBase(MegaBase):
    def __init__(self, tree, outfile, wrapper, **kwargs):
        super(WHAnalyzerBase, self).__init__(tree, outfile, **kwargs)
        # Cython wrapper class must be passed
        self.tree = wrapper(tree)
        self.out = outfile
        self.histograms = {}
        self.histo_locations = {} #just a mapping of the histograms we have to avoid changing self.histograms indexing an screw other files
        self.hfunc   = { #maps the name of non-trivial histograms to a function to get the proper value, the function MUST have two args (evt and weight). Used in fill_histos later
            'nTruePU' : lambda row, weight: (row.nTruePU,None),
            'weight'  : lambda row, weight: (weight,None) if weight is not None else (1.,None),
            'Event_ID': lambda row, weight: (array.array("f", [row.run,row.lumi,int(row.evt)/10**5,int(row.evt)%10**5] ), None),
            }

    @staticmethod
    def build_wh_folder_structure():
        # Build list of folders, and a mapping of
        # (sign, ob1, obj2, ...) => ('the/path', (weights))
        #folders = []
        flag_map = {}
        for sign in ['ss', 'os']:
            for failing_objs in [(), (1,), (2,), (3,), (1,3), (2, 3), (1,2), (1,2,3)]:
                cut_key = [sign == 'ss']
                region_label = ''
                for i in range(1,4):
                    if i in failing_objs:
                        region_label += 'f' + str(i)
                        cut_key.append(False)
                    else:
                        region_label += 'p' + str(i)
                        cut_key.append(True)
                # Figure out which objects to weight for FR
                weights_to_apply = []
                # Single fake
                if len(failing_objs) == 1:
                    weights_to_apply.append(
                        (failing_objs, "w%i" % failing_objs))
                    # A version using the QCD fake rate
                    weights_to_apply.append(
                        (failing_objs, "q%i" % failing_objs))
                if len(failing_objs) == 2:
                    # in the 1-2 case, apply both.  Otherwise, just apply the
                    # first (a light lepton)
                    if 3 not in failing_objs:
                        weights_to_apply.append(
                            (failing_objs, "w%i%i" % failing_objs))
                        # Using QCD rate
                        weights_to_apply.append(
                            (failing_objs, "q%i%i" % failing_objs))
                    else:
                        weights_to_apply.append(
                            (failing_objs, "w%i" % failing_objs[0]))
                        # Using QCD rate
                        weights_to_apply.append(
                            (failing_objs, "q%i" % failing_objs[0]))

                if len(failing_objs) == 3:
                    weights_to_apply.append( ((3,), "w3") )
                    weights_to_apply.append( ((1,3,), "w13") )
                    weights_to_apply.append( ((2,3,), "w23") )
                    # QCD weight versions
                    weights_to_apply.append( ((1,3,), "q13") )
                    weights_to_apply.append( ((2,3,), "q23") )
                    # Needed for f3 CR
                    weights_to_apply.append( ((1,2), "w12"))
                    weights_to_apply.append( ((1,), "w1"))
                    weights_to_apply.append( ((2,), "w2"))
                    weights_to_apply.append( ((1,2), "q12"))
                    weights_to_apply.append( ((1,), "q1"))
                    weights_to_apply.append( ((2,), "q2"))

                #folders_to_add = [ (sign, region_label) ]
                # Which objects to weight for each region
                weights_to_add = []
                for failing_objs, weight_to_apply in weights_to_apply:
                    #folders_to_add.append( (sign, region_label, weight_to_apply) )
                    weights_to_add.append(weight_to_apply)

                flag_map[tuple(cut_key)] = ((sign, region_label), tuple(weights_to_add))
                #folders.extend(folders_to_add)

        return flag_map

    def fill_histos(self, histos, folder, row, weight, filter_label = ''):
        '''fills histograms'''
        #find all keys matching
        folder_str = '/'.join(folder + ('',))
        for attr in self.histo_locations[folder_str]:
            name = attr
            if filter_label:
                if not attr.startswith(filter_label+'$'):
                    continue
                attr = attr.replace(filter_label+'$', '')
            value = self.histograms[folder_str+name]
            if value.InheritsFrom('TH2'):
                if attr in self.hfunc:
                    result, out_weight = self.hfunc[attr](row, weight)
                    r1, r2 = result
                    if out_weight is None:
                        value.Fill( r1, r2 ) #saves you when filling NTuples!
                    else:
                        value.Fill( r1, r2, out_weight )
                else:
                    attr1, attr2 = tuple(attr.split('#'))
                    v1 = getattr(row,attr1)
                    v2 = getattr(row,attr2)
                    value.Fill( v1, v2, weight ) if weight is not None else value.Fill( v1, v2 )
            else:
                if attr in self.hfunc:
                    result, out_weight = self.hfunc[attr](row, weight)
                    if out_weight is None:
                        value.Fill( result ) #saves you when filling NTuples!
                    else:
                        value.Fill( result, out_weight )
                else:
                    value.Fill( getattr(row,attr), weight ) if weight is not None else value.Fill( getattr(row,attr) )
        return None

    def begin(self):
        # Loop over regions, book histograms
        for _, folders in self.build_wh_folder_structure().iteritems():
            base_folder, weight_folders = folders
            folder = "/".join(base_folder)
            self.book_histos(folder) # in subclass
            # Each of the weight subfolders
            for weight_folder in weight_folders:
                self.book_histos("/".join(base_folder + (weight_folder,)))

        # Add WZ control region
        self.book_histos('ss/p1p2p3_enhance_wz')
        # Where second light lepton fails
        self.book_histos('ss/p1f2p3_enhance_wz')
        self.book_histos('ss/p1f2p3_enhance_wz/w2')

        # Add charge-fake control region - probability that obj1 will flip into
        # ss/p1p2p3
        self.book_histos('os/p1p2p3/c1')
        self.book_histos('os/p1p2f3/c1')
        self.book_histos('os/p1p2p3/c2')
        self.book_histos('os/p1p2f3/c2')
        self.book_histos('os/p1p2p3/c1_sysup')
        self.book_histos('os/p1p2f3/c1_sysup')
        self.book_histos('os/p1p2p3/c2_sysup')
        self.book_histos('os/p1p2f3/c2_sysup')
        for key in self.histograms:
            charpos  = key.rfind('/')
            location = key[ : charpos]+'/'
            name     = key[ charpos + 1 :]
            if location in self.histo_locations:
                self.histo_locations[location].append(name)
            else:
                self.histo_locations[location] = [name]

        #Makes the cut flow histogram
        cut_flow_step = ['bare', 'WH Event',
                         #'obj1 GenMatching', 'obj2 GenMatching', 'obj3 GenMatching',
                         'trigger',
                         'obj1 Presel', 'obj2 Presel',
                         ## 'pt requirements', 'eta requirements',
                         ## 'MissingHits', 'HasConversion', 'JetBtag', 'ChargeIdTight', 'DZ',
                         'obj3 Presel',
                         ## 'obj1 IDIso', 'obj2 IDIso', 'obj3 IDIso',
                         'LT', 'vetos', 'charge_fakes', 'sign cut', 'anti WZ', 
                         ]
        self.book('ss', "CUT_FLOW", "Cut Flow", len(cut_flow_step), 0, len(cut_flow_step))
        xaxis = self.histograms['ss/CUT_FLOW'].GetXaxis()
        self.cut_flow_histo = self.histograms['ss/CUT_FLOW']
        self.cut_flow_map   = {}
        for i, name in enumerate(cut_flow_step):
            xaxis.SetBinLabel(i+1, name)
            self.cut_flow_map[name] = i+0.5

    def process(self):
        # For speed, map the result of the region cuts to a folder path
        # string using a dictionary
        # key = (sign, obj1, obj2, obj3)
        cut_region_map = self.build_wh_folder_structure()

        # Reduce number of self lookups and get the derived functions here
        histos = self.histograms
        preselection = self.preselection
        sign_cut = self.sign_cut
        obj1_id = self.obj1_id
        obj2_id = self.obj2_id
        obj3_id = self.obj3_id
        fill_histos = self.fill_histos
        anti_wz_cut = self.anti_wz
        
        weight_func = self.event_weight

        cut_flow_histo = self.cut_flow_histo
        cut_flow_trk   = cut_flow_tracker(cut_flow_histo)

        # Which weight folders correspond to which weight functions
        weight_map = {
            'w1' : (self.obj1_weight, ),
            'w2' : (self.obj2_weight, ),
            'w3' : (self.obj3_weight, ),
            'w12' : (self.obj1_weight, self.obj2_weight),
            'w13' : (self.obj1_weight, self.obj3_weight),
            'w23' : (self.obj2_weight, self.obj3_weight),

            'q1' : (self.obj1_qcd_weight, ),
            'q2' : (self.obj2_qcd_weight, ),
            'q3' : (self.obj3_qcd_weight, ),
            'q12' : (self.obj1_qcd_weight, self.obj2_qcd_weight),
            'q13' : (self.obj1_qcd_weight, self.obj3_qcd_weight),
            'q23' : (self.obj2_qcd_weight, self.obj3_qcd_weight),
        }

        for i, row in enumerate(self.tree):
            ## if i > 0:
            ##     raise

            for cut_label, cut_settings in optimizer.grid_search.iteritems():
                if 'WZJetsTo3LNu' in os.environ['megatarget']:
                    if 'ZToTauTau' in os.environ['megatarget']:
                        if not (row.isGtautau or row.isZtautau):
                            continue
                    else:
                        if (row.isGtautau or row.isZtautau):
                            continue
                        
                cut_flow_trk.new_row(row.run,row.lumi,row.evt)
                cut_flow_trk.Fill('bare')
                # Apply basic preselection
                if not cut_flow_trk.disabled:
                    if row.processID != 26:
                        continue

                cut_flow_trk.Fill('WH Event')
                #print cut_label, cut_settings
                if not preselection(row, cut_flow_trk, cut_settings['LT']):
                    continue
                # Get the generic event weight
                event_weight = weight_func(row)

                # Get the cuts that define the region
                sign_result = sign_cut(row)
                obj1_id_result = obj1_id(row, cut_settings['leading_iso'], cut_settings['subleading_iso'])
                obj2_id_result = obj2_id(row, cut_settings['leading_iso'], cut_settings['subleading_iso'])
                obj3_id_result = obj3_id(row)
                anti_wz = anti_wz_cut(row)

                if not cut_flow_trk.disabled:
                    if obj1_id_result:
                        cut_flow_trk.Fill('obj1 IDIso')
                        if obj2_id_result:
                            cut_flow_trk.Fill('obj2 IDIso')
                            if obj3_id_result:
                                cut_flow_trk.Fill('obj3 IDIso')
                                if sign_result:
                                    cut_flow_trk.Fill('sign cut')
                                    if anti_wz :
                                        cut_flow_trk.Fill('anti WZ')

                # Figure out which folder/region we are in
                region_result = cut_region_map.get(
                    (sign_result, obj1_id_result, obj2_id_result, obj3_id_result))

                # Ignore stupid regions we don't care about
                if region_result is None:
                    continue

                if anti_wz:
                    base_folder, weights = region_result
                    # Fill the un-fr-weighted histograms
                    fill_histos(histos, base_folder, row, event_weight, cut_label)

                    # Now loop over all necessary weighted regions and compute & fill
                    for weight_folder in weights:
                        # Compute product of all weights
                        fr_weight = event_weight
                        # Figure out which object weight functions to apply
                        for subweight in weight_map[weight_folder]:
                            fr = subweight(row, cut_settings['leading_iso'], cut_settings['subleading_iso'])
                            fr_weight *= fr/(1.-fr)

                        # Now fill the histos for this weight folder
                        fill_histos(histos, base_folder + (weight_folder,), row, fr_weight, cut_label)

                    if not sign_result and obj1_id_result and obj2_id_result:
                        # Object 1 can only flip if it is OS with the tau
                        obj1_obj3_SS     = self.obj1_obj3_SS(row)
                        if (obj1_obj3_SS and hasattr(self,'obj1_charge_flip')) \
                            or ( not obj1_obj3_SS and hasattr(self,'obj2_charge_flip')): #there is the function --> we have to compute it, otherwise skip and save some precious filling time!
                            charge_flip_prob = self.obj1_charge_flip(row)       if obj1_obj3_SS else self.obj2_charge_flip(row)
                            charge_flip_sysu = self.obj1_charge_flip_sysup(row) if obj1_obj3_SS else self.obj2_charge_flip_sysup(row)
                            directory        = 'os/p1p2%s3/c1' if obj1_obj3_SS else 'os/p1p2%s3/c2'
                            directory        = directory % ('p' if obj3_id_result else 'f')
                            charge_flip_prob = charge_flip_prob/(1. - charge_flip_prob)
                            fill_histos(histos, (directory,), row, event_weight*charge_flip_prob, cut_label)
                            fill_histos(histos, (directory+'_sysup',), row, event_weight*charge_flip_sysu, cut_label)

                elif sign_result and obj1_id_result and obj3_id_result:
                    # WZ object topology fails. Check if we are in signal region.
                    if self.enhance_wz(row):
                        # Signal region
                        if obj2_id_result:
                            fill_histos(histos, ('ss/p1p2p3_enhance_wz',), row, event_weight, cut_label)
                        else:
                            fill_histos(histos, ('ss/p1f2p3_enhance_wz',), row, event_weight, cut_label)
                            fake_rate_obj2 = self.obj2_weight(row)
                            fake_weight = fake_rate_obj2/(1.-fake_rate_obj2)
                            fill_histos(histos, ('ss/p1f2p3_enhance_wz/w2',), row, event_weight*fake_weight, cut_label)
        cut_flow_trk.flush()

    def finish(self):
        self.write_histos()

if __name__ == "__main__":
    import pprint
    pprint.pprint(WHAnalyzerBase.build_wh_folder_structure())
