'''

Generic base class for ZH leptonic analysis.

ZH has three objects:
    obj1
    obj2
    obj3

Objects 1 & 2 must be SS in the signal region.   The OS region is kept as a
control.  The subclasses must define the following functions:

    self.preselection - includes loose preselection on objects and Z identification (MM or EE)
    slef.sign_cut - returns true for OS (signal), false for SS amo
    self.probe1_id
    self.probe2_id

    self.event_weight(row) - general corrections

    self.book_histos(folder) # books histograms in a given folder (region)
    self.fill_histos(histos, row, weight) # fill histograms in a given region

'''

from FinalStateAnalysis.PlotTools.MegaBase import MegaBase
#import pprint

class ZHAnalyzerBase(MegaBase):
    def __init__(self, tree, outfile, wrapper, channel, **kwargs):
        #print '__init__ called'
        super(ZHAnalyzerBase, self).__init__(tree, outfile, **kwargs)
        # Cython wrapper class must be passed
        self.tree = wrapper(tree)
        self.out = outfile
        self.histograms = {}
        self.channel = channel
        #print '__init__->self.channel %s' % self.channel
        
    ## def get_channel(self):
    ##     return 'LL'

    #@classmethod #staticmethod
    def build_zh_folder_structure(self):
        # Build list of folders, and a mapping of
        # (sign, ob1, obj2, ...) => ('the/path', (weights))
        #folders = []
        #print dir(self)
        #print 'folder->self.channel %s'% self.channel
        channel  = self.channel
        #print channel
        flag_map = {}
        for sign in ['ss', 'os']:
            for failing_objs in [(), (1,2), (1,), (2,),]:
                region_label = '_'.join(['%sIsoFailed' % self.H_decay_products()[obj-1] for obj in failing_objs]) if len(failing_objs) else 'All_Passed'
                weightsAvail = [None, self.obj1_weight, self.obj2_weight]
                flag_map[(sign,region_label)] = {
                    'selection' : {
                        self.sign_cut  : (sign == 'os'),
                        self.probe1_id : (not (1 in failing_objs) ),
                        self.probe2_id : (not (2 in failing_objs) ),
                        },
                    'weights'   : [weightsAvail[i] for i in failing_objs], 
                    }
        #print flag_map.keys()
        return flag_map

    def begin(self):
        # Loop over regions, book histograms
        for folders, regionInfo in self.build_zh_folder_structure().iteritems():
            folder = "/".join(folders)
            self.book_histos(folder) # in subclass
            # Each of the weight subfolders
            wToApply = regionInfo['weights']
            for w in wToApply:
                subf = "/".join([folder, w.__name__])
                self.book_histos(subf)
            if len(wToApply) > 1:
                subf = "/".join([folder, 'all_weights_applied'] )
                self.book_histos(subf)
            
    def process(self):
        # For speed, map the result of the region cuts to a folder path
        # string using a dictionary
        cut_region_map = self.build_zh_folder_structure()

        # Reduce number of self lookups and get the derived functions here
        histos       = self.histograms
        preselection = self.preselection
        id_functions = cut_region_map[ cut_region_map.keys()[0] ]['selection'].keys() 
        fill_histos  = self.fill_histos
        weight_func  = self.event_weight

        counter = 0
        for row in self.tree:
            # Apply basic preselection
            if not preselection(row):
                continue
            counter += 1
            #if counter < 200: print "passed preselection!"
            #map id results
            row_id_map = dict( [ (f, f( row) )  for f in id_functions] )

            # Get the generic event weight
            event_weight = weight_func(row)

            # Figure out which folder/region we are in, multiple regions allowed
            for folder, region_info in cut_region_map.iteritems():
                selection = region_info['selection']
                #if counter < 200:
                #print [ (row_id_map[f] == res) for f, res in selection.iteritems() if res is not None], all( [ (row_id_map[f] == res) for f, res in selection.iteritems() if res is not None] )
                if all( [ (row_id_map[f] == res) for f, res in selection.iteritems() if res is not None] ): #all cuts match the one of the region, None means the cut is not needed
                                                                                                            #if counter < 200: print "region found!: ",folder                   
                    fill_histos(histos, folder, row, event_weight)
                    wToApply = [ (w, w(row) )  for w in region_info['weights'] ]
                    for w_fcn, w_val in wToApply:
                        fill_histos(histos, folder+(w_fcn.__name__,), row, event_weight*w_val)
                    if len(wToApply) > 1:
                        w_prod = reduce(lambda x, y: x*y, [x for y,x in wToApply])
                        fill_histos(histos, folder+('all_weights_applied',), row, event_weight*w_prod)

    def book_general_histos(self, folder):
        '''Book general histogram, valid for each analyzer'''
        self.book(folder, "nTruePU", "NPU", 62, -1.5, 60.5)
        self.book(folder, "weight", "Event weight", 100, 0, 5)
        #self.book(folder, "weight_nopu", "Event weight without PU", 100, 0, 5)
        self.book(folder, "rho", "Fastjet #rho", 100, 0, 25)
        self.book(folder, "nvtx", "Number of vertices", 31, -0.5, 30.5)
        return None

    def book_kin_histos(self, folder, Id):
        '''books pt/jetPt/absEta for any object'''
        IdToName = {'m' : 'Muon', 'e' : 'Electron', 't' : 'Tau'}
        number   = Id[1] if len(Id) == 2 else ''
        self.book(folder, "%sPt" % Id,     "%s %s Pt" % (IdToName[Id[0]], number),     100, 0, 100)
        self.book(folder, "%sJetPt" % Id,  "%s %s Jet Pt" % (IdToName[Id[0]], number), 100, 0, 200)
        self.book(folder, "%sAbsEta" % Id, "%s %s AbsEta" % (IdToName[Id[0]], number), 100, 0, 2.4)
        return None        

    def book_mass_histos(self, folder, *args):
        IdToName = {'m' : 'Muon', 'e' : 'Electron', 't' : 'Tau'}
        def get_name(Id):
            number  = Id[1] if len(Id) == 2 else ''
            name = IdToName[Id[0]]
            return (name, number)
        for pos, obj1 in enumerate(args):
            for obj2 in args[pos+1:]:
                self.book(folder, "%s_%s_Mass" % (obj1, obj2), "%s %s - %s %s Mass" % (get_name(obj1) + get_name(obj2) ), 150, 0, 150)

    def book_resonance_histos(self, folder, products, name):
        self.book(folder, "%s_%s_Pt"     % products, "%s candidate Pt"              % name, 100, 0, 100)
        #self.book(folder, "%s_%s_AbsEta" % products, "%s candidate AbsEta"          % name, 100, 0, 2.4)
        self.book(folder, "%s_%s_Mass"   % products, "%s candidate Mass"            % name, 150, 0, 150)
        self.book(folder, "%s_%s_DR"     % products, "%s decay products #DeltaR"    % name, 100, 0, 100)
        self.book(folder, "%s_%s_DPhi"   % products, "%s decay products #Delta#phi" % name, 100, 0, 100)

    def book_Z_histos(self, folder):
        self.book_resonance_histos(folder, self.Z_decay_products(), 'Z')

    def book_H_histos(self, folder):
        self.book_resonance_histos(folder, self.H_decay_products(), 'H')
            
    def fill_histos(self, histos, folder, row, weight):
        '''fills histograms'''
        #find all keys mathing
        folder_str = '/'.join(folder + ('',))
        for key, value in histos.iteritems():
            if folder_str not in key:
                continue
            attr = key[ key.rfind('/') + 1 :]
            if attr == 'nTruePU':
                value.Fill(row.nTruePU)
            elif attr == 'weight':
                value.Fill(weight)
            else:
                value.Fill( getattr(row,attr), weight )
        return None

    def finish(self):
        self.write_histos()

if __name__ == "__main__":
    import pprint
    pprint.pprint(ZHAnalyzerBase.build_zh_folder_structure())
