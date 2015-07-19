import rootpy.plotting.views as views
from FinalStateAnalysis.PlotTools.Plotter        import Plotter
from FinalStateAnalysis.PlotTools.BlindView      import BlindView
from FinalStateAnalysis.PlotTools.PoissonView    import PoissonView
from FinalStateAnalysis.PlotTools.MedianView     import MedianView
from FinalStateAnalysis.PlotTools.ProjectionView import ProjectionView
#from FinalStateAnalysis.PlotTools.FixedIntegralView import FixedIntegralView
from FinalStateAnalysis.PlotTools.RebinView  import RebinView
from FinalStateAnalysis.MetaData.data_styles import data_styles, colors
from FinalStateAnalysis.PlotTools.decorators import memo
from FinalStateAnalysis.MetaData.datacommon  import br_w_leptons, br_z_leptons
from FinalStateAnalysis.PlotTools.SubtractionView      import SubtractionView, PositiveView
from optparse import OptionParser
import os
import ROOT
import glob
import math
import logging
import pdb
import array
from fnmatch import fnmatch
from yellowhiggs import xs, br, xsbr

from BasePlotter import BasePlotter


ROOT.gROOT.SetBatch()
ROOT.gStyle.SetOptStat(0)
jobid = os.environ['jobid']

print jobid
mc_samples = [
    'ggHiggsToETau',
    'vbfHiggsToETau',
    'GluGluToHToTauTau_M-125_8TeV-powheg-pythia6',
    'VBF_HToTauTau_M-125_8TeV-powheg-pythia6',
    #'Zjets_M50', 
    'Zjets_M50_skimmedLL',
    'Z1jets_M50_skimmedLL',
    'Z2jets_M50_skimmedLL',
    'Z3jets_M50_skimmedLL',
    'Z4jets_M50_skimmedLL',
    'Zjets_M50_skimmedTT',
    'Z1jets_M50_skimmedTT',
    'Z2jets_M50_skimmedTT',
    'Z3jets_M50_skimmedTT',
    'Z4jets_M50_skimmedTT',
    'TTJets*',
    'T_t*',
    'Tbar_t*', 
    'WplusJets_madgraph_skimmed',
    'Wplus1Jets_madgraph',
   # 'Wplus1Jets_madgraph_tapas',
    'Wplus2Jets_madgraph',
   # 'Wplus2Jets_madgraph_tapas',
    'Wplus3Jets_madgraph',
    'Wplus4Jets_madgraph',
    'WWJets*',
    'WZJets*',
    'ZZJets*',
    'data*'
]

files = []
lumifiles = []
channel = 'ee'
for x in mc_samples:
#    files.extend(glob.glob('results/%s/EEAnalyzerMVA/%s.root' % (jobid, x)))
    files.extend(glob.glob('results/%s/EEAnalyzerMVA_noTrCorr/%s.root' % (jobid, x)))
    lumifiles.extend(glob.glob('inputs/%s/%s.lumicalc.sum' % (jobid, x)))
    


period = '8TeV'
sqrts = 7 if '7TeV' in jobid else 8

def remove_name_entry(dictionary):
    return dict( [ i for i in dictionary.iteritems() if i[0] != 'name'] )


sign = ['os']
jets = [0, 1, 2, 3]
#outputdir = 'plots/%s/EEAnalyzerMVA/%s/' % (jobid, channel)
outputdir = 'plots/%s/EEAnalyzerMVA_noTrCorr/%s/' % (jobid, channel)
if not os.path.exists(outputdir):
    os.makedirs(outputdir)

plotter = BasePlotter(files, lumifiles, outputdir) 




EWKDiboson = views.StyleView(
    views.SumView( 
        *[ plotter.get_view(regex) for regex in \
          filter(lambda x : x.startswith('WW') or x.startswith('WZ') or x.startswith('ZZ') or x.startswith('WG'), mc_samples )]
    ), **remove_name_entry(data_styles['WW*'#,'WZ*', 'WG*', 'ZZ*'
])
)
Wplus = views.StyleView(views.SumView(  *[ plotter.get_view(regex) for regex in filter(lambda x :  x.startswith('Wplus'), mc_samples )]), **remove_name_entry(data_styles['Wplus*Jets*']))
DYLL = views.StyleView(views.SumView( *[ plotter.get_view(regex) for regex in filter(lambda x :  x.endswith('skimmedLL'), mc_samples )]), **remove_name_entry(data_styles['Z*jets*LL']))
DYTT = views.StyleView(views.SumView( *[ plotter.get_view(regex) for regex in filter(lambda x :  x.endswith('jets_M50_skimmedTT'), mc_samples )]), **remove_name_entry(data_styles['Z*jets*TT']))
TT = views.StyleView(views.SumView(  *[ plotter.get_view(regex) for regex in  filter(lambda x : x.startswith('TT') , mc_samples)]), **remove_name_entry(data_styles['TTJets*']))
singleT = views.StyleView(views.SumView(  *[ plotter.get_view(regex) for regex in  filter(lambda x : x.startswith('T_') or x.startswith('Tbar_'), mc_samples)]), **remove_name_entry(data_styles['T*_t*']))
SMH = views.StyleView(views.SumView( *[ plotter.get_view(regex) for regex in filter(lambda x : 'HToTauTau' in x , mc_samples)]), **remove_name_entry(data_styles['GluGluToHToTauTau*']))

plotter.views['EWKDiboson']={'view' : EWKDiboson }
plotter.views['Wplus']={'view' : Wplus }
plotter.views['DYLL']={'view' : DYLL }
plotter.views['DYTT']={'view' : DYTT }
plotter.views['singleT']={'view' : singleT }
plotter.views['SMH']={'view' : SMH }
plotter.views['TT']={'view' : TT }
#plotter.views['DY']={'view' : DY }
new_mc_samples =[]

new_mc_samples.extend(['EWKDiboson', 'SMH', 'singleT','Wplus',  'TT',
    'DYLL', 'DYTT'
])
print new_mc_samples

binx = array.array('l', [5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,100, 120, 140, 160,200 ])

histoname = [('e1Pt', 'e1 p_{T} (GeV)', 5), ('e1Phi','e1 #phi',5) ,('e1Eta','e1 #eta',2),
             ('e2Pt','e2 p_{T} (GeV)',5) ,('e2Phi','e2 #phi',5), ('e2Eta','e2 #eta',2),
             ('e1e2_DeltaPhi','e1-e2 #Delta#phi',2), ('e1e2_DeltaR','e1-e2 #DeltaR',2), ('e1e2Mass', 'e1-e2 Inv Mass (GeV)',5),
             ('type1_pfMetEt', 'type1PF MET (GeV)', 2) , ('pfMetEt', 'PF MET (GeV)', 2), 
             ('type1_pfMetPhi','type1 PF MET #phi', 5), ('pfMetPhi','PF MET #phi',5),  
             ('pfMetEt_par','PF MET parallel (GeV)',2), ('pfMetEt_perp','PF MET perpendicular (GeV)',2),
             ('type1_pfMetEt_par', 'type1PF MET parallel (GeV)', 2),('type1_pfMetEt_perp', 'type1PF MET perpendicular (GeV)',2),
             ('e1PFMET_DeltaPhi','e1-type1PFMET #Delta#phi',2), ('e1PFMET_Mt','e1-type1PFMET M_{T} (GeV) ',5),
             ('e2PFMET_DeltaPhi','e2-type1PFMET #Delta#phi',2),('e2PFMET_Mt','e2-type1PFMET M_{T} (GeV)',5), 
             ('nPV_unweighted', 'unweighted N of vertices', 1), ('nPV', 'number of vertices', 1)
]


plotter.mc_samples = new_mc_samples
foldernames = []
for i in sign :
    foldernames.append(i)
    for j in jets :
        foldernames.append(i+'/'+str(int(j)))
        
def create_mapper(mapping):
    def _f(path):
        for key, out in mapping.iteritems():
            if key == path:
                path = path.replace(key,out)
                print 'path', path
        return path
    return _f

def get_ss(x):
    return x.replace('os/', 'ss/')

print foldernames
mymapper = {"os/": "ss/",
            "os/0/": "ss/0/",
            "os/1": "ss/1",
            "os/2": "ss/2",
            "os/3": "ss/3"
}
QCD =  views.TitleView(views.StyleView(SubtractionView(
    views.PathModifierView( plotter.data,  get_ss),
    views.PathModifierView( TT, get_ss),
    views.PathModifierView( singleT,  get_ss),
    views.PathModifierView( SMH,   get_ss ),
    views.PathModifierView( DYTT,  get_ss ),
    views.PathModifierView( DYLL,   get_ss ),
    views.PathModifierView( Wplus,  get_ss ),
    views.PathModifierView( EWKDiboson,  get_ss ))
                                       ,**data_styles['QCD*']), 'QCD')


plotter.views['QCD']= {'view': QCD}
plotter.mc_samples.extend(['QCD'])

                                  
for foldername in foldernames:
    if foldername.startswith("ss") and  bool('QCD' in plotter.mc_samples)==True:
        plotter.mc_samples.remove('QCD')
    if foldername.startswith("os") and  bool('QCD' in plotter.mc_samples)==False:
        plotter.mc_samples.extend(['QCD'])

    for n,h in enumerate(histoname) :
        
        
        plotter.pad.SetLogy(True)
        #plotter.plot('QCD', foldername+'/'+h[0], 'hist')
        print foldername
        #plotter.plot_mc_vs_data(foldername, h[0], rebin=h[2], xaxis=h[1], leftside=False, show_ratio=True, ratio_range=0.2, sorted=True)
        
        plotter.plot_with_bkg_uncert(foldername,h[0], rebin=h[2], xaxis=h[1], leftside=False, show_ratio=True, ratio_range=0.5, sort=True, obj1='e1', obj2='e2')




        if not os.path.exists(outputdir+foldername):
            os.makedirs(outputdir+foldername)
            
        plotter.save(foldername+'/'+h[0])

            
 
