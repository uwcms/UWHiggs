import rootpy.plotting.views as views
from FinalStateAnalysis.PlotTools.SimplePlotter        import SimplePlotter
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
from optparse import OptionParser
import os
import ROOT
import glob
import math
import logging
from fnmatch import fnmatch
from yellowhiggs import xs, br, xsbr

ROOT.gROOT.SetBatch()
ROOT.gStyle.SetOptStat(0)
jobid = os.environ['jobid']

print jobid
mc_samples = [
    'ggHiggsToETau',
    'vbfHiggsToETau',
    'GluGluToHToTauTau_M-125_8TeV-powheg-pythia6',
    'VBF_HToTauTau_M-125_8TeV-powheg-pythia6',
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
    'Wplus1Jets_madgraph_tapas',
    'Wplus2Jets_madgraph',
    'Wplus2Jets_madgraph_tapas',
    'Wplus3Jets_madgraph',
    'Wplus4Jets_madgraph',
    'WWJets*',
    'WZJets*',
    'ZZJets*'
]

files = []
lumifiles = []
channel = 'et'
for x in mc_samples:
    files.extend(glob.glob('results/%s/LFVHETauAnalyzer/%s.root' % (jobid, x)))
    lumifiles.extend(glob.glob('inputs/%s/%s.lumicalc.sum' % (jobid, x)))



period = '8TeV'
sqrts = 7 if '7TeV' in jobid else 8

def remove_name_entry(dictionary):
    return dict( [ i for i in dictionary.iteritems() if i[0] != 'name'] )


#sign = ['os', 'ss']
process = ['gg']
#ptcut = [0, 40]
njets= [0,1,2,3]
sign = ['os']

ptcut = [0] #was 0 



outputdir = 'plots/%s/LFVHETauAnalyzer/%s/' % (jobid, channel)
if not os.path.exists(outputdir):
    os.makedirs(outputdir)

plotter = Plotter(files, lumifiles, outputdir, None, forceLumi=19800) 
#plotter = SimplePlotter(files, lumifiles, outputdir)




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
singleT = views.StyleView(views.SumView(  *[ plotter.get_view(regex) for regex in  filter(lambda x : x.startswith('T_') or x.startswith('Tbar_'), mc_samples)]), **remove_name_entry(data_styles['T*_t*']))

SMH = views.StyleView(views.SumView( *[ plotter.get_view(regex) for regex in filter(lambda x : 'HToTauTau' in x , mc_samples)]), **remove_name_entry(data_styles['GluGluToHToTauTau*']))


plotter.views['EWKDiboson']={'view' : EWKDiboson }
plotter.views['Wplus']={'view' : Wplus }
plotter.views['DYLL']={'view' : DYLL }
plotter.views['DYTT']={'view' : DYTT }
plotter.views['singleT']={'view' : singleT }
plotter.views['SMH']={'view' : SMH }


new_mc_samples = filter( lambda x : not x.startswith('T_') and not x.startswith('Tbar_') and not  x.endswith('jets_M50_skimmedLL') and not  x.endswith('jets_M50_skimmedTT') and not x.startswith('Wplus') and not  x.startswith('WW') and not x.startswith('WZ') and not  x.startswith('ZZ') and not x.startswith('WG') and not x.endswith('HiggsToETau') and not 'HToTauTau' in x, mc_samples)

new_sigsamples= filter(lambda x: x.endswith('HiggsToETau'), mc_samples)

print new_sigsamples 
new_mc_samples.extend(['EWKDiboson', 'SMH', 'singleT','Wplus', 'DYLL', 'DYTT' 
])
print new_mc_samples

histoname = [('tPt','#tau p_{T} (GeV)',5), ('tPhi','#tau #phi',5), ('tEta','#tau #eta',2), 
             ('ePt', 'e p_{T} (GeV)', 5), ('ePhi','e #phi', 5), ('eEta','e #eta',2), 
             ('et_DeltaPhi','e-#tau #Delta#phi',1), ('et_DeltaR','e-#tau #DeltaR',1), ('tPFMET_DeltaPhi','#tau-PFMET #Delta#phi',2) ,
             ('tPFMET_Mt','#tau-PFMET M_{T} (GeV)',5),  ('tMVAMET_DeltaPhi','#tau-MVAMET #Delta#phi',2), ('tMVAMET_Mt','#tau-MVAMET M_{T} (GeV)',5),
             ('ePFMET_DeltaPhi','e-PFMET #Delta#phi',2), ('ePFMET_Mt','e-PFMET M_{T} (GeV)',5), ('eMVAMET_DeltaPhi','e-MVAMET #Delta#phi',2),
             ('eMVAMET_Mt','e-MVAMET #M_{T} (GeV)',5), ('jetN_20','Number of jets',1), ('jetN_30','Number of jets',1), 
             ('h_collmass_pfmet','M_{e#tau}coll (GeV)',1), ('h_collmass_mvamet','M_{e#tau}coll (GeV)',1), ('h_vismass','M_{e#tau} vis (GeV)',1) ]

plotter.mc_samples = new_mc_samples
#plotter.mc_samples = mc_samples
for i in sign :
    for j in process:
        for k in ptcut : 
            for nj in  njets:


                for n,h in enumerate(histoname) :
                    foldername = i+'/'+j+'/ept'+str(int(k))+'/'+str(int(nj))

                    #plotter.canvas.SetLogy(True)
                    plotter.plot_mc(foldername, ['ggHiggsToETau','vbfHiggToETau'],h[0], rebin=h[2], xaxis=h[1], leftside=False, show_ratio=False, ratio_range=1.5,  rescale=10)
                    #plotter.simpleplot_mc(foldername,h[0], rebin=h[2], xaxis= h[1], leftside=False)
                    if not os.path.exists(outputdir+foldername):
                        os.makedirs(outputdir+foldername)

                    plotter.save(foldername+'/mc_'+h)

                    foldername = i+'/'+j+'/ept'+str(k)+'/'+str(nj)+'/selected'
                    #plotter.canvas.SetLogy(True)
                    plotter.plot_mc(foldername, ['ggHiggsToETau','vbfHiggsToETau'],h[0], rebin=h[2], xaxis= h[1], leftside=False, show_ratio=False,ratio_range=3,  rescale=10)
                    if not os.path.exists(outputdir+foldername):
                        os.makedirs(outputdir+foldername)

                    plotter.save(foldername+'/mc_'+h[0])
