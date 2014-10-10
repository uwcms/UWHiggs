#from mauro plotters
import os
from sys import argv, stdout, stderr
import ROOT
import sys
from FinalStateAnalysis.PlotTools.MegaBase import make_dirs
from FinalStateAnalysis.MetaData.data_styles import data_styles
from FinalStateAnalysis.PlotTools.BlindView import BlindView,  blind_in_range
from FinalStateAnalysis.PlotTools.SubtractionView      import SubtractionView, PositiveView

import glob
import logging
import sys
from BasePlotter import BasePlotter

jobid = os.environ['jobid']
#jobid = 'MCntuples_3March' 
channel = 'et'
import rootpy.plotting.views as views
        
ROOT.gROOT.SetStyle("Plain")
ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)

print "\nPlotting %s for %s\n" % (channel, jobid)

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
    #'WplusJets_madgraph_skimmed',
    'Wplus1Jets_madgraph*',
    #'Wplus1Jets_madgraph_tapas',
    'Wplus2Jets_madgraph*',
    #'Wplus2Jets_madgraph_tapas',
    'Wplus3Jets_madgraph',
    'Wplus4Jets_madgraph',
    'WWJets*',
    'WZJets*',
    'ZZJets*',
    'Fake*',
    'data*'
]

        
files = []
lumifiles = []

for x in mc_samples:
        files.extend(glob.glob('results/%s/LFVHETauAnalyzerMVA_noOverlap/%s.root' % (jobid, x)))
        lumifiles.extend(glob.glob('inputs/%s/%s.lumicalc.sum' % (jobid, x)))

sign = ['os', 'ss']
jets = [0, 1, 2]
processtype=['gg']
threshold=['ept30']

outputdir = 'plots/%s/LFVHETauAnalyzerMVA/%s/' % (jobid, channel)
if not os.path.exists(outputdir):
        os.makedirs(outputdir)


        
def remove_name_entry(dictionary):
    return dict( [ i for i in dictionary.iteritems() if i[0] != 'name'] )




plotter = BasePlotter(channel,files, lumifiles, outputdir, forceLumi=19800) 

ggMCSignal =views.ScaleView( views.TitleView(views.StyleView(plotter.get_view('ggHiggsToETau'), **remove_name_entry(data_styles['ggHiggsToETau*'])), 'LFV gg Higgs, BR = 1.'),1.)
vbfMCSignal = views.ScaleView(views.TitleView(views.StyleView(plotter.get_view('vbfHiggsToETau'), **remove_name_entry(data_styles['vbfHiggsToETau*'])), 'LFV vbf Higgs, BR = 1.'),1.)

plotter.views['ggMCSignal']= {'view': ggMCSignal}
plotter.views['vbfMCSignal']= {'view': vbfMCSignal}

EWKDiboson = views.TitleView(views.StyleView(
    views.SumView( 
        *[ plotter.get_view(regex) for regex in \
          filter(lambda x : x.startswith('WW') or x.startswith('WZ') or x.startswith('ZZ') or x.startswith('WG'), mc_samples )]
    ), **remove_name_entry(data_styles['WW*'#,'WZ*', 'WG*', 'ZZ*'
])
), 'EWK Dibosons')
Wplus = views.TitleView(views.StyleView(views.SumView(  *[ plotter.get_view(regex) for regex in filter(lambda x :  x.startswith('Wplus'), mc_samples )]), **remove_name_entry(data_styles['Wplus*Jets*'])), 'W + jets') 
DYLL = views.TitleView(views.StyleView(views.SumView( *[ plotter.get_view(regex) for regex in filter(lambda x :  x.endswith('skimmedLL'), mc_samples )]), **remove_name_entry(data_styles['Z*jets*LL'])), 'DY (#rightarrow ll)  + jets') 
DYLLfake = views.SubdirectoryView(views.SumView( *[ plotter.get_view(regex) for regex in filter(lambda x :  x.endswith('jets_M50_skimmedLL'), mc_samples )]), 'tLoose')

DYTT = views.TitleView(views.StyleView(views.SumView( *[ plotter.get_view(regex) for regex in filter(lambda x :  x.endswith('jets_M50_skimmedTT'), mc_samples )]), **remove_name_entry(data_styles['Z*jets*TT'])), 
'DY (#rightarrow #tau#tau)  + jets')
TT = views.TitleView(views.StyleView(views.SumView(  *[ plotter.get_view(regex) for regex in  filter(lambda x : x.startswith('TT') , mc_samples)]), **remove_name_entry(data_styles['TTJets*'])), 't#bar{t}')
singleT = views.TitleView(views.StyleView(views.SumView(  *[ plotter.get_view(regex) for regex in  filter(lambda x : x.startswith('T_') or x.startswith('Tbar_'), mc_samples)]), **remove_name_entry(data_styles['T*_t*'])), 'Single Top')
SMH = views.TitleView(views.StyleView(views.SumView( *[ plotter.get_view(regex) for regex in filter(lambda x : 'HToTauTau' in x , mc_samples)]), **remove_name_entry(data_styles['GluGluToHToTauTau*'])), 'SM Higgs')


plotter.views['EWKDiboson']={'view' : EWKDiboson }
plotter.views['Wplus']={'view' : Wplus }
#plotter.views['DYLL']={'view' : DYLL }
plotter.views['DYTT']={'view' : DYTT }
plotter.views['singleT']={'view' : singleT }
plotter.views['SMH']={'view' : SMH }
plotter.views['TT']={'view' : TT }
#plotter.views['DY']={'view' : DY }
new_mc_samples =[]

#DYLL  = views.TitleView(views.StyleView(views.SubtractionView(views.SumView( *[ plotter.get_view(regex) for regex in filter(lambda x :  x.endswith('jets_M50_skimmedTT'), mc_samples )]), views.SubdirectoryView(views.SumView( *[ plotter.get_view(regex) for regex in filter(lambda x :  x.endswith('jets_M50_skimmedLL'), mc_samples )]), 'tLoose')), **remove_name_entry(data_styles['Z*jets*LL'])), 'DY (#rightarrow ll)  + jets')
finalDYLL = SubtractionView(DYLL,DYLLfake)
plotter.views['finalDYLL']={'view' : finalDYLL }

new_mc_samples.extend(['EWKDiboson', 
                       'SMH', 'singleT',#'Wplus',  
                       'TT',
                       'finalDYLL', 'DYTT'#, 'ggMCSignal', 'vbfMCSignal'
])
def get_fakeTaus(x):
        y=x
        if x.startswith('os') or x.startswith('ss'):
                y = x.replace('.*s/', 'tLoose/*s/')
        
        return y

#myFake = views.TitleView(views.StyleView(views.SumView(  *[ plotter.get_view(regex) for regex in  filter(lambda x : x.startswith('Fake'), mc_samples)]),**data_styles['Fakes*']), 'Fakes')

Fakes =  views.SubdirectoryView(views.TitleView(views.StyleView(views.SumView( *[ plotter.get_view(regex) for regex in filter(lambda x :  'Fake*' in x, mc_samples)]), **remove_name_entry(data_styles['Fakes*'])),'Fakes'), 'tLoose/')

#SumView(views.PathModifierView( views.SumView(  *[ plotter.get_view(regex) for regex in  filter(lambda x : x.startswith('Fake'), mc_samples)]),  get_fakeTaus)),**remove_name_entry(data_styles['Fakes*'])), 'Fakes')
plotter.views['Fakes']= {'view': Fakes}
new_mc_samples.extend(['Fakes'])

print new_mc_samples

histoname = [('tPt', 'p_{T}(#tau) (GeV)', 5), ('tEta', '#eta(#tau)', 2),  ('tPhi', '#phi(#tau)', 5), 
             ('ePt', 'p_{T}(e) (GeV)', 5), ('eEta', '#eta(e)', 2),  ('ePhi', '#phi(e)', 5), 
             ('et_DeltaPhi', 'e#tau #Delta#phi', 1.), ('et_DeltaR', 'e#tau #Delta R', 1.),
             ('h_collmass_pfmet', 'M_{coll}(e#tau) (GeV)', 1.), ('h_vismass', 'M_{vis} (GeV)', 1.),
             ('jetN_30', 'number of jets (p_{T} > 30 GeV)', 1.) , ('ePFMET_Mt', 'M_{T} e-PFMET', 5), 
             ('tPFMET_Mt', 'M_{T} #tau-PFMET', 5) ]


plotter.mc_samples = new_mc_samples
foldernames = []
for i in sign:
        for j in processtype:
                for k in threshold:
                        for jn in jets: 

                                foldernames.append(i+'/'+j+'/'+k +'/'+str(jn))
                                foldernames.append(i+'/'+j+'/'+k +'/'+str(jn)+'/selected')



for foldername in foldernames:
        for n,h in enumerate(histoname) :
        
        
                plotter.pad.SetLogy(False)
                #print foldername
                
                #plotter.simpleplot_mc(foldername,['ggMCSignal','vbfMCSignal'],  h[0], rebin=int(h[2]), xaxis=h[1], leftside=False, forceLumi=19800)
                #plotter.plot_with_bkg_uncert(foldername,h[0], rebin=int(h[2]), xaxis=h[1], leftside=False, show_ratio=True, ratio_range=0.5, sort=True, obj=['e'])
                plotter.plot_without_uncert(foldername,h[0], rebin=int(h[2]), xaxis=h[1], leftside=False, show_ratio=True, ratio_range=0.5, sort=True, obj=['e'])
                #plotter.plot_mc_vs_data(foldername, h[0], rebin=int(h[2]), xaxis=h[1], leftside=False, show_ratio=False, ratio_range=0.2)



                if not os.path.exists(outputdir+foldername):
                        os.makedirs(outputdir+foldername)
            
                plotter.save(foldername+'/'+h[0])


