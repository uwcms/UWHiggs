import rootpy.plotting.views as views
#from FinalStateAnalysis.PlotTools.SimplePlotter        import SimplePlotter
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
    'ZZJets*',
    'data*'
]

files = []
lumifiles = []
channel = 'eet'
for x in mc_samples:
    #print x
    files.extend(glob.glob('results/%s/TauFakeRateAnalyzerMVA/%s.root' % (jobid, x)))
    lumifiles.extend(glob.glob('inputs/%s/%s.lumicalc.sum' % (jobid, x)))
    


period = '8TeV'
sqrts = 7 if '7TeV' in jobid else 8

def remove_name_entry(dictionary):
    return dict( [ i for i in dictionary.iteritems() if i[0] != 'name'] )


#sign = ['ss','os']
sign = ['os']
tauiso = ['tNoCuts', 'tSuperSuperLoose', 'tSuperLoose', 'tLoose', 'tTigh']
#sign = ['os']
#tauiso=['tTigh']
outputdir = 'plots/%s/TauFakeRateAnalyzerMVA/%s/' % (jobid, channel)
if not os.path.exists(outputdir):
    os.makedirs(outputdir)

plotter = Plotter(files, lumifiles, outputdir) 




EWKDiboson = views.StyleView(
    views.SumView( 
        *[ plotter.get_view(regex) for regex in \
          filter(lambda x : x.startswith('WW') or x.startswith('WZ') or x.startswith('ZZ') or x.startswith('WG'), mc_samples )]
#          filter(lambda x : x.startswith('WZ') , mc_samples )]
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


new_mc_samples = filter( lambda x : x.startswith('TTJets') , mc_samples)
#new_mc_samples = []
#new_sigsamples= filter(lambda x: x.endswith('HiggsToETau'), mc_samples)

#print new_sigsamples 
new_mc_samples.extend(['EWKDiboson', 'SMH', 'singleT','Wplus', 'DYLL', 'DYTT' , 
])
#new_mc_samples.extend(['EWKDiboson','DYLL', 'DYTT'])
print new_mc_samples


#histoname = ['e1Pt','e1Phi','e1Eta','e2Pt','e2Phi','e2Eta',
#             'e1e2Mass',  'tPt','tPhi','tEta']
#axistitle = ['e1 p_{T} (GeV)','e1 #phi','e1 #eta', 'e2 p_{T} (GeV)','e2 #phi','e2 #eta', 'e1-e2 Inv Mass (GeV)','#tau p_{T} (GeV)','#tau #phi','#tau #eta']

histoname = [('e1Pt','e1 p_{T} (GeV)', 2),('e1Phi','e1 #phi',4),('e1Eta','e1 #eta',2),('e2Pt','e2 p_{T} (GeV)',2),('e2Phi','e2 #phi',4),('e2Eta','e2 #eta',2),
             ('e1e2Mass', 'e1-e2 Inv Mass (GeV)',1), ('tPt','#tau p_{T} (GeV)',2), ('tPhi','#tau #phi',4),('tEta','#tau #eta',2),('tAbsEta','#tau |#eta|',2),
             ('etDR', 'e #tau dR', 2),  ('etDPhi', 'e #tau #Delta#phi', 2), ('ztDR', 'Z #tau dR', 2),  ('ztDPhi', 'Z #tau #Delta#phi', 2), ('type1_pfMetEt', 'type1_pfMet', 2), ( 'jetN_30', 'Number of jets, p_{T}>30', 1), ('bjetCSVVeto30', 'Number of b-jets',1) , ('tRawIso3Hits', 'tRawIso3Hits', 1)
]

#rebins = [5, 5, 2, 5, 5, 2, 1, 5, 5, 2, 1]
#rebins = []
#for n in histoname :
#    rebins.append(1)


plotter.mc_samples = new_mc_samples

print plotter.mc_samples

for i in tauiso :
    for s in sign :        
        #for histo, axis in zip(histoname, axistitle):
#        for n,h in enumerate(histoname) :
        foldername = s+'/'+i
        if not os.path.exists(outputdir+foldername):
            os.makedirs(outputdir+foldername)
 
        if i == 'tNoCuts':
            plotter.plot_mc_vs_data(foldername,'CUT_FLOW', rebin=1, xaxis='CUT_FLOW', leftside=False, show_ratio=False, ratio_range=1.5, sorted=True)
            plotter.save(foldername+'/CUT_FLOW')
            
        foldername = s+'/'+i
        for h in histoname:
            #plotter.canvas.SetLogy(True)
            #plotter.plot_data(foldername, h, rebin=rebins[n],  xaxis= axistitle[n], leftside=False)
            #plotter.plot_mc_vs_data(foldername,h, rebin=rebins[n], xaxis= axistitle[n], leftside=False, show_ratio=True, ratio_range=1.5, sort=True)
            
            plotter.plot_mc_vs_data(foldername,h[0], rebin=h[2], xaxis=h[1], leftside=False, show_ratio=True, ratio_range=1.5, sorted=True)
            plotter.save(foldername+'/'+h[0])
            
            
#        foldername = s+'/'+i+'/tptregion'
#        if not os.path.exists(outputdir+foldername):
#            os.makedirs(outputdir+foldername)
       
#        for h in histoname:
#            plotter.plot_mc_vs_data(foldername,h[0], rebin=h[2], xaxis=h[1], leftside=False, show_ratio=True, ratio_range=1.5, sorted=True)
#            plotter.save(foldername+'/'+h[0])


        jets=0
        while jets <  4 :
            foldername = s+'/'+i+'/'+str(int(jets))
            if not os.path.exists(outputdir+foldername):
                os.makedirs(outputdir+foldername)
            for h in histoname:
                plotter.plot_mc_vs_data(foldername,h[0], rebin=h[2], xaxis=h[1], leftside=False, show_ratio=True, ratio_range=1.5, sorted=True)
                plotter.save(foldername+'/'+h[0])
 #           foldername = s+'/'+i+'/'+str(int(jets))+'/tptregion'
 #           if not os.path.exists(outputdir+foldername):
 #               os.makedirs(outputdir+foldername)
 #           for h in histoname:
 #               plotter.plot_mc_vs_data(foldername,h[0], rebin=h[2], xaxis=h[1], leftside=False, show_ratio=True, ratio_range=1.5, sorted=True)
 #               plotter.save(foldername+'/'+h[0])
            jets+=1
                
