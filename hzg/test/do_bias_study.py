#! /usr/bin/env python

import os, sys
from math import sqrt

analysis_path = os.path.join(os.environ['hzganalysisroot'],
                             os.environ['hzganalysisname'])

channels = ['muon','electron']

input_mc_samples = ['ZGToLLG-8TeV/ZGToLLG_8TeV-madgraph.Summer12_DR53X-PU_S10_START53_V7A-v1.AODSIM.root',
                    'DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball-8TeV/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball.Summer12_DR53X-PU_S10_START53_V7A-v1.AODSIM.root']

study_inputs = {}

for channel in channels:
    base_dir = os.path.join(analysis_path,channel)
    study_inputs[channel] = \
             [os.path.join(base_dir,sample) for sample in input_mc_samples]

import ROOT
from ROOT import RooWorkspace, RooMCStudy, RooFit, TFile, gDirectory,\
     RooDataSet, RooFormulaVar, RooConstVar, RooArgList, RooArgSet,\
     RooMinuit
ROOT.gSystem.Load("libUWHiggshzg.so")
pwd = gDirectory.GetPath()

int_lumi = 19.6*1000 #19.6/fb in /pb

def prepare_truth_models(ws,truth,cat,mass):    
    for channel in study_inputs:
        bkg_data = RooDataSet("bkgdata_%s"%(channel),
                              "M_{ll#gamma} with Errors",
                              ws.set("observables_weight"),
                              "weight")
        
        for k,file in enumerate(study_inputs[channel]):
            f_in = TFile.Open(file,"READ")
            gDirectory.cd(pwd)
            n_events = f_in.Get("eventCount").GetBinContent(1)
            tree = f_in.Get("selected_zg").CloneTree()
            tree.SetName("tree_%s_%i"%(channel,k))
            f_in.Close()            
            
            d_in = RooDataSet("temp",
                              "M_{ll#gamma} with Errors",
                              tree,
                              ws.set("observables"))
            norm = RooConstVar("norm","norm",n_events)
            weight = RooFormulaVar("weight",
                                   "weight",
                                   "1.3*@0*@1/@2", #1.3 gives 19.6/fb
                                   RooArgList(ws.var("puWeight"),
                                              ws.var("procWeight"),
                                              norm)
                                   )
            d_in.addColumn(weight)
            d_in_weight = RooDataSet("temp_weight",
                                     "M_{ll#gamma} with Errors",
                                     d_in,
                                     ws.set("observables_weight"),
                                     '','weight')
            bkg_data.append(d_in_weight)
            
        # split off data for each category, create the
        # toy dataset truth models
        data = bkg_data.reduce("Mz + Mzg > 185 && r94cat == %i"%cat)
        getattr(ws,'import')(data,
                             RooFit.Rename('bkgdata_%s_%i'%(channel,
                                                            cat)))
        #make RooDecay 'truth' model
        ws.factory(
            'RooGaussModel::MzgResoShape_%s_cat%i(Mzg,'\
            'bias_%s_cat%i[120,90,150],sigma_%s_cat%i[1,0.01,10])'%(
            channel,cat,
            channel,cat,
            channel,cat)
            )
        ws.factory(
            'RooDecay::MzgTruthModelBase_%s_cat%i(Mzg,'\
            'tau_%s_cat%i[5,0,50],MzgResoShape_%s_cat%i,'
            'RooDecay::SingleSided)'%(channel,cat,
                                      channel,cat,
                                      channel,cat)
            )
        nevts = data.sumEntries('Mzg > %f && Mzg < %f'%(mass-1.5,mass+1.5))
        ws.factory(
            'RooExtendPdf::MzgTruthModel_%s_cat%i('\
            'MzgTruthModelBase_%s_cat%i,'\
            'norm_truth_%s_cat%i[%f,%f,%f],"ROI")'%(channel,cat,
                                                    channel,cat,
                                                    channel,cat,
                                                    nevts,
                                                    0.75*nevts,1.25*nevts)
            )
        ws.pdf('MzgTruthModel_%s_cat%i'%(channel,cat)).fitTo(
            ws.data('bkgdata_%s_%i'%(channel,cat)),
            RooFit.SumW2Error(True)
            )
        
poly_order ={1:4,2:5,3:5,4:5}
def build_fitting_models(ws,cat,mass):
    ws.var('Mzg').setBins(20000,'cache')

    cs=['c%i_cat%i[-1e-6,1.01]'%(k+1,cat) for k in range(poly_order[cat])]
    config = ['Mzg',
              'stepVal_cat%i[0.1,0,1]'%cat]
    config.append('{'+','.join(['1.0']+cs)+'}')
    #build standard candle model (RooStepBerntein(x)Gaus)        
    ws.factory(
        'RooStepBernstein::RSBFitModelTruth_cat%i(%s)'%(cat,
                                                        ','.join(config))
        )
    ws.factory(
        "RooGaussian::RSBFitModelReso_cat%i(%s)"%(
        cat,
        ','.join(['Mzg',
                  'rsb_bias_cat%i[0]'%cat,
                  'rsb_sigma_cat%i[5,0.01,20]'%cat])
        )
        )
    ws.factory('FCONV::RSBFitModelBase_cat%i(Mzg,'\
               'RSBFitModelTruth_cat%i,'\
               'RSBFitModelReso_cat%i)'%(cat,cat,cat))
    # setup extended pdf with dummy values for norm
    # set proper values during toys loop
    ws.factory(
        'RooExtendPdf::RSBFitModel_cat%i(RSBFitModelBase_cat%i,'\
        'norm_rsb_cat%i[1000,500,1500],"ROI")'%(cat,cat,cat))
    #build alternate fitting model (Bernstein(x)logistics)
    ws.factory(
        'EXPR::RSBFitModelAltReso_cat%i("'\
        'exp(-(@0-@1)/@2)/(@2*(1.0+exp(-(@0-@1)/@2))**2)'\
        '",%s)'%(
        cat,
        ','.join(['Mzg',
                  'rsb_altbias_cat%i[0]'%cat,
                  'rsb_altsigma_cat%i[5,0.01,20]'%cat])
        )
        )
    ws.factory('FCONV::RSBFitModelAltBase_cat%i(Mzg,'\
               'RSBFitModelTruth_cat%i,'\
               'RSBFitModelAltReso_cat%i)'%(cat,cat,cat))
    # setup extended pdf with dummy values for norm
    # set proper values during toys loop
    ws.factory(
        'RooExtendPdf::RSBFitModelAlt_cat%i(RSBFitModelAltBase_cat%i,'\
        'norm_altrsb_cat%i[1000,500,1500],"ROI")'%(cat,cat,cat))
    
    
    #ws.data('bkgdata_electron_1').Print("V")
    sumEntries = ws.data('bkgdata_electron_%i'%cat).sumEntries('Mzg > %f && Mzg < %f'%(mass-1.5,mass+1.5))
    

def gen_data_and_fit(ws, iterations,cat, mass):
    
    for channel in channels:        
        
        n_sample = int(ws.data('bkgdata_%s_%i'%(channel,cat)).sumEntries())
        ws.factory('N_ROI_true_%s_cat%i[0]'%(channel,cat))
        ws.factory('N_ROI_erf_%s_cat%i[0]'%(channel,cat))
        ws.factory('N_ROI_sigm_%s_cat%i[0]'%(channel,cat))
        ws.factory('N_err_ROI_true_%s_cat%i[0]'%(channel,cat))
        ws.factory('N_err_ROI_erf_%s_cat%i[0]'%(channel,cat))
        ws.factory('N_err_ROI_sigm_%s_cat%i[0]'%(channel,cat))
        ws.factory('pull_ROI_sigm_%s_cat%i[0]'%(channel,cat))
        ws.factory('pull_ROI_erf_%s_cat%i[0]'%(channel,cat))
        ws.defineSet('biasVars_%s_cat%i'%(channel,cat),
                     'N_ROI_true_%s_cat%i,N_ROI_erf_%s_cat%i,N_ROI_sigm_%s_cat%i,'\
                     'N_err_ROI_true_%s_cat%i,N_err_ROI_erf_%s_cat%i,'\
                     'N_err_ROI_sigm_%s_cat%i,pull_ROI_sigm_%s_cat%i,'\
                     'pull_ROI_erf_%s_cat%i'%(channel,cat,channel,cat,
                                              channel,cat,channel,cat,
                                              channel,cat,channel,cat,
                                              channel,cat,channel,cat)
                     )
        biasData = RooDataSet('biasData_%s_cat%i'%(channel,cat),
                              'bias data',
                              ws.set('biasVars_%s_cat%i'%(channel,cat)))
        
        for i in xrange(iterations):
            
            #build data from our truth model
            truth = ws.pdf('MzgTruthModel_%s_cat%i'%(channel,cat))
            toy_data = truth.generate(                    
                RooArgSet(ws.var('Mzg')),
                n_sample,
                RooFit.Extended(),
                RooFit.Name('toy_%s_cat%i_%i'%(channel,cat,i))
                )
            
            true_ROI_yield = toy_data.sumEntries('Mzg > %f && Mzg < %f'%(mass-1.5,mass+1.5))
            ws.var('N_ROI_true_%s_cat%i'%(channel,cat)).setVal(
                true_ROI_yield
                )
            ws.var('N_err_ROI_true_%s_cat%i'%(channel,cat)).setVal(
                sqrt(true_ROI_yield)
                )
            
            #setup the normalizations
            sumEntries = toy_data.sumEntries('Mzg > %f && Mzg < %f'%(mass-1.5,mass+1.5))
            ws.var('norm_altrsb_cat%i'%cat).setMin(sumEntries*0.70)
            ws.var('norm_altrsb_cat%i'%cat).setMax(sumEntries*1.30)
            ws.var('norm_altrsb_cat%i'%cat).setVal(sumEntries)
            
            ws.var('norm_rsb_cat%i'%cat).setMin(sumEntries*0.70)
            ws.var('norm_rsb_cat%i'%cat).setMax(sumEntries*1.30)
            ws.var('norm_rsb_cat%i'%cat).setVal(sumEntries)  
            
            minos_var = RooArgSet(ws.var('norm_rsb_cat%i'%cat))
            sigm_nll = ws.pdf('RSBFitModelAlt_cat%i'%cat).createNLL(
                toy_data
                )
            sigm_min = RooMinuit(sigm_nll)
            sigm_min.simplex()
            sigm_min.migrad()
            sigm_min.hesse()
            #sigm_min.minos(minos_var)
            
            fit_sigm_norm = ws.var('norm_altrsb_cat%i'%cat).getVal()
            fit_sigm_err  = sqrt(fit_sigm_norm)#ws.var('norm_altrsb_cat%i'%cat).getError()
            fit_sigm_pull = (fit_sigm_norm - true_ROI_yield)/fit_sigm_err
            
            ws.var('N_ROI_sigm_%s_cat%i'%(channel,cat)).setVal(
                fit_sigm_norm
                )
            ws.var('N_err_ROI_sigm_%s_cat%i'%(channel,cat)).setVal(
                fit_sigm_err
                )                
            ws.var('pull_ROI_sigm_%s_cat%i'%(channel,cat)).setVal(
                fit_sigm_pull
                )
            
            minos_var = RooArgSet(ws.var('norm_rsb_cat%i'%cat))
            gaus_nll = ws.pdf('RSBFitModel_cat%i'%cat).createNLL(
                toy_data                    
                )
            gaus_min= RooMinuit(gaus_nll)
            gaus_min.simplex()
            gaus_min.migrad()
            gaus_min.hesse()
            #gaus_min.minos(minos_var)
            
            fit_erf_norm = ws.var('norm_rsb_cat%i'%cat).getVal()
            fit_erf_err  = sqrt(fit_erf_norm)#ws.var('norm_rsb_cat%i'%cat).getError()
            fit_erf_pull = (fit_erf_norm - true_ROI_yield)/fit_erf_err
            
            ws.var('N_ROI_erf_%s_cat%i'%(channel,cat)).setVal(
                fit_erf_norm
                )
            ws.var('N_err_ROI_erf_%s_cat%i'%(channel,cat)).setVal(
                fit_erf_err
                )
            ws.var('pull_ROI_erf_%s_cat%i'%(channel,cat)).setVal(
                fit_erf_pull
                )
            
            biasData.add(ws.set('biasVars_%s_cat%i'%(channel,cat)))
            getattr(ws,'import')(toy_data)
            getattr(ws,'import')(biasData)
            


# if being executed run bias study
if __name__ == '__main__':
    ntoys = int(sys.argv[1])
    truth = sys.argv[2]
    category = int(sys.argv[3])
    mass = float(sys.argv[4])
    
    
    bs = RooWorkspace('bias_study')

    bs.factory("procWeight[0]")
    bs.factory("puWeight[0]")
    bs.factory("weight[0]")
    bs.factory("Mzg[100,180]")
    bs.var("Mzg").setRange("ROI",mass-1.5,mass+1.5)
    bs.factory("Mz[0]")
    bs.factory("dMzg[0,25]")
    bs.factory("dMz[0,25]")
    bs.factory("r94cat[cat1=1,cat2=2,cat3=3,cat4=4]")
    bs.defineSet("observables",
                 "Mzg,Mz,dMzg,dMz,r94cat,procWeight,puWeight")
    bs.defineSet("observables_weight",
                 "Mzg,Mz,dMzg,dMz,r94cat,procWeight,puWeight,weight")

    prepare_truth_models(bs,truth,category,mass)

    build_fitting_models(bs,category,mass)

    gen_data_and_fit(bs, ntoys, category,mass)

    out_f = TFile.Open("bias_study_ntoys%i_true%s_cat%i_m%s.root"%(ntoys,truth,
                                                                   category,str(mass).replace('.','p')),
                       "recreate")
    bs.Write()
    out_f.Close()
