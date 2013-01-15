#!/usr/bin/env python

from UWHiggs.hzg.datacard.directory_prep import directory_prep
from UWHiggs.hzg.datacard.metadata_association import metadata_association
from UWHiggs.hzg.datacard.categories_map import categories_map,\
     make_background_for_cat,make_signal_for_cat
from UWHiggs.hzg.plotting.plotting_groups import make_plot_groups
from UWHiggs.hzg.plotting.hzg_plots import process_plot_groups

import os

from RecoLuminosity.LumiDB import argparse

parser = argparse.ArgumentParser(description="build hzg datacards and plots")
parser.add_argument('--unblind', default=False,action='store_true',
                    help='run without blinding cuts')

args = parser.parse_args()

import ROOT
from ROOT import TFile,TH1F,TTree,RooWorkspace,gDirectory,RooDataSet,RooFit,\
     RooArgSet, RooFormulaVar, RooArgList, RooConstVar, TGraph

def initialize_workspace(ws):
    #make the category definition
    cattype   = os.environ['hzgcategorytype']
    catstring = categories_map[cattype]['catstring']
    catname   = categories_map[cattype]['leafname']
    ws.factory(catstring)

    #fitting variables (init with no implicit range)
    ws.factory("Mzg[0]")
    
    #(possible) weights
    ws.factory('procWeight[0]')
    ws.factory('puWeight[0]')
    ws.factory('weight[0]')

    ws.defineSet('vars_with_weights','Mzg,procWeight,puWeight,%s'%catname)
    ws.defineSet('vars_with_weights_final',
                 'Mzg,procWeight,puWeight,weight,%s'%catname)
    ws.defineSet('vars','Mzg,%s'%catname)

def make_weighted_dataset(subproc,ws,tree,mc_events):
    data = RooDataSet('%s_shape_data'%subproc,
                      'M_{ll#gamma} Shape Data for %s'%subproc,
                      tree,
                      ws.set('vars_with_weights')                      
                      )

    mc_yield_var = RooConstVar('temp','temp',mc_events)
    weighter = RooFormulaVar('weight','weight','@0*@1/@2',
                             RooArgList( ws.var('procWeight'),
                                         ws.var('puWeight'),
                                         mc_yield_var )
                             )
    data.addColumn(weighter)

    data_total_weight = RooDataSet('%s_shape_data'%subproc,
                                   'M_{ll#gamma} Shape Data for %s'%subproc,
                                   data,
                                   ws.set('vars_with_weights_final'),
                                   '','weight')

    data_pu_weight = RooDataSet('%s_shape_data_puonly'%subproc,
                                'M_{ll#gamma} Shape Data for %s'%subproc,
                                data,
                                ws.set('vars_with_weights_final'),
                                '','puWeight')
    
    return data_total_weight, data_pu_weight

def extract_higgs_data_in_categories(subproc,input_file,ws):
    pwd = gDirectory.GetPath()    
    fin = TFile.Open(input_file,'read')
    gDirectory.cd(pwd)

    tree = fin.Get('selected_zg')    

    mc_tot_events = float(fin.Get('eventCount').GetBinContent(1))    

    data_total_weight,data_pu_weight = make_weighted_dataset(subproc,
                                                             ws,tree,
                                                             mc_tot_events)

    ws.factory('%s_acceff[%f]'%(subproc,
                                data_pu_weight.sumEntries()/mc_tot_events))
    
    fin.Close()
    
    #save with and without process cross section normalization
    getattr(ws,'import')(data_total_weight)
    getattr(ws,'import')(data_pu_weight)
    return mc_tot_events

def extract_bkg_data_in_categories(subproc,input_file,ws):
    pwd = gDirectory.GetPath()    
    fin = TFile.Open(input_file,'read')
    gDirectory.cd(pwd)

    tree = fin.Get('selected_zg')    

    mc_tot_events = float(fin.Get('eventCount').GetBinContent(1))

    data_proper_weight,data_pu_weight = make_weighted_dataset(subproc,
                                                              ws,tree,
                                                              mc_tot_events)
    
    fin.Close()    
    
    data_in_ws = ws.data('mc_background_shape_data')
    if not not data_in_ws:
        getattr(ws,'import')(data_proper_weight)
        data_in_ws.append(data_proper_weight)        
    else:
        getattr(ws,'import')(data_proper_weight,
                             RooFit.Rename('mc_background_shape_data'))
        getattr(ws,'import')(data_proper_weight)
    return mc_tot_events
                             

def extract_data_in_categories(channel,input_file,ws):
    pwd = gDirectory.GetPath()    
    fin = TFile.Open(input_file,'read')
    gDirectory.cd(pwd)

    tree = fin.Get('selected_zg')

    data = RooDataSet('%s_data'%channel,
                      'real data %s channel'%channel,
                      tree,
                      ws.set('vars') )    


    data_in_ws = ws.data('%s_data'%channel)
    if not not data_in_ws:
        data_in_ws.append(data)
    else:
        getattr(ws,'import')(data)

def create_master_workspaces(meta_data):
    pwd = gDirectory.GetPath()
    ws_list = {}
    for (sample,chanlist) in meta_data.getAssociation().iteritems():
        ws_list[sample] = {}
        for (channel,proclist) in chanlist.iteritems():            
            #make the workspace we're going to use
            this_ws = RooWorkspace('%s-%s'%(channel,sample))
            initialize_workspace(this_ws)
            ws_list[sample][channel] = {}
            for (process,subproclist) in proclist.iteritems():
                print sample, channel, process
                if 'HToZG'in process:
                    ws_list[sample][channel][process] = []
                for (subproc,info) in subproclist.iteritems():
                    print '\tprocessing: %s'%subproc
                    input_file = info['input_file']
                    info['num_mc_events'] = -1
                    if input_file == '':
                        print '\t no input file found! Skipping!'
                        continue
                    if 'data' not in process:                        
                        print '\t mc input = %s'%input_file.split('/')[-1]
                        if 'HToZG' in subproc:                            
                            info['num_mc_events'] = \
                                  extract_higgs_data_in_categories(subproc,
                                                                   input_file,
                                                                   this_ws)
                            ws_list[sample][channel][process].append(subproc)
                        else:
                            info['num_mc_events'] = \
                                  extract_bkg_data_in_categories(subproc,
                                                                 input_file,
                                                                 this_ws)
                    else:
                        print '\t data input = %s'%input_file.split('/')[-1]
                        extract_data_in_categories(channel,input_file,this_ws)
            #end loop over processes and data
            fout_name = '%s_%s_master_workspace.root'%(channel,sample)
            fout = TFile.Open(fout_name,'recreate')
            fout.cd()
            this_ws.Write()
            fout.Close()
            gDirectory.cd(pwd)
            ws_list[sample][channel]['filename'] = fout_name        
    return ws_list

def build_category_workspaces(ws_list,metadata):
    pwd = gDirectory.GetPath()
    assoc = metadata.getAssociation()
    category_type = os.environ['hzgcategorytype']
    signal_type = os.environ['hzgsignalmodel']
    background_type = os.environ['hzgbkgmodel']
    categories = categories_map[category_type]['categories']
    catname = categories_map[category_type]['leafname']
    for (sample,chanlist) in ws_list.iteritems():
        for (channel,chaninfo) in chanlist.iteritems():
            master_file = TFile.Open(chaninfo['filename'],'read')
            gDirectory.cd(pwd)
            
            master_ws = master_file.Get('%s-%s'%(channel,sample))
            processes = chaninfo.keys()
            processes.remove('filename')

            #make workspaces for each category
            cat_workspaces = {}
            if not len(processes): continue
            
            for category in categories:                        
                cat_workspaces[category] =\
                                         RooWorkspace('%s-%s-cat%i'%(channel,
                                                                     sample,
                                                                     category),
                                                      '%s-%s-cat%i'%(channel,
                                                                     sample,
                                                                     category))
                cat_ws = cat_workspaces[category]
                
                cat_data = master_ws.data('%s_data'%channel)\
                           .reduce('%s == %i'%(catname,category))
                getattr(cat_ws,'import')(
                    cat_data,
                    RooFit.Rename('%s_data_cat%i'%(channel,category))
                    )

                cat_bkg = master_ws.data('mc_background_shape_data')\
                          .reduce('%s == %i'%(catname,category))
                getattr(cat_ws,'import')(cat_bkg)
                #make background model
                bkg_mdl = make_background_for_cat(category_type,
                                                  category,
                                                  background_type)
                for line in bkg_mdl:
                    cat_workspaces[category].factory(line)
                    
            #build background model + data workspaces
            for category in categories:                
                #write out category information
                cat_file = TFile.Open(
                    '%s_%s_data_category_%i_workspace.root'\
                    %(channel,sample,category),
                    'recreate')
                cat_workspaces[category].Write()
                cat_file.Close()  

            #generate signal workspaces for each mass point
            #and category
            acceptances = {}
            proc_workspaces = {}
            for process in processes:
                acceptance = TGraph(0)
                acceptance.SetName('g%s_acceff'%process)
                proc_workspaces[process] = RooWorkspace('%s-%s-%s'%(channel,
                                                                    sample,
                                                                    process))
                proc_name_short = process.split('To')[0]                
                subproc_list = chaninfo[process]
                print sample, channel, process
                masses = {}
                for subproc in subproc_list:
                    masses[assoc[sample][channel][process]\
                           [subproc]['mass']] = subproc
                for k,mass in enumerate(sorted(masses.keys())):
                    friendly_mass_name = str(mass).replace('.','p')
                    sig_fit_range = 'fit_%s'%friendly_mass_name
                    subproc = masses[mass]
                    #get acc*eff (calculated in master ws)
                    acceptance.SetPoint(
                        k,mass,
                        master_ws.var('%s_acceff'%subproc).getVal()
                        )
                    
                    for category in categories:
                        #get and import signal dataset for this cat
                        cat_sig_nopu = master_ws.data(
                            '%s_shape_data_puonly'%subproc
                            ).reduce('%s == %i'%(catname,category))
                        cat_sig = master_ws.data(
                            '%s_shape_data'%subproc
                            ).reduce('%s == %i'%(catname,category))
                        getattr(proc_workspaces[process],'import')(
                            cat_sig_nopu,
                            RooFit.Rename(
                            '%s_shape_data_puonly_cat%i'%(subproc,category)
                            )
                            )
                        getattr(proc_workspaces[process],'import')(
                            cat_sig,
                            RooFit.Rename('%s_shape_data_cat%i'%(subproc,
                                                                 category))
                            )
                        #build signal model
                        sig_model = make_signal_for_cat(category_type,
                                                        category,
                                                        signal_type,
                                                        mass,
                                                        proc_name_short)
                        for line in sig_model:
                            proc_workspaces[process].factory(line)
                            
                        #end dataset import and model defintions
                        #setup fitting range
                    proc_workspaces[process].var("Mzg").setRange(
                        sig_fit_range,
                        mass-30,mass+30
                        )
                #for each mass and category in this process
                #fit the signal datasets
                for k,mass in enumerate(sorted(masses.keys())):
                    subproc = masses[mass]
                    mname = str(mass).replace('.','p')
                    for category in categories:
                        #fit the signal model to pu-only weighted data
                        sig_model = 'signal_model_m%s_%s_cat%i'%(
                            mname,
                            proc_name_short,
                            category
                            )
                        sig_data  = proc_workspaces[process].data(
                            '%s_shape_data_puonly_cat%i'%(subproc,category)
                            )                        
                        sig_pdf   = proc_workspaces[process].pdf(sig_model)
                        # fit in range around peak to remove outliers that
                        # cause pdf to be zero
                        # get errors proportional to number of actual events
                        #sig_pdf.fitTo(sig_data,
                        #              RooFit.Range(sig_fit_range),
                        #              RooFit.SumW2Error(True)
                        #              )
                        #set all fitted parameters constant
                        params = sig_pdf.getParameters(sig_data)
                        params_it = params.iterator()
                        while ( not not params_it.Next() ):
                            params_it.setConstant(True)
                        del params
                
                #embed acceptance data in workspace
                getattr(proc_workspaces[process],'import')(acceptance)
            
            for process in processes:                
                #write out category information
                proc_file = TFile.Open(
                    '%s_%s_process_%s_workspace.root'\
                    %(channel,sample,process),
                    'recreate')
                proc_workspaces[process].Write()
                proc_file.Close()                          

            #close channel
            master_file.Close()

if __name__ == '__main__':
    analysis_root = os.environ['hzganalysisroot']
    analysis_name = os.environ['hzganalysisname']

    unblind = args.unblind
    
    dp = directory_prep(analysis_root,analysis_name)
    dp.build_inputs()

    mda = metadata_association(dp)
    
    master_ws_list = create_master_workspaces(mda)

    build_category_workspaces(master_ws_list,mda)

    #make the associated plots
    pgs = make_plot_groups(mda, higgs_masses = [125], higgs_prod = 'ggH')

    plot_path = os.path.join(analysis_root,analysis_name)
    plot_path = os.path.join(plot_path,'plots')

    if unblind:
        process_plot_groups(pgs,'test',[])
    else:
        process_plot_groups(pgs,plot_path,
                            ['(zg.M() > 150 || zg.M() < 120)']
    

