import os
import re
import glob
import ROOT
from math import sqrt


def fakerate_central_histogram(nbin, xmin, xmax ):
    
    
    
    frfit_dir = os.path.join('results', os.environ['jobid'], 'fakerate_fits')+'/'
    
    
    myfile = ROOT.TFile (frfit_dir+'t_os_tLoose_tTigh_tAbsEta.root')
    myWS = myfile.Get('fit_efficiency')
    efficiency_x = myWS.var('x')
    efficiency_func = myWS.function("efficiency") 
    frhisto = efficiency_func.createHistogram("efficiency_x", efficiency_x)
    results = myWS.genobj("fitresult_chi2_efficiency_xy_data")
    parameters=results.floatParsFinal() 
    weight = []

    for n in range(1, nbin) :
        
        eta = (n -0.5)* (xmax - xmin)/nbin  
        etamin = (n-1)* (xmax - xmin)/nbin  
    
        weight.append((parameters[0].getVal()+parameters[1].getVal()*eta+parameters[2].getVal()*eta*eta, etamin))
    
    return weight

def fakerate_p1s_histogram(nbin, xmin, xmax) :
    
    frfit_dir = os.path.join('results', os.environ['jobid'], 'fakerate_fits')+'/'
    
    myfile = ROOT.TFile (frfit_dir+'t_os_tLoose_tTigh_tAbsEta.root')
    myWS = myfile.Get('fit_efficiency')
    efficiency_x = myWS.var('x')
    efficiency_func = myWS.function("efficiency") 
    frhisto = efficiency_func.createHistogram("efficiency_x", efficiency_x)
    results = myWS.genobj("fitresult_chi2_efficiency_xy_data")
    parameters=results.floatParsFinal() 
    covMatrix = results.covarianceMatrix()
    weight = []

    for n in range(1, nbin) :
        
        eta = (n -0.5)* (xmax - xmin)/nbin  
        etamin = (n-1)* (xmax - xmin)/nbin  
        
        err2= covMatrix(0,0) + covMatrix(1,1) * eta*eta + covMatrix(2,2) *eta*eta*eta*eta + 2*eta*covMatrix(0,1) + 2*eta*eta*covMatrix(0,2) + 2*eta*eta*eta*covMatrix(1,2) 
        myweight = (parameters[0].getVal()+parameters[1].getVal()*eta+parameters[2].getVal()*eta*eta)+ sqrt(err2)
        
        weight.append((myweight, etamin))
 
    return weight

def fakerate_m1s_histogram(nbin, xmin, xmax) :
    
    frfit_dir = os.path.join('results', os.environ['jobid'], 'fakerate_fits')+'/'
    
    myfile = ROOT.TFile (frfit_dir+'t_os_tLoose_tTigh_tAbsEta.root')
    myWS = myfile.Get('fit_efficiency')
    efficiency_x = myWS.var('x')
    efficiency_func = myWS.function("efficiency") 
    frhisto = efficiency_func.createHistogram("efficiency_x", efficiency_x)
    results = myWS.genobj("fitresult_chi2_efficiency_xy_data")
    parameters=results.floatParsFinal() 
    covMatrix = results.covarianceMatrix()

    weight = []

    for n in range(1, nbin) :
        
        eta = (n -0.5)* (xmax - xmin)/nbin  
        etamin = (n-1)* (xmax - xmin)/nbin  
    
        err2= covMatrix(0,0) + covMatrix(1,1) * eta*eta + covMatrix(2,2) *eta*eta*eta*eta + 2*eta*covMatrix(0,1) + 2*eta*eta*covMatrix(0,2) + 2*eta*eta*eta*covMatrix(1,2) 
        myweight = (parameters[0].getVal()+parameters[1].getVal()*eta+parameters[2].getVal()*eta*eta) - sqrt(err2)
        
        weight.append((myweight, etamin))

    return weight

