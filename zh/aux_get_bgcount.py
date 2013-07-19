#!/usr/bin/python
import ROOT
import math 

 
def print_table(channel, weighted=True):
# print a tex table for a given channel
    my_file = ROOT.TFile("results/2013-06-29-8TeV-v1-ZH_light/ZHAnalyze%s/data.root" % channel, "READ")
    column_values = []
    background_estimates = []
    errors = []
    for column, folder, weight_folder in [('PPPP','All_Passed',''), ('PPFF', 'Leg3Failed_Leg4Failed', 'all_weights_applied') , ('PPFP', 'Leg3Failed', 'leg3_weight'), ('PPPF', 'Leg4Failed', 'leg4_weight')]:
        path_to_histo = 'os/' + folder
        histo_unweighted = my_file.Get(path_to_histo + '/kinematicDiscriminant1')
        if weighted:
            path_to_histo += '/' + weight_folder
        # we can get any existing histo here, we just need the integral
        #print path_to_histo
        
        histo = my_file.Get(path_to_histo + '/kinematicDiscriminant1')
        if weight_folder:
            s = str(round(histo.Integral(),3)) + ' (' + str(round(histo_unweighted.Integral(),3)) + ')'
            #s = str(round(histo.Integral(),2)) + ' $\\pm$ ' 
            errors.append(round(1/math.sqrt(histo_unweighted.Integral()) * histo.Integral(), 2))
            #errors.append(0.1)
            #s += str(round(error,2))
            background_estimates.append(str(round(histo.Integral(),2)))
            column_values.append(s)
        else: column_values.append(str(round(histo.Integral(),2)))
        #print str(histo.Integral())
        # now print them out, TeX style
        #line = " & ".join(str(column_values)) + "\\\\"
    #return line;
    line = channel 
    for count in column_values:
        line += ' & ' + count
    
    # estimate background contribution
    fr_error = round(math.sqrt(errors[0] + errors[1] + errors[2]), 2)
    line += ' & ' + str(float(background_estimates[1]) + float(background_estimates[2]) - float(background_estimates[0]))
    line += ' $\\pm$ ' + str(fr_error)

    return line
print '''
\documentclass[final,letterpaper,twoside,12pt]{article}
\usepackage{chngpage}
\usepackage{tabularx}
\\begin{document}
\\begin{table}[htbp]
\\begin{adjustwidth}{-4em}{-4em}
\centering
\\begin{tabular}{|c|c|c|c|c|c|} \\hline
channel & PPPP & PPFF & PPFP & PPPF & Estimated Background \\\\ 
\\hline \\hline '''

print print_table('MMTT') + '  \\\\ \\hline'
print print_table('EETT') + '  \\\\ \\hline'
print print_table('MMMT') + '  \\\\ \\hline'
print print_table('EEMT') + '  \\\\ \\hline'
print print_table('MMET') + '  \\\\ \\hline'
print print_table('EEET') + '  \\\\ \\hline'
print print_table('MMEM') + '  \\\\ \\hline'
print print_table('EEEM') + '  \\\\ \\hline'
print print_table('COMB') + '  \\\\ \\hline'

print '''
\end{tabular}
\caption {Background counts in each channel, }
\label{tab:cqdata0}
\end{adjustwidth}\end{table}

\end{document}
'''
