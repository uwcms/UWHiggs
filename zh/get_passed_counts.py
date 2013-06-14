#!/usr/bin/python
import ROOT
 

 
def print_table(channel, weighted=True):
# print a tex table for a given channel
    data = ROOT.TFile("results/2013-04-30-8TeV-v1-ZH_light/ZHAnalyze%s/data.root" % channel, "READ")
    zz = ROOT.TFile("results/2013-04-30-8TeV-v1-ZH_light/ZHAnalyze%s/ZZJetsTo4L_pythia.root" % channel, "READ")
    zjets = ROOT.TFile("results/2013-04-30-8TeV-v1-ZH_light/ZHAnalyze%s/Zjets_M50.root" % channel, "READ")
    wz = ROOT.TFile("results/2013-04-30-8TeV-v1-ZH_light/ZHAnalyze%s/WZJetsTo3LNu_pythia.root" % channel, "READ")
    column_values = []
    background_estimates = []
    for my_file, folder, weight_folder, lumi_weight in [(zz,'All_Passed','', 19.250/25866.365), (data,'Leg3Failed_Leg4Failed', 'all_weights_applied', 1.0), (data,'Leg3Failed','leg3_weight', 1.0), (data, 'Leg4Failed', 'leg4_weight', 1.0), (zjets,'All_Passed','', 19.250/4.320), (wz,'All_Passed','', 19.250/1899.212), (data, 'All_Passed', '', 1.0)]:
        path_to_histo = 'os/' + folder
        histo_unweighted = my_file.Get(path_to_histo + '/kinematicDiscriminant1')
        if weighted:
            path_to_histo += '/' + weight_folder
        # we can get any existing histo here, we just need the integral
        
        histo = my_file.Get(path_to_histo + '/kinematicDiscriminant1')
        histo.Scale(lumi_weight)
        if weight_folder:
            s = str(round(histo.Integral(),3)) + ' (' + str(round(histo_unweighted.Integral(),3)) + ')'
            background_estimates.append(str(round(histo.Integral(),3)))
            column_values.append(s)
        else: column_values.append(str(round(histo.Integral(),3)))
        #print str(histo.Integral())
        # now print them out, TeX style
        #line = " & ".join(str(column_values)) + "\\\\"
    #return line;
    line = channel 
    line += ' & ' + column_values[0]
        
    # estimate background contribution
    fr_contribution = float(background_estimates[1]) + float(background_estimates[2]) - float(background_estimates[0])
    line += ' & ' + str(fr_contribution)

    for i in [3,4, 5]:
        line += ' & ' + column_values[i]
    
    line += ' & ' + str(float(column_values[0]) + fr_contribution)
    line += ' & ' + column_values[6]

    return line
print '''
\documentclass[final,letterpaper,twoside,12pt]{article}
\usepackage{chngpage}
\usepackage{tabularx}
\\begin{document}
\\begin{table}[htbp]
\\begin{adjustwidth}{-4em}{-4em}
\centering
\\begin{tabular}{|c|c|c|c|c|c|c|c|c|} \\hline
channel & ZZ & fakerate & PPPF (fr) & Zjets & WZ & ZZ + fakerate & observed data \\\\ 
\\hline \\hline '''

print print_table('MMTT') + '  \\\\ \\hline'
print print_table('EETT') + '  \\\\ \\hline'
print print_table('EEET') + '  \\\\ \\hline'
print print_table('MMET') + '  \\\\ \\hline'
print print_table('MMMT') + '  \\\\ \\hline'
print print_table('EEMT') + '  \\\\ \\hline'
print print_table('EEEM') + '  \\\\ \\hline'
print print_table('MMEM') + '  \\\\ \\hline'
print print_table('COMB') + '  \\\\ \\hline'

print '''
\end{tabular}
\caption {Background counts in each channel, }
\label{tab:cqdata0}
\end{adjustwidth}\end{table}

\end{document}
'''
