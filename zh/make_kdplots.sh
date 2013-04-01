## Generates both individual and combined kinematic discriminator plots

cd results/2012-10-24-8TeV-v1-Higgs/kd_plots/root_files/
rm output-COMB.root
hadd -f output-COMB.root output-*.root
cd ../../../../ 
root -l make_kdplots.cxx
