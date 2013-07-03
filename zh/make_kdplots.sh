## Generates both individual and combined kinematic discriminator plots

source jobid.sh
export jobid=$jobid8
mkdir -p results/$jobid/kd_plots/root_files/
root -lx make_ind_kdplots.cxx
cd results/$jobid/kd_plots/root_files/
rm output-COMB.root
hadd -f output-COMB.root output-*.root
cd ../../../../ 
root -lx kd_plot_combined.cxx
