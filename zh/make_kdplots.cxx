// Generates both individual and combined kinematic discriminator plots

.x kdplots.cxx("MMMT","m3","t");
.x kdplots.cxx("EEMT","m","t");
.x kdplots.cxx("MMET","e","t");
.x kdplots.cxx("EEET","e3","t");
.x kdplots.cxx("MMEM","e","m3");
.x kdplots.cxx("EEEM","e3","m");
.x kdplots.cxx("MMTT","t1","t2");
.x kdplots.cxx("EETT","t1","t2");

.x kd_plots_combined.cxx()
