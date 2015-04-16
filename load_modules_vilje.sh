#! /bin/sh
source /etc/profile.d/modules.sh
module purge
module load intelcomp/13.0.1

module load gcc/4.8.2
module load proj/4.7.0

# some modules loads older versions of icc, we want icc 14
module unload intelcomp/13.0.1
module unload intelcomp/12.0.5.220
module load intelcomp/14.0.1

module load boost/1.55.0 #was 1.48.0
module load mpt/2.09
module load netcdf/4.3.1
module load udunits/2.2.4 #was 2.1.24
module load openjpeg/1.5.0

# some modules loads older versions of icc, we want icc 14
module unload intelcomp/12.0.5.220
module load intelcomp/14.0.1

