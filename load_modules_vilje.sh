#! /bin/sh
source /etc/profile.d/modules.sh

# proj/4.8.0 requires intelcomp/13.0.1
module load intelcomp/13.0.1

module load gcc/4.8.2
module load proj/4.8.0
module load boost/1.48.0

# boost/1.48.0 loads mpt/2.06, but netcdf/4.3.1 requires mpt/2.09
module unload mpt/2.06
module load mpt/2.09

# some modules loads older versions of icc, we want icc 14
module unload intelcomp/13.0.1
module unload intelcomp/12.0.5.220
module load intelcomp/14.0.1

module load netcdf/4.3.1
module load udunits/2.1.24
module load openjpeg/1.5.0

# some modules loads older versions of icc, we want icc 14
module unload intelcomp/12.0.5.220
module load intelcomp/14.0.1
