#!/bin/sh
#==============================================================================
# Master script for ppalgs (ducting/icing/contrails) runs
# Expects one argument: UTC
# Path: ~forecast/atmos/harmonie/job/
#==============================================================================
#
# >> met.no/FoU	15.11.2011  Ole Vignes				... first version
# >> MET/IT  	21.04.2015  Martin Lilleeng Sætra	Adapted for ppalgs


# Trace commands
set -x

# Get required arguments
usage="Usage: $0 UTC"
utc=${1?$usage}

# Load necessary modules
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
# END Load necessary modules

# Set jobname. This is used in the start and end messages
jobname="am_ppalgs${utc} = Arome 2.5 ducting/icing/contrails ${utc}z"
maxnodes=1

# Include common HIRLAM shell functions
source $HOME/status/Status_Funcs.sh
###source $HOME/atmos/hirlam/job/Hirlam_Funcs.sh

Runstatus "start ...." "$jobname" $maxnodes
trap 'Runstatus "FAILED ..." "$jobname" $maxnodes; exit 159' 0

# Source local setup
[ -f /etc/profile ] && . /etc/profile

inpdir=$HOME/atmos/arome/input/
wrkdir=$HOME/run/ppalgs/
[ -d $wrkdir ] || mkdir $wrkdir
cd $wrkdir || { echo "Could not chdir $wrkdir"; exit; }
# Cleanup
rm -f ppalgs-*.nc
rm -f *msg

# Create start message
echo "$jobname start" > arome2_5_ppalgs.start_msg

# Get times
DTG=$( ls /prod/cooper/run/AM25_oper/fc*+000grib )
yyyy=$( echo $DTG | cut -c30-33 )
mm=$( echo $DTG | cut -c34-35 )
dd=$( echo $DTG | cut -c36-37 )
hh=$( echo $DTG | cut -c38-39 )

## Location of program
ppalgs=$HOME/src/ppalgs/bin/ppalgs
fimex=$HOME/bin/fimex-0.58

## Output file names
output_ml=arome2_5km_ppalgs${utc}.nc
output_pl=arome2_5km_ppalgs_plevels${utc}.nc

## Set times at which to do ppalgs
times=00,03,06,09,12,15,18,21,24,30,36,42,48,54,60,66
ntimes=$( perl -e 'split(",",shift);print $#_+1' $times )

# Fimex needs these
cp $inpdir/ducting_pl.ncml .
cp $inpdir/icing_index_pl.ncml .
cp $inpdir/contrails_pl.ncml .

for time in ${times//,/ }; do
  ## The following block is executed in parallel over all tasks
  { 
    gribfile=/prod/cooper/run/AM25_oper/fc$yyyy$mm$dd$hh+0${time}grib
    
    ## Copy nc template files
    tmp_nc_file_ml=ppalgs-$yyyy$mm$dd$hh-$time.nc
    cp $inpdir/ppalgs_template.nc $tmp_nc_file_ml

    ## Perform the various ducting/icing/contrails computations
    $ppalgs $gribfile $tmp_nc_file_ml --input-type=grib --output-type=nc -d

    $ppalgs $gribfile $tmp_nc_file_ml --input-type=grib --output-type=nc -i

    $ppalgs $gribfile $tmp_nc_file_ml --input-type=grib --output-type=nc -c
    
    ## Convert ml -> pl using fimex (FIXME: use linear_weak_extra when available)
    tmp_nc_file_pl=ppalgs-$yyyy$mm$dd$hh-${time}_plevels.nc
    
    $fimex --ncml.config=ducting_pl.ncml \
     --verticalInterpolate.type pressure --verticalInterpolate.method linear \
     --verticalInterpolate.level1 1000,925,850,800,700,500,400,300,250,200,150,100,70,50,30,10 \
     --input.file $tmp_nc_file_ml --input.type nc4 --extract.selectVariables ducting_ml \
     --output.file $tmp_nc_file_pl --output.type nc4
     
    $fimex --ncml.config=icing_index_pl.ncml \
     --verticalInterpolate.type pressure --verticalInterpolate.method linear \
     --verticalInterpolate.level1 1000,925,850,800,700,500,400,300,250,200,150,100,70,50,30,10 \
     --input.file $tmp_nc_file_ml --input.type nc4 --extract.selectVariables icing_index_ml \
     --output.file $tmp_nc_file_pl --output.type nc4
     
    $fimex --ncml.config=contrails_pl.ncml \
     --verticalInterpolate.type pressure --verticalInterpolate.method linear \
     --verticalInterpolate.level1 1000,925,850,800,700,500,400,300,250,200,150,100,70,50,30,10 \
     --input.file $tmp_nc_file_ml --input.type nc4 --extract.selectVariables contrails_ml \
     --output.file $tmp_nc_file_pl --output.type nc4
  } 1>task$task.log 2>&1 &

  task=$(( $task + 1 ))

done

## Wait for all tasks to finish
wait
ls -ltr

## Aggregate timesteps with fimex
agrdir_ml=$HOME/run/ppalgs/aggregate_ml/
rm -f $agrdir_ml/*
[ -d $agrdir_ml ] || mkdir $agrdir_ml
cp $inpdir/aggregate.ncml $agrdir_ml

agrdir_pl=$HOME/run/ppalgs/aggregate_pl/
rm -f $agrdir_pl/*
[ -d $agrdir_pl ] || mkdir $agrdir_pl
cp $inpdir/aggregate.ncml $agrdir_pl

## Copy individual timesteps to aggregation folder
for time in ${times//,/ }; do
  tmp_nc_file_ml=ppalgs-$yyyy$mm$dd$hh-$time.nc
  cp $tmp_nc_file_ml $agrdir_ml
  
  tmp_nc_file_pl=ppalgs-$yyyy$mm$dd$hh-${time}_plevels.nc
  cp $tmp_nc_file_pl $agrdir_pl
done

cd $agrdir_ml || { echo "Could not chdir $agrdir_ml"; exit; }
$fimex --input.file=aggregate.ncml --output.file=$output_ml --output.type=nc4
cp $output_ml ..

cd $agrdir_pl || { echo "Could not chdir $agrdir_pl"; exit; }
$fimex --input.file=aggregate.ncml --output.file=$output_pl --output.type=nc4
cp $output_pl ..

## Create msg-files for get-job
cd $wrkdir || { echo "Could not chdir $wrkdir"; exit; }
echo 1 > ${output_ml}msg
echo 1 > ${output_pl}msg
echo "$jobname end" > arome2_5_ppalgs.end_msg

Runstatus "end ......" "$jobname" $maxnodes
trap 0
exit
