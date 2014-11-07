#! /bin/env bash

#$1 --> category (0, 1, 2) --> 0 jets, 1 jet, vbf
#$2 --> files extensions produced by plotting (used inoptimization)
category=$1
if [[ -z "$2" ]]; 
then
    file_extension=$2
else
    file_extension=$1
fi
input_dir=plots/$jobid/lfvet/
output_dir=limits/$jobid/$category/
mkdir -p output_dir

#copy basic config
cp card_config/cgs.$category.conf $output_dir/cgs.config
cp card_config/unc.$category.conf $output_dir/unc.config
cp card_config/unc.$category.vals $output_dir/unc.vals

#copy shape file
cp $input_dir/shapes.$file_extension.root $output_dir/shapes.root

#add dynamically assigned variables
cat $input_dir/unc.$file_extension.conf >> $output_dir/unc.config
cat $input_dir/unc.$file_extension.vals >> $output_dir/unc.vals

#add bbb errors 
if [ $category == '0' ]; then
    command='gg0etau:fakes'
elif [ $category == '1' ]; then
    command='boostetau:fakes'    
else
    command='vbfetau:fakes,ttbar,singlet,ztautau,SMVBF126,WWVBF126'
fi

#TODO: normalize or not?
#will it work with input and output in the same place?
add_bbb_errors.py --normalize -i $output_dir -o $output_dir $command

#build datacard
pushd $output_dir
mkdir 126
create-datacard.py -i shapes.root -o 126/datacard.txt
popd

exit 0
