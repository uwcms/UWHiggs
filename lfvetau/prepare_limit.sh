#! /bin/env bash

#set -o nounset
set -o errexit

#$1 --> category (0, 1, 2) --> 0 jets, 1 jet, vbf
#$2 --> files extensions produced by plotting (used inoptimization)
category=$1
file_extension=$2
if [ -z "$file_extension" ]; 
then
    file_extension=$category
fi
input_dir=plots/$jobid/lfvet/
output_dir=limits/$jobid/$category/
mkdir -p $output_dir

#copy basic config
cp card_config/cgs.$category.conf $output_dir/cgs.conf
cp card_config/unc.$category.conf $output_dir/unc.conf
cp card_config/unc.$category.vals $output_dir/unc.vals

#copy shape file
cp $input_dir/shapes.$file_extension.root $output_dir/shapes.root

#add dynamically assigned variables
cat $input_dir/unc.$file_extension.conf >> $output_dir/unc.conf
cat $input_dir/unc.$file_extension.vals >> $output_dir/unc.vals

#add bbb errors 
if [ $category == '0' ]; then
    command='gg0etau:fakes'    
    category_name='gg0etau'
    samples='fakes'
elif [ $category == '1' ]; then
    command='boostetau:fakes'    
    category_name='boostetau'
    samples='fakes'
else
    command='vbfetau:fakes,ttbar,singlet,ztautau,SMVBF126,WWVBF126'
    category_name='vbfetau'
    samples='fakes ttbar singlet ztautau SMVBF126'
fi

#TODO: normalize or not?
#will it work with input and output in the same place?
echo >> $output_dir/unc.conf
echo >> $output_dir/unc.vals    
echo >> $output_dir/unc.conf
echo >> $output_dir/unc.vals    
echo "#bbb uncertainties" >> $output_dir/unc.conf
echo "#bbb uncertainties" >> $output_dir/unc.vals    
for sample in $samples; do
    echo "adding bbb errors for $sample"
    echo add_stat_shapes.py $output_dir/shapes.root $output_dir/shapes.root --normalize --filter $category_name/$sample --prefix $category_name
    bbb_added=$(add_stat_shapes.py $output_dir/shapes.root $output_dir/shapes.root --normalize --filter $category_name/$sample --prefix $category_name | grep _bin_)
    echo >> $output_dir/unc.conf
    echo >> $output_dir/unc.vals    
    for unc in $bbb_added; do
	echo $unc shape >> $output_dir/unc.conf
	echo $category_name $sample $unc 1.0 >>  $output_dir/unc.vals    
    done
done

#build datacard
pushd $output_dir
mkdir -p 126
pushd 126
rm -rf shapes.root
ln -s ../shapes.root
popd
create-datacard.py -i shapes.root -o 126/datacard.txt
popd

exit 0
