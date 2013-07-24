#! /bin/bash

# ARGS
# 1 --> analyzer
# 2 --> sample name
# 3 --> input file list
# 4 --> output file name
echo $@

analyzer=$1
sample=$2
inputs=$3
output=$4

pushd $CMSSW_BASE
source $CMSSW_BASE/src/UWHiggs/environment.sh
popd

echo 'Setting RUN_OPTIMIZATION'
export RUN_OPTIMIZATION=True
echo 'Setting megaworkers to 1'
export megaworkers=1
export megachain=`echo $inputs | sed 's|,|\n|g' | wc -l`
echo 'chain setted to '$megachain
shaSum=`echo $inputs | sha1sum | awk '{print $1}'`
input_txt=$shaSum.txt
echo 'input txt file will be: '$input_txt
echo $inputs | sed 's|,|\n|g' | awk '{print "root://cmsxrootd.hep.wisc.edu/"$1}'> $input_txt
target=$sample-$shaSum.root
echo 'Setting target to '$target
export megatarget=$target

echo
echo 'inputs'
echo
echo $inputs
echo
echo 'input file'
echo
cat $input_txt
echo
echo 'current PWD'
echo
echo $PWD

sandbox=$PWD

pushd $CMSSW_BASE/src/UWHiggs/wh/
echo time mega $analyzer $sandbox/$input_txt $sandbox/$target --workers 1 --chain $megachain
time mega $analyzer $sandbox/$input_txt $sandbox/$target --workers 1 --chain $megachain
popd

mv $target $output
rm $input_txt
