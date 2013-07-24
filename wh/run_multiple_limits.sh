#! /bin/env bash

set -o nounset
set -o errexit

export RUN_OPTIMIZATION=True
#export megaworkers=9
#python optimizer.py
#touch WHAnalyzerBase.py
#./run.sh

#make ALL the shape files
#mmt_working_points=`python optimizer.py | grep eid13Tight | grep -E ':80$' | tr '\n' ','`
#mmt_working_points=`echo "${mmt_working_points%?}"` #remove trailing comma
#echo python WHPlotterMMT.py --dry-run --prefixes=$mmt_working_points
#echo
#echo
#python WHPlotterMMT.py --dry-run --prefixes=$mmt_working_points

#emt_working_points=`python optimizer.py | grep -E ':80$' | tr '\n' ','`
#emt_working_points=`echo "${emt_working_points%?}"` #remove trailing comma
#echo python WHPlotterEMT.py --dry-run --prefixes=$emt_working_points
#echo
#echo
#python WHPlotterEMT.py --dry-run --prefixes=$emt_working_points

#eet_working_points=`python optimizer.py | tr '\n' ','`
#eet_working_points=`echo "${eet_working_points%?}"` #remove trailing comma
#echo python WHPlotterEET.py --dry-run --prefixes=$eet_working_points
#echo
#echo
#python WHPlotterEET.py --dry-run --prefixes=$eet_working_points

#for each available WP compute the limit
plotdir=$plots
carddir=$results/cards
cwd=$PWD
for working_point in $(python optimizer.py); do
    echo $working_point
    wp_label=`echo $working_point | sed 's|:|_|g'`
    #lead_id=`echo $working_point | awk -F':' '{print $1}'`
    #sublead_id=`echo $working_point | awk -F':' '{print $2}'`
    #LT_thr=`echo $working_point | awk -F':' '{print $3}'`
    #tauId=`echo $working_point | awk -F':' '{print $4}'`
    #tauPt=`echo $working_point | awk -F':' '{print $5}'`
    #charge_fakes=`echo $working_point | awk -F':' '{print $6}'`
    for channel in 'mmt' 'emt'; do
	shapefile=$plotdir/$channel/$wp_label'_'$channel'_shapes_8TeV.root'
	if [ -f $shapefile ]; then
	    if [ ! -f automatic_optimization/$wp_label'_'$channel'_'$channel'_ProfileLikelihood_significance_limit.json' ]; then
	    	limitdir=$carddir/$channel
	    	mkdir -p $limitdir
	    	rm -f $limitdir/shapes.root
	    	echo -e "\n\nStat shapes + Merphing\n\n"
	    	add_stat_shapes.py $shapefile $limitdir/shapes.root --filter '#{prefix}/fakes' --prefix CMS_vhtt_8TeV_$channel'_fakeshape' > $limitdir/$channel'_shapes_8TeV_statshapes.txt'
	    	horizontal-morphing.py --categories=$channel --samples='WH_hww{MASS}' --uncerts='' --masses='110,120,130,140' --step-size=5 --extrapolate=145,150,155,160 $limitdir/shapes.root

	    	echo -e "\n\nConfiguring input files\n\n"
	    	#copy cgs.conf
	    	cp card_config/cgs.conf.$channel.8TeV $limitdir/cgs.conf
	    	#copy and modify unc.conf
	    	cp card_config/unc.conf.$channel.8TeV $limitdir/unc.conf

	    	cd $limitdir
	    	echo '' >> unc.conf
	    	echo '# Stat shape uncertainties' >> unc.conf
	    	cat $channel'_shapes_8TeV_statshapes.txt' | xargs -n 1 -I {} echo '{} shape' >> unc.conf
	    	cd $cwd

	    	#copy and modify unc.vals
	    	cp card_config/unc.vals.$channel.8TeV $limitdir/unc.vals
	    	# Append all the stat shape types

	    	cd $limitdir
	    	echo '' >> unc.vals
	    	echo '# Stat shape uncertainties' >> unc.vals
	    	cat $channel'_shapes_8TeV_statshapes.txt' | xargs -n 1 -I {} echo '#{channel} fakes {} 1.0' >> unc.vals
	    	cd $cwd

	    	echo -e "\n\nGet Fake Sys\n\n"
	    	python get_fake_systematic.py $limitdir/shapes.root $channel CMS_vhtt_$channel'_fakes_8TeV' >> $limitdir/unc.vals

	    	cd $limitdir
	    	mkdir -p 125/
	    	cd 125/
	    	echo -e "\n\nLinking\n\n"
	    	rm -f shapes.root
	    	ln -s ../shapes.root
	    	cd ..
	    	echo -e "\n\nCreate datacard\n\n"
	    	create-datacard.py -i shapes.root -o 125/vhtt_$channel.txt 125 --categories $channel
	    	echo -e "\n\nCompute significance\n\n"
	    	compute_significance.sh 125 'YES'
	    	cd $cwd

	    	harvest_limits.py $limitdir
	    	#add_tag_to_json.py $limitdir/*.json -l jobid -t $jobid
	    	##add_tag_to_json.py $limitdir/*.json -l lumi -t #{get_lumi(samples_map[channel])}"
    	    	#add_tag_to_json.py $limitdir/*.json -l leading_lepton_id -t $lead_id
    	    	#add_tag_to_json.py $limitdir/*.json -l subleading_lepton_id -t $sublead_id
    	    	#add_tag_to_json.py $limitdir/*.json -l LT_threshold -t $LT_thr
    	    	#add_tag_to_json.py $limitdir/*.json -l hiFakeTauID -t $tauId
    	    	#add_tag_to_json.py $limitdir/*.json -l tauPt -t $tauPt
    	    	#add_tag_to_json.py $limitdir/*.json -l charge_fakes_rejection_power -t $charge_fakes
	    	for json_file in $(ls $limitdir/*.json); do
	    		name=`basename $json_file`
	    		cp -v $json_file automatic_optimization/$wp_label'_'$channel'_'$name
	    	done
	    else
		echo "skipping $channel $working_point, it is already done"
	    fi
	else
	    echo "$working_point not available for channel $channel"
	fi
    done
done
