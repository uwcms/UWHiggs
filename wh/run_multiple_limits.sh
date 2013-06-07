#! /bin/env bash

set -o nounset
set -o errexit

export RUN_OPTIMIZATION=True
#export megaworkers=9
#python optimizer.py
#touch WHAnalyzerBase.py
#./run.sh

for working_point in $(python optimizer.py); do
    echo $working_point
    wp_label=`echo $working_point | sed 's|:|_|g'`
    lead_id=`echo $working_point | awk -F':' '{print $1}'`
    sublead_id=`echo $working_point | awk -F':' '{print $2}'`
    LT_thr=`echo $working_point | awk -F':' '{print $3}'`
    if [ ! -f automatic_optimization/$wp_label'_'mmt_Asymptotic_expected_limit.json ]; then
	python WHPlotterMMT.py --dry-run --prefix=$working_point
	python WHPlotterEMT.py --dry-run --prefix=$working_point
	python WHPlotterEET.py --dry-run --prefix=$working_point
	./limits.sh
	echo
	cp results/$jobid/cards/shapes.root automatic_optimization/$wp_label'_'shapes.root
	for json_file in $(ls results/$jobid/cards/*/*.json); do #$jobid
    	    base=`basename $json_file`
    	    dest=automatic_optimization/$wp_label'_'$base
    	    cp $json_file $dest
    	    add_tag_to_json.py $dest -l leading_lepton_id -t $lead_id
    	    add_tag_to_json.py $dest -l subleading_lepton_id -t $sublead_id
    	    add_tag_to_json.py $dest -l LT_threshold -t $LT_thr
	    echo
	done
    else
	echo "working point already done, skipping"
    fi
done


#make ALL the shape files
#for each available WP compute the limit
plotdir= $plots
carddir= $results/cards
for working_point in $(python optimizer.py); do
    echo $working_point
    wp_label=`echo $working_point | sed 's|:|_|g'`
    lead_id=`echo $working_point | awk -F':' '{print $1}'`
    sublead_id=`echo $working_point | awk -F':' '{print $2}'`
    LT_thr=`echo $working_point | awk -F':' '{print $3}'`
    for channel in 'mmt' 'emt' 'eet'; do
	shapefile=$plotdir/$channel/$wp_label'_'$channel'_shapes_8TeV.root'
	if [ -f $shapefile]; then
	    limitdir=$carddir/$channel
	    mkdir -p $limitdir
	    rm -f $limitdir/shapes.root
	    add_stat_shapes.py $shapefile $limitdir/shapes.root --filter '#{prefix}/fakes' --prefix CMS_vhtt_8TeV_$channel'_fakeshape' > $limitdir/$channel'_shapes_8TeV_statshapes.txt'
	    horizontal-morphing.py --categories=$channel --samples='WH_hww{MASS}' --uncerts='' --masses='110,120,130,140' --step-size=5 --extrapolate=145,150,155,160 $limitdir/shapes.root
	    #copy cgs.conf
	    cp card_config/cgs.conf.$channel.8TeV $limitdir/cgs.conf
	    #copy and modify unc.conf
	    cp card_config/unc.conf.$channel.8TeV $limitdir/unc.conf

	    cd $limitdir
	    echo '' >> unc.conf
	    echo '# Stat shape uncertainties' >> unc.conf
	    cat $channel'_shapes_8TeV_statshapes.txt' | xargs -n 1 -I {} echo '{} shape' >> unc.conf
	    cd -

	    #copy and modify unc.vals
	    cp card_config/unc.vals.$channel.8TeV $limitdir/unc.vals
	    # Append all the stat shape types

	    cd $limitdir
	    echo '' >> unc.vals
	    echo '# Stat shape uncertainties' >> unc.vals
	    cat $channel'_shapes_8TeV_statshapes.txt' | xargs -n 1 -I {} echo '#{channel} fakes {} 1.0' >> unc.vals
	    cd -

	    python get_fake_systematic.py $limitdir/shapes.root $channel CMS_vhtt_$channel'_fakes_8TeV' >> $limitdir/unc.vals

	    cd $limitdir
	    mkdir -p 125/
	    create-datacard.py -i shapes.root -o 125/vhtt_$channel.txt 125 --categories $channel
	    compute_significance.sh 125 'YES'
	    cd -

	    harvest_limits.py $limitdir
	    add_tag_to_json.py $limitdir/*.json -l jobid -t $jobid
	    #add_tag_to_json.py $limitdir/*.json -l lumi -t #{get_lumi(samples_map[channel])}"
    	    add_tag_to_json.py $limitdir/*.json -l leading_lepton_id -t $lead_id
    	    add_tag_to_json.py $limitdir/*.json -l subleading_lepton_id -t $sublead_id
    	    add_tag_to_json.py $limitdir/*.json -l LT_threshold -t $LT_thr
	    ls $limitdir/*.json | xargs -n 1 -I {} cp {} automatic_optimization/$wp_label'_'$channel'_'{}
	else
	    echo "$working_point not available for channel $channel"
	fi
    done
done
