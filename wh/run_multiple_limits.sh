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