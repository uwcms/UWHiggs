'''This is garbage, but I have short time'''

import itertools
lep_id = [
    'h2taucuts',
    'h2taucuts020',
    #'idiso02'
    ]

LT_cut = [60, 70, 80, 90]

loop = '''
#SET cuts
export leadLeptonIsoTag=%s
export subleadLeptonIsoTag=%s
export ltThreshold=%s
#Runs analysis
./run.sh
./limits.sh

#copies limits
cp results/$jobid8/cards/??t/*.json results/$jobid8/cards/.
ls results/$jobid8/cards/ | grep json | xargs -n1 -I{} cp results/$jobid8/cards/{} automatic_optimization/$leadLeptonIsoTag'_'$subleadLeptonIsoTag'_'LT$ltThreshold'_'{}
cp results/$jobid8/cards/shapes.root automatic_optimization/$leadLeptonIsoTag'_'$subleadLeptonIsoTag'_'LT$ltThreshold'_'shapes.root
rm -rf results/$jobid8/WHAnalyze*
rm -rf results/$jobid8/plots


'''

header = '''
#!/bin/bash
# Run all optimization

set -o nounset
set -o errexit

source jobid.sh
mkdir -p automatic_optimization

'''

print header

for cuts in itertools.product(lep_id, lep_id, LT_cut):
    print loop % cuts
    
