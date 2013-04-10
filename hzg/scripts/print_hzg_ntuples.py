from UWHiggs.hzg.hzg_pattuples import allTuples

for tupleSet in allTuples.itervalues():
    for tupleName in tupleSet:
        if( tupleName != 'tupleName' and
            tupleName != 'tupleDate' and
            tupleName != 'tupleRoot' ):
            for tuple in tupleSet[tupleName]:
                print '%s%s/%s-%s'%(tupleSet['tupleRoot'],
                                    tuple,
                                    tupleSet['tupleDate'],
                                    tupleSet['tupleName'])
