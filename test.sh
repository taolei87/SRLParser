#!/bin/sh

export LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:."

type="lab"
args=$1
runid=$2
shift
shift

 ln -s $args.train data/$args.train.$runid
 ln -s $args.test data/$args.test.$runid

	java -classpath "bin:lib/trove.jar" -Xmx20000m parser.DependencyParser model-file:runs/$args.model.$runid unimap-file:unimap/$args.uni.map test test-file:data/$args.test.$runid output-file:runs/$args.$runid.test.out $@ 

perl eval09.pl -p -q -g data/$args.test.$runid -s runs/$args.$runid.test.out 2> /dev/null


    rm data/$args.train.$runid
    rm data/$args.test.$runid
	#rm data/$args/$args.$type.model

