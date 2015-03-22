#!/bin/sh

export LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:."

type="lab"
args=$1
runid=$2
shift
shift

 ln -s $args.train data/$args.train.$runid
 ln -s $args.dev data/$args.dev.$runid

	java -classpath "bin:lib/trove.jar" -Xmx30000m parser.DependencyParser model-file:runs/$args.model.$runid train train-file:data/$args.train.$runid unimap-file:unimap/$args.uni.map test test-file:data/$args.dev.$runid output-file:runs/$args.$runid.out $@ > runs/$args.$runid.log

    rm data/$args.train.$runid
    rm data/$args.dev.$runid
	#rm data/$args/$args.$type.model

