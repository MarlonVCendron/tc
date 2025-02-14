#!/bin/sh

. ./globalvars.sh

if [ "$SLEEPS" == "false" ]; then
    ARGS=""
else
    ARGS="--sleeps"
fi

make -C $DIR -j8 $BIN && mpirun -n $NP $DIR/$BIN \
	--dir $OUTDIR \
	--prefix rf1 --size 4096 --save \
	--monf ./data/rf1.pat \
	--recfile $RECFILE --xi $XI \
	--stim ./data/shapes1.pat \
	--wie 0.2 --wee 0.1 --wext 0.05 \
	--simtime $SIMTIME --tauf $TAUF --taud $TAUD \
	--intsparse $INTSPARSENESS \
	--extsparse 0.10 \
	--off 2.0 --on 1.0 \
	--beta $BETA --eta $ETA --bgrate $BGRATE --scale $SCALE --weight_a $WEIGHTA --alpha $ALPHA --delta 0.02 \
	$ARGS

