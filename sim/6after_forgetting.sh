#!/bin/sh

. ./globalvars.sh

if [ "$SLEEPS" == "false" ]; then
    ARGS=""
else
    ARGS="--sleeps"
fi

make -C $DIR -j8 $BIN && mpirun -n $NP $DIR/$BIN \
	--dir $OUTDIR \
	--load $OUTDIR/rf5 \
	--prefix rf6 --size 4096 --save \
	--monf ./data/rf1.pat \
	--stim ./data/shapes1.pat \
	--wie 0.2 --wee 0.1 --wext 0.1 \
	--simtime 2400 --tauf $TAUF --taud $TAUD \
	--intsparse $INTSPARSENESS \
	--extsparse 0.10 \
	--off 10.0 --on 0.1 \
	--beta $BETA --eta $ETA --bgrate $BGRATE --scale $SCALE --weight_a $WEIGHTA --alpha $ALPHA --delta 0.02  \
	$ARGS
