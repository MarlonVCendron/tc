#!/bin/sh

OUTDIR="./data/sim_sleep_base"
RECFILE="./data/rf_discsR8.mtx"

# Create directory if it does not exist
mkdir -p $OUTDIR

DIR="."
BIN="sim_sleep"
NP=4
BETA="0.05"
ETA="1e-3"
ALPHA=4
WEIGHTA="0.0"
XI="0.5"
BGRATE="5"
SCALE="25"
INTSPARSENESS=0.05
SIMTIME=1800
# SIMTIME=2400
SLEEPS=true
TAUF="0.6" # Was 0.6 (600ms) here we are trying something shorter
TAUD="0.15" # Was 0.2 (200ms) in original sim, but this seemingly works better

