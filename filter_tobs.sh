#!/bin/sh

GRAPHFILE=$1
GRAPH=$(basename $1)
OUTPUTPREFIX="tobs_"
OUTPUTDIR="/toshiba/pinchedplanar/tobs"

BLOCKSIZE=${2:-512k}

SPINDLE="/toshiba/pinchedplanar/nauty27r1/spindle"
SPINDLEFLAGS="-t"

PARALLEL="/usr/local/bin/parallel"

#with --pipepart
#${PARALLEL} --progress --pipepart --blocksize ${BLOCKSIZE} -a ${GRAPHFILE} ${SPINDLE} ${SPINDLEFLAGS} > ${OUTPUTDIR}/${OUTPUTPREFIX}${GRAPH}

#without --pipepart, which seems to choke on large files
cat ${GRAPHFILE} | ${PARALLEL} --pipe --blocksize ${BLOCKSIZE} ${SPINDLE} ${SPINDLEFLAGS} > ${OUTPUTDIR}/${OUTPUTPREFIX}${GRAPH}
