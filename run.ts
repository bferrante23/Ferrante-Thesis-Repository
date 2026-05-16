#!/bin/tcsh -f

rm -f log.*
set NI="/usr/bin/nice --adjustment=19 "

set dx="4"
set dy="4"

$NI $PY $CA/get_ts_diff.py -inp dates.run.1 -ptsfile ts_ra.txt1 -ftype unwrap_ll_flat.grd -dx $dx -dy $dy -tsdir . >& log.cand.tst 

