#!/bin/bash
bsub -Is -q geyser -W 6:00 -n 16 -P P48500028 -J "rcmd" $SHELL < Rcmd_reads.sh; bkill -J "rcmd"

