#!/bin/bash
bsub -Is -q geyser -W 1:00 -n 16 -P P48500028 -J "rcmd" $SHELL < Rcmd_procSNO.sh; bkill -J "rcmd"
