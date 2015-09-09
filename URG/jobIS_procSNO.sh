#!/bin/bash
bsub -Is -q geyser -W 1:00 -n 16 -P ######## -J "rcmd" $SHELL < Rcmd_procSNO.sh; bkill -J "rcmd"

