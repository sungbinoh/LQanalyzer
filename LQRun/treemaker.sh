#!/bin/sh
### sets all configurable variables to defaul values
source functions.sh
###
###########################################################
## CONFIGURE JOB ####
###########################################################
#
#  THIS SECTION IS FOR USER:
#
############################################################

cycle="SKTreeMaker"

stream="egamma"
njobs=1
data_lumi="AtoD"
loglevel="INFO"
logstep=1000
#declare -a input_samples=("A" "B" "C" "D" "DY10to50" "DY50plus" "ttbar" "Wjets" "WZ" "ZZ" "W-W-" "W+W+" "Wgamma")
#declare -a input_samples=("A" "B" "C" "D")

#declare -a input_samples=("DY10to50" "DY50plus" "ttbar" "Wjets" "WZ" "ZZ" "W-W-" "W+W+" "Wgamma")    
declare -a input_samples=("W-W-")

### submit this configured job (uses bin/submit.sh)
source submit.sh
