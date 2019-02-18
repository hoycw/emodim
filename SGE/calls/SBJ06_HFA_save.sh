#!/bin/sh
# The SGE_TASK variable specifies the data sets, as in $HOME/test.qsub

cd /home/knight/emodim/scripts/

# don't change this variable
# used by the submit script to define which data sets to analyze
SBJ="${SGE_TASK}"

# define function
FUNCTION='SBJ06_HFA_save'

# set up matlab function call
func_call="${FUNCTION}('${SBJ}', '${pipeline_id}', '${an_id}')"

# define commands to execute via SGE
echo ${SBJ}
echo ${func_call}
echo $$
echo ${func_call} > NotBackedUp/tmpSGE/${FUNCTION}_${SBJ}.m
echo "wrote tmp.m inside SBJ06.sh call"
time matlab -nodesktop -nosplash -nodisplay < NotBackedUp/tmpSGE/${FUNCTION}_${SBJ}.m
echo "sent to MATLAB tmp.m inside SBJ06.sh call"
rm NotBackedUp/tmpSGE/${FUNCTION}_${SBJ}.m
