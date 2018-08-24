#!/bin/sh
# The SGE_TASK variable specifies the data sets, as in $HOME/test.qsub

cd /home/knight/emodim/scripts/

# don't change this variable
# used by the submit script to define which data sets to analyze
DATASET="${SGE_TASK}"

# define function
FUNCTION='SBJ01_import_data'

# set up matlab function call
func_call="${FUNCTION}('${DATASET}', '${pipeline_id}')"

# define commands to execute via SGE
echo ${DATASET}
echo ${func_call}
echo $$
echo ${func_call} > NotBackedUp/tmpSGE/${FUNCTION}_${DATASET}.m
time matlab -nodesktop -nosplash -nodisplay < NotBackedUp/tmpSGE/${FUNCTION}_${DATASET}.m
rm NotBackedUp/tmpSGE/${FUNCTION}_${DATASET}.m
