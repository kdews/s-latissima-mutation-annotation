## Add these lines to your script before activating a conda env
# cond=~/.conda_for_sbatch.sh
# if [[ -f "$cond" ]]
# then
#   source "$cond"
# else
#   echo "Error on source of $cond"
#   exit 1
# fi
test="$('/home1/kdeweese/miniforge3/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
eval "$test"
