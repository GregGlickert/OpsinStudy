#!/bin/bash
#SBATCH --output=/dev/null # dont want an output from this guy so bye bye
#SBATCH --error=/dev/null

sbatch submit_stable.sh
sbatch submit_unstable.sh

