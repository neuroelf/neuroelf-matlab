#!/bin/sh

# Directives
#PBS -q batch0
#PBS -N PT_glm_searchlight
#PBS -W group_list=yetipsych
#PBS -l nodes=1,walltime=01:55:00,mem=6144mb
#PBS -M jw2661@columbia.edu
#PBS -m a
#PBS -V

# array spec
#PBS -t 2000-2001

# Set output and error directories
#PBS -o localhost:/vega/psych/users/jw2661/searchlight_test/
#PBS -e localhost:/vega/psych/users/jw2661/searchlight_test/
#PBS -j oe

# run MATLAB
matlab -r "cd /vega/psych/users/jw2661/searchlight_test ; part=$PBS_ARRAYID; outof=4000; PT_glm_searchlight_run; exit;" >> /vega/psych/users/jw2661/searchlight_test/searchlight_output

