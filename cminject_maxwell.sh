#!/bin/bash

#SBATCH --partition=cfel
#SBATCH --time=1:00:00
#SBATCH --nodes=1
#SBATCH --job-name cminject
#If you want to receive mails about your jobs, add the following
#(add "SBATCH" as above and replace "someemail" by your email address):
# --mail-type ALL
# --mail-user someemail

source /etc/profile.d/modules.sh
module load maxwell julia

CMINJECT_PATH="/gpfs/cfel/group/cmi/labs/CMInject.jl/"
CMINJECT_SH="${CMINJECT_PATH}/cminject.sh"

$CMINJECT_SH $@
