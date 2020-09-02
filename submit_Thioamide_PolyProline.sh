#!/bin/bash
#$ -q all.q
#$ -cwd
#$ -N Thioamide_PolyProline
#$ -S /bin/bash
#$ -o Thioamide_PolyProline.log
#$ -e Thioamide_PolyProline.error
#$ -l h_rt=900:00:00
#$ -t 2-6:1

# Enable Additional Software
. /etc/profile.d/modules.sh
module unload cmsub
module load shared Anaconda/2019.10
source activate lion

# Setup Output dir
targetdir="Thioamide_PolyProline/GEN_${SGE_TASK_ID}"
mkdir -p $targetdir
cp ${JOB_NAME}.py $targetdir/${JOB_NAME}.py
cp phe_cyanated.txt $targetdir
cp thioamideN.txt $targetdir
cp TBL.params $targetdir
cp lys_thioacetyl.txt $targetdir
chmod +x $targetdir/*
cd $targetdir

# Execute Script
python ${JOB_NAME}.py -T SC -F CNF -PNUM ${SGE_TASK_ID} 
python ${JOB_NAME}.py -T SC -F CNF -PNUM ${SGE_TASK_ID} -Rnd_Chi True 
python ${JOB_NAME}.py -T SC -F CNF -PNUM ${SGE_TASK_ID} -Rnd_Chi True -Inter_CIS True 
python ${JOB_NAME}.py -T SC -F CNF -PNUM ${SGE_TASK_ID} -Rnd_Chi True -Inter_CIS True -Term_CIS True 
python ${JOB_NAME}.py -T SC -F CNF -PNUM ${SGE_TASK_ID} -Rnd_Chi True -Inter_CIS True -Term_CIS True -CHPI_CIS True