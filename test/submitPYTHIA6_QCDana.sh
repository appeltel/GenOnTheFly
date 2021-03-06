#!/bin/bash

export USER_ACCOUNT=cms_stage3
export CMSRUN_TUNE=Z2
export CMSRUN_ENERGY=5020
export CMSRUN_MAXEVENTS=400000

export CMSRUN_PROCESSTYPE=MB
export CMSRUN_PTHATLOW=0
export CMSRUN_PTHATHIGH=20
export CMSRUN_OUTPUT=AnaQCD_${CMSRUN_TUNE}_${CMSRUN_ENERGY}TeV_${CMSRUN_PROCESSTYPE}_${CMSRUN_PTHATLOW}to${CMSRUN_PTHATHIGH}.root
sbatch runPYTHIA6_Z2.slurm --account=$USER_ACCOUNT --export=CMSRUN_TUNE,CMSRUN_PROCESSTYPE,CMSRUN_PTHATLOW,CMSRUN_PTHATHIGH,CMSRUN_ENERGY,CMSRUN_MAXEVENTS,CMSRUN_OUTPUT

export CMSRUN_PROCESSTYPE=NSD
export CMSRUN_PTHATLOW=20
export CMSRUN_PTHATHIGH=30
export CMSRUN_OUTPUT=AnaQCD_${CMSRUN_TUNE}_${CMSRUN_ENERGY}TeV_${CMSRUN_PROCESSTYPE}_${CMSRUN_PTHATLOW}to${CMSRUN_PTHATHIGH}.root
sbatch runPYTHIA6_Z2.slurm --account=$USER_ACCOUNT --export=CMSRUN_TUNE,CMSRUN_PROCESSTYPE,CMSRUN_PTHATLOW,CMSRUN_PTHATHIGH,CMSRUN_ENERGY,CMSRUN_MAXEVENTS,CMSRUN_OUTPUT

export CMSRUN_PROCESSTYPE=NSD
export CMSRUN_PTHATLOW=30
export CMSRUN_PTHATHIGH=50
export CMSRUN_OUTPUT=AnaQCD_${CMSRUN_TUNE}_${CMSRUN_ENERGY}TeV_${CMSRUN_PROCESSTYPE}_${CMSRUN_PTHATLOW}to${CMSRUN_PTHATHIGH}.root
sbatch runPYTHIA6_Z2.slurm --account=$USER_ACCOUNT --export=CMSRUN_TUNE,CMSRUN_PROCESSTYPE,CMSRUN_PTHATLOW,CMSRUN_PTHATHIGH,CMSRUN_ENERGY,CMSRUN_MAXEVENTS,CMSRUN_OUTPUT

export CMSRUN_PROCESSTYPE=NSD
export CMSRUN_PTHATLOW=50
export CMSRUN_PTHATHIGH=80
export CMSRUN_OUTPUT=AnaQCD_${CMSRUN_TUNE}_${CMSRUN_ENERGY}TeV_${CMSRUN_PROCESSTYPE}_${CMSRUN_PTHATLOW}to${CMSRUN_PTHATHIGH}.root
sbatch runPYTHIA6_Z2.slurm --account=$USER_ACCOUNT --export=CMSRUN_TUNE,CMSRUN_PROCESSTYPE,CMSRUN_PTHATLOW,CMSRUN_PTHATHIGH,CMSRUN_ENERGY,CMSRUN_MAXEVENTS,CMSRUN_OUTPUT

export CMSRUN_PROCESSTYPE=NSD
export CMSRUN_PTHATLOW=80
export CMSRUN_PTHATHIGH=120
export CMSRUN_OUTPUT=AnaQCD_${CMSRUN_TUNE}_${CMSRUN_ENERGY}TeV_${CMSRUN_PROCESSTYPE}_${CMSRUN_PTHATLOW}to${CMSRUN_PTHATHIGH}.root
sbatch runPYTHIA6_Z2.slurm --account=$USER_ACCOUNT --export=CMSRUN_TUNE,CMSRUN_PROCESSTYPE,CMSRUN_PTHATLOW,CMSRUN_PTHATHIGH,CMSRUN_ENERGY,CMSRUN_MAXEVENTS,CMSRUN_OUTPUT

export CMSRUN_PROCESSTYPE=NSD
export CMSRUN_PTHATLOW=120
export CMSRUN_PTHATHIGH=170
export CMSRUN_OUTPUT=AnaQCD_${CMSRUN_TUNE}_${CMSRUN_ENERGY}TeV_${CMSRUN_PROCESSTYPE}_${CMSRUN_PTHATLOW}to${CMSRUN_PTHATHIGH}.root
sbatch runPYTHIA6_Z2.slurm --account=$USER_ACCOUNT --export=CMSRUN_TUNE,CMSRUN_PROCESSTYPE,CMSRUN_PTHATLOW,CMSRUN_PTHATHIGH,CMSRUN_ENERGY,CMSRUN_MAXEVENTS,CMSRUN_OUTPUT

export CMSRUN_PROCESSTYPE=NSD
export CMSRUN_PTHATLOW=170
export CMSRUN_PTHATHIGH=230
export CMSRUN_OUTPUT=AnaQCD_${CMSRUN_TUNE}_${CMSRUN_ENERGY}TeV_${CMSRUN_PROCESSTYPE}_${CMSRUN_PTHATLOW}to${CMSRUN_PTHATHIGH}.root
sbatch runPYTHIA6_Z2.slurm --account=$USER_ACCOUNT --export=CMSRUN_TUNE,CMSRUN_PROCESSTYPE,CMSRUN_PTHATLOW,CMSRUN_PTHATHIGH,CMSRUN_ENERGY,CMSRUN_MAXEVENTS,CMSRUN_OUTPUT

export CMSRUN_PROCESSTYPE=NSD
export CMSRUN_PTHATLOW=230
export CMSRUN_PTHATHIGH=300
export CMSRUN_OUTPUT=AnaQCD_${CMSRUN_TUNE}_${CMSRUN_ENERGY}TeV_${CMSRUN_PROCESSTYPE}_${CMSRUN_PTHATLOW}to${CMSRUN_PTHATHIGH}.root
sbatch runPYTHIA6_Z2.slurm --account=$USER_ACCOUNT --export=CMSRUN_TUNE,CMSRUN_PROCESSTYPE,CMSRUN_PTHATLOW,CMSRUN_PTHATHIGH,CMSRUN_ENERGY,CMSRUN_MAXEVENTS,CMSRUN_OUTPUT

export CMSRUN_PROCESSTYPE=NSD
export CMSRUN_PTHATLOW=300
export CMSRUN_PTHATHIGH=380
export CMSRUN_OUTPUT=AnaQCD_${CMSRUN_TUNE}_${CMSRUN_ENERGY}TeV_${CMSRUN_PROCESSTYPE}_${CMSRUN_PTHATLOW}to${CMSRUN_PTHATHIGH}.root
sbatch runPYTHIA6_Z2.slurm --account=$USER_ACCOUNT --export=CMSRUN_TUNE,CMSRUN_PROCESSTYPE,CMSRUN_PTHATLOW,CMSRUN_PTHATHIGH,CMSRUN_ENERGY,CMSRUN_MAXEVENTS,CMSRUN_OUTPUT


export CMSRUN_PROCESSTYPE=NSD
export CMSRUN_PTHATLOW=380
export CMSRUN_PTHATHIGH=10000
export CMSRUN_OUTPUT=AnaQCD_${CMSRUN_TUNE}_${CMSRUN_ENERGY}TeV_${CMSRUN_PROCESSTYPE}_${CMSRUN_PTHATLOW}to${CMSRUN_PTHATHIGH}.root
sbatch runPYTHIA6_Z2.slurm --account=$USER_ACCOUNT --export=CMSRUN_TUNE,CMSRUN_PROCESSTYPE,CMSRUN_PTHATLOW,CMSRUN_PTHATHIGH,CMSRUN_ENERGY,CMSRUN_MAXEVENTS,CMSRUN_OUTPUT

