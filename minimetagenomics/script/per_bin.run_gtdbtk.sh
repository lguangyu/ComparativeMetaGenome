#!/bin/bash
#SBATCH -J gtdbtk
#SBATCH -o gtdbtk.log
#SBATCH -p short -N 1 -c 6
#SBATCH --time 1-00:00:00 --mem 224GB

export PATH="$HOME/opt/prodigal/2.6.3-2016-02-11:$PATH" # prodigal
export PATH="$HOME/opt/hmmer/3.3.2-2020-11-26/bin:$PATH" # hmmer
export PATH="$HOME/opt/pplacer/1.1.alpha19-2016-12-30:$PATH" # pplacer
export PATH="$HOME/opt/fastani/1.32-2020-09-30/bin:$PATH" # fastANI
export PATH="$HOME/opt/fasttree/2.1.11-2019-05-17:$PATH" # fasttree
export PATH="$HOME/opt/mash/2.3-2021-02-26:$PATH" # mash

export GTDBTK_DATA_PATH="$HOME/DATABASE/GTDBTK_DB/release202"

gtdbtk classify_wf \
	--genome_dir bins -x fasta \
	--out_dir gtdbtk \
	--cpus $SLURM_CPUS_PER_TASK \
	--pplacer_cpus $SLURM_CPUS_PER_TASK \
	--force # for test&dev, force overwriting

