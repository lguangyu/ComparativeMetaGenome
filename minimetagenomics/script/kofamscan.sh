#!/bin/bash
# run kofamscan in conda (and snake make) result in broken pipe
# apparently caused by ruby version/library conflict
# therefore run manually

#SBATCH -J kofamscan
#SBATCH -o .log/kofamscan.log
#SBATCH -p short -N1 -c16
#SBATCH --mem 4GB --time 12:00:00

mkdir -p annotation
finish_mark="annotation/kofamscan.txt.end" # generated upon finishing
tmpdir="/dev/shm/$(whoami)_kofamscan_$RANDOM"

rm -f $finish_mark
mkdir -p $tmpdir
~/opt/kofamscan/1.3.0-2020-05-17/exec_annotation \
	-p ~/DATABASE/KOFAM_DB/2021-01-24/profiles \
	-k ~/DATABASE/KOFAM_DB/2021-01-24/ko_list \
	--cpu $SLURM_CPUS_PER_TASK \
	-f detail-tsv --no-report-unannotated \
	-o annotation/kofamscan.txt \
	--tmp-dir $tmpdir \
	gene_calling/resolved.faa \
&& touch $finish_mark
rm -rf $tmpdir

