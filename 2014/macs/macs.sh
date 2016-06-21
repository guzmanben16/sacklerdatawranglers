#!/bin/bash
#$ -S /bin/bash
#$ -cwd

module load macs/2.0.10.20131216

export REF=/phoenix/iGenomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta

macs2 callpeak -t myAlign.bam -c myControl.bam
