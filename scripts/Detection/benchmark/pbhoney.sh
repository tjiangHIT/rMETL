# !/usr/bin/shell

# reads="/data/tjiang/workspace/ALU_scripts/hg19_sim/reads/20X/_0001.fastq"
reads="/data/tjiang/workspace/ALU_scripts/hg19_sim/reads/5X/_0001.fastq"
ref="/data/tjiang/workspace/ALU_scripts/hg19_sim/fake_ref/fake_chr1.fa"
threads=32
out=$1_0001
workspace=$1

# Blasr
blasr ${reads} ${ref} -sam -nproc ${threads} -clipping soft -sa ${ref}.sa -out /dev/stdout | samtools view -hb - | samtools sort -m 4G -@ 4 -T ${workspace} -o ${out}.bam -
samtools index ${out}.bam

# detect
Honey.py pie -o ${out}.final -n 8 ${out}.bam ${ref} --temp ${workspace} | samtools sort -m 4G -@ 4 -T ${workspace} -o ${out}.final.bam -
samtools index ${out}.final.bam

Honey.py tails -o ${out}.tails ${out}.final.bam

Honey.py spots --reference ${ref} -n 8 -o ${out}.spots ${out}.final.bam