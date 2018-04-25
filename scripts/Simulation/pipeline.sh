#!/bin/bash 


# SCRIPTS
PYTHON_SCRIPT=/home/tjiang/Tools/SV_install/ALU_simulation/sim_alu.py
PBSIM_PATH=/home/tjiang/Tools/pbsim-1.0.3-Linux-amd64/
PBSIM_SCRIPT=${PBSIM_PATH}Linux-amd64/bin/pbsim
#DATA
REF_GENOME=/home/tjiang/Ref/hg19.new.fa
# REF_GENOME=${3}
ALU_LOCUS=/data/tjiang/workspace/ALU_scripts/Alu_hg19.bed
ALU_SEQUENCE=/home/tjiang/dbRIP/Alu_hg19_v2h.txt

WORK_PATH=${1}
DEEPTH=${2}

starttime=`date +'%Y-%m-%d %H:%M:%S'`
echo ">>> Fake genome construction."
python2 ${PYTHON_SCRIPT} ${REF_GENOME} ${ALU_LOCUS} ${ALU_SEQUENCE} ${WORK_PATH}

echo ">>> Build sim_reads."
${PBSIM_SCRIPT} ${WORK_PATH}simulation.fa --depth ${DEEPTH} --model_qc ${PBSIM_PATH}data/model_qc_clr --length-mean 8000 --accuracy-mean 0.84 --accuracy-min 0.8 --accuracy-max 1 --difference-ratio 9:107:43 --prefix ${WORK_PATH}

echo ">>> Map sim_reads."
ngmlr -t 32 -r ${REF_GENOME} -q ${WORK_PATH}_0001.fastq -o ${WORK_PATH}_0001.sam

echo ">>> Alignments processing."
samtools view -Sb ${WORK_PATH}_0001.sam | samtools sort -O bam -T /data/tjiang/ - > ${WORK_PATH}_0001.bam
samtools index ${WORK_PATH}_0001.bam

echo ">>> Variant calling."
sniffles -m ${WORK_PATH}_0001.bam -b ${WORK_PATH}_0001.bed -t 32 -s 5

endtime=`date +'%Y-%m-%d %H:%M:%S'`
start_seconds=$(date --date="$starttime" +%s);
end_seconds=$(date --date="$endtime" +%s);
echo "Time:"$((end_seconds-start_seconds))"s"