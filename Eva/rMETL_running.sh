# !/usr/bin/shell

threads=$1
DataType=$2

# rMETL running
if [ "${DataType}" == "pac" ]; then
	# echo "pac"
	TYPE="pacbio"
	BAM="/data/tjiang/workspace/NA12878/na12878_pacbio_mts_ngmlr-0.2.3_mapped.bam"
	TEMP_DIR="/data/tjiang/workspace/NA12878/lab/eva_pac/${threads}X/"
else
	# echo "ont"
	TYPE="ont"
	BAM="/data/tjiang/workspace/NA12878/ngm_Nanopore_human_ngmlr-0.2.3_mapped.bam"
	TEMP_DIR="/data/tjiang/workspace/NA12878/lab/eva_ont/${threads}X/"
fi

REF="/home/tjiang/Ref/hs37d5/hs37d5.fa"
FASTA="${TEMP_DIR}potential_ME.fa"
ME_REF="/data/tjiang/workspace/database/super_TE.fa"
SAM="${TEMP_DIR}cluster.sam"

# echo ${TYPE}
# echo ${TEMP_DIR}
# echo ${FASTA}
# echo ${SAM}

# extract
rMETL extract ${BAM} ${REF} ${TEMP_DIR} ${TEMP_DIR} -t ${threads} -x ${TYPE}
# map
# rMETL map ${FASTA} ${ME_REF} ${TEMP_DIR} -t ${threads} -x ${TYPE}
# # call
# rMETL call ${SAM} ${REF} vcf ${TEMP_DIR} --sample NA12878
