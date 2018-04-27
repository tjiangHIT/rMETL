# !/usr/bin/shell

threads=$1
DataType=$2

# rMETL running
if [ "${DataType}" == "pac" ]; then
	# echo "pac"
	TYPE="pacbio"
else
	# echo "ont"
	TYPE="ont"
fi

# extract
rMETL extract ${BAM} ${REF} ${TEMP_DIR} ${TEMP_DIR} -t ${threads} -x ${TYPE}
# map
rMETL map ${FASTA} ${ME_REF} ${TEMP_DIR} -t ${threads} -x ${TYPE}
# call
rMETL call ${SAM} ${REF} vcf ${TEMP_DIR} --sample NA12878