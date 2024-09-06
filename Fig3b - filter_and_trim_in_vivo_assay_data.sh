# Bash script that uses cutadapt v1.16 to trim CS1 and CS2 off all initial TdT reads.
# It also filters out reads of low quality or if they are missing CS1 or CS2.
# Finally it filters out reads of length 0.
#
# Author: Jonathan Strutz
# Affiliation: PhD Student, Northwestern University
# Contact: jonathanstrutz2021@u.northwestern.edu

set -x
starting_dir=$PWD

module load fastqc

# Each NGS run is stored in a separate folder in run_dir
run_dir="/Users/marijamilisavljevic/Desktop/GMCF_1832_QC"


# This is the only part that changes between runs
data_dir="${run_dir}/FASTQ_CKC767/"

for dir in "${data_dir}/"*/ ; do

    # .fastq files can also be called .fq
    # Assumes lane is always 001
	read1zip=`find ${dir} -type f -name '*R1_001.f*q.gz'`
	read2zip=`find ${dir} -type f -name '*R2_001.f*q.gz'`

	read1="${read1zip%.*}"
	read2="${read2zip%.*}"

	gzip -cfd "$read1zip" > "$read1" 
	gzip -cfd "$read2zip" > "$read2"

	trim1="${read1%.*}_trimmed.fq"
	trim2="${read2%.*}_trimmed.fq"

	cutadapt \
		-m 1 \
		-e 0.05 \
		-q 30 \
		-O 10 \
		--discard-untrimmed \
		-a AGACCAAGTCTCTGCTACCGTA \
		-A TGTAGAACCATGTCGTCAGTGT \
		-o "${trim1}" \
		-p "${trim2}" \
		"${read1}" "${read2}"

	gzip -cf "$trim1" > "${trim1}.gz"
	gzip -cf "$trim2" > "${trim2}.gz"

	fastqc "$trim1"
	fastqc "$trim2"

done

cd "$starting_dir"
