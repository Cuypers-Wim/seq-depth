#!/bin/bash

##	Wim Cuypers (cuyperswim001@gmail.com)
##
## Script to count the number of sequenced bases in fastq files
##
## Where I got my inspiration:
##   Argument parsing in bash: https://pretzelhands.com/posts/command-line-flags
##   Seqtk  by Heng Li: https://lh3.github.io/2018/11/12/seqtk-code-walkthrough
##   Nullarbor by Torsten Seemann: https://github.com/tseemann/nullarbor
##   Gnu parallel usage: 
##       https://stackoverflow.com/questions/60527000/how-to-execute-a-script-in-parallel-in-bash
##       https://davetang.org/muse/2013/11/18/using-gnu-parallel/


# TODO

# add samtools depth processing as optional 


# Time used to calculate sequencing depth of 282 fastq.gz files (141 forward and 141 reverse) using 40 threads

# Time count_fastq_bases_seqtk = 
#real    0m39.168s
#user    22m0.571s
#sys     0m10.886s

# Time count_fastq_bases without seqtk = 
# real    1m4.799s
# user    34m51.587s
# sys     5m26.791s

# the output for both options was identical

wrong_input_message()
{
echo -e "The following argument was not recognised: \n"
echo -e "${OTHER_ARGUMENTS[@]} \n"
}

usage()
{
	echo -e "Calculate the sequencing depth of a bunch of fastq files" \
    "\nExample:" \
	"\n $0 -i path/to/fastq_file_directory -p prefix -g 500000 -o path/to/output_directory \n"
}


print_help()
{

usage 

cat <<EOF
MANDATORY ARGUMENTS:
-i, --input  Path to a folder containg fastq files
-g, --genome_size  Size of the reference genome
OPTIONS:
-p, --prefix  Prefix used to name output files (default 'coolproject')
-o, --outdir  Output folder (default '.')
-h, --help  Prints this help
	
	
EOF

}


# Default values of arguments

INPUT_FASTQFILES="$PWD"
PREFIX="coolproject"
GENOME_SIZE=500000
OUTPUT_DIRECTORY="$HOME"
OTHER_ARGUMENTS=()

#
# FUNCTIONS
#

debug()
{

echo "-----these are the set variables after argument parsing-------"

echo "--input fastqfiles--"
echo $INPUT_FASTQFILES

echo "--prefix--"
echo $PREFIX

echo "--genome size--"
echo $GENOME_SIZE

echo "--output directory--"
echo $OUTPUT_DIRECTORY

echo "--other arguments--"
echo $OTHER_ARGUMENTS

echo "-------------------------------------------------------------"

}


# file_info
#
# Check whether all fastq files are compressed, uncompressed, or a mixed
# Only if no files with suffix '.fastq' or '.fastq.gz' are found, the script willl exit
#
# Input arguments: 
#    * gzipped fastq file via gnu parallel
#    * outputfile prefix
#

file_info()
{
	
	echo "file information:"
	
	# number of files
	NF="$(ls *.fastq* | wc -l)"
	# number of compressed files
	NFC="$(file $(ls *.fastq*) | grep "compressed" | wc -l)"
	
	# let's see what we have here
	
	# all compressed?
	if [ "$NF" -eq "$NFC" ]
	then
	echo "all files are compressed"
	
	# all uncompressed?
	elif [ "$NF" > 0 ] && [ "$NFC" == 0]
	then
	echo "Files are uncompressed"
	
	# mixed compressed - not compressed?
	elif [ "$NF" > 0 ] && [ "$NFC" > 0 ] && [ "$NF" -neq "$NFC" ]
	then
	echo "Warning, not all files are compressed. That's ok for this program, but maybe you deserve to know"
	
	# have you not read the usage?
	elif [ "$NF" == 0 ]
	then
	echo "ERROR, no files are present with suffix 'fastq' or 'fastq.gz' "
	usage
	exit 1
	
	else
	echo "there are $NF files, of which $NFC are compressed"
	
	fi
	
}

# count_fastq_bases
#
# Unzips and reads gzipped fastq files
#
# Input arguments: 
#    * gzipped fastq file via gnu parallel
#    * outputfile prefix
#

count_fastq_bases()
{

FASTQ="$1"
FASTQfileName="$(basename ${FASTQ%.*})"
outputfileName="$2"

# both gzipped and uncompressed fastq can be handled

if (file $FASTQ | grep "compressed" ) ; 
then
	zcat $FASTQ | paste - - - - \
	| cut -f2 | tr -d '\n' | wc -c \
	| awk -v OFS='\t' -v name="$FASTQfileName" '{print name,$1}' >> "$outputfileName"
else
	cat $FASTQ | paste - - - - \
	| cut -f2 | tr -d '\n' | wc -c \
	| awk -v OFS='\t' -v name="$FASTQfileName" '{print name,$1}' >> "$outputfileName"
fi

}

count_fastq_bases_seqtk()
{

FASTQ="$1"
FASTQfileName="$(basename ${FASTQ%.*})"
outputfileName="$2"

seqtk fqchk "$FASTQ" | grep "ALL" | awk -v OFS='\t' -v name="$FASTQfileName" '{print name, $2}' >> "$outputfileName"

}

# export  functions for use in next function

export -f count_fastq_bases
export -f count_fastq_bases_seqtk


# get_fastq_basecounts
#
# Uses gnu parallel to excute 'count_fastq_bases' on each fastq file 
#  (one will be processed per thread)
# If Heng Li's seqtk is installed on the system, 'seqtk fqchk' is used to extract a basecount
#  which is faster than unzipping and counting bases after 'zcat', extracting sequence lines and 'wc -l'
#
# Input arguments: 
#    * fastq.gz files
#    * outputfile prefix
#


get_fastq_basecounts()
{

out_file="$PREFIX"_basecounts.txt

touch "$out_file"

# if seqtk is installed, use it
# seqtk handles both fastq and fastq.gz

if [ -x "$(command -v seqtk)" ]
then
	echo "using seqtk"
	parallel count_fastq_bases_seqtk ::: "$(ls *.fastq*)" ::: "$out_file"
else 
	parallel count_fastq_bases ::: "$(ls *.fastq*)" ::: "$out_file"
fi

echo "$PREFIX processing done"

}

# $(ls *.@(fastq|fastq.gz))


# basecountTotal_coverage
#
# Processes text files with basecounts per fastq file (previous steps)
#  and sums the basecounts for R1 and R2 of a sequencing run
#
# Input arguments: 
#    * outputfile prefix
#

basecountTotal_coverage()
{

# sum basecounts for forward and reversed reads for the same sample

sort "$PREFIX"_basecounts.txt \
| paste - - \
| awk -v OFS='\t' -v size="$GENOME_SIZE" '{print $1, $2+$4, ($2+$4)/size}' \
| sed 's/_R1.fastq//g' > "$OUTPUT_DIRECTORY"/"$PREFIX"_summed_basecounts.tab

rm "$PREFIX"_basecounts.txt

}

# set variables

# data="/mnt/data/wcuypers/groupC1/coverage"
# outDir="/mnt/data/wcuypers/groupC1/coverage/test2"
# prefix="groupC1"

#######
# GenomeSize Typi CT18
# genomeSize=4809037
#######

#######
# GenomeSize Concord
# genomeSize=4811936
#######
# run 

#
# MAIN
#

# Parse arguments

if [[ "$#" -eq 0 ]] ; then
    echo 'please provide arguments'
    exit 0
fi

while [ "$#" -gt 0 ]
do
	case "$1" in
	-i|--input)
        INPUT_FASTQFILES="$2"
        shift
        shift
        ;;
		
	-p|--prefix)
        PREFIX="$2"
        shift
        shift
        ;;
		
	-g|--genome_size)
        GENOME_SIZE="$2"
        shift
        shift
        ;;
		
	-o|--outdir)
        OUTPUT_DIRECTORY="$2"
        shift
        shift
        ;;
		
	-h|--help)
        print_help
        exit
        ;;
		
    *)
		OTHER_ARGUMENTS+=("$1")
		shift
		wrong_input_message
        usage 
		exit
        ;;
		
    esac
done # End of for loop parsing command line options.


# main commands

cd "$INPUT_FASTQFILES"
file_info
get_fastq_basecounts "$PREFIX"
basecountTotal_coverage

# debug
