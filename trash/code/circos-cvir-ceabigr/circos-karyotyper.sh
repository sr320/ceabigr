#!/bin/bash

# Script to create karyotype file for Circos

# Expects two inputs:

# 1. FastA index file
# 2. Abbreviation for species of interest

# Copy script to desired output directory and use like this:

# ./circos_pgen_karyotype.sh -f fasta.fai -s pg

# Creates output file named: karyotype.${species_abbreviation}.txt;
# formatted like:
# chr - pg1 1 0 89643857 chr1


fflag=false
sflag=false

usage() { echo "How to use:"; 
          echo "./circos-karyotyper.sh -f <path/to/fasta-index-file> -s <species_abbreviation>"
        }


# Set -f and -s flags for passing argument to script.
while getopts f:s:h flag
do
	case "${flag}" in
		f) fflag=true; fasta_index=${OPTARG};;
		s) sflag=true; species_abbreviation=${OPTARG};;
                h) usage; exit;;
                :) echo "Missing argument for -$OPTARG" >&2; exit 1;;
		*) echo "Unused flag: -${flag}";
		   echo "Use -f to specify path to FastA index. Use -s to specify species abbreviation."
		   exit;;
	esac
done

# Handle missing options.
if ((OPTIND == 1))
then
    echo "No options specified"
fi

shift "$((OPTIND - 1))"

if ! $fflag
then
    echo "Need FastA index path with -f flag" >&2
    exit 1
fi

if ! $sflag
then
    echo "Need species abbreviation with -s flag" >&2
    exit 1
fi


# Set output filename
karyo_file=karyotype.${species_abbreviation}.txt

# Make sure karytoype file doesn't already exist.
if [[ -f "${karyo_file}" ]]; then

  echo "Found existing file named ${karyo_file}."
  echo "Please rename or backup file, then re-run this script."
  exit 1
fi

# Parse FastA index file and format output file.
while IFS=$'\t' read -r name length offset linebases linewidth
do
	printf "chr - %s\n" "${species_abbreviation}${name} ${name} 0 ${length} chr${name}"
done < "${fasta_index}" > "${karyo_file}"




