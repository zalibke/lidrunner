#!/bin/bash

# bash script to gel together entire pipeline to run LDhat locally on specified pop
# gcloud auth application-default login #MUST AUTHENTICATE FIRST TO ACCESS MALARIAGEN API
# arguments

POP_NAME=$1
WINDOW_WIDTH=$2
SITES_PER_WINDOW=$3
CHUNK=$4
WINDOWS_PER_CHUNK=$5
N=$6
RANDSAMP=$7

echo "lidrunner v1.2"
echo "Zane Libke - July 2024"
echo "zane.libke23@imperial.ac.uk"
#echo "INPUT PARAMETERS: SAMPLE_QUERY = "$SAMPLE_QUERY""
echo "POP_NAME = "$POP_NAME""
echo "WINDOW_WIDTH = "$WINDOW_WIDTH""
echo "SITES_PER_WINDOW = "$SITES_PER_WINDOW""
echo "CHUNK = "$CHUNK""
echo "WINDOWS_PER_CHUNK = "$WINDOWS_PER_CHUNK""
echo "N = "$N""
echo "RANDSAMP = "$RANDSAMP""

mkdir -p ../data/input
mkdir -p ../data/output

echo "converting "$POP_NAME" hdf5 to Rdata..."
Rscript 2-convert_hdf5.R "$POP_NAME"
echo "${POP_NAME} RData files saved to ../data/"

echo "running LDhat "
Rscript 3-run_and_extract_LDhat.R "$POP_NAME" "$WINDOW_WIDTH" "$SITES_PER_WINDOW" "$CHUNK" "$WINDOWS_PER_CHUNK" "$N" "$RANDSAMP"

echo "plotting LDhat results..."
Rscript 4-plot_LDhat.R "$POP_NAME"

echo "raw data and plots saved"

# clean up
mv ../data/output "../data/output_$POP_NAME"
mv ../data/input "../data/input_$POP_NAME"

