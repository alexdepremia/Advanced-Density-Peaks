#!/bin/bash
if [ $# != 2 ]
then
  echo "Usage: ./Mat2input.sh matrix_file_name input_file_for_DPA"
  exit
fi
MATFILE=${1}
INPFILE=${2}
awk '{for (i=NR+1;i<=NF;i++) print NR,i,$i}' ${MATFILE} > ${INPFILE}
