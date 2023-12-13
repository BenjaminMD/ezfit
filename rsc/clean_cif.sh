#/bin/bash

# used to clean up a cif file using vesta

VESTA_PATH=$1
CIF_PATH=$2

if [ $# -eq 0 ]
  then
    echo "No arguments supplied"
    exit 1
fi

if [ ! -f $VESTA_PATH ]
  then
    echo $VESTA_PATH "not found!"
    exit 1
fi

if [ ! -f $CIF_PATH ]
  then
    echo $CIF_PATH "not found!"
    exit 1
fi

CIF_Out_Path=${CIF_PATH%.*}_clean.cif

$VESTA_PATH -nogui -i $CIF_PATH -o -format=cif $CIF_Out_Path
