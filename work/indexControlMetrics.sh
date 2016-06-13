#!/bin/bash

FILE=$1
HTSLIB=$2

sort -n -k2 $FILE | $HTSLIB/bgzip -c >$FILE.gz
$HTSLIB/tabix -b 2 -e 2 $FILE.gz
