#!/bin/bash

# Setting to strict mode
set -euo pipefail

# Filter percent table
cat $1 | awk -F'\t' '{ if (!seen[$1 "-" $4]) \ { print $0; seen[$1 "-" $4]=1 } }' > $2

