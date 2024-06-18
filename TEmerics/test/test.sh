#!/bin/bash

cleanup () {
	rc=$?
	cd $user_dir
	echo "Exit status: $rc"
	rm -rf .snakemake
}
trap cleanup EXIT

set -eo pipefail
set -u
set -x
user_dir=$PWD
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" > /dev/null 2>&1 && pwd)"
cd $script_dir

snakemake \
	--snakefile="../workflow/Snakefile" \
	--cores 4 \
	--configfile="../config/config.yaml" \
	--use-conda \
	--printshellcmds \
	--rerun-incomplete \
	--verbose
