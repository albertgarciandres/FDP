#!/bin/bash

cpat.py -x ${DIR}/flyPrebuild/fly_Hexamer.tsv -d ${DIR}/flyPrebuild/Fly_logitModel.RData \
			--top-orf=100 \
			--antisense \
			-g ${DIR}/minimap2/${Condition}/tmp/$transcript.$transcript_trinity.fa \
			-o ${DIR}/minimap2/${Condition}/tmp/CPAT/$transcript_trinity
