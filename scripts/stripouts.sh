#!/bin/bash
cd $1
find . -type f -name '*-impout.*' -delete
find . -type f -name '*-seqout.*' -delete
find . -type f -name '*-rev.*' -delete
find . -type f -name '*-indivout.*' -delete
find . -type f -name '*.newick' -delete
find . -type f -name '*.nwk' -delete
find . -type f -name '*outtree*' -delete
cd ..