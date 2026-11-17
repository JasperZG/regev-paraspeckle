#!/bin/zsh
set -e
cd /Users/jasper/Desktop/regev-paraspeckle
export PATH=$PWD/tools/env/bin:$PATH
echo "gunzip"; gunzip -c ref/GRCh38.p13.genome.fa.gz > ref/GRCh38.p13.genome.fa
samtools faidx ref/GRCh38.p13.genome.fa
echo "features"; python analysis/gene_features_genome.py
echo "probe scan"; python analysis/probe_scan_genome.py > work/probe_scan_genome.log 2>&1 &
echo "deseq already run separately"
echo "master"; python analysis/build_master.py
echo "cqn"; python analysis/cqn.py > work/cqn.log 2>&1
wait
echo "chain done"
