#!/bin/bash
Q=4
N=2
time ../cse_indel_simple.py ../data/bsubtilis/bsubtilis_NC_000964.fasta ../data/bsubtilis/GAIIx-bs.bam $Q $N -d 0 -s ../serialized/genome_annotate_bsubtilis_indels.pkl &> ../output/bsubtilis/indel_q$Q\_n$N\_o0.out.corrected

Q=8
N=2
time ../cse_indel_simple.py ../data/bsubtilis/bsubtilis_NC_000964.fasta ../data/bsubtilis/GAIIx-bs.bam $Q $N -d 0 -s ../serialized/genome_annotate_bsubtilis_indels.pkl &> ../output/bsubtilis/indel_q$Q\_n$N\_o0.out.corrected

Q=10
N=2
time ../cse_indel_simple.py ../data/bsubtilis/bsubtilis_NC_000964.fasta ../data/bsubtilis/GAIIx-bs.bam $Q $N -d 0 -s ../serialized/genome_annotate_bsubtilis_indels.pkl &> ../output/bsubtilis/indel_q$Q\_n$N\_o0.out.corrected
