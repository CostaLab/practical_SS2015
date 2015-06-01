#!/bin/bash
#Q=4
#N=2
#time ../../discovering_cse.py ../data/bsubtilis/bsubtilis_NC_000964.fasta ../data/bsubtilis/GAIIx-bs.bam $Q $N -d 0  &> ../benchmarks/bsubtilis/ref_q$Q\_n$N\_o0.out
#
#Q=8
#N=2
#time ../../discovering_cse.py ../data/bsubtilis/bsubtilis_NC_000964.fasta ../data/bsubtilis/GAIIx-bs.bam $Q $N -d 0  &> ../benchmarks/bsubtilis/ref_q$Q\_n$N\_o0.out

Q=4
N=2
time ../reference_manuel_scipy_chi.py ../data/bsubtilis/bsubtilis_NC_000964.fasta ../data/bsubtilis/GAIIx-bs.bam $Q $N -d 0  &> ../benchmarks/bsubtilis/new_q$Q\_n$N\_o0.out

Q=8
N=2
time ../reference_manuel_scipy_chi.py ../data/bsubtilis/bsubtilis_NC_000964.fasta ../data/bsubtilis/GAIIx-bs.bam $Q $N -d 0  &> ../benchmarks/bsubtilis/new_q$Q\_n$N\_o0.out
