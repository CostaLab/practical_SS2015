#!/bin/bash
Q=6
N=4
time ../../discovering_cse.py /mnt/code/rwth/lab-bioinf/costalab/barni/data/bordetella_pertussis/18323/bordetella_pertussis_18323.fasta /mnt/code/rwth/lab-bioinf/costalab/barni/data/bordetella_pertussis/18323/IonTorrentPGM/ERX137391/ERR161541/ERR161541.bam  $Q $N -d 0  &> ../benchmarks/bordetella_pertussis/ref_q$Q\_n$N\_o0.out

Q=6
N=4
time ../reference_manuel_scipy_chi.py /mnt/code/rwth/lab-bioinf/costalab/barni/data/bordetella_pertussis/18323/bordetella_pertussis_18323.fasta /mnt/code/rwth/lab-bioinf/costalab/barni/data/bordetella_pertussis/18323/IonTorrentPGM/ERX137391/ERR161541/ERR161541.bam $Q $N -d 0  &> ../benchmarks/bordetella_pertussis/new_q$Q\_n$N\_o0.out

