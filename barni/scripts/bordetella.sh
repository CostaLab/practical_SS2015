#!/bin/bash
Q=8
N=2
time ../cse_indel_simple.py /mnt/code/rwth/lab-bioinf/costalab/barni/data/bordetella_pertussis/18323/bordetella_pertussis_18323.fasta /mnt/code/rwth/lab-bioinf/costalab/barni/data/bordetella_pertussis/18323/IonTorrentPGM/ERX137391/ERR161541/ERR161541.bam $Q $N -d 0 -s ../serialized/genome_annotate_bordetella_pertussis_indels.pkl &> ../output/bordetella_pertussis/ionTorrent/ERX137391_indel_q$Q\_n$N\_o0.out

Q=10
N=2
time ../cse_indel_simple.py /mnt/code/rwth/lab-bioinf/costalab/barni/data/bordetella_pertussis/18323/bordetella_pertussis_18323.fasta /mnt/code/rwth/lab-bioinf/costalab/barni/data/bordetella_pertussis/18323/IonTorrentPGM/ERX137391/ERR161541/ERR161541.bam $Q $N -d 0 -s ../serialized/genome_annotate_bordetella_pertussis_indels.pkl &> ../output/bordetella_pertussis/ionTorrent/ERX137391_indel_q$Q\_n$N\_o0.out

