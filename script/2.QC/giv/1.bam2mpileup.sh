perl /data/person/huz/my_software/Damage-estimator/split_mapped_reads.pl \
	-bam /data/person/huz/DNAseq/UCEC_final/02_mapping/03_bqsr/FD_L_037_C.bqsr.bam \
	-genome /data/database/homo_sapiens/GRCh37/Sequence/WholeGenomeFasta/human_g1k_v37_decoy.fasta \
	-mpileup1 FD_L_037_C.bqsr.R1.mpileup -mpileup2 FD_L_037_C.bqsr.R2.mpileup \
	-Q 20 -q 20
