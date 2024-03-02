## somatic variants calling
1. mutect2
2. Pindel -w 10
3. Strelka  isSkipDepthFilters = 0
4. VarScan --min-coverage 3 --min-var-freq 0.08 --p-value 0.10 --somatic-p-value 0.05 --strand-filter 1
5. SomaticSniper (bam-somaticsniper -F vcf -q 1 -Q 15
6. Manta


msisensor
CNVkit
ASCAT
