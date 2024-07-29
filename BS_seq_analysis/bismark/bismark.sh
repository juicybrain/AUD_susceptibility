


# get the target region and extend hundreds of NNNNNN in both up- and down-stream
# prepare the genome
bismark_genome_preparation  ~/workspace/alcohol_study/alcohol_sanger_sequencing/ab1/all_fasta/rgs9_f_rc_bismark


# trim the fastq file

trim_galore  --clip_R1 50 --three_prime_clip_R1 50 all.trim.f.rc.fastq

# alignment

bismark  --bam --bowtie2 --non_directional  ~/workspace/alcohol_study/alcohol_sanger_sequencing/ab1/all_fasta/rgs9_f_rc_bismark  all.trim.f.rc_trimmed.fq  -o ./all.rc_bismark

# extract the bam files by group



# extract mC

cat <(samtools view -H all.trim.f.rc_trimmed_bismark_bt2.bam) <(samtools view all.trim.f.rc_trimmed_bismark_bt2.bam | grep "^1-") | samtools view -bS -o all.trim.f.rc_trimmed_bismark_bt2.Rep1.bam
cat <(samtools view -H all.trim.f.rc_trimmed_bismark_bt2.bam) <(samtools view all.trim.f.rc_trimmed_bismark_bt2.bam | grep "^2-") | samtools view -bS -o all.trim.f.rc_trimmed_bismark_bt2.Rep2.bam
cat <(samtools view -H all.trim.f.rc_trimmed_bismark_bt2.bam) <(samtools view all.trim.f.rc_trimmed_bismark_bt2.bam | grep "^3-") | samtools view -bS -o all.trim.f.rc_trimmed_bismark_bt2.Rep3.bam
cat <(samtools view -H all.trim.f.rc_trimmed_bismark_bt2.bam) <(samtools view all.trim.f.rc_trimmed_bismark_bt2.bam | grep "^4-") | samtools view -bS -o all.trim.f.rc_trimmed_bismark_bt2.Rep4.bam