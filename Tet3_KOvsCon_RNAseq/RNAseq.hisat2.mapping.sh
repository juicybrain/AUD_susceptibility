#!/bin/bash

        wd="/home/..."
        keywd=*R1.fastq.gz
        i=0
        name=()
        rRNA="/home/yli/ref/institute/rRNA/rRNA.mm10"




        for line in `ls $wd$keywd`
        do
        name[${i}]=${line%%".R1"*}
        let i=${i}+1
        done

        for spl in ${name[@]}
            do
            echo $spl
#      trim_galore --paired  $spl\_R1_001.fastq.gz $spl\_R2_001.fastq.gz

            sample=${spl##*"/"}
            echo ${sample}
           bowtie2 -p 8 -x $rRNA -U ${sample}.R1.fastq.gz 2>${sample}.rRNA.mapping.log | samtools view -S -F 4 - | cut -f 1 >${sample}_R1.rRNA.list &
           bowtie2 -p 8 -x $rRNA -U ${sample}.R2.fastq.gz 2>>${sample}.rRNA.mapping.log | samtools view -S -F 4 - | cut -f 1 >${sample}_R2.rRNA.list

           wait
           perl /home/yli/bin/filter_rRNA.PE.pl  ${sample}_R1.rRNA.list ${sample}_R2.rRNA.list ${sample}.R1.fastq.gz  ${sample}.R2.fastq.gz ${sample}.R1.clean.fq.gz ${sample}.R2.clean.fq.gz
           wait


#          "( hisat2 -p 36 --dta -x /home/yli/ref/index/hisat/mm10/genome -1 ${sample}\_R1.clean.fq.gz -2 ${sample}\_R2.clean.fq.gz --dta-cufflinks 2>${sample}\.log | samtools view -Sb -F 4  -q 20 -@ 10 - | samtools sort -n -@ 10 - -o ${sample}\.u.s.bam)"
           hisat2 -p 2 --dta -x /home/yli/ref/index/hisat/mm10/genome --rna-strandness RF -1 ${sample}\.R1.clean.fq.gz -2 ${sample}\.R1.clean.fq.gz --summary-file ${sample}\.map.log  2>${sample}\.log | samtools view -Sb -F 4  -q 20 -@ 4 - | samtools sort -n -@ 4 - -o ${sample}\.u.s.bam
           samtools sort  -n $sample\.u.s.bam --threads 4 -o $sample\.u.s.n.bam

#            featureCounts -s 2 -t exon  -g gene_id -a /home/yli/ref/institute/Encode/gencodevM7/ENCFF871VGR.gtf -o ${sample}\.FC.reads.counts.txt ${sample}\.u.s.n.bam 2>${sample}\.reads.FC.log &
#            featureCounts -p -s 2 -B -t exon -g gene_id -a /home/yli/ref/institute/Encode/gencodevM7/ENCFF871VGR.gtf -o ${sample}\.FC.pe_fragments.counts.txt ${sample}\.u.s.n.bam 2>${sample}\.pe_freagemet.FC.log&
#            featureCounts -T 6 -p -t exon -g gene_id -a /home/yli/ref/institute/gencode/gencode.vM24.basic.annotation.gtf -s 2 -o ${sample}\.basic.readscount.s.txt ${sample}\.u.s.n.bam 2>${sample}\.basicRC.s.log &


             featureCounts -T 4 -p -t exon  -g gene_id -a /home/yli/ref/institute/gencode/gencode.vM24.annotation.gtf -s 2 -o ${sample}\.main.readscounts.s.txt ${sample}\.u.s.n.bam 2>${sample}\.mainRC.s.log &
