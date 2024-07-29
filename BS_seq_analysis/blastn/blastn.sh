
makeblastdb -in rgs9_f_rc.fa -dbtype nucl -out rgs9_database

blastn -query all.trim.f.rc.fasta  -db rgs9_database -outfmt 1 >  all.trim.f.rc.fasta.res1.txt