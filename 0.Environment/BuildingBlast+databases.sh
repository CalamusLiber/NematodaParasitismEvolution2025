# 配置blast+和DIAMOND的nr和UniVec比对数据库
/To/Your/Directory/blast/bin/makeblastdb -in nr.fa -input_type fasta -dbtype prot -title nr -parse_seqids -out nr -logfile build_blast+br.log
/To/Your/Directory/blast/bin/makeblastdb -in UniVec.fasta -input_type fasta -dbtype nucl -title UniVec -parse_seqids -out UniVec -logfile build_blast+UniVec.log
/To/Your/Directory/miniconda3/bin/diamond makedb --in nr.fa -d nr --taxonmap prot.accession2taxid.FULL --taxonnodes nodes.dmp --taxonnames names.dmp  
/To/Your/Directory/miniconda3/bin/diamond makedb --in uniref100.fasta -d uniref90 --taxonmap uniref90.taxlist