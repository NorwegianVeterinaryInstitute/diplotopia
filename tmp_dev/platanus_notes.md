https://bioinformaticsworkbook.org/dataAnalysis/GenomeAssembly/Arabidopsis/AT_platanus-genome-assembly.html#gsc.tab=0


1. Trim 

fofn file of names : fastq (not gz) 
platanus_trim -i pe.fofn -t 16

> using trimmed from Arvind 


2. Assemble
platanus assemble \
   -o platanus \
   -f {SRR3157034,SRR3166543}_?.fastq.trimmed \
   -t 12 \
   -m 115 \
   -tmp $TMPDIR

3. Scaffold
platanus scaffold \
    -o platanus \
    -c platanus_contig.fa \
    -b platanus_contigBubble.fa \
    -IP1 {SRR3157034,SRR3166543}_?.fastq.trimmed \
    -OP2 SRR3156163_?.fastq.int_trimmed \
    -n2 7000 \
    -a2 8000 \
    -d2 1000 \
    -OP3 SRR3156596_?.fastq.int_trimmed \
    -n3 19000 \
    -a3 20000 \
    -d3 2000 \
    -t 12 \
    -tmp $TMPDIR

    -n minimum iq
    - a average insert size
    -d standard deviation
    
    -OP is insert sizes can have several

    out : platanus_scaffold.fa: assembled sequences with gaps
platanus_scaffoldBubble.fa: removed bubble sequences
platanus_scaffoldComponent.tsv: table describing contig joins


4. gap closing
platanus gap_close \
    -o platanus \
    -c platanus_scaffold.fa \
    -IP1 {SRR3157034,SRR3166543}_?.fastq.trimmed \
    -OP2 SRR3156163_?.fastq.int_trimmed \
    -OP3 SRR3156596_?.fastq.int_trimmed \
    -t 16 \
    -tmp $TMPDIR