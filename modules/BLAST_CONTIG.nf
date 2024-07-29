//BLAST sequence for each contig for each assembly

process BLAST_CONTIG { 

    //conda (params.enable_conda ? 'bioconda::bwa=0.7.8' : null)
	container 'quay.io/biocontainers/blast:2.14.1--pl5321h6f7f691_0'

    tag "${ID}_${contigID}"

    label 'process_high'


    input:
    tuple val(ID), path(contig)

    output: 
    tuple val(ID), val(contigID), path("*taxo.tsv"), emit: contigstaxo_ch

    script: 
    contigID = contig.getBaseName(1) - ".part_"  - ID
    

    """
    blastn -db $params.blastDB/nt -query $contig -out ${ID}_${contigID}_taxo.tsv -max_target_seqs 5 \
    -outfmt "7 qseqid sseqid  qlen slen pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids sscinames scomnames" \
    &> ${ID}_${contigID}_blast.log
    """
}
// Normally it should get the correspondence - But maybe the db was not in the path ? ASK Thomas 
// taxdb database is in the path defined by BLASTDB environment variable. After that you will be able to ask for several additional taxonomic information in the tabular output. For example:
// https://www.ncbi.nlm.nih.gov/books/NBK279690/#CmdLineAppsManual.Quick_start


// get the taxonomy db -> seems we need the new taxonomy 
// https://www.ncbi.nlm.nih.gov/books/NBK569841/


      