Opportunistic development of the Diplotopia nextflow pipeline is going on … and it is starting to look good ! 

 

Diplotopia is a pipeline aiming at processing data for diploid (and some parts can be used for polyploid or even haploid organisms). 

 

The pipeline is composed of several tracks that can be used either individually (eg similar purpose) or sequentially. 

 

Not all planned tracks are developed yet. But we now start to have some usable tracks (some still with a bit of pain which would require minor modifications) which you might find useful for your projects. So should you have interest, please let us know, so we know what you would need most … and maybe contribute to the opportunistic development that is ongoing, and would certainly force us to improve current documentation!.  

 

Right now the pipeline tracks that are advanced enough to be used are and can do: 

Filtering out contaminant contigs from genomes/assemblies. The idea is that you can be flexible to filter out organisms, you do not need to filter by what is not your species of interest, instead you can filter out contigs that do not belong to the same level of taxonomical classification (genus, family, order). This is quite useful when you are working with a newly sequenced organism with sparse prior genetic information in public databases. This track allows you to filter out contaminants by keeping sequences similar to organisms your subject species is likely related to. (FILTER_CONTIGS track)
Trial of different assemblers (newly implemented, TRYSSEMBLY track). We now have Masurca pipeline (short reads / hybrid assembly, based on SOAP - Flye and masurca own assembler), Dipspades (Spades for diploids). Platanus (will try to get platanus-allee) and Redundans (which contains an early version of platanus) are currently being  implemented (under testing). 
Assembly quality and comparison of assemblies. The idea of this track (COMPASS) is to help you evaluate which assembly from several assemblers are best suited for your research question of interest. It provides many metrics for assembly quality evaluation (eg, total genome size, fragmentation, synteny according to a reference) and completeness evaluation. You can then choose the assembly that may have the most single copy genes, or the most genes in duplicated families,  or the one that is most contiguous … After all, assemblies are only a model of the truth, and you must use the best model to answer the scientific question at hand!
Cleaning assembly from non previously recognized orthologous sequences (HAPLOPURGE track). Haploids and polyploids have multiple alleles at the same locus in the genome. However, when the levels of divergence is a bit too high, the assemblers can fail to recognize the alleles as belonging to the same locus. In consequence, the assemblies obtained can be much larger that the actual haploid genome (haploid genome is supposed to represent the genome of an organism) because the contigs of the alleles at those locus can be duplicated and can eventually produce chimeric sequence also. Haplopure is a track that aims at identifying alleles, repeats and junk that should not be included in an “as good as possible assembly”. 
Variant calling using freebayes, filtering and normalization, both at the individual level and population level (VARWRRUM  : Variants At Regions Where Reference Reads Uniquely Map, sorry for the name … but I started to get nuts because of the complexity during the development stage). We use a haplotype representation of the genome (I call that haplosembly) as reference to call variants. It is possible to call variants for each sample according to a specific reference, but also variants that are unique to a group of samples (eg. resistant, vs non resistant samples to a certain drug). 
 

Feature development (hopefully one day)

Merging different assemblies into one : as different assemblers use different algorithms which make different assumptions regarding genome structure and composition, all assemblies contain part of the truth and part of errors. The idea is that combining different assemblies should take us a bit closer to the true representation of the genome of the organism of interest. 
Implementing more assemblers for testing in TRYSSEMBLY (planned: NECAT for nanopore reads, and maybe others as tested in this article)
Improving documentation (ongoing but now published on website)
Improving options implementation to increase usability at a larger scale
Reads decontamination prior to assembly (reducing data complexity, facilitating assemblies)
FunGen : was supposed to be another type of variant calling, but right now, its on suspend… unless someone wants to use it !. 
Adding own scripts to improve visualisation of genome assembly quality (eg. as done in Sapro poster) and likely improvement depending on needs !
 

The development of this pipeline has been done essentially to analyze Saprolegnia data  (which is a difficult genome to analyze as it is composed of ca. 40% repeats) for the Saprolegnia pilot project), and need for analysis of flatworms (TapeResist/PAHW) triggered the further development of the track TRYSSEMBLY. 

 

Hoping that you will find some of the tracks reusable for your projects. Let us now ! 

All the best.

see less
2  DIPLOTOPIA pipeline and program descriptions – Diplotopia pipeline documentation