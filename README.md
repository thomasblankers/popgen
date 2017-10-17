# popgen
GCcontent.sh simply calculates the GC content of a multifasta file. 

MKTtest is a pipeline to do the Mcdonald Kreitman test assuming you have a VCF file, a genome/transcriptome assembly in fasta format, and a gff or gtf file with the annotation of the open reading frames and coding sequence regions.

TajimaD.r is an R function to calculate Tajima's D based on the the nucleotide diversity, pi, the number of segregating sites, S, and the number of sequences. 

VCF2FASTA and VCF2indFASTA convert vcf files and a reference fasta into individual sequences, either a single sequence per individual using the IUPAC coding for heterozygote ambiguities (VCF2FASTA_IUPAC) or both haplotypes for diploid individuals (VCF2indFASTA). 
