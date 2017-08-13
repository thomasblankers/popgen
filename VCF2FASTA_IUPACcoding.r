## I wrote this script to generate a large number of 2000 bp sequences randomly drawn from the reference genome for RAxML
## note that the scripts prints out one sequence per diploid sample using the IUPAC coding for heterozygote genotypes.
## If you simply want to output all haplotypes in a population (2 haploytpes per individual) for a given set of sequences please refer to the script VCF2indFASTA.r
## use VCFtools to generate --012 output of VCF file, create an .Allele file by copying the CHROM,REF,ALT columns from the VCF file
## I advise to subset reference fasta using samtools faidx to contain only those loci for which alignment is required. Large fasta files will make this a painfully slow process.
## faidx <ref_fasta> <regions_file> > <fastafile_subset>

# read in SNPs
library(seqinr)
geno<-t(read.delim("VCF.012", header=FALSE))[-1,]
allele<-read.delim("VCF.Allele", sep="/", header=FALSE)
snps<-read.delim("VCF.pos", header=FALSE)
ind<-read.delim("VCF.indv", header=FALSE)
genotypes<-data.frame(cbind(snps,allele,geno))
colnames(genotypes)<-c("contig","pos","ref","alt",as.character(ind[,1]))

## sample all SNPs to retain one SNP per contig at random and select nucleotide in the range of 2 kb around the SNP (plus 1000 bp and minus 1000 bp).

 one_SNP_per_contig<-data.frame(unique(snps[,1]))
 one_SNP_per_contig$SNP<-NA
 for(i in 1:nrow(one_SNP_per_contig)) {
	contig_snps<-snps[snps[,1]==as.character(one_SNP_per_contig[i,1]),]
	if(nrow(contig_snps)!=1){
		one_SNP_per_contig$SNP[i]<-sample(contig_snps[,2],1)
		} else {
		one_SNP_per_contig$SNP[i]<-contig_snps[,2]}
	}
one_SNP_per_contig$plus1000<-one_SNP_per_contig$SNP+1000
one_SNP_per_contig$min1000<-one_SNP_per_contig$SNP-1000

# update original SNP positions relative to fragment start and end
one_SNP_per_contig_adjusted<-merge(one_SNP_per_contig_adjusted,genotypes,"contig",sort=FALSE)
one_SNP_per_contig_adjusted$pos_adjusted<-one_SNP_per_contig_adjusted$pos-one_SNP_per_contig_adjusted$min1000_adjusted
one_SNP_per_contig_adjusted<-one_SNP_per_contig_adjusted[one_SNP_per_contig_adjusted[,"pos_adjusted"] > 0,]
one_SNP_per_contig_adjusted<-one_SNP_per_contig_adjusted[one_SNP_per_contig_adjusted[,"pos_adjusted"] < 2002,]


fa<-read.fasta(fastafilename, seqtype="DNA", as.string=TRUE, forceDNAtolower = FALSE)




alignment<-list()
# alignment is going to be you output list of fasta sequences with one entry for each individual, which will be all sequences for that individual concatenated
for(i in 1:nrow(ind)) {
	temp.loci<-NULL

	for(j in 1:length(fa)) {

		contig_name<-strsplit(attributes(fa[[j]])$name,":")[[1]][1] # this depends on how the sequences are named in the fasta file - in my case: <contig_name>:<region>. So I only retain the contig name by splitting based on ":" and taking the first element
		temp.genotypes<-one_SNP_per_contig_adjusted[one_SNP_per_contig_adjusted[,"contig"]==contig_name,c("ref","alt","pos_adjusted",paste(ind[i,1]))]
		temp.seq<-fa[[j]]

		## replace nucleotides in reference fasta ("temp.seq") with genotypes in VCF file
		## 0 -> REF, 2 -> ALT, 1 -> ambiguous
		## IUPAC Ambiguities: M (A or C), R (A or G), W (A or T), S (C or G), Y (C or T), K (G or T)
		for(k in 1:nrow(temp.genotypes)) {
			if(temp.genotypes[k,paste(ind[i,1])]==0) { temp.seq[temp.genotypes[k,"pos_adjusted"]]<-as.character(temp.genotypes[k,"ref"])}
			if(temp.genotypes[k,paste(ind[i,1])]==2) { temp.seq[temp.genotypes[k,"pos_adjusted"]]<-as.character(temp.genotypes[k,"alt"])}
			if(temp.genotypes[k,paste(ind[i,1])]==1 && paste0(temp.genotypes[k,"ref"],temp.genotypes[k,"alt"])=="AC" ) { temp.seq[temp.genotypes[k,"pos_adjusted"]]<-"M"}
			if(temp.genotypes[k,paste(ind[i,1])]==1 && paste0(temp.genotypes[k,"ref"],temp.genotypes[k,"alt"])=="CA" ) { temp.seq[temp.genotypes[k,"pos_adjusted"]]<-"M"}
			if(temp.genotypes[k,paste(ind[i,1])]==1 && paste0(temp.genotypes[k,"ref"],temp.genotypes[k,"alt"])=="AG" ) { temp.seq[temp.genotypes[k,"pos_adjusted"]]<-"R"}
			if(temp.genotypes[k,paste(ind[i,1])]==1 && paste0(temp.genotypes[k,"ref"],temp.genotypes[k,"alt"])=="GA" ) { temp.seq[temp.genotypes[k,"pos_adjusted"]]<-"R"}
			if(temp.genotypes[k,paste(ind[iwrite.fasta(alignment, names=attributes(alignment)$name, file.out="filename.fa")
,1])]==1 && paste0(temp.genotypes[k,"ref"],temp.genotypes[k,"alt"])=="AT" ) { temp.seq[temp.genotypes[k,"pos_adjusted"]]<-"W"}
			if(temp.genotypes[k,paste(ind[i,1])]==1 && paste0(temp.genotypes[k,"ref"],temp.genotypes[k,"alt"])=="TA" ) { temp.seq[temp.genotypes[k,"pos_adjusted"]]<-"W"}
			if(temp.genotypes[k,paste(ind[i,1])]==1 && paste0(temp.genotypes[k,"ref"],temp.genotypes[k,"alt"])=="CG" ) { temp.seq[temp.genotypes[k,"pos_adjusted"]]<-"S"}
			if(temp.genotypes[k,paste(ind[i,1])]==1 && paste0(temp.genotypes[k,"ref"],temp.genotypes[k,"alt"])=="GC" ) { temp.seq[temp.genotypes[k,"pos_adjusted"]]<-"S"}
			if(temp.genotypes[k,paste(ind[i,1])]==1 && paste0(temp.genotypes[k,"ref"],temp.genotypes[k,"alt"])=="CT") { temp.seq[temp.genotypes[k,"pos_adjusted"]]<-"Y"}
			if(temp.genotypes[k,paste(ind[i,1])]==1 && paste0(temp.genotypes[k,"ref"],temp.genotypes[k,"alt"])=="TC") { temp.seq[temp.genotypes[k,"pos_adjusted"]]<-"Y"}
			if(temp.genotypes[k,paste(ind[i,1])]==1 && paste0(temp.genotypes[k,"ref"],temp.genotypes[k,"alt"])=="GT" ) { temp.seq[temp.genotypes[k,"pos_adjusted"]]<-"K"}
			if(temp.genotypes[k,paste(ind[i,1])]==1 && paste0(temp.genotypes[k,"ref"],temp.genotypes[k,"alt"])=="TG" ) { temp.seq[temp.genotypes[k,"pos_adjusted"]]<-"K"}

			}

		temp.loci<-c(temp.loci,unlist(temp.seq)) ## this creates a list with each element being all the 2000 bp regions of a single contig pasted together for a given individual

	}


	alignment[[i]]<-temp.loci
	attributes(alignment)$name[i]<-paste0(ind[i,1])
}

# write.fasta(alignment, names=attributes(alignment)$name, file.out="filename.fa")
# this write the output to file.out
