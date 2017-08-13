## use this functions to output the two haploytpes for each diploid individual in your population for a set of sequences


# read in VCf file and rename column 1 to "contig"
vcf_path="path/to/vcf/your_SNPs.annotated.MAF10.recode.vcf"
vcf_file<-readLines(vcf_path)
vcf_SNPs<-data.frame(vcf_file[grep(pattern="#CHROM",vcf_file):length(vcf_file)])
vcf_SNPs <- data.frame(do.call('rbind', strsplit(as.character(vcf_SNPs[,1]),'\t',fixed=TRUE)))
colnames(vcf_SNPs) <- as.character(unlist(vcf_SNPs[1,]))
vcf_SNPs<-vcf_SNPs[-1,]
colnames(vcf_SNPs)[1]="contig"

# make individual fasta files
library(seqinr)
fastafilename="/path/to/your/genome.fa"
ref_fa<-read.fasta(fastafilename, seqtype="DNA", as.string=TRUE, forceDNAtolower = FALSE)
# subset fasta to include only contigs in the vcf file
# alternatively, you can use samtools faidx to create a subset of your fasta file and read that in using read.fasta, that would be substantially faster
keep_fasta<-which(attr(ref_fa,"name") %in% vcf_SNPs$contig)
new_ref_fa<-list()
for(i in 1:length(keep_fasta)) {
	new_ref_fa[[i]]<-ref_fa[[keep_fasta[i]]]
	}

N=N # 2 times the number of diploid individuals in VCF file
A=A # column in VCF file with the first individual
ind_fa<-list() # create empty list
# ind_fa will be a list with two entries per contig per individual.
for(i in 0:(length(new_ref_fa)-1)) { # start at zero for counting purposes, iterate over contigs
	temp_SNPs<-vcf_SNPs[which(vcf_SNPs$contig == attr(new_ref_fa[[i+1]],"name")),] # subset VCF file for the current contig
		for(j in 1:nrow(temp_SNPs)) {
			for(k in seq(from=1, to=N, by=2)) { # for each entry in row j and in column 1+A
				ind_fa[[((i*N)+k)]]<-new_ref_fa[[i+1]] # create two individual copies of the reference sequence of the current contig, one for each haplotype / fake haploid individual. This is copy 1
				attr(ind_fa[[((i*N)+k)]],"name")<-colnames(temp_SNPs)[(((k+1)*0.5)+A-1)] # copy the individuals' name into the fasta sequence attributes
				ind_fa[[((i*N)+(k+1))]]<-new_ref_fa[[i+1]] # this is copy 2
				attr(ind_fa[[((i*N)+(k+1))]],"name")<-colnames(temp_SNPs)[(((k+1)*0.5)+A-1)] # same individual name
				# now replace the reference allele in the two copies of the reference fasta depending on the genotype in the vcf file: 0/0 keep two ref copies, 0/1: keep one ref copy and rplace the the other with the alternative allele, 1/1: replace with alternative allele in both copies
				if(length(grep("0/0", temp_SNPs[j,(((k+1)*0.5)+A-1)]))==1) { substr(ind_fa[[((i*N)+k)]][1],temp_SNPs$POS[j],temp_SNPs$POS[j])<-as.character(temp_SNPs$REF[j]);substr(ind_fa[[((i*N)+(k+1))]][1],temp_SNPs$POS[j],temp_SNPs$POS[j])<-as.character(temp_SNPs$REF[j]) }
				if(length(grep("0/1", temp_SNPs[j,(((k+1)*0.5)+A-1)]))==1) { substr(ind_fa[[((i*N)+k)]][1],temp_SNPs$POS[j],temp_SNPs$POS[j])<-as.character(temp_SNPs$REF[j]);substr(ind_fa[[((i*N)+(k+1))]][1],temp_SNPs$POS[j],temp_SNPs$POS[j])<-as.character(temp_SNPs$ALT[j]) }
				if(length(grep("1/1", temp_SNPs[j,(((k+1)*0.5)+A-1)]))==1) { substr(ind_fa[[((i*N)+k)]][1],temp_SNPs$POS[j],temp_SNPs$POS[j])<-as.character(temp_SNPs$ALT[j]);substr(ind_fa[[((i*N)+(k+1))]][1],temp_SNPs$POS[j],temp_SNPs$POS[j])<-as.character(temp_SNPs$ALT[j]) }
				}
			}
		}

# write.fasta(ind_fa, names=attributes(ind_fa)$name, file.out="filename.fa")
# this write the output to file.out
