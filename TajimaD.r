TajimaD<-function(pi,S,n) {
	# a simple function to calculate Tajima's D from nucleotide diversity pi, the number of segregating sites, and the number of sequences
	harmonic<-function(x) { out<-c(); for(i in 1:(x-1)) { out<-c(out,1/i)}; sum_out=sum(out); sum_out} #calculate the (n-1)th harmonic number
	harmonic_sq<-function(x) { out<-c(); for(i in 1:(x-1)) { out<-c(out,1/(i^2))}; sum_out=sum(out); sum_out} #calculate the squared (n-1)th harmonic number
	a1=harmonic(n)
	a2=harmonic_sq(n)
	b1=(n+1)/(3*(n-1))
	b2=(2*((n^2)+n+3))/((9*n)*(n-1))
	c1=b1-(1/a1)
	c2=b2-((n+2)/(a1*n))+(a2/(a1^2))
	e1=c1/a1
	e2=c2/((a1^2)+a2)
	d=(pi-(S/a1)) # calculate difference between observed nucleotide diversity and expected nucleotide diversity
	D=d/sqrt(e1*S+(e2*S)*(S-1)) # calculate Tajima's D
	D
	}
