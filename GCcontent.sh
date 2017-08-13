# this calculates the GC content of a (multi)fasta file
grep "^[^>;]" genome.fa > noheaderFasta #no fasta sequence headers
GC="$(tr -d -C CG < noheaderFasta | wc -c)" #count all "G" and "C" characters
ALL="$(tr -d -C CGATN < noheaderFasta | wc -c)" #count all characters, including Ns
echo "scale=4 ; $GC / $ALL" | bc # calculate fraction
