# generate vector of BLAST output files
# used rpf and lysozyme sequences in "lennon2blastdb.fasta" as database
# blasted translated gene sequences from prokka annotated soil genomes
hitFiles=list.files(pattern="hits")

# create vector to append rpf hits to
allRPFhits="A"

# loop through genomes and filter BLAST results
###### Casting broad net
for(i in 1:length(hitFiles)){
# columns:
# 1 - query id
# 2 - subject id
# 3 - % identity
# 4 - alignment length
# 5 - mismatches
# 6 - gap opens
# 7 - query start
# 8 - query end
# 9 - subject start
# 10 - subject end
# 11 - evalue
# 12 - bit score
# 13 - sequence

	# load current BLAST results file
	hits=read.table(hitFiles[i],header=FALSE,stringsAsFactors=FALSE)
	# remove hits with an evalue>1e-10
	hits=hits[hits[,11]<1e-10,]

	# vector of unique genome hits
	putativeRPFgenes=sort(unique(hits[,1]))

	# generate dataframe to store the best match for each putativeRPF gene
	bestHits=hits[1,]
	for(j in 1:length(putativeRPFgenes)){
		cur=hits[hits[,1]==putativeRPFgenes[j],]
		if(sum(cur[,11]==min(cur[,11]))==1){
			bestHits=rbind(bestHits,cur[cur[,11]==min(cur[,11]),])
		}else if(sum(cur[,11]==min(cur[,11]))>1){
			cur=cur[order(cur[,11]),]
			bestHits=rbind(bestHits,cur[1,])
		}
	}
	bestHits=bestHits[-1,]

	rpfHits=bestHits[grep("rpf",bestHits[,2]),]
	
	allRPFhits=c(allRPFhits,rpfHits[,1])
}

allRPFhits=allRPFhits[-1]

# tabulate hits per genome
table(substr(allRPFhits,1,7))



###### finding best few hits
# create vector to append rpf hits to
allRPFhits="A"

for(i in 1:length(hitFiles)){
# columns:
# 1 - query id
# 2 - subject id
# 3 - % identity
# 4 - alignment length
# 5 - mismatches
# 6 - gap opens
# 7 - query start
# 8 - query end
# 9 - subject start
# 10 - subject end
# 11 - evalue
# 12 - bit score
# 13 - sequence

	# load current BLAST results file
	hits=read.table(hitFiles[i],header=FALSE,stringsAsFactors=FALSE)
	# remove hits with an evalue>1e-10
	hits=hits[hits[,11]<1e-200,]
	
	if(nrow(hits)>0){
		# vector of unique genome hits
		putativeRPFgenes=sort(unique(hits[,1]))

		# generate dataframe to store the best match for each putativeRPF gene
		bestHits=hits[1,]
		for(j in 1:length(putativeRPFgenes)){
			cur=hits[hits[,1]==putativeRPFgenes[j],]
			if(sum(cur[,11]==min(cur[,11]))==1){
				bestHits=rbind(bestHits,cur[cur[,11]==min(cur[,11]),])
			}else if(sum(cur[,11]==min(cur[,11]))>1){
				cur=cur[order(cur[,11]),]
				bestHits=rbind(bestHits,cur[1,])
			}
		}
		bestHits=bestHits[-1,]

		rpfHits=bestHits[grep("rpf",bestHits[,2]),]
	
		allRPFhits=c(allRPFhits,rpfHits[,1])
	}
}

allRPFhits=allRPFhits[-1]

# tabulate hits per genome
table(substr(allRPFhits,1,7))

#### pulling these hits from open reading frames from prokka annotated genomes
print(allRPFhits)
#### this part not on github because of size of genome files
