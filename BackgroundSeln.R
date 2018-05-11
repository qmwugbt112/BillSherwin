# Simulation of background selection without recombination
# to see if in the short term the most polymorphic loci
# show associative overdominance

numNeutral <- 200	# number of neutral loci
numSelectd <- 200 	# number under selection

uRateN <- 10^-4		# neutral mutation rate
uRateS <- 10^-4		# seleted mutation rate

ne2 <- 100			# 2 x Ne

maxGens <- 30*ne2	# number of gens to simulate

# Distribution of selected effects
# parameters 2,10^-3 give mode s=10^-3 few beyond 6x that
sDistn <- function(n)rgamma(n, shape=2,scale=10^-3)

# Each col is the selection vs a hom at the corresponding
# locus.  effects0 has no selection effects does have it 
effects0<- matrix(0,#rep(sDistn(numSelectd),each=ne2/2),
						nrow=ne2/2,
						ncol=numSelectd)

effects<- matrix(rep(sDistn(numSelectd),each=ne2/2),
						nrow=ne2/2,
						ncol=numSelectd)						

# matrix of neutral and selected alleles
nMat <- matrix(	sample(c(T,F),ne2*numNeutral,replace=T), 
				nrow=ne2, 
				ncol=numNeutral)
sMat <- matrix(	sample(c(T,F),ne2*numSelectd,replace=T), 
				nrow=ne2, 
				ncol=numSelectd)
				
ancestor <- 1:ne2

nN 	<- length(nMat)	
nS	<- length(sMat)
nIndex <- 1:nN
sIndex <- 1:nS
sHalfI <- 1:(ne2/2)

pat <- seq(1,ne2,2)

totS <-NULL
totN <-NULL
totA <-NULL


# burn-in run to the point when only the descendants of
# one ancestral haplotype are alive
for (g in 1:maxGens){
	
	# mutate both types of locus
	mutnTargets <- sample(nIndex,rpois(1,uRateN*nN))
	nMat[mutnTargets] <- !nMat[mutnTargets]
	mutnTargets <- sample(nIndex,rpois(1,uRateS*nS))
	sMat[mutnTargets] <- !sMat[mutnTargets]

	# apply selection on homozygotes
	s <- (sMat[pat,] & sMat[-pat,])*effects
	w <- apply(s,1,function(x) prod(1-x))
	breeders <- sample(	1:ne2,
						ne2,
						replace=T,
						prob=rep(w,each=2)
						)
	nMat <- nMat[breeders,]
	sMat <- sMat[breeders,]	
	ancestor <- ancestor[breeders]
	
	totS <- c(totS,dim(unique(sMat))[1])
	totN <- c(totN,dim(unique(nMat))[1])
	totA <- c(totA,length(unique(ancestor)))
	#if (length(unique(ancestor))==1) {
	#	print('generation ', g)
	#	finalGen <- g+1
	#	break					
	#	}
	finalGen <- g+1
}

for (g in 1:maxGens){
	
	# mutate both types of locus
	mutnTargets <- sample(nIndex,rpois(1,uRateN*nN))
	nMat[mutnTargets] <- !nMat[mutnTargets]
	mutnTargets <- sample(nIndex,rpois(1,uRateS*nS))
	sMat[mutnTargets] <- !sMat[mutnTargets]

	# apply selection on homozygotes
	s <- (sMat[pat,] & sMat[-pat,])*effects
	w <- apply(s,1,function(x) prod(1-x))
	breeders <- sample(	1:ne2,
						ne2,
						replace=T,
						prob=rep(w,each=2)
						)
	nMat <- nMat[breeders,]
	sMat <- sMat[breeders,]	
	ancestor <- ancestor[breeders]
	
	totS <- c(totS,dim(unique(sMat))[1])
	totN <- c(totN,dim(unique(nMat))[1])
	totA <- c(totA,length(unique(ancestor)))
	if (length(unique(ancestor))==1) {
		print(paste('burn-in for ', g, 'generations'))
		finalGen <- g+1
		break					
		}
	finalGen <- g+1
}

if (finalGen >= maxgens) stop('burn-in unsuccessful')
ancestor <- 1:ne2

for (g in 1:maxGens){
	
	# mutate both types of locus
	mutnTargets <- sample(nIndex,rpois(1,uRateN*nN))
	nMat[mutnTargets] <- !nMat[mutnTargets]
	mutnTargets <- sample(nIndex,rpois(1,uRateS*nS))
	sMat[mutnTargets] <- !sMat[mutnTargets]

	# apply selection on homozygotes
	s <- (sMat[pat,] & sMat[-pat,])*effects
	w <- apply(s,1,function(x) prod(1-x))
	breeders <- sample(	1:ne2,
						ne2,
						replace=T,
						prob=rep(w,each=2)
						)
	nMat <- nMat[breeders,]
	sMat <- sMat[breeders,]	
	ancestor <- ancestor[breeders]
	
	totS <- c(totS,dim(unique(sMat))[1])
	totN <- c(totN,dim(unique(nMat))[1])
	totA <- c(totA,length(unique(ancestor)))
	if (length(unique(ancestor))==1) {
		print(paste('burn-in for ', g, 'generations'))
		finalGen <- g+1
		break					
		}
	finalGen <- g+1
}


plot(1:finalGen,totN, type='l', ylim=c(0,max(c(totS,totN))/2))
lines(1:finalGen,totS, col='red')
lines(1:finalGen,totA, col='blue')
abline(h=7.5)
abline(h=1,col='blue')

# calculate the frequency distribution of haplotypes
hapDistn <-function(mat){
	# find the unique haplotypes
	uHap <- unique(mat)
	countVec <- rep(0,nrow(uHap))
	for (r in 1:nrow(uHap)){
		countVec[r] <- sum(
			apply(mat,1,function(x) identical(x,uHap[r,]))
			)
						
	}
	return(countVec)
}

hist(hapDistn(nMat),
		main=paste('Haplotype frequencies (',totA[finalGen],' ancestral )'),
		xlab='Haplotype count',
		ylab='Number of haplotypes',
		breaks=0:100)
hist(colSums(nMat),
		main='Allele frequencies',
		xlab='Allele Count',
		ylab='Number of loci',
		breaks=0:100)


