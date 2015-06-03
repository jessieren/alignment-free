


countKmerComputeStats <- function(phageFileName, phageNCBIName, hostFileName, hostNCBIName, phageKmerCountPath, hostKmerCountPath, w, order, multiStatsNames)
{
	## count kmer
	for ( k in 1:w )
	{
		
		if( !file.exists( file.path(phageKmerCountPath, paste(phageNCBIName, "_k", k, "_singleStrand_wordcount", sep="") ) ) )
		{
			countKmerCommand <- paste("/panfs/cmb-panasas2/renj/NCBI_phage_host/code/countKmer.out -l -k ", k, " -i ", phageFileName, " -s ", phageNCBIName, " -o ", phageKmerCountPath, sep="")
			print("===== count kmer phage ======")
			#print(countKmerCommand)
			system(countKmerCommand, ignore.stdout = TRUE, ignore.stderr = TRUE)
			#Sys.sleep(10)
			## the command will be excuted one by one system(countKemrCommand, wait=TRUE)
			## you can test with system("sleep 1", wait=TRUE); system("echo test")
		}
	}
	
	for ( k in 1:w )
	{
		if( !file.exists( file.path(hostKmerCountPath, paste(hostNCBIName, "_k", k, "_singleStrand_wordcount", sep="") ) ) )
		{
			countKmerCommand <- paste("/panfs/cmb-panasas2/renj/NCBI_phage_host/code/countKmer.out -l -k ", k, " -i ", hostFileName, " -s ", hostNCBIName, " -o ", hostKmerCountPath, sep="")
			print("===== count kmer host ======")
			#print(countKmerCommand)
			system(countKmerCommand, ignore.stdout = TRUE, ignore.stderr = TRUE)
			#Sys.sleep(10)
		}
	}
	
	
	## compute dissimilarity measures ##
	## d2, N2 is similarity
	## phage - host

	print("===== compute multiStats ======")
	computeSimilarityCommand <- paste("/panfs/cmb-panasas2/renj/NCBI_phage_host/code/computeD2MC_multiStats_stdout.out -a ", phageNCBIName, " -b ", order, " -c ", hostNCBIName, " -d ", order, " -k ", w, " -i ", phageKmerCountPath, " -j ", hostKmerCountPath, sep="")
	print(computeSimilarityCommand)
	stdout <- system(computeSimilarityCommand, intern = TRUE)
	statScores <- rep(NA, length(multiStatsNames))
	names(statScores) <- multiStatsNames
	for( l in stdout)
	{
		stat <- strsplit(l, ", ")[[1]][1]
		value <- as.numeric(strsplit(l, ", ")[[1]][2])
		if( stat == multiStatsNames[which(l==stdout)] & !is.na(value) )
		{
			statScores[stat] <- value
		}
	}
	as.matrix(statScores)

	
}






## To use, for example,


w <- 6
order <- 0

setwd("/panfs/cmb-panasas2/renj/NCBI_phage_host/results")
hostFilePath <- "/panfs/cmb-panasas2/renj/NCBI_phage_host/host"
phageFilePath <- "/panfs/cmb-panasas2/renj/NCBI_phage_host/phage"
kmerCountPath <- "/panfs/cmb-panasas2/renj/NCBI_phage_host/kmerCount"

multiStatsCandi <- c("D2", "d2", "D2star", "d2star", "D2shepp", "d2shepp", "D2NGS", "d2NGS", "D2starNGS", "d2starNGS", "D2sheppNGS", "d2sheppNGS", "S2", "hao", "JS", "Eu", "Ma", "Ch", "EuF", "willner")
if ( w < 3 )
{
	## no HAO
	multiStatsNames <- multiStatsCandi[-which("hao"==multiStatsCandi)]
	
} else if (w > 4) {
	## no willner
	multiStatsNames <- multiStatsCandi[-which("willner"==multiStatsCandi)]
	
} else {
	## k=3, 4
	multiStatsNames <- multiStatsCandi
}


phageNCBIName <- "JN698997"
phageFileName <- file.path(phageFilePath, paste(phageNCBIName, ".fasta", sep=""))
phageKmerCountPath <- file.path(kmerCountPath, phageNCBIName)

hostNCBIName <- "NZ_CP009494.1"
hostFileName <- file.path(hostFilePath, paste(hostNCBIName, ".fasta", sep=""))
hostKmerCountPath <- file.path(kmerCountPath, hostNCBIName)

countKmerComputeStats(phageFileName, phageNCBIName, hostFileName, hostNCBIName, phageKmerCountPath, hostKmerCountPath, w, order, multiStatsNames)



