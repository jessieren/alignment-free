

countKmerComputeStats <- function(phageFileName, phageNCBIName, hostFileName, hostNCBIName, phageKmerCountPath, hostKmerCountPath, w, order, multiStatsNames, countKmerCodePath, computeD2MCCodePath)
{
	## count kmer
	for ( k in 1:w )
	{
		
		if( !file.exists( file.path(phageKmerCountPath, paste(phageNCBIName, "_k", k, "_singleStrand_wordcount", sep="") ) ) )
		{
			countKmerCommand <- file.path(countKmerCodePath, paste("countKmer.out -l -k ", k, " -i ", phageFileName, " -s ", phageNCBIName, " -o ", phageKmerCountPath, sep=""))
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
			countKmerCommand <- file.path(computeD2MCCodePath, paste("countKmer.out -l -k ", k, " -i ", hostFileName, " -s ", hostNCBIName, " -o ", hostKmerCountPath, sep=""))
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
	computeSimilarityCommand <- file.path(computeD2MCCodePath, paste("computeD2MC_multiStats_stdout.out -a ", phageNCBIName, " -b ", order, " -c ", hostNCBIName, " -d ", order, " -k ", w, " -i ", phageKmerCountPath, " -j ", hostKmerCountPath, sep=""))
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




## ListA is fatty, ListB is skinny
computeStatsMulti2All <- function(NCBINameListA, NCBINameListB, kmerCountPath, w, order, multiStatsNames, infoListDir, countKmerCodePath, computeD2MCCodePath)
{
	## check whether kmercount files are EXISTING
	#	for( NCBIName in c(as.vector(NCBINameListA), as.vector(NCBINameListB)) )
	#	{
	#		for ( k in 1:w )
	#		{
	#			kmerCountFile <- file.path(kmerCountPath, NCBIName, paste(NCBIName, "_k", k, "_singleStrand_wordcount", sep="") )
	#			if( !file.exists( kmerCountFile ) )
	#			{
	#				stop(paste("ERROR: no existing kmerFile ", kmerCountFile, sep=""))
	#				#countKmerCommand <- paste("/panfs/cmb-panasas2/renj/NCBI_phage_host/code/countKmer.out -l -k ", k, " -i ", phageFileName, " -s ", phageNCBIName, " -o ", phageKmerCountPath, sep="")
	#				#print("===== count kmer phage ======")
	#				#print(countKmerCommand)
	#				#system(countKmerCommand, ignore.stdout = TRUE, ignore.stderr = TRUE)
	#				#Sys.sleep(10)
	#				## the command will be excuted one by one system(countKemrCommand, wait=TRUE)
	#				## you can test with system("sleep 1", wait=TRUE); system("echo test")
	#			}
	#		}
	#	}
	
	## generate the list files
	## infoListDir <- "/auto/cmb-12/fs3/renj/NCBI_phage_host/info/multi2AllList"
	listFilePathA <- generateListFilePath(NCBINameListA, infoListDir, order)
	listFilePathB <- generateListFilePath(NCBINameListB, infoListDir, order)
	
	## compute dissimilarity measures ##
	## d2, N2 is similarity
	## phage - host
	## 20150524 paste("k", w, "_order", order, sep="")
	#similarityDirPhage <- file.path(similarityDir, phageNCBIName, paste("k", w, "_order", order, sep=""))
	#multiStatsFileName <- paste("k", w, "_", phageNCBIName, "_order", order, "_", hostNCBIName, "_order", order, "_statScores", sep="")
	#multiStatsPathFileName <- file.path(similarityDirPhage, multiStatsFileName)
	
	#if( !file.exists(multiStatsPathFileName) )
	#{
	print("===== compute multiStats ======")
	#computeSimilarityCommand <- paste("/auto/cmb-12/fs3/renj/NCBI_phage_host/code/computeD2MC_multiStats_multi2All_stdout.out -k ", w, " -i ", listFilePathA, " -j ", listFilePathB, sep="")
	#print(computeSimilarityCommand)
	#system(computeSimilarityCommand, ignore.stdout = TRUE, ignore.stderr = TRUE)
	#stdout <- system(computeSimilarityCommand, intern = TRUE, ignore.stdout = FALSE)
	cprogram <- file.path(computeD2MCCodePath, "computeD2MC_multiStatsCombine_multi2All_stdout.out")
	computeSimilarityCommand <- paste(" -k ", w, " -i ", listFilePathA, " -j ", listFilePathB, sep="")
	stdout <- system2(cprogram, computeSimilarityCommand, stdout=TRUE)
	statScores <- array(NA, dim=c(length(NCBINameListA), length(NCBINameListB), length(multiStatsNames)), dimnames=list(NCBINameListA, NCBINameListB, multiStatsNames) )
	for( l in stdout)
	{
		speciesLine <- strsplit(l, ":")[[1]]
		if( speciesLine[1]=="speciesA" )
		{
			IDA <- as.numeric(speciesLine[2])+1
		}else if( speciesLine[1]=="speciesB" )
		{
			IDB <- as.numeric(speciesLine[2])+1
		}else{
			scoreLine <- strsplit(l, ", ")[[1]]
			stat <- scoreLine[1]
			value <- as.numeric(scoreLine[2])
			if( stat %in% multiStatsNames & !is.na(value) )
			{
				statScores[IDA, IDB, stat] <- value
			}
		}
	}
	statScores
}






## To use, for example,


w <- 6
order <- 0


setwd("/panfs/cmb-panasas2/renj/NCBI_phage_host/results")
hostFilePath <- "/panfs/cmb-panasas2/renj/NCBI_phage_host/host"
phageFilePath <- "/panfs/cmb-panasas2/renj/NCBI_phage_host/phage"
kmerCountPath <- "/panfs/cmb-panasas2/renj/NCBI_phage_host/kmerCount"
countKmerCodePath <- "/panfs/cmb-panasas2/renj/NCBI_phage_host/code"
computeD2MCCodePath <- "/panfs/cmb-panasas2/renj/NCBI_phage_host/code"


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

countKmerComputeStats(phageFileName, phageNCBIName, hostFileName, hostNCBIName, phageKmerCountPath, hostKmerCountPath, w, order, multiStatsNames, countKmerCodePath, computeD2MCCodePath)



