
###############################################################################################################
#	I- Filtrage duplicats, selection des reads de junction et identification des bornes des reads

#	/software/bin/run_cluster.pl -cpu=40 -cmd="Rscript  /data-isilon/Cagnard/RNAseqSEA/RNAseqSEA_AllTnjs/RNAseqSEA_AllTnjs.r idprojet=NGS2023_6585  align=hisat2 esp=HG19_MT ReName=/data-isilon/sequencing/ngs/NGS2023_6585/HG19_MT/Echs.txt"

  options(scipen=999)

	args = commandArgs(trailingOnly=TRUE)
	saveRDS(args, "./ArgsRNAseqSEAallTnjs.rds")
	#	args = readRDS("./ArgsRNAseqSEAallTnjs.rds")

	suprTmpFiles=TRUE
	source("/data-isilon/Cagnard/RNAseqSEA/RNAseqSEA_AllTnjs/RNAseqSEA_AllTnjs_biblio.r")
	tmpArgs = t(rbind(apply(as.matrix(unlist(args)), 1, function(x) unlist(strsplit(x, "=")))))
	if(sum(grepl("idprojet", tmpArgs[,1]))>0)	idprojet = tmpArgs[grep("idprojet", tmpArgs[,1], ignore.case=TRUE),2]
	esp = GenVer(idprojet)
	align="star"; ReName=NULL; 	nCPU=40; limDecal=5; limCountJun=2; RMdup=FALSE

	if(sum(grepl("^esp$", tmpArgs[,1]))>0)	esp = tmpArgs[grep("^esp$", tmpArgs[,1], ignore.case=TRUE),2]
	if(sum(grepl("nCPU", tmpArgs[,1]))>0)	nCPU = as.numeric(tmpArgs[grep("nCPU", tmpArgs[,1], ignore.case=TRUE),2])
	if(sum(grepl("align", tmpArgs[,1]))>0)	align = tmpArgs[grep("align", tmpArgs[,1], ignore.case=TRUE),2]
	if(sum(grepl("limDecal", tmpArgs[,1]))>0)	limDecal = as.numeric(tmpArgs[grep("limDecal", tmpArgs[,1], ignore.case=TRUE),2])
	if(sum(grepl("limCountJun", tmpArgs[,1]))>0)	limCountJun = as.numeric(tmpArgs[grep("limCountJun", tmpArgs[,1], ignore.case=TRUE),2])
	if(sum(grepl("ReName", tmpArgs[,1]))>0)	ReName = tmpArgs[grep("ReName", tmpArgs[,1], ignore.case=TRUE),2]
	if(sum(grepl("RMdup", tmpArgs[,1]))>0) RMdup = as.logical(tmpArgs[grep("RMdup", tmpArgs[,1], ignore.case=TRUE),2])
	if(sum(grepl("suprTmpFiles", tmpArgs[,1]))>0) suprTmpFiles = as.logical(tmpArgs[grep("suprTmpFiles", tmpArgs[,1], ignore.case=TRUE),2])
	nfos = c(idprojet, esp, align, Sys.info()[c("nodename", "login", "user", "effective_user")], getwd())
	slogRNAseqSEA_AllTnjs(nfos)
	
	RNAseqSEApath = "/data-isilon/Cagnard/RNAseqSEA/"
	if(!is.null(ReName))	BamsReNames = as.matrix(read.table(ReName, sep="\t", header=TRUE))
	
	nfo = c("in total", "duplicates", "mapped", "properly paired")
	stats = matrix(0, ncol=length(nfo), nrow=0)
	colnames(stats) = nfo
	pathRes = paste("/data-isilon/sequencing/ngs/", idprojet, "/", esp, "/RNAseqSEA/", sep="")
	if(!file.exists(pathRes))	dir.create(pathRes)
	pathBamJuncs = paste(pathRes, "/BamJuncs/", sep="")
	if(!file.exists(pathBamJuncs))	dir.create(pathBamJuncs)
	pathBed = paste(pathRes, "/BedJuncs/", sep="")
	if(!file.exists(pathBed))	dir.create(pathBed)

###############################################################################################################
#	I Splits les bams
	
write(paste("\n\n#######################################\n#\n#\t\t I Splits les bams \n#\n########################################", sep=""), file="")

	write(paste("\tSelection des reads de junction et filtrage des duplicates", sep=""), file="")
	sourceBamsPath = paste("/data-isilon/sequencing/ngs/", idprojet, "/", esp, "/align/", align, "/", sep="")
	bams = list.files(sourceBamsPath, full.names=TRUE, pattern="\\.bam$")
	if(length(bams)==0){stop(paste("\t\t!!! Aucun bam dans ", sourceBamsPath, " !!!", sep=""))}

	addChr=""
	if((length(bams)>0)&(sum(grepl("chr", system(paste("samtools idxstats -@ ", nCPU, " ", bams[1], sep=""), intern=TRUE)))>0)) addChr = "chr"
		
	if(RMdup)
	{
		rmdupBamPath = paste(pathRes, "/rmDupsBams/", sep="")
		if(!file.exists(rmdupBamPath))	dir.create(rmdupBamPath)
		
		runif(1)
		library(parallelMap)
		parallelStart(mode = "multicore", cpus=5, show.info=TRUE) 
		f = function(I) rmdup(I)
		y = parallelMap(f, c(1:length(bams)))
		parallelStop()
		
		bams = list.files(rmdupBamPath, full.names=TRUE, pattern="\\.bam$")
	}

	# testBamChr
	addChr = ""
	testChrs = system(paste("samtools idxstats ", bams[1], " | cut -f 1", sep=""), intern=TRUE)
	if(sum(grepl("^chr", testChrs))>=1)  addChr = "chr"
	
	splitBamPath = paste(pathRes, "/SplitBams/", sep="")
	if(!file.exists(splitBamPath))	dir.create(splitBamPath)
	
	runif(1)
	library(parallelMap)
	parallelStart(mode = "multicore", cpus=20, show.info=TRUE) 
	f = function(B) splitBamSamTools(B)
	y = parallelMap(f, c(1:length(bams)))
	parallelStop()
	
	bamsSplit = list.files(splitBamPath, recursive=TRUE, full.names=TRUE, pattern="\\.bam$")
	#	samples = list.files(splitBamPath, recursive=FALSE, full.names=FALSE)
	samples = gsub(".bam$", "", basename(bams))
	
###############################################################################################################
#	II Selection des reads de jonctions => bed

write(paste("\n\n#######################################\n#\n#\t\t II Selection des reads de jonctions => bed \n#\n########################################", sep=""), file="")

	runif(1)
	library(parallelMap)
	parallelStart(mode = "multicore", cpus=10, show.info=TRUE) 
	f = function(B) selectJunc(B)
	y = parallelMap(f, c(1:length(bamsSplit)))
	parallelStop()
	
	filenames = list.files(pathBed, full.names=TRUE)
	empty = filenames[file.size(filenames) == 0L]
	unlink(empty)
	
###############################################################################################################
#	III bed en RDS (=> divise en plusieurs parties ? => MT)

write(paste("\n\n#######################################\n#\n#\t\t III bed en RDS \n#\n########################################", sep=""), file="")

	listBed = list.files(pathBed, full.names=TRUE, pattern="\\.bed$")

	runif(1)
	library(parallelMap)
	parallelStart(mode = "multicore", cpus=25, show.info=TRUE) 
	f = function(S) bed2rds(S)
	y = parallelMap(f, 1:length(listBed))
	parallelStop() 
	
	#	if(suprTmpFiles) unlink(pathBamJuncs, recursive = TRUE)
	
###############################################################################################################
#	IV calcul des positions des junctions	=> fichier a conserver avec les coords de juncs brutes

write(paste("\n\n#######################################\n#\n#\t\t IV calcul des positions des junctions et comptage \n#\n########################################", sep=""), file="")

	listBedRds = list.files(pathBed, full.names=TRUE, pattern="\\.rds$")
	AllJuncPath = paste(pathRes, "/AllJunctionsCounts/", sep="")
	if(!file.exists(AllJuncPath))	dir.create(AllJuncPath)
	
	runif(1)
	library(parallelMap)
	parallelStart(mode = "multicore", cpus=40, show.info=TRUE) 
	f = function(P) coordsJunc(P)
	y = parallelMap(f, c(1:length(listBedRds)))
	parallelStop()

	if(suprTmpFiles) unlink(pathBed, recursive = TRUE)
	
###############################################################################################################
#	V Selection de la liste de genes

write(paste("\n\n#######################################\n#\n#\t\t V Selection de la liste de genes \n#\n########################################", sep=""), file="")

	if(sum(grepl("HG38", esp))==1)
	{
		GCv = "43"
		#	Old_juncBED = readRDS("/data-isilon/Cagnard/RNAseqSEA/rds/Junc_HG38_gencode40.rds")
		juncBED = readRDS(paste(RNAseqSEApath, "/Refs/Junc_HG38_gencode", GCv, ".rds", sep=""))
		#	old_refBED  = readRDS("/data-isilon/Cagnard/RNAseqSEA/rds/HG38_gencode40.rds")
		refBED = readRDS(paste(RNAseqSEApath, "/Refs/HG38_gencode", GCv, ".rds", sep=""))
		nomStruct = paste( esp, "_gencode", GCv, sep="")
	}
	if(sum(grepl("HG19", esp))==1)
	{
		GCv = "43"
		juncBED = readRDS(paste(RNAseqSEApath, "/Refs/Junc_HG19_gencode", GCv, ".rds", sep=""))
		refBED = readRDS(paste(RNAseqSEApath, "/Refs/HG19_gencode", GCv, ".rds", sep=""))
		nomStruct = paste( esp, "_gencode", GCv, sep="")
	}
	if(sum(grepl("MM38", esp))==1)
	{
		GCv = "25"
		juncBED = readRDS(paste(RNAseqSEApath, "/Refs/Junc_MM38_gencode", GCv, ".rds", sep=""))
		refBED = readRDS(paste(RNAseqSEApath, "/Refs/MM38_gencode", GCv, ".rds", sep=""))
		nomStruct = paste( esp, "_gencode", GCv, sep="")
	}
	if(sum(grepl("MM39", esp))==1)
	{
		GCv = "32"
		juncBED = readRDS(paste(RNAseqSEApath, "/Refs/Junc_MM39_gencode", GCv, ".rds", sep=""))
		refBED = readRDS(paste(RNAseqSEApath, "/Refs/MM39_gencode", GCv, ".rds", sep=""))
		nomStruct = paste( esp, "_gencode", GCv, sep="")
	}

	geneListsPath = paste(pathRes, "/GenesLists/", sep="")
	if(!file.exists(geneListsPath))	dir.create(geneListsPath)

	starts = refBED[order(refBED[, "exon_chrom_start"]),]
	starts = starts[!duplicated(starts[,"ensembl_gene_id"]),]
	starts = starts[order(starts[,"ensembl_gene_id"]),]
	ends = refBED[order(refBED[, "exon_chrom_end"], decreasing=TRUE),]
	ends = ends[!duplicated(ends[,"ensembl_gene_id"]),]
	ends = ends[order(ends[,"ensembl_gene_id"]),]
	bedRef = cbind(starts[,c("chromosome_name", "exon_chrom_start")], ends[,c("exon_chrom_end", "ensembl_gene_id")])
	bedRef = bedRef[order(bedRef[,"exon_chrom_start"]),]
	bedRef = bedRef[order(bedRef[,"chromosome_name"]),]
	
#	selection des ENSG ID des gènes qui ont des junctions
	listCoordsBed = list.files(AllJuncPath, full.names=TRUE, pattern="\\.rds$")

	runif(1)
	library(parallelMap)
	parallelStart(mode = "multicore", cpus=nCPU, show.info=TRUE) 
	f = function(C) caracterizeJuncs(C)
	y = parallelMap(f, c(1:length(listCoordsBed)))
	parallelStop()

#	Concatene les listes de gènes

write(paste("\n\n#######################################\n#\n#\t\t\t Concatene les listes de gènes \n#\n########################################", sep=""), file="")

	AllGenesLists = list.files(geneListsPath, full.names=TRUE, pattern="\\.rds$")

	for(S in 1:length(samples))
	{
		write(paste("\t ", S, "/", length(samples), " - ", samples[S], " : Concatene genes lists ", sep=""), file="")
		
		tmpGenesLists = AllGenesLists[grepl(paste("^", samples[S], "_", sep=""), basename(AllGenesLists))]
		
		PatGenesList = NULL
		for(tGL in 1:length(tmpGenesLists))	PatGenesList = c(PatGenesList, readRDS(tmpGenesLists[tGL]))
		
		PatGenesList = unique(PatGenesList)
		
		saveRDS(PatGenesList, file=paste(geneListsPath, "/", samples[S], "_GenesLists.rds", sep=""))
		unlink(tmpGenesLists)
	}

#	concat les bedsRDS

write(paste("\n\n#######################################\n#\n#\t\t\t Concatene les bedsRDS \n#\n########################################", sep=""), file="")

	for(S in 1:length(samples))
	{
		write(paste("\t ", S, "/", length(samples), " - ", samples[S], " : Concatene junctions counts ", sep=""), file="")
		tmpBeds = listCoordsBed[grepl(paste("^", samples[S], "_", sep=""), basename(listCoordsBed))]
		PatRDS = readRDS(tmpBeds[1])
		if(length(tmpBeds)>1)	for(B in 2:length(tmpBeds))
		{
			PatRDS = rbind(PatRDS, readRDS(tmpBeds[B]))
		}	
		saveRDS(PatRDS, file=paste(AllJuncPath, "/", samples[S], "_junctions.rds", sep=""))
		unlink(tmpBeds)
	}

###############################################################################################################
#	VI	Caracterisation des jonctions

write(paste("\n\n#######################################\n#\n#\t\t VI	Caracterisation des jonctions \n#\n########################################", sep=""), file="")

	AllresPath = paste(pathRes, "/AllRes/", sep="")
	if(!file.exists(AllresPath))	dir.create(AllresPath)
	ErrorPath =  paste(AllresPath, "/Errors/", sep="")
	if(!file.exists(ErrorPath))	dir.create(ErrorPath)
	
	# Annots = readRDS("/data-isilon/Cagnard/MagicMorgan/Methode_RNAseq_DevL/Annots/AllEns.Rds")

	pathResRI = paste(pathRes, "/AllresRI/", sep="")
	if(!file.exists(pathResRI))	dir.create(pathResRI)
	for(S in 1:length(samples))	dir.create(paste(pathResRI, "/", samples[S], "/", sep=""))
	
	pathResSE = paste(pathRes, "/AllresSE/", sep="")
	if(!file.exists(pathResSE))	dir.create(pathResSE)
	for(S in 1:length(samples))	dir.create(paste(pathResSE, "/", samples[S], "/", sep=""))
	
#	table Pat - Gene
	ToDoTab = matrix(0, ncol=2, nrow=0)
	geneListsPath = paste(pathRes, "/GenesLists/", sep="")
	AllGenesLists = list.files(geneListsPath, full.names=TRUE, pattern="\\.rds$")
	for(GL in 1:length(AllGenesLists))
	{
		genesList = readRDS(AllGenesLists[GL])
		nomPat = gsub("_GenesLists.rds", "", basename(AllGenesLists[GL]))
		ToDoTab = rbind(ToDoTab, cbind(rep(nomPat, length(genesList)), genesList))
	}

	#	source("/data-isilon/Cagnard/RNAseqSEA/RNAseqSEA_AllTnjs/RNAseqSEA_AllTnjs_biblio.r")
	runif(1)
	library(parallelMap)
	parallelStart(mode = "multicore", cpus=nCPU, show.info=TRUE) 
	f = function(G) analysePatGene(G) #	testanalysePatGene(G) 
	y = parallelMap(f, c(1:nrow(ToDoTab)))
	parallelStop()
	
	if(suprTmpFiles) unlink(geneListsPath, recursive=TRUE)
	
###############################################################################################################
#	VII	Condense les résultats par patient et type
	
write(paste("\n\n#######################################\n#\n#\t\t VII	Condense les résultats par patient et type \n#\n########################################", sep=""), file="")

	if(FALSE)
	{
	  args = readRDS("./ArgsRNAseqSEAallTnjs.rds")
  	suprTmpFiles=TRUE
  	source("/data-isilon/Cagnard/RNAseqSEA/RNAseqSEA_AllTnjs/RNAseqSEA_AllTnjs_biblio.r")
  	tmpArgs = t(rbind(apply(as.matrix(unlist(args)), 1, function(x) unlist(strsplit(x, "=")))))
  	if(sum(grepl("idprojet", tmpArgs[,1]))>0)	idprojet = tmpArgs[grep("idprojet", tmpArgs[,1], ignore.case=TRUE),2]
  	if(sum(grepl("^esp$", tmpArgs[,1]))>0)	esp = tmpArgs[grep("^esp$", tmpArgs[,1], ignore.case=TRUE),2]
  	if(sum(grepl("nCPU", tmpArgs[,1]))>0)	nCPU = as.numeric(tmpArgs[grep("nCPU", tmpArgs[,1], ignore.case=TRUE),2])
  	if(sum(grepl("align", tmpArgs[,1]))>0)	align = tmpArgs[grep("align", tmpArgs[,1], ignore.case=TRUE),2]
  	if(sum(grepl("suprTmpFiles", tmpArgs[,1]))>0) suprTmpFiles = as.logical(tmpArgs[grep("suprTmpFiles", tmpArgs[,1], ignore.case=TRUE),2])
  	sourceBamsPath = paste("/data-isilon/sequencing/ngs/", idprojet, "/", esp, "/align/", align, "/", sep="")
  	bams = list.files(sourceBamsPath, full.names=TRUE, pattern="\\.bam$")
  	samples = gsub(".bam$", "", basename(bams))
  	
  	pathRes = paste("/data-isilon/sequencing/ngs/", idprojet, "/", esp, "/RNAseqSEA/", sep="")
  	AllresPath = paste(pathRes, "/AllRes/", sep="")
	}
	
	tabResCondens = rbind(cbind(samples, rep("SE", length(samples))), cbind(samples, rep("RI", length(samples))))
	pathResAnalysis = paste("/data-isilon/sequencing/ngs/", idprojet, "/", esp, "/analysis/", sep="")
	if(!file.exists(pathResAnalysis))	dir.create(pathResAnalysis)
	pathAllResAnalysis = paste(pathResAnalysis, "/AllRes/", sep="")
	if(!file.exists(pathAllResAnalysis))	dir.create(pathAllResAnalysis)
	
	runif(1)
	library(parallelMap)
	parallelStart(mode = "multicore", cpus=40, show.info=TRUE) 
	f = function(R) condenseResPat(R)
	y = parallelMap(f, c(1:nrow(tabResCondens)))
	parallelStop()

	#unlink(pathResSE, recursive=TRUE)
	#unlink(pathResRI, recursive=TRUE)
	
#	Condense les résultats par type
	AllResRDSpatRIfiles = list.files(AllresPath, full.names=TRUE, pattern="^allRes_RI*.")
	allResRI = NULL
	for(N in 1:length(AllResRDSpatRIfiles)){allResRI[[length(allResRI)+1]] = readRDS(AllResRDSpatRIfiles[N]); write(paste(gsub(".rds", "", basename(AllResRDSpatRIfiles[N])), " (", N, "/", length(AllResRDSpatRIfiles), ")", sep=""), file="")}
	allResRI = do.call("rbind",allResRI)
	allResRI = allResRI[order(allResRI[,min(grep("_Start", colnames(allResRI)))]),,drop=FALSE]
	allResRI = allResRI[order(allResRI[,grep("^Chr$", colnames(allResRI))]),,drop=FALSE]
	#evs = apply(allResRI, 1, function(x) paste(x[1], x[length(x)], sep="_", collapse="_"))	
	#allResRI = allResRI[!duplicated(evs),,drop=FALSE]
	saveRDS(allResRI, file=paste(pathAllResAnalysis, "allResRI.rds", sep=""))
	
	AllResRDSpatSEfiles = list.files(AllresPath, full.names=TRUE, pattern="^allRes_SE*.")
	allResSE = NULL
	for(N in 1:length(AllResRDSpatSEfiles)){allResSE[[length(allResSE)+1]] = readRDS(AllResRDSpatSEfiles[N]); write(paste(gsub(".rds", "", basename(AllResRDSpatSEfiles[N])), " (", N, "/", length(AllResRDSpatSEfiles), ")", sep=""), file="")}
	allResSE = do.call("rbind",allResSE)
	allResSE = allResSE[order(allResSE[,min(grep("_Start", colnames(allResSE)))]),,drop=FALSE]
	allResSE = allResSE[order(allResSE[,grep("^Chr$", colnames(allResSE))]),,drop=FALSE]
	#evs = apply(allResSE, 1, function(x) paste(x[1], x[length(x)], sep="_", collapse="_"))	
	#allResSE = allResSE[!duplicated(evs),,drop=FALSE]
	saveRDS(allResSE, file=paste(pathAllResAnalysis, "allResSE.rds", sep=""))

	#	Stats
	#	Stats(idprojet, esp, nCPU, align)	

	formatResAll(pathAllResAnalysis)

#	RDS BDD, [chr, start, stop, ev, count, score, proj, pat]

	