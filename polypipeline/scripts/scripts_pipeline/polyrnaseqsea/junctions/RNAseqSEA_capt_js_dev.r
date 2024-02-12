
options(bitmapType='cairo')

args = commandArgs(trailingOnly=TRUE)
# saveRDS(args, paste(getwd(), "/ArgsRNAseqSEAcaptJS_gene.rds", sep=""))
#	args = readRDS(paste(getwd(), "/ArgsRNAseqSEAcaptJS.rds", sep=""));

tmpArgs = t(rbind(apply(as.matrix(unlist(args)), 1, function(x) unlist(strsplit(x, "=")))))
if(sum(grepl("biblioPath", tmpArgs[,1]))>0)	biblioPath = tmpArgs[grep("biblioPath", tmpArgs[,1], ignore.case=TRUE),2]
source(biblioPath)

if(sum(grepl("resPath", tmpArgs[,1]))>0)	resPath = tmpArgs[grep("resPath", tmpArgs[,1], ignore.case=TRUE),2]
basePath = dirname(dirname(resPath))
if(sum(grepl("gcrds", tmpArgs[,1]))>0)	gcrds = tmpArgs[grep("gcrds", tmpArgs[,1], ignore.case=TRUE),2]
if(sum(grepl("AnnotsPath", tmpArgs[,1]))>0)	AnnotsPath = tmpArgs[grep("AnnotsPath", tmpArgs[,1], ignore.case=TRUE),2]

if(sum(grepl("gene", tmpArgs[,1]))>0)	gene = tmpArgs[grep("gene", tmpArgs[,1], ignore.case=TRUE),2]
if(sum(grepl("esp", tmpArgs[,1]))>0)	esp = tmpArgs[grep("esp", tmpArgs[,1], ignore.case=FALSE),2]
if(sum(grepl("gcvers", tmpArgs[,1]))>0)	gcvers = tmpArgs[grep("gcvers", tmpArgs[,1], ignore.case=FALSE),2]
if(sum(grepl("bamsProj", tmpArgs[,1]))>0)
{
	bamsProj = tmpArgs[grep("bamsProj", tmpArgs[,1], ignore.case=TRUE),2]
	bamsProj = unlist(strsplit(bamsProj, ","))
}
if(sum(grepl("idprojet", tmpArgs[,1]))>0)	idprojet = tmpArgs[grep("idprojet", tmpArgs[,1], ignore.case=TRUE),2]

if(sum(grepl("titre", tmpArgs[,1]))>0)
{
	titre = tmpArgs[grep("titre", tmpArgs[,1], ignore.case=TRUE),2]
}else{
	titre = date
}

align="hisat2"
if(sum(grepl("^align$", tmpArgs[,1])))	align = tmpArgs[grep("align", tmpArgs[,1], ignore.case=TRUE),2]

nCPUmax = 1
if(sum(grepl("nCPUmax", tmpArgs[,1]))>0)	nCPUmax = as.numeric(tmpArgs[grep("nCPUmax", tmpArgs[,1], ignore.case=TRUE),2])

limDecal = 1	# 5	Marge acceptable d'approximation des bornes des jonctions a concatener
if(sum(grepl("limDecal", tmpArgs[,1]))>0)	limDecal = as.numeric(tmpArgs[grep("limDecal", tmpArgs[,1], ignore.case=TRUE),2])

if(sum(grepl("sambambaPath", tmpArgs[,1]))>0)	sambambaPath = tmpArgs[grep("sambambaPath", tmpArgs[,1], ignore.case=TRUE),2]
if(sum(grepl("samtoolsPath", tmpArgs[,1]))>0)	samtoolsPath = tmpArgs[grep("samtoolsPath", tmpArgs[,1], ignore.case=TRUE),2]
if(sum(grepl("picardPath", tmpArgs[,1]))>0)	picardPath = tmpArgs[grep("picardPath", tmpArgs[,1], ignore.case=TRUE),2]

#	remplace les arguments par le SampleType
pathtFileSamplesTypes = paste(resPath, "/RNAseqSEA/SamplesTypes", if(!is.null(titre)){paste("_", titre, sep="")}, ".txt", sep="")

pathtFileENSg = list.files(paste(resPath, "/RNAseqSEA/", sep=""), full.names=TRUE)
pathtFileENSg = pathtFileENSg[grepl("^ENSg", basename(pathtFileENSg))]

#nfos = c(idprojet, esp, align, Sys.info()[c("nodename", "login", "user", "effective_user")], getwd())
#slogRNAseqSEA_captjs(nfos)

if(file.exists(pathtFileENSg))
{
	#	mac2unix ENSg.txt; dos2unix ENSg.txt
	cmd = paste("mac2unix ", pathtFileENSg, "; dos2unix ", pathtFileENSg, sep="")
	system(paste(cmd, sep=""))
}

if(file.exists(pathtFileSamplesTypes))
{
	SamplesTypes = as.matrix(read.table(pathtFileSamplesTypes, sep="\t", quote="", header=TRUE))
	
	for(B in 1:nrow(SamplesTypes))
	{
		if(file.exists(paste(basePath, "/", SamplesTypes[B,"Proj"], "/", esp, "/align/", align, "/", SamplesTypes[B,"Sample"],".bam", sep="")))
		{
			write(paste("\tOK\t", SamplesTypes[B,"Sample"], "/", SamplesTypes[B,"Proj"], sep=""), file="")
		}else{
			write(paste("\tNot found !!!\t", SamplesTypes[B,"Sample"], "/", SamplesTypes[B,"Proj"], sep=""), file="")
		}
	}
	bamsProj = unique(SamplesTypes[,"Proj"])
	
	if(sum(!grepl("t", SamplesTypes[,"Type"]))==0)
	{
		Comps=TRUE	#	si ctrl ou pat renseigne pour tous les echs
	}else{
		write(paste("\t Pas de type => pas de comparaisons", sep=""), file="")
		Comps = FALSE
	}
}else{
	Comps = FALSE
}

projPath = paste(resPath, "/analysis/", sep="")
if(!file.exists(projPath))	dir.create(projPath)
baseSamplePath = paste(projPath, "RNAseqSEA", if(!is.null(titre)){paste("_", titre, sep="")}, "/", sep="")
if(!file.exists(baseSamplePath))	dir.create(baseSamplePath)
pathResGenes = paste(baseSamplePath, "resGenes/", sep="")
if(!file.exists(pathResGenes))	dir.create(pathResGenes)
pathRes = paste(pathResGenes, gene, "/", sep="")
if(!file.exists(pathRes))	dir.create(pathRes)

pathPosGene = paste(pathRes, "juncPairPos/", sep="")
if(!file.exists(pathPosGene))	dir.create(pathPosGene)
pathBamsGene = paste(pathRes, "bams/", sep="")
if(!file.exists(pathBamsGene))	dir.create(pathBamsGene)
rmdupPath = paste(pathRes, "rmdup/", sep="")
if(!file.exists(rmdupPath))	dir.create(rmdupPath)
bedpath = paste(pathRes, "mbed/", sep="")
if(!file.exists(bedpath))	dir.create(bedpath)
juncBamPath = paste(pathRes, "juncBam/", sep="")
if(!file.exists(juncBamPath))	dir.create(juncBamPath)

db_name = paste(gsub("_.*", "", esp), "_gencode", gcvers, sep="")
db_structure = readRDS(gcrds)
nomGene = unique(db_structure[db_structure[,"ensembl_gene_id"]%in%gene,"external_gene_name"])

fStructGene = db_structure[,"ensembl_gene_id"]%in%gene
if(sum(fStructGene)>0)
{
	exons = db_structure[db_structure[,"ensembl_gene_id"]%in%gene,]
	exons = exons[as.numeric(exons[,"level"])<=2,]	#	filtrage level Get level 1 & 2 annotation (manually annotated) only
}else{
	write(paste("#\n#\t", gene, " est absent de la ref ", db_name, "\n#", sep=""), file="")
}

reverse = (sum((exons[,"exon_chrom_end"]-exons[,"exon_chrom_start"])<0)>0)
#	pb si le gene est reverse, test si start < end de chaque exon
if(reverse)
{
	exons = exons[,c("chromosome_name", "ensembl_gene_id", "ensembl_transcript_id", "ensembl_exon_id", "exon_chrom_start", "exon_chrom_end", "external_gene_name", "level")]
	colnames(exons) = c("chromosome_name", "ensembl_gene_id", "ensembl_transcript_id", "ensembl_exon_id", "exon_chrom_end", "exon_chrom_start", "external_gene_name", "level")
}

if(nrow(exons)>0)
{
	Chr = exons[1,"chromosome_name"]
	Start = min(exons[,c("exon_chrom_start", "exon_chrom_end")])
	End = max(exons[,c("exon_chrom_start", "exon_chrom_end")])
	
	########################################################################################
#	Junctions possibles, junctions SPE de transcrit
	exons = exons[order(as.numeric(exons[,"exon_chrom_start"])),,drop=FALSE]
	exons = exons[order(exons[,"ensembl_transcript_id"]),,drop=FALSE]
	JuncT = list()	#	Liste des jonctions par transcrit
	for(Tid in unique(exons[,"ensembl_transcript_id"]))	#	Tid = unique(exons[,"ensembl_transcript_id"])[1]
	{
		Texons = exons[exons[,"ensembl_transcript_id"]==Tid,,drop=FALSE]
		Texons = Texons[order(as.numeric(Texons[,"exon_chrom_start"])),,drop=FALSE]
		Tjunc = cbind(Texons[1:(nrow(Texons)-1),"exon_chrom_end"], Texons[2:(nrow(Texons)),"exon_chrom_start"])
		
		Tjunc = Tjunc[!is.na(apply(Tjunc, 1, function(x) sum(x))),,drop=FALSE]
		
		JuncT[[length(JuncT)+1]] = apply(Tjunc, 1, function(x) paste(x[1], "_", x[2], sep=""))
		names(JuncT)[length(JuncT)] = Tid
	}
	
	AllJunc = unique(unlist(JuncT))	#	Toutes les jonctions possibles
	if(length(JuncT)>1)
	{
		TabAllJunc = t(apply(as.matrix(AllJunc), 1, function(x) as.numeric(grepl(x, JuncT))))
		colnames(TabAllJunc) = names(JuncT)
		colnames(TabAllJunc) = c(paste("Junc\t", colnames(TabAllJunc)[1], sep=""), colnames(TabAllJunc)[2:length(colnames(TabAllJunc))])
	}else{
		TabAllJunc = matrix(1, ncol=length(Tid), nrow=length(AllJunc))
		colnames(TabAllJunc) = names(JuncT)
		colnames(TabAllJunc) = c(paste("Junc\t", colnames(TabAllJunc)[1], sep=""))
	}
	rownames(TabAllJunc) = AllJunc
	write.table(TabAllJunc, file=paste(pathRes, "/TableAll_Junctions_Transcrits.txt", sep=""), sep="\t", quote=FALSE)
	#	write.table(TabAllJunc, file=paste(pathRes, "/TableAll_Junctions_All_Transcrits.txt", sep=""), sep="\t", quote=FALSE)
	
#	Jonctions specifique de transcrit (parmis ceux pris en compte)
	JuncCount = apply(TabAllJunc, 1, function(x) sum(x))
	speJ = names(JuncCount)[JuncCount==1]
	speJ = cbind(speJ, apply(as.matrix(speJ), 1, function(x) names(JuncT)[grepl(x, JuncT)]))
	colnames(speJ) = c("JunctionSPE", "ensembl_transcript_id")
	
#	Bed des junctions
	bamfiles = NULL
	for(B in bamsProj)
		bamfiles = c(bamfiles, list.files(paste("/data-isilon/sequencing/ngs/", B, "/", esp, "/align/", align, "/", sep=""), full.name=TRUE))
	bamfiles = bamfiles[!grepl(".bai", bamfiles)]	
	urls_bams = gsub("/data-isilon/sequencing/ngs/", "http://www.polyweb.fr/NGS/", bamfiles)
	write.table(urls_bams, file=paste(basePath, "/urls_bams.txt", sep=""), sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
	
	runif(1)
	
	if(file.exists("./Echs.txt"))	Echs = as.matrix(read.table("./Echs.txt", header=TRUE, sep="\t", quote=""))
	
	#	type de nom chr
	CHRbam = ""
	chrStats = system(paste("samtools idxstats ", bamfiles[1], " | cut -f 1 | head -n 10", sep=""), intern=TRUE)
	if(sum(grepl("chr", chrStats))>1) CHRbam = "chr"
	
	nCPU = nCPUmax
	library(parallelMap)
	parallelStart(mode = "multicore", cpus=nCPU, show.info=TRUE) 
	f = function(B) geneBams(B)
	y = parallelMap(f, c(1:length(bamfiles)))
	parallelStop()
	
	bamsGene = list.files(pathBamsGene, full.name=TRUE)
	bamsGene = bamsGene[!grepl(".bai", bamsGene)]
	
	CountsMat = matrix(0, ncol=4, length(bamsGene))
	colnames(CountsMat) = c("bam", "TotReadsGene", "rmdup", "JuncReadsGene")
	CountsMat[,"bam"] = gsub(".bam", "", basename(bamsGene))
	for(cB in 1:length(bamsGene))
		CountsMat[cB, "TotReadsGene"] = system(paste("samtools view ", bamsGene[cB], " | wc -l", sep=""), intern=TRUE)
	
	library(parallelMap)
	parallelStart(mode = "multicore", cpus=nCPUmax, show.info=TRUE) 
	f = function(I) rmdup(I)
	y = parallelMap(f, c(1:length(bamsGene)))
	parallelStop()  
	#	unlink(pathBamsGene, recursive=TRUE)
	
	bamsRMdups = list.files(rmdupPath, full.name=TRUE, pattern="\\.bam$")
	nCPU = nCPUmax
	library(parallelMap)
	parallelStart(mode = "multicore", cpus=nCPU, show.info=TRUE) 
	f = function(J) exJunBed(J)
	y = parallelMap(f, c(1:length(bamsRMdups)))
	parallelStop()
	
	bamsRMdups = list.files(rmdupPath, full.name=TRUE, pattern="\\.bam$")
	nCPU = nCPUmax
	library(parallelMap)
	parallelStart(mode = "multicore", cpus=nCPU, show.info=TRUE) 
	f = function(J) JuncBam(J)
	y = parallelMap(f, c(1:length(bamsRMdups)))
	parallelStop()
		
	for(cB in 1:length(bamsRMdups))
		CountsMat[cB, "rmdup"] = system(paste("samtools view ", bamsRMdups[cB], " | wc -l", sep=""), intern=TRUE)
	
#	Bed des junctions
	listBed = list.files(bedpath, full.names=TRUE)
	beds = listBed[grepl(Chr, basename(listBed))]
	
	nCPU = nCPUmax
	library(parallelMap)
	parallelStart(mode = "multicore", cpus=nCPU, show.info=TRUE) 
	f = function(S) SelectJun(S)
	y = parallelMap(f, beds)
	parallelStop() 
	
	samples = list.files(pathPosGene, full.name=TRUE, pattern="\\.rds$")
	if(length(samples)>0)	for(cB in 1:length(samples))
		{
			tmp = readRDS(samples[cB])
			NcB = grep(gsub(".rds", "", gsub(paste("juncPairPos_", gene, "_", Chr, "_", sep=""), "", basename(samples[cB]))), CountsMat[,"bam"])
			CountsMat[NcB, "JuncReadsGene"] = nrow(tmp)
		}
	#	unlink(bedpath, recursive=TRUE)
	
	write.table(CountsMat, file=paste(pathRes, "/Counts_", gene, ".txt", sep=""), sep="\t", quote=FALSE, row.names=FALSE)
	
	#######################################################################################
#	 graphique
	
	if(length(samples)>0)
	{
		pathEv = paste(pathRes, "SpliceRes/", sep="")
		if(!file.exists(pathEv))	dir.create(pathEv)
		pathFigs = paste(pathEv, "/Figs/", sep="")
		if(!file.exists(pathFigs))	dir.create(pathFigs)
		pathResSample = paste(pathEv, "/ResSample/", sep="")
		if(!file.exists(pathResSample))	dir.create(pathResSample)
		pathResDataJunc = paste(pathEv, "/DataJunc/", sep="")
		if(!file.exists(pathResDataJunc))	dir.create(pathResDataJunc)
		pathResNFOJunc = paste(pathEv, "/nfoJunc/", sep="")
		if(!file.exists(pathResNFOJunc))	dir.create(pathResNFOJunc)

#	Concat de toutes les junctions des Samples, meme les anormales associe un nom unique
		AllgeneJun = matrix(0, ncol=2, nrow=0)
		for(S in 1:length(samples))
		{
			geneJun = readRDS(samples[S])
			AllgeneJun = rbind(AllgeneJun, geneJun[,1:2])
			
			for(C in nrow(AllgeneJun):1)
			{
				minStartDec = min(abs(AllgeneJun[C,"StartJun"]-as.numeric(exons[,"exon_chrom_end"])))
				minEndDec = min(abs(AllgeneJun[C,"EndJun"]-as.numeric(exons[,"exon_chrom_start"])))
				if((minStartDec>0)&(minStartDec<=limDecal))	AllgeneJun[C,"StartJun"] = min(as.numeric(exons[abs(AllgeneJun[C,"StartJun"]-as.numeric(exons[,"exon_chrom_end"]))<=limDecal,"exon_chrom_end"]))
				if((minEndDec>0)&(minEndDec<=limDecal))	AllgeneJun[C,"EndJun"] = min(as.numeric(exons[abs(AllgeneJun[C,"EndJun"]-as.numeric(exons[,"exon_chrom_start"]))<=limDecal,"exon_chrom_start"]))
			}
		}
		AllgeneJun = cbind(AllgeneJun, apply(AllgeneJun, 1, function(x) paste(x[1], "_", x[2], sep="")))	#	Ajoute le Start-End
		AllgeneJun = AllgeneJun[!duplicated(AllgeneJun),,drop=FALSE]	#	Supprime les duplicats
		AllgeneJun = AllgeneJun[order(as.numeric(AllgeneJun[,1])),,drop=FALSE]	#	Classe par Start Jun
		rownames(AllgeneJun) = AllgeneJun[,3]
		AllgeneJun = cbind(AllgeneJun,paste("J", c(1:nrow(AllgeneJun)), sep=""))
		colnames(AllgeneJun) = c("StartJun", "EndJun", "Coords", "JunctionUniqueName")
		
		nomsJ = NULL
		for(E in 1:length(AllJunc))	if(sum(grepl(AllJunc[E], rownames(AllgeneJun)))>0){nomsJ = c(nomsJ, AllgeneJun[rownames(AllgeneJun)==AllJunc[E], "JunctionUniqueName"])}else{nomsJ = c(nomsJ, "")}
		AllJunc = cbind(nomsJ, AllJunc)
		write.table(AllJunc, file=paste(pathEv, "/JonctionsNormales.txt", sep=""), sep="\t", quote=FALSE, col.names=c(paste("JunID\t", colnames(AllJunc)[1], sep=""), colnames(AllJunc)[2:length(colnames(AllJunc))]))
		
		TabNFO = matrix(0, ncol=6, nrow=0)	#	Taille des annots
		colnames(TabNFO) = c("ENSID", "nJ", "Chr", "Start", "End", "Interpretation")
		
		#	table de comptage de toutes les junctions
		countAllJunc = cbind(rep(gene, nrow(AllgeneJun)), rep(nomGene, nrow(AllgeneJun)), rep(Chr, nrow(AllgeneJun)), AllgeneJun[,c("StartJun", "EndJun", "JunctionUniqueName"),drop=FALSE])
		rownames(countAllJunc) = countAllJunc[,"JunctionUniqueName"]
		colnames(countAllJunc)[1:3] = c("ENSID", "Gene", "Chr")  
		countAllJunc = cbind(countAllJunc, matrix("", ncol=1, nrow=nrow(countAllJunc)))
		colnames(countAllJunc)[ncol(countAllJunc)] = "Type"
		countAllJunc[AllJunc[grep("J", AllJunc[,"nomsJ"]),"nomsJ"],"Type"] = "Canonique"
	
		for(S in 1:length(samples))
		{	
			nom = gsub(".rds", "", basename(samples[S]))
			nom = gsub("juncPairPos_", "", nom)
			write(paste("\n\t#\t", S, "/", length(samples), " ", nom, sep=""), file="")
			
			#	nfo file
			nfoname = paste(pathResNFOJunc, "/", nom, "_nfo_junc.txt", sep="")
			write(paste("nJ\tStart\tEnd\t", nom, "_CountJ\tInterpretation", sep=""), file=nfoname, append=FALSE)	#	Titres
			
			geneJun = readRDS(samples[S])	
			geneJun = cbind(geneJun, apply(geneJun, 1, function(x) x[2]-x[1]))
			colnames(geneJun)[ncol(geneJun)] = "JuncLength" 
			write.table(geneJun, file=paste(pathResDataJunc, "/", nom, "_geneJun.txt", sep=""), sep="\t", quote=FALSE, col.names=c(paste("JunID\t", colnames(geneJun)[1], sep=""), colnames(geneJun)[2:length(colnames(geneJun))]))
			
			#	=> correction des pos a 15 bases pres
			geneJunCorr = geneJun
			geneJunCorr = geneJunCorr[order(geneJunCorr[,"EndJun"]),,drop=FALSE]
			geneJunCorr = geneJunCorr[order(geneJunCorr[,"StartJun"]),,drop=FALSE]
			for(C in nrow(geneJunCorr):1)
			{
				minStartDec = min(abs(geneJunCorr[C,"StartJun"]-as.numeric(exons[,"exon_chrom_end"])))
				minEndDec = min(abs(geneJunCorr[C,"EndJun"]-as.numeric(exons[,"exon_chrom_start"])))
				if((minStartDec>0)&(minStartDec<=limDecal))	geneJunCorr[C,"StartJun"] = min(as.numeric(exons[abs(geneJunCorr[C,"StartJun"]-as.numeric(exons[,"exon_chrom_end"]))<=limDecal,"exon_chrom_end"]))
				if((minEndDec>0)&(minEndDec<=limDecal))	geneJunCorr[C,"EndJun"] = min(as.numeric(exons[abs(geneJunCorr[C,"EndJun"]-as.numeric(exons[,"exon_chrom_start"]))<=limDecal,"exon_chrom_start"]))
			}
			
			if(nrow(geneJunCorr)>=2)
			{
				geneJunCorr = geneJunCorr[order(geneJunCorr[,"EndJun"]),,drop=FALSE]
				geneJunCorr = geneJunCorr[order(geneJunCorr[,"StartJun"]),,drop=FALSE]
				geneJunNew = matrix(0, ncol=ncol(geneJun), nrow=0)
				geneJunNew = rbind(geneJunNew, geneJunCorr[1,,drop=FALSE])
				for(G in 2:nrow(geneJunCorr))
				{
					if((geneJunNew[nrow(geneJunNew),"StartJun"]==geneJunCorr[G,"StartJun"])&(geneJunNew[nrow(geneJunNew),"EndJun"]==geneJunCorr[G,"EndJun"]))
					{
						geneJunNew[nrow(geneJunNew),"CountJun"] = geneJunNew[nrow(geneJunNew),"CountJun"] + geneJunCorr[G,"CountJun"]
					}else{
						geneJunNew = rbind(geneJunNew, geneJunCorr[G,,drop=FALSE])
					}
				}
				geneJun = geneJunNew
			}
			
			#	filtre si les start ET end sont en dehors du gene... pb des marges
			geneJun = geneJun[!(((geneJun[,"StartJun"]>=End)&(geneJun[,"EndJun"]>=Start))|((geneJun[,"StartJun"]<=End)&(geneJun[,"EndJun"]<=Start))),,drop=FALSE]

			######################################################################################
			rownames(geneJun) = AllgeneJun[apply(geneJun, 1, function(x) paste(x[1], "_", x[2], sep="")),ncol(AllgeneJun)]	# Nomme chaque junction
			write.table(geneJun, file=paste(pathResDataJunc, "/", nom, "_geneJunNew.txt", sep=""), sep="\t", quote=FALSE, col.names=c(paste("JunID\t", colnames(geneJun)[1], sep=""), colnames(geneJun)[2:length(colnames(geneJun))]))
			
			#	maj table counts all
			countAllJunc = cbind(countAllJunc, matrix(0, ncol=1, nrow=nrow(countAllJunc)))
			countAllJunc[as.character(rownames(geneJun)),ncol(countAllJunc)] = as.numeric(geneJun[,"CountJun"])
			colnames(countAllJunc)[ncol(countAllJunc)] = gsub(paste(gene, "_", Chr, "_", sep=""), "", nom)
			
			############################################################################################
			#	Tableau Jonctions SPE
			tmpJ = apply(geneJun, 1, function(x) paste(x[1], "_", x[2], sep=""))
			speJ = cbind(speJ, matrix(0, ncol=1, nrow=nrow(speJ)))
			colnames(speJ)[ncol(speJ)] = nom
			speJ[,nom] = apply(speJ, 1, function(x) if(sum(tmpJ==x[1])>0){geneJun[tmpJ==x[1],"CountJun"]}else{return(0)})
			
			############################################################################################
			#	Selection des junctions
#			JunH = 4	#hauteur de decalage entre chaque jonction
#			Tr = 1
#			Trs = names(table(exons[,"ensembl_transcript_id"]))	#	Noms transcrits alt
#			#bitmap(file = paste(pathFigs, "/", gsub(".rds", "", basename(samples[S])), ".png", sep=""), height=50, width=80, pointsize=30, bg="white", res=200)	#	TSEN54
#			
#			# SVG file ?
#			
#			svg(file = paste(pathFigs, "/", gsub(".rds", "", basename(samples[S])), ".svg", sep=""), height=30, width=70, pointsize=30, bg="white")
#						
#			op <- par(bg = "white")
#			plot(x=c(min(as.numeric(exons[,"exon_chrom_start"]))-2000, max(as.numeric(exons[,"exon_chrom_end"]))), y=c(-10, length(Trs)*10+nrow(geneJun)*JunH+100), type = "n", xlab = "", ylab = "", main = gene, axes = FALSE)
#			
#			#	Jonctions
#			if(nrow(geneJun)>0)	for(J in 1:nrow(geneJun))
#			{
#				segments(x0=geneJun[J,"StartJun"], y0=Tr*10+4+J*JunH, x1=geneJun[J,"EndJun"], col="green")
#				segments(x0=geneJun[J,"StartJun"], y0=Tr*10+4+J*JunH, y1=-5, col="green", lty="dotted")
#				text(x=(geneJun[J,"EndJun"]-geneJun[J,"StartJun"])/2+geneJun[J,"StartJun"], y=(Tr*10+J*JunH)+3.3, label=geneJun[J,"CountJun"], cex=0.5)
#				segments(x0=geneJun[J,"EndJun"], y0=Tr*10+4+J*JunH, y1=-5, col="green", lty="dotted")
#				text(x=(geneJun[J,"EndJun"]-geneJun[J,"StartJun"])/2+geneJun[J,"StartJun"], y=(Tr*10+J*JunH)+4.7, label=rownames(geneJun)[J], cex=0.5, col="green")
#			}
#			Trs = Trs[Trs%in%unique(exons[,"ensembl_transcript_id"])]
#			
#			#	Exons
#			for(M in 1:length(Trs))	#	Tr = 1
#			{
#				tmpTr = exons[exons[,"ensembl_transcript_id"]==Trs[M], c("ensembl_exon_id", "exon_chrom_start", "exon_chrom_end")]
#				segments(x0=min(as.numeric(tmpTr[,"exon_chrom_start"])), y0=M*10, x1=max(as.numeric(tmpTr[,"exon_chrom_end"])), col="grey")
#				
#				for(E in 1:nrow(tmpTr))
#					rect(as.numeric(tmpTr[E,"exon_chrom_start"]), M*10-4, as.numeric(tmpTr[E,"exon_chrom_end"]), M*10+4, col = "black")
#				
#				text(x=min(as.numeric(tmpTr[,"exon_chrom_start"]))-3500, y=M*10, labels=Trs[M], col="black", cex=0.5)
#			}
			
			##############################################################################
			#	Traitement des evenements de splicing
			##############################################################################
			
			cols = c("Junc_RI", "ENSID", "Gene", "Chr", "Junc_RI_Start", "Junc_RI_End", "Junc_RI_Count", "Junc_Normale", "Junc_Normale_Start", "Junc_Normale_End", "Junc_Normale_Count", "Type")
			matEV_RI = matrix(0, nrow=0, ncol=length(cols))
			colnames(matEV_RI) = cols
			
			#	Si start et end de la junctions sont connus
			test_RI = ((!(geneJun[,"StartJun"]%in%exons[,"exon_chrom_end"])|!(geneJun[,"EndJun"]%in%exons[,"exon_chrom_start"]))&(!(geneJun[,"StartJun"]%in%exons[,"exon_chrom_start"])|!(geneJun[,"EndJun"]%in%exons[,"exon_chrom_end"])))
			
			tmpExons = exons[order(as.numeric(exons[,"exon_chrom_start"])),]
			tmpExons = tmpExons[!duplicated(tmpExons[,"ensembl_exon_id"]),]
			Map = rbind(cbind(paste("Start_E", c(1:nrow(tmpExons)), sep="_"), tmpExons[,"exon_chrom_start"], tmpExons[,"ensembl_exon_id"]),
					cbind(paste("End_E", c(1:nrow(tmpExons)), sep="_"), tmpExons[,"exon_chrom_end"], tmpExons[,"ensembl_exon_id"]))
			colnames(Map) = c("N", "Pos", "EnsID_JuncType")
			
			if(sum(test_RI)>0)
			{
				J_anorm = geneJun[test_RI, c("StartJun", "EndJun", "CountJun"),drop=FALSE]	#	Junctions anormales
				resNFOjun = matrix(0, ncol=28, nrow=0)
				colnames(resNFOjun) = c("ENSID", "Junction", "Chr", colnames(J_anorm), "StartConnu", "EndConnu", "Start_in_exon", "End_in_exon", "Start_in_intron", "End_in_intron", "Start_in_3p",
						"End_in_3p", "Start_in_5p", "End_in_5p", "Start_in_exon1", "End_in_exon1", "Start_in_DerExon", "End_in_DerExon", "End_EQ_StartExon1", "Etart_EQ_EndDerExon",
						"EtartEnd_Meme_exon", "StartEnd_Meme_intron", "IE_contigus", "EI_contigus", "FinJrecouvreExon_SE", "DebutJrecouvreExon_SE")
				
				TabNFO = cbind(TabNFO, rep(0, nrow(TabNFO)))	#	Ajoute une colonne / ech
				colnames(TabNFO)[ncol(TabNFO)] = nom
				
				for(J in 1:nrow(J_anorm))
				{
					nJ = rownames(J_anorm)[J]
					rangJ = grep(paste(nJ, "$", sep=""), rownames(geneJun))
					
					#	Borne connue / inconnue
					test_StartConnu = J_anorm[J,"StartJun"]%in%exons[,"exon_chrom_end"]
					test_EndConnu = J_anorm[J,"EndJun"]%in%exons[,"exon_chrom_start"]
					
					test_EndJ_EndX = test_StartJ_StartX = FALSE
					if(!test_EndConnu)	test_EndJ_EndX = J_anorm[J,"EndJun"]%in%exons[,"exon_chrom_end"]	#	La junction se termine avec un exon sauté => SE
					if(!test_StartConnu)	test_StartJ_StartX = J_anorm[J,"StartJun"]%in%exons[,"exon_chrom_start"]	#	La junction debute avec un exon connu
					
					#	Jonction au moins partiellement dans un exon
					Jstart_in_exon = (J_anorm[J, "StartJun"]>exons[,"exon_chrom_start"])&(J_anorm[J, "StartJun"]<exons[,"exon_chrom_end"])
					test_Jstart_in_exon = sum(Jstart_in_exon)>0
					Jend_in_exon = (J_anorm[J, "EndJun"]>exons[,"exon_chrom_start"])&(J_anorm[J, "EndJun"]<exons[,"exon_chrom_end"])
					test_Jend_in_exon = sum(Jend_in_exon)>0
					
					#	Jonction au moins partiellement dans un intron
					test_Jstart_in_intron = !(test_Jstart_in_exon|test_StartConnu)&(J_anorm[J, "StartJun"]>min(exons[,"exon_chrom_start"]))&(J_anorm[J, "StartJun"]<max(exons[,"exon_chrom_end"]))
					test_Jend_in_intron = !(test_Jend_in_exon|test_EndConnu)&(J_anorm[J, "EndJun"]>min(exons[,"exon_chrom_start"]))&(J_anorm[J, "EndJun"]<max(exons[,"exon_chrom_end"]))
					
					#	Hors gene
					test_Jstart_in_3p = (J_anorm[J, "StartJun"]>max(exons[,"exon_chrom_end"]))
					test_Jend_in_3p = (J_anorm[J, "EndJun"]>max(exons[,"exon_chrom_end"]))
					
					test_Jend_in_5p = (J_anorm[J, "EndJun"]<min(exons[,"exon_chrom_start"]))
					test_Jstart_in_5p = (J_anorm[J, "StartJun"]<min(exons[,"exon_chrom_start"]))
					
					#	1 er exon
					test_Jstart_in_exon1 = (J_anorm[J, "StartJun"]>min(exons[,"exon_chrom_start"]))&(J_anorm[J, "StartJun"]<min(exons[,"exon_chrom_end"]))
					test_Jend_in_exon1 = (J_anorm[J, "EndJun"]>min(exons[,"exon_chrom_start"]))&(J_anorm[J, "EndJun"]<min(exons[,"exon_chrom_end"]))
					
					#	Der exon
					test_Jstart_in_Derexon = (J_anorm[J, "StartJun"]>max(exons[,"exon_chrom_start"]))&(J_anorm[J, "StartJun"]<max(exons[,"exon_chrom_end"]))
					test_Jend_in_Derexon = (J_anorm[J, "EndJun"]>max(exons[,"exon_chrom_start"]))&(J_anorm[J, "EndJun"]<max(exons[,"exon_chrom_end"]))
					
					test_Jend_EQ_Startexon1 = J_anorm[J, "EndJun"]==min(exons[,"exon_chrom_start"])	#	Fin Junc en debut exon1
					test_Jstart_EQ_EndDerexon = J_anorm[J, "StartJun"]== max(exons[,"exon_chrom_end"]) #	Debut Junc en fin dernier exon
					
					#	coords exon
					if(test_Jstart_in_exon)	# coords exon inclurant Jstart
					{
						Deb_Ex_Incl_Jstart = max(exons[exons[,"exon_chrom_start"]<J_anorm[J, "StartJun"], "exon_chrom_start"])
						Fin_Ex_Incl_Jstart = max(exons[exons[,"exon_chrom_end"]>J_anorm[J, "StartJun"], "exon_chrom_end"])
					}
					if(test_Jend_in_exon)	# coords exon inclurant Jend
					{
						Deb_Ex_Incl_Jend = max(exons[exons[,"exon_chrom_start"]<J_anorm[J, "EndJun"], "exon_chrom_start"])
						Fin_Ex_Incl_Jend = max(exons[exons[,"exon_chrom_end"]>J_anorm[J, "EndJun"], "exon_chrom_end"])
					}
					#	Start et End dans le meme exon
					test_JstartJend_Meme_exon = FALSE
					if(test_Jstart_in_exon&test_Jend_in_exon)	test_JstartJend_Meme_exon = (Deb_Ex_Incl_Jstart==Deb_Ex_Incl_Jend)&(Fin_Ex_Incl_Jstart==Deb_Ex_Incl_Jend)
					
					#	coords intron
					if(test_Jstart_in_intron)	# coords intron inclurant Jstart
					{
						Deb_Int_Incl_Jstart = max(exons[exons[,"exon_chrom_end"]<J_anorm[J, "StartJun"], "exon_chrom_end"])
						Fin_Int_Incl_Jstart = min(exons[exons[,"exon_chrom_start"]>J_anorm[J, "StartJun"], "exon_chrom_start"])
					}
					if(test_Jend_in_intron)	# coords intron inclurant Jend
					{
						Deb_Int_Incl_Jend = max(exons[exons[,"exon_chrom_end"]<J_anorm[J, "EndJun"], "exon_chrom_end"])
						Fin_Int_Incl_Jend = min(exons[exons[,"exon_chrom_start"]>J_anorm[J, "EndJun"], "exon_chrom_start"])
					}
					
					#	Start et End dans le meme intron
					test_JstartJend_Meme_intron = FALSE
					if(test_Jstart_in_intron&test_Jend_in_intron)	test_JstartJend_Meme_intron = (Deb_Int_Incl_Jstart==Deb_Int_Incl_Jend)&(Fin_Int_Incl_Jstart==Fin_Int_Incl_Jend)
					
					#	exon intron contigus
					IE_contig = FALSE
					if(test_Jstart_in_intron&test_Jend_in_exon)	IE_contig = Fin_Int_Incl_Jstart==Deb_Ex_Incl_Jend #	Start in intron ET End in exon => Fin intron == debut exon		
					
					EI_contig = FALSE
					if(test_Jstart_in_exon&test_Jend_in_intron)	EI_contig = Fin_Ex_Incl_Jstart==Deb_Int_Incl_Jend #	Start in exon ET end in intron => Fin exon == debut intron
					
					#	Tableau res de la junction
					resNFOjun = rbind(resNFOjun, c(gene, rownames(J_anorm)[J], Chr, J_anorm[J,], sum(test_StartConnu), sum(test_EndConnu), sum(test_Jstart_in_exon), sum(test_Jend_in_exon),
									sum(test_Jstart_in_intron), sum(test_Jend_in_intron), sum(test_Jstart_in_3p), sum(test_Jend_in_3p), sum(test_Jstart_in_5p), sum(test_Jend_in_5p),
									sum(test_Jstart_in_exon1), sum(test_Jend_in_exon1), sum(test_Jstart_in_Derexon), sum(test_Jend_in_Derexon), sum(test_Jend_EQ_Startexon1),
									sum(test_Jstart_EQ_EndDerexon), sum(test_JstartJend_Meme_exon), sum(test_JstartJend_Meme_intron), sum(IE_contig), sum(EI_contig), sum(test_EndJ_EndX), sum(test_StartJ_StartX)))
					
					#	Recherhe des jonctions incluant OU INCLUES dans la jonction anormale
					#	Si la jonc Anorm est plus courte et est inclue dans une autre
					canonJunc = AllgeneJun[rownames(AllgeneJun)%in%rownames(TabAllJunc),,drop=FALSE]
					canonGeneJun = geneJun[rownames(geneJun)%in%canonJunc[,"JunctionUniqueName"],,drop=FALSE]
					
					test_J_courte = (J_anorm[J, "StartJun"]==canonGeneJun[,"StartJun"])&(J_anorm[J, "EndJun"]<canonGeneJun[,"EndJun"])|
							(J_anorm[J, "StartJun"]>canonGeneJun[,"StartJun"])&(J_anorm[J, "EndJun"]==canonGeneJun[,"EndJun"])
					
					#	Si la jonction Anorm est plus logue et inclue une autre
					test_J_longue = (J_anorm[J, "StartJun"]==canonGeneJun[,"StartJun"])&(J_anorm[J, "EndJun"]>canonGeneJun[,"EndJun"])|
							(J_anorm[J, "StartJun"]<canonGeneJun[,"StartJun"])&(J_anorm[J, "EndJun"]==canonGeneJun[,"EndJun"])
					
					#	Selectionne la J normale avec le max de comp
					test_J_norm = test_J_longue|test_J_courte
					if(sum(test_J_norm)==0) # pas de J_norm ?
					{
						J_norm = matrix(c("---", "---", 0), ncol=3, nrow=1)
						rownames(J_norm) = "---"
					}
					if(sum(test_J_norm)>0)	J_norm = canonGeneJun[test_J_norm,c("StartJun", "EndJun", "CountJun"),drop=FALSE]
					if(sum(test_J_norm)>1)	J_norm = J_norm[J_norm[,"CountJun"]==max(as.numeric(J_norm[,"CountJun"])),,drop=FALSE]
					if(nrow(J_norm)>1)	J_norm = J_norm[1,,drop=FALSE]	#	Si plusieurs junc canoniques a comptages egaux...
					
					newRI = c(nJ, gene, nomGene, Chr, J_anorm[J,,drop=FALSE], rownames(J_norm), J_norm)					
					newRI = matrix(newRI, ncol=length(newRI), byrow=TRUE)
					
					#	Type evenement
					#	correction des priorites
					if(test_Jstart_in_exon & test_StartConnu) test_Jstart_in_exon = FALSE
					if(test_Jstart_in_intron & test_StartConnu) test_Jstart_in_intron = FALSE
					if(test_Jend_in_exon & test_EndConnu) test_Jend_in_exon = FALSE
					if(test_Jend_in_intron & test_EndConnu) test_Jend_in_intron = FALSE
					
					##################################################################################
					if(test_Jstart_in_3p&test_Jend_in_3p)	#	Start in 3p ET end in 3p
					{
						colRI = "red"
						expl = "Start in 3p ET end in 3p => RI en 3p"
						nomJ = paste(nJ, "=> RI_3p", sep="")
						line = paste(nomJ, "\t", gene, "\t", nomGene, "\t", "\t", geneJun[J,"StartJun"], "\t", geneJun[J,"EndJun"], "\t", geneJun[nJ,"CountJun"], "\t", expl, sep="")
						write(line, file=nfoname, append=TRUE)
						write(line, file="")
						matEV_RI = rbind(matEV_RI, cbind(newRI, as.matrix(rep("SE", nrow(newRI)))))
						
						if(sum(TabNFO[,"nJ"]%in%nomJ)==0)	TabNFO = rbind(TabNFO, c(gene, nomJ, Chr, geneJun[J,"StartJun"], geneJun[J,"EndJun"], expl, rep(0, ncol(TabNFO)-6)))
						TabNFO[TabNFO[,"nJ"]%in%nomJ,nom] = geneJun[nJ,"CountJun"]
					}
					
					if(test_Jstart_in_exon&test_Jend_in_3p)	#	Start in exon ET end in 3p
					{
						colRI = "red"
						expl = "Start in exon ET end in 3p => raccourci exon en 3 et RI"
						nomJ = paste(nJ, "=> A3SS+RI", sep="")
						line = paste(nomJ, "\t", gene, "\t", nomGene, "\t", "\t", geneJun[J,"StartJun"], "\t", geneJun[J,"EndJun"], "\t", geneJun[nJ,"CountJun"], "\t", expl, sep="")
						write(line, file=nfoname, append=TRUE)
						write(line, file="")
						matEV_RI = rbind(matEV_RI, cbind(newRI, as.matrix(rep("SE", nrow(newRI)))))
						
						if(sum(TabNFO[,"nJ"]%in%nomJ)==0)	TabNFO = rbind(TabNFO, c(gene, nomJ, Chr, geneJun[J,"StartJun"], geneJun[J,"EndJun"], expl, rep(0, ncol(TabNFO)-6)))
						TabNFO[TabNFO[,"nJ"]%in%nomJ,nom] = geneJun[nJ,"CountJun"]
					}
					
					if(test_StartJ_StartX&test_Jend_in_exon1)	#	Start avec exon 1 ET end in exon 1
					{
						colRI = "blue"
						expl = "Start avec exon 1 et fin dans exon 1 => raccourci exon 1 en 5"
						nomJ = paste(nJ, "=> SE", sep="")
						line = paste(nomJ, "\t", gene, "\t", nomGene, "\t", "\t", geneJun[J,"StartJun"], "\t", geneJun[J,"EndJun"], "\t", geneJun[nJ,"CountJun"], "\t", expl, sep="")
						write(line, file=nfoname, append=TRUE)
						write(line, file="")
						matEV_RI = rbind(matEV_RI, cbind(newRI, as.matrix(rep("SE", nrow(newRI)))))
						
						if(sum(TabNFO[,"nJ"]%in%nomJ)==0)	TabNFO = rbind(TabNFO, c(gene, nomJ, Chr, geneJun[J,"StartJun"], geneJun[J,"EndJun"], expl, rep(0, ncol(TabNFO)-6)))
						TabNFO[TabNFO[,"nJ"]%in%nomJ,nom] = geneJun[nJ,"CountJun"]
					}
					
					if(test_StartJ_StartX&test_EndConnu)	#	(Start exon connu) ET (End connu)
					{
						colRI = "blue"
						expl = "Start avec exon et fin avec intron"
						nomJ = paste(nJ, "=> SE", sep="")
						line = paste(nomJ, "\t", gene, "\t", nomGene, "\t", "\t", geneJun[J,"StartJun"], "\t", geneJun[J,"EndJun"], "\t", geneJun[nJ,"CountJun"], "\t", expl, sep="")
						write(line, file=nfoname, append=TRUE)
						write(line, file="")
						matEV_RI = rbind(matEV_RI, cbind(newRI, as.matrix(rep("SE", nrow(newRI)))))
						
						if(sum(TabNFO[,"nJ"]%in%nomJ)==0)	TabNFO = rbind(TabNFO, c(gene, nomJ, Chr, geneJun[J,"StartJun"], geneJun[J,"EndJun"], expl, rep(0, ncol(TabNFO)-6)))
						TabNFO[TabNFO[,"nJ"]%in%nomJ,nom] = geneJun[nJ,"CountJun"]
					}
					
					#	A3SS+SE+RI_3p
					if(test_Jstart_in_exon1&test_Jend_in_3p)	#	(Start in premier exon) ET (End in 3p)
					{
						colRI = "red"
						expl = "Start in premier exon et End in 3p"
						nomJ = paste(nJ, "=>A3SS+SE+RI_3p", sep="")
						line = paste(nomJ, "\t", gene, "\t", nomGene, "\t", "\t", geneJun[J,"StartJun"], "\t", geneJun[J,"EndJun"], "\t", geneJun[nJ,"CountJun"], "\t", expl, sep="")
						write(line, file=nfoname, append=TRUE)
						write(line, file="")
						matEV_RI = rbind(matEV_RI, cbind(newRI, as.matrix(rep("A3SS+SE+RI_3p", nrow(newRI)))))
						
						if(sum(TabNFO[,"nJ"]%in%nomJ)==0)	TabNFO = rbind(TabNFO, c(gene, nomJ, Chr, geneJun[J,"StartJun"], geneJun[J,"EndJun"], expl, rep(0, ncol(TabNFO)-6)))
						TabNFO[TabNFO[,"nJ"]%in%nomJ,nom] = geneJun[nJ,"CountJun"]
					}
					
					#	SE+RI_3p
					if(test_StartJ_StartX&test_Jend_in_3p)	#	(Start exon connu) ET (End in 3p)
					{
						colRI = "red"
						expl = "Start exon connu et End in 3p"
						nomJ = paste(nJ, "=>SE+RI_3p", sep="")
						line = paste(nomJ, "\t", gene, "\t", nomGene, "\t", "\t", geneJun[J,"StartJun"], "\t", geneJun[J,"EndJun"], "\t", geneJun[nJ,"CountJun"], "\t", expl, sep="")
						write(line, file=nfoname, append=TRUE)
						write(line, file="")
						matEV_RI = rbind(matEV_RI, cbind(newRI, as.matrix(rep("RI_5p+SE", nrow(newRI)))))
						
						if(sum(TabNFO[,"nJ"]%in%nomJ)==0)	TabNFO = rbind(TabNFO, c(gene, nomJ, Chr, geneJun[J,"StartJun"], geneJun[J,"EndJun"], expl, rep(0, ncol(TabNFO)-6)))
						TabNFO[TabNFO[,"nJ"]%in%nomJ,nom] = geneJun[nJ,"CountJun"]
					}
					
					#	RI_5p+RI_3p
					if(test_Jstart_in_5p&test_Jend_in_3p)	#	(Start in 5p) ET (End in 5p)
					{
						colRI = "red"
						expl = "Start en 5p et End in 3p"
						nomJ = paste(nJ, "=>RI_5p+RI_3p", sep="")
						line = paste(nomJ, "\t", gene, "\t", nomGene, "\t", "\t", geneJun[J,"StartJun"], "\t", geneJun[J,"EndJun"], "\t", geneJun[nJ,"CountJun"], "\t", expl, sep="")
						write(line, file=nfoname, append=TRUE)
						write(line, file="")
						matEV_RI = rbind(matEV_RI, cbind(newRI, as.matrix(rep("RI_5p+SE", nrow(newRI)))))
						
						if(sum(TabNFO[,"nJ"]%in%nomJ)==0)	TabNFO = rbind(TabNFO, c(gene, nomJ, Chr, geneJun[J,"StartJun"], geneJun[J,"EndJun"], expl, rep(0, ncol(TabNFO)-6)))
						TabNFO[TabNFO[,"nJ"]%in%nomJ,nom] = geneJun[nJ,"CountJun"]
					}
					
					#	RI_amont+SE
					if(test_Jstart_in_intron&test_EndJ_EndX)	#	(Start in intron) ET (End exon connu)
					{
						colRI = "red"
						expl = "Start en 5p et end en fin d un exon connu"
						nomJ = paste(nJ, "=>SE+RI", sep="")
						line = paste(nomJ, "\t", gene, "\t", nomGene, "\t", "\t", geneJun[J,"StartJun"], "\t", geneJun[J,"EndJun"], "\t", geneJun[nJ,"CountJun"], "\t", expl, sep="")
						write(line, file=nfoname, append=TRUE)
						write(line, file="")
						matEV_RI = rbind(matEV_RI, cbind(newRI, as.matrix(rep("RI_5p+SE", nrow(newRI)))))
						
						if(sum(TabNFO[,"nJ"]%in%nomJ)==0)	TabNFO = rbind(TabNFO, c(gene, nomJ, Chr, geneJun[J,"StartJun"], geneJun[J,"EndJun"], expl, rep(0, ncol(TabNFO)-6)))
						TabNFO[TabNFO[,"nJ"]%in%nomJ,nom] = geneJun[nJ,"CountJun"]
					}
					
					#	RI_5p+SE
					if(test_Jstart_in_5p&test_EndJ_EndX)	#	(Start en 5p) ET (End exon connu)
					{
						colRI = "red"
						expl = "Start en 5p et end en fin d un exon connu"
						nomJ = paste(nJ, "=>SE+RI", sep="")
						line = paste(nomJ, "\t", gene, "\t", nomGene, "\t", "\t", geneJun[J,"StartJun"], "\t", geneJun[J,"EndJun"], "\t", geneJun[nJ,"CountJun"], "\t", expl, sep="")
						write(line, file=nfoname, append=TRUE)
						write(line, file="")
						matEV_RI = rbind(matEV_RI, cbind(newRI, as.matrix(rep("RI_5p+SE", nrow(newRI)))))
						
						if(sum(TabNFO[,"nJ"]%in%nomJ)==0)	TabNFO = rbind(TabNFO, c(gene, nomJ, Chr, geneJun[J,"StartJun"], geneJun[J,"EndJun"], expl, rep(0, ncol(TabNFO)-6)))
						TabNFO[TabNFO[,"nJ"]%in%nomJ,nom] = geneJun[nJ,"CountJun"]
					}
					
					#	SE+RI
					if(test_StartJ_StartX&test_Jend_in_intron)	#	(Start exon connu) ET (End in intron)
					{
						colRI = "red"
						expl = "Start au debut d un exon et Retention intergenique"
						nomJ = paste(nJ, "=>SE+RI", sep="")
						line = paste(nomJ, "\t", gene, "\t", nomGene, "\t", "\t", geneJun[J,"StartJun"], "\t", geneJun[J,"EndJun"], "\t", geneJun[nJ,"CountJun"], "\t", expl, sep="")
						write(line, file=nfoname, append=TRUE)
						write(line, file="")
						matEV_RI = rbind(matEV_RI, cbind(newRI, as.matrix(rep("RI_5p", nrow(newRI)))))
						
						if(sum(TabNFO[,"nJ"]%in%nomJ)==0)	TabNFO = rbind(TabNFO, c(gene, nomJ, Chr, geneJun[J,"StartJun"], geneJun[J,"EndJun"], expl, rep(0, ncol(TabNFO)-6)))
						TabNFO[TabNFO[,"nJ"]%in%nomJ,nom] = geneJun[nJ,"CountJun"]
					}
					
					#	RI_5p
					if(test_Jstart_in_5p&test_Jend_in_5p&!test_EndConnu&!test_StartConnu)	#	(Start in 5p) ET (End in 5p)
					{
						colRI = "red"
						expl = "Retention intergenique en 5p"
						nomJ = paste(nJ, "=>RI_5p", sep="")
						line = paste(nomJ, "\t", gene, "\t", nomGene, "\t", "\t", geneJun[J,"StartJun"], "\t", geneJun[J,"EndJun"], "\t", geneJun[nJ,"CountJun"], "\t", expl, sep="")
						write(line, file=nfoname, append=TRUE)
						write(line, file="")
						matEV_RI = rbind(matEV_RI, cbind(newRI, as.matrix(rep("RI_5p", nrow(newRI)))))
						
						if(sum(TabNFO[,"nJ"]%in%nomJ)==0)	TabNFO = rbind(TabNFO, c(gene, nomJ, Chr, geneJun[J,"StartJun"], geneJun[J,"EndJun"], expl, rep(0, ncol(TabNFO)-6)))
						TabNFO[TabNFO[,"nJ"]%in%nomJ,nom] = geneJun[nJ,"CountJun"]
					}
					
					#	RIG_5p
					if(test_Jstart_in_5p&test_EndConnu&test_Jend_EQ_Startexon1)	#	(Start in 5p) ET (End connu) ET (End = start 1er exon)
					{
						colRI = "red"
						expl = "Retention intergenique en 5p"
						nomJ = paste(nJ, "=>RIG_5p", sep="")
						line = paste(nomJ, "\t", gene, "\t", nomGene, "\t", "\t", geneJun[J,"StartJun"], "\t", geneJun[J,"EndJun"], "\t", geneJun[nJ,"CountJun"], "\t", expl, sep="")
						write(line, file=nfoname, append=TRUE)
						write(line, file="")
						matEV_RI = rbind(matEV_RI, cbind(newRI, as.matrix(rep("RIG_5p", nrow(newRI)))))
						
						if(sum(TabNFO[,"nJ"]%in%nomJ)==0)	TabNFO = rbind(TabNFO, c(gene, nomJ, Chr, geneJun[J,"StartJun"], geneJun[J,"EndJun"], expl, rep(0, ncol(TabNFO)-6)))
						TabNFO[TabNFO[,"nJ"]%in%nomJ,nom] = geneJun[nJ,"CountJun"]
					}
					
					#	EL
					if(test_Jstart_in_exon&test_Jend_in_exon&test_JstartJend_Meme_exon)	#	(Start in exon) ET (End in exon) ET (Start et End dans le meme exon)
					{
						colRI = "black"
						expl = "Exon Loss"
						nomJ = paste(nJ, "=>EL", sep="")
						line = paste(nomJ, "\t", gene, "\t", nomGene, "\t", ChrgeneJun[J,"StartJun"], "\t", geneJun[J,"EndJun"], "\t", geneJun[nJ,"CountJun"], "\t", expl, sep="")
						write(line, file=nfoname, append=TRUE)
						write(line, file="")
						matEV_RI = rbind(matEV_RI, cbind(newRI, as.matrix(rep("EL", nrow(newRI)))))
						
						if(sum(TabNFO[,"nJ"]%in%nomJ)==0)	TabNFO = rbind(TabNFO, c(gene, nomJ, Chr, geneJun[J,"StartJun"], geneJun[J,"EndJun"], expl, rep(0, ncol(TabNFO)-6)))
						TabNFO[TabNFO[,"nJ"]%in%nomJ,nom] = geneJun[nJ,"CountJun"]
					}
					
					#	A5SS+A3SS
					if(test_Jstart_in_exon&test_Jend_in_exon&(!test_JstartJend_Meme_exon))	#	(Start in exon) ET (End in exon) ET (Start et end dans 2 exons differents)
					{
						colRI = "black"
						expl = "Exon1 raccourci en aval + Exon2 raccourci en amont"
						nomJ = paste(nJ, "=>A5SS+A3SS", sep="")
						line = paste(nomJ, "\t", gene, "\t", nomGene, "\t", "\t", geneJun[J,"StartJun"], "\t", geneJun[J,"EndJun"], "\t", geneJun[nJ,"CountJun"], "\t", expl, sep="")
						write(line, file=nfoname, append=TRUE)
						write(line, file="")
						matEV_RI = rbind(matEV_RI, cbind(newRI, as.matrix(rep("A5SS+A3SS", nrow(newRI)))))
						
						if(sum(TabNFO[,"nJ"]%in%nomJ)==0)	TabNFO = rbind(TabNFO, c(gene, nomJ, Chr, geneJun[J,"StartJun"], geneJun[J,"EndJun"], expl, rep(0, ncol(TabNFO)-6)))
						TabNFO[TabNFO[,"nJ"]%in%nomJ,nom] = geneJun[nJ,"CountJun"]
					}
					
					#	A5SS
					if(test_StartConnu&test_Jend_in_exon)	#	(start connu) ET (End in exon)
					{
						colRI = "black"
						expl = "Exon2 raccourci en amont"
						nomJ = paste(nJ, "=>A5SS", sep="")
						line = paste(nomJ, "\t", gene, "\t", nomGene, "\t", "\t", geneJun[J,"StartJun"], "\t", geneJun[J,"EndJun"], "\t", geneJun[nJ,"CountJun"], "\t", expl, sep="")
						write(line, file=nfoname, append=TRUE)
						write(line, file="")
						matEV_RI = rbind(matEV_RI, cbind(newRI, as.matrix(rep("A5SS", nrow(newRI)))))
						
						if(sum(TabNFO[,"nJ"]%in%nomJ)==0)	TabNFO = rbind(TabNFO, c(gene, nomJ, Chr, geneJun[J,"StartJun"], geneJun[J,"EndJun"], expl, rep(0, ncol(TabNFO)-6)))
						TabNFO[TabNFO[,"nJ"]%in%nomJ,nom] = geneJun[nJ,"CountJun"]
					}
					
					#	A3SS
					if(test_Jstart_in_exon&test_EndConnu)	#	(start in exon) ET (End connu)
					{
						colRI = "black"
						expl = "Exon1 raccourci en aval"
						nomJ = paste(nJ, "=>A3SS", sep="")
						line = paste(nomJ, "\t", gene, "\t", nomGene, "\t", "\t", geneJun[J,"StartJun"], "\t", geneJun[J,"EndJun"], "\t", geneJun[nJ,"CountJun"], "\t", expl, sep="")
						write(line, file=nfoname, append=TRUE)
						write(line, file="")
						matEV_RI = rbind(matEV_RI, cbind(newRI, as.matrix(rep("A3SS", nrow(newRI)))))
						
						if(sum(TabNFO[,"nJ"]%in%nomJ)==0)	TabNFO = rbind(TabNFO, c(gene, nomJ, Chr, geneJun[J,"StartJun"], geneJun[J,"EndJun"], expl, rep(0, ncol(TabNFO)-6)))
						TabNFO[TabNFO[,"nJ"]%in%nomJ,nom] = geneJun[nJ,"CountJun"]
					}
					
					#	RI_aval
					if(test_Jstart_in_intron&test_EndConnu)	#	(Start in intron) ET (End connu)
					{
						colRI = "red"
						expl = "Jonction en aval de RI"
						nomJ = paste(nJ, "=>RI_aval", sep="")
						line = paste(nomJ, "\t", gene, "\t", nomGene, "\t", "\t", geneJun[J,"StartJun"], "\t", geneJun[J,"EndJun"], "\t", geneJun[nJ,"CountJun"], "\t", expl, sep="")
						write(line, file=nfoname, append=TRUE)
						write(line, file="")
						matEV_RI = rbind(matEV_RI, cbind(newRI, as.matrix(rep("RI_aval", nrow(newRI)))))
						
						if(sum(TabNFO[,"nJ"]%in%nomJ)==0)	TabNFO = rbind(TabNFO, c(gene, nomJ, Chr, geneJun[J,"StartJun"], geneJun[J,"EndJun"], expl, rep(0, ncol(TabNFO)-6)))
						TabNFO[TabNFO[,"nJ"]%in%nomJ,nom] = geneJun[nJ,"CountJun"]
					}
					
					#	RI_amont
					if(test_StartConnu&test_Jend_in_intron)	#	(Start connu) ET (End in intron)
					{
						colRI = "red"
						expl = "Jonction en amont de RI"
						nomJ = paste(nJ, "=>RI_amont", sep="")
						line = paste(nomJ, "\t", gene, "\t", nomGene, "\t", "\t", geneJun[J,"StartJun"], "\t", geneJun[J,"EndJun"], "\t", geneJun[nJ,"CountJun"], "\t", expl, sep="")
						write(line, file=nfoname, append=TRUE)
						write(line, file="")
						matEV_RI = rbind(matEV_RI, cbind(newRI, as.matrix(rep("RI_amont", nrow(newRI)))))
						
						if(sum(TabNFO[,"nJ"]%in%nomJ)==0)	TabNFO = rbind(TabNFO, c(gene, nomJ, Chr, geneJun[J,"StartJun"], geneJun[J,"EndJun"], expl, rep(0, ncol(TabNFO)-6)))
						TabNFO[TabNFO[,"nJ"]%in%nomJ,nom] = geneJun[nJ,"CountJun"]
					}
					
					#	RIs
					if(test_Jstart_in_intron&test_Jend_in_intron&test_JstartJend_Meme_intron)	#	(Start in intron) ET (End in intron) ET (Start et End dans le meme intron)
					{
						colRI = "red"
						expl = "Junction entre RI"
						nomJ = paste(nJ, "=>RIs", sep="")
						line = paste(nomJ, "\t", gene, "\t", nomGene, "\t", "\t", geneJun[J,"StartJun"], "\t", geneJun[J,"EndJun"], "\t", geneJun[nJ,"CountJun"], "\t", expl, sep="")
						write(line, file=nfoname, append=TRUE)
						write(line, file="")
						matEV_RI = rbind(matEV_RI, cbind(newRI, as.matrix(rep("RIs", nrow(newRI)))))
						
						if(sum(TabNFO[,"nJ"]%in%nomJ)==0)	TabNFO = rbind(TabNFO, c(gene, nomJ, Chr, geneJun[J,"StartJun"], geneJun[J,"EndJun"], expl, rep(0, ncol(TabNFO)-6)))
						TabNFO[TabNFO[,"nJ"]%in%nomJ,nom] = geneJun[nJ,"CountJun"]
					}
					
					#	SE+RIs
					if(test_Jstart_in_intron&test_Jend_in_intron&(!test_JstartJend_Meme_intron))	#	(Start in intron) ET (End in intron) ET (Start et End dans 2 introns differents) 
					{
						colRI = "purple"
						expl = "Junction introns / RIs avec saut exon1"
						nomJ = paste(nJ, "=>SE+RIs", sep="")
						line = paste(nomJ, "\t", gene, "\t", nomGene, "\t", "\t", geneJun[J,"StartJun"], "\t", geneJun[J,"EndJun"], "\t", geneJun[nJ,"CountJun"], "\t", expl, sep="")
						write(line, file=nfoname, append=TRUE)
						write(line, file="")
						matEV_RI = rbind(matEV_RI, cbind(newRI, as.matrix(rep("SE+RIs", nrow(newRI)))))
						
						if(sum(TabNFO[,"nJ"]%in%nomJ)==0)	TabNFO = rbind(TabNFO, c(gene, nomJ, Chr, geneJun[J,"StartJun"], geneJun[J,"EndJun"], expl, rep(0, ncol(TabNFO)-6)))
						TabNFO[TabNFO[,"nJ"]%in%nomJ,nom] = geneJun[nJ,"CountJun"]
					}
					
					#	A5SS+RI_amont
					if(test_Jstart_in_intron&test_Jend_in_exon&IE_contig)	#	(start intron) ET (End in exon) ET (IE contigus)
					{
						colRI = "brown"
						expl = "Jonction entre RI et exon raccourci en amont + start dans exon suivant end dans intron"
						nomJ = paste(nJ, "=>A5SS+RI_amont", sep="")
						line = paste(nomJ, "\t", gene, "\t", nomGene, "\t", "\t", geneJun[J,"StartJun"], "\t", geneJun[J,"EndJun"], "\t", geneJun[nJ,"CountJun"], "\t", expl, sep="")
						write(line, file=nfoname, append=TRUE)
						write(line, file="")
						matEV_RI = rbind(matEV_RI, cbind(newRI, as.matrix(rep("A5SS+RI_amont", nrow(newRI)))))
						
						if(sum(TabNFO[,"nJ"]%in%nomJ)==0)	TabNFO = rbind(TabNFO, c(gene, nomJ, Chr, geneJun[J,"StartJun"], geneJun[J,"EndJun"], expl, rep(0, ncol(TabNFO)-6)))
						TabNFO[TabNFO[,"nJ"]%in%nomJ,nom] = geneJun[nJ,"CountJun"]
					}
					
					#	A5SS+SE+RI_amont
					if(test_Jstart_in_intron&test_Jend_in_exon&(!IE_contig))	#	(start intron) ET (End in exon) ET (IE NON contigus)
					{
						colRI = "purple"
						expl = "Jonction entre RI et exon raccourci en amont avec start dans exon, end dans intron NON contigus"
						nomJ = paste(nJ, "=>A5SS+SE+RI_amont", sep="")
						line = paste(nomJ, "\t", gene, "\t", nomGene, "\t", "\t", geneJun[J,"StartJun"], "\t", geneJun[J,"EndJun"], "\t", geneJun[nJ,"CountJun"], "\t", expl, sep="")
						write(line, file=nfoname, append=TRUE)
						write(line, file="")
						matEV_RI = rbind(matEV_RI, cbind(newRI, as.matrix(rep("A5SS+SE+RI_amont", nrow(newRI)))))
						
						if(sum(TabNFO[,"nJ"]%in%nomJ)==0)	TabNFO = rbind(TabNFO, c(gene, nomJ, Chr, geneJun[J,"StartJun"], geneJun[J,"EndJun"], expl, rep(0, ncol(TabNFO)-6)))
						TabNFO[TabNFO[,"nJ"]%in%nomJ,nom] = geneJun[nJ,"CountJun"]
					}
					
					#	RI_5p+A5SS
					if(test_Jstart_in_5p&test_Jend_in_exon1)	#	(start dans 5p) ET (End dans 1er exon)
					{
						colRI = "brown"
						expl = "Junction entre 5p et Exon raccourci en amont"
						nomJ = paste(nJ, "=>RI_5p+A5SS", sep="")
						line = paste(nomJ, "\t", gene, "\t", nomGene, "\t", "\t", geneJun[J,"StartJun"], "\t", geneJun[J,"EndJun"], "\t", geneJun[nJ,"CountJun"], "\t", expl, sep="")
						write(line, file=nfoname, append=TRUE)
						write(line, file="")
						matEV_RI = rbind(matEV_RI, cbind(newRI, as.matrix(rep("RI_5p+A5SS", nrow(newRI)))))
						
						if(sum(TabNFO[,"nJ"]%in%nomJ)==0)	TabNFO = rbind(TabNFO, c(gene, nomJ, Chr, geneJun[J,"StartJun"], geneJun[J,"EndJun"], expl, rep(0, ncol(TabNFO)-6)))
						TabNFO[TabNFO[,"nJ"]%in%nomJ,nom] = geneJun[nJ,"CountJun"]
					}
					
					#	RI_5p+SE+A5SS
					if(test_Jstart_in_5p&test_Jend_in_exon&(!test_Jend_in_exon1))	#	(start dans 5p) ET (End dans exon) ET (End pas in exon1)
					{
						colRI = "brown"
						expl = "Junction entre 5p, saut exon et Exon raccourci en amont"
						nomJ = paste(nJ, "=>RI_5p+SE+A5SS", sep="")
						line = paste(nomJ, "\t", gene, "\t", nomGene, "\t", "\t", geneJun[J,"StartJun"], "\t", geneJun[J,"EndJun"], "\t", geneJun[nJ,"CountJun"], "\t", expl, sep="")
						write(line, file=nfoname, append=TRUE)
						write(line, file="")
						matEV_RI = rbind(matEV_RI, cbind(newRI, as.matrix(rep("RI_5p+SE+A5SS", nrow(newRI)))))
						
						if(sum(TabNFO[,"nJ"]%in%nomJ)==0)	TabNFO = rbind(TabNFO, c(gene, nomJ, Chr, geneJun[J,"StartJun"], geneJun[J,"EndJun"], expl, rep(0, ncol(TabNFO)-6)))
						TabNFO[TabNFO[,"nJ"]%in%nomJ,nom] = geneJun[nJ,"CountJun"]
					}
					
					#	A3SS+RI_aval
					if(test_Jstart_in_exon&test_Jend_in_intron&EI_contig)	#	(Start in exon) ET (End in intron) ET (EI contigus)
					{
						colRI = "brown"
						expl = "Exon raccourci en aval + RI + start dans exon suivant end dans intron"
						nomJ = paste(nJ, "=>A3SS+RI_aval", sep="")
						line = paste(nomJ, "\t", gene, "\t", nomGene, "\t", "\t", geneJun[J,"StartJun"], "\t", geneJun[J,"EndJun"], "\t", geneJun[nJ,"CountJun"], "\t", expl, sep="")
						write(line, file=nfoname, append=TRUE)
						write(line, file="")
						matEV_RI = rbind(matEV_RI, cbind(newRI, as.matrix(rep("A3SS+RI_aval", nrow(newRI)))))
						
						if(sum(TabNFO[,"nJ"]%in%nomJ)==0)	TabNFO = rbind(TabNFO, c(gene, nomJ, Chr, geneJun[J,"StartJun"], geneJun[J,"EndJun"], expl, rep(0, ncol(TabNFO)-6)))
						TabNFO[TabNFO[,"nJ"]%in%nomJ,nom] = geneJun[nJ,"CountJun"]
					}
					
					#	A3SS+SE+RI_aval
					if(test_Jstart_in_exon&test_Jend_in_intron&(!EI_contig))	#	(Start in exon) ET (End in intron) ET (EI NON contigus)
					{
						colRI = "purple"
						expl = "Exon raccourci en aval + RI + start dans exon, end dans intron NON contigus"
						nomJ = paste(nJ, "=>A3SS+SE+RI_aval", sep="")
						line = paste(nomJ, "\t", gene, "\t", nomGene, "\t", "\t", geneJun[J,"StartJun"], "\t", geneJun[J,"EndJun"], "\t", geneJun[nJ,"CountJun"], "\t", expl, sep="")
						write(line, file=nfoname, append=TRUE)
						write(line, file="")
						matEV_RI = rbind(matEV_RI, cbind(newRI, as.matrix(rep("A3SS+SE+RI_aval", nrow(newRI)))))
						
						if(sum(TabNFO[,"nJ"]%in%nomJ)==0)	TabNFO = rbind(TabNFO, c(gene, nomJ, Chr, geneJun[J,"StartJun"], geneJun[J,"EndJun"], expl, rep(0, ncol(TabNFO)-6)))
						TabNFO[TabNFO[,"nJ"]%in%nomJ,nom] = geneJun[nJ,"CountJun"]
					}
					
					#	A3SS+RI_3p
					if(test_Jstart_in_Derexon&test_Jend_in_3p)	#	(Start in Dernier exon) ET (End in 3p)
					{
						colRI = "brown"
						expl = "Dernier Exon raccourci en aval + fin de junction en 5p"
						nomJ = paste(nJ, "=>A3SS+RI_3p", sep="")
						line = paste(nomJ, "\t", gene, "\t", nomGene, "\t", "\t", geneJun[J,"StartJun"], "\t", geneJun[J,"EndJun"], "\t", geneJun[nJ,"CountJun"], "\t", expl, sep="")
						write(line, file=nfoname, append=TRUE)
						write(line, file="")
						matEV_RI = rbind(matEV_RI, cbind(newRI, as.matrix(rep("A3SS+RI_3p", nrow(newRI)))))
						
						if(sum(TabNFO[,"nJ"]%in%nomJ)==0)	TabNFO = rbind(TabNFO, c(gene, nomJ, Chr, geneJun[J,"StartJun"], geneJun[J,"EndJun"], expl, rep(0, ncol(TabNFO)-6)))
						TabNFO[TabNFO[,"nJ"]%in%nomJ,nom] = geneJun[nJ,"CountJun"]
					}
					
					#	SE+RI_3p
					if(test_StartConnu&test_Jend_in_3p)	#	(Start connu) ET (End in 3p)
					{
						colRI = "purple"
						expl = "Jonction prolongee en 3p avec saut exon"
						nomJ = paste(nJ, "=>SE+RI_3p", sep="")
						line = paste(nomJ, "\t", gene, "\t", nomGene, "\t", "\t", geneJun[J,"StartJun"], "\t", geneJun[J,"EndJun"], "\t", geneJun[nJ,"CountJun"], "\t", expl, sep="")
						write(line, file=nfoname, append=TRUE)
						write(line, file="")
						matEV_RI = rbind(matEV_RI, cbind(newRI, as.matrix(rep("SE+RI_3p", nrow(newRI)))))
						
						if(sum(TabNFO[,"nJ"]%in%nomJ)==0)	TabNFO = rbind(TabNFO, c(gene, nomJ, Chr, geneJun[J,"StartJun"], geneJun[J,"EndJun"], expl, rep(0, ncol(TabNFO)-6)))
						TabNFO[TabNFO[,"nJ"]%in%nomJ,nom] = geneJun[nJ,"CountJun"]
					}
					
					#	RI_5p+SE
					if(test_Jstart_in_5p&test_EndConnu)	#	(Start in 5p) ET (End connu)
					{
						colRI = "purple"
						expl = "Jonction prolongee en 5p avec saut exon"
						nomJ = paste(nJ, "=>RI_5p+SE", sep="")
						line = paste(nomJ, "\t", gene, "\t", nomGene, "\t", "\t", geneJun[J,"StartJun"], "\t", geneJun[J,"EndJun"], "\t", geneJun[nJ,"CountJun"], "\t", expl, sep="")
						write(line, file=nfoname, append=TRUE)
						write(line, file="")
						matEV_RI = rbind(matEV_RI, cbind(newRI, as.matrix(rep("RI_5p+SE", nrow(newRI)))))
						
						if(sum(TabNFO[,"nJ"]%in%nomJ)==0)	TabNFO = rbind(TabNFO, c(gene, nomJ, Chr, geneJun[J,"StartJun"], geneJun[J,"EndJun"], expl, rep(0, ncol(TabNFO)-6)))
						TabNFO[TabNFO[,"nJ"]%in%nomJ,nom] = geneJun[nJ,"CountJun"]
					}
					
					#	RI_5p+SE+RI
					if(test_Jstart_in_5p&test_Jend_in_intron)	#	(Start in 5p) ET (End in intron)
					{
						colRI = "purple"
						expl = "Jonction prolongée en 5p avec saut d’exon et RI"
						nomJ = paste(nJ, "=>RI_5p+SE+RI", sep="")
						line = paste(nomJ, "\t", gene, "\t", nomGene, "\t", "\t", geneJun[J,"StartJun"], "\t", geneJun[J,"EndJun"], "\t", geneJun[nJ,"CountJun"], "\t", expl, sep="")
						write(line, file=nfoname, append=TRUE)
						write(line, file="")
						matEV_RI = rbind(matEV_RI, cbind(newRI, as.matrix(rep("RI_5p+SE+RI", nrow(newRI)))))
						
						if(sum(TabNFO[,"nJ"]%in%nomJ)==0)	TabNFO = rbind(TabNFO, c(gene, nomJ, Chr, geneJun[J,"StartJun"], geneJun[J,"EndJun"], expl, rep(0, ncol(TabNFO)-6)))
						TabNFO[TabNFO[,"nJ"]%in%nomJ,nom] = geneJun[nJ,"CountJun"]
					}
					
					#	RIG_3p
					if(test_Jstart_EQ_EndDerexon&test_Jend_in_3p)	#	(Start = End dernier exon)) ET (End in 3p)
					{
						colRI = "red"
						expl = "Retention intergenique 3p"
						nomJ = paste(nJ, "=>RIG_3p", sep="")
						line = paste(nomJ, "\t", gene, "\t", nomGene, "\t", "\t", geneJun[J,"StartJun"], "\t", geneJun[J,"EndJun"], "\t", geneJun[nJ,"CountJun"], "\t", expl, sep="")
						write(line, file=nfoname, append=TRUE)
						write(line, file="")
						matEV_RI = rbind(matEV_RI, cbind(newRI, as.matrix(rep("RIG_3p", nrow(newRI)))))
						
						if(sum(TabNFO[,"nJ"]%in%nomJ)==0)	TabNFO = rbind(TabNFO, c(gene, nomJ, Chr, geneJun[J,"StartJun"], geneJun[J,"EndJun"], expl, rep(0, ncol(TabNFO)-6)))
						TabNFO[TabNFO[,"nJ"]%in%nomJ,nom] = geneJun[nJ,"CountJun"]
					}
					
					#	RI_aval+SE+RIG_3p
					if(test_Jstart_in_intron&test_Jend_in_3p)	#	(Start in intron) ET (End in 3p)
					{
						colRI = "purple"
						expl = "Retention intron en amont, saut exon et Retention intergenique 3p"
						nomJ = paste(nJ, "=>RI_amont+SE+RIG_3p", sep="")
						line = paste(nomJ, "\t", gene, "\t", nomGene, "\t", "\t", geneJun[J,"StartJun"], "\t", geneJun[J,"EndJun"], "\t", geneJun[nJ,"CountJun"], "\t", expl, sep="")
						write(line, file=nfoname, append=TRUE)
						write(line, file="")
						matEV_RI = rbind(matEV_RI, cbind(newRI, as.matrix(rep("RI_aval+SE+RIG_3p", nrow(newRI)))))
						
						if(sum(TabNFO[,"nJ"]%in%nomJ)==0)	TabNFO = rbind(TabNFO, c(gene, nomJ, Chr, geneJun[J,"StartJun"], geneJun[J,"EndJun"], expl, rep(0, ncol(TabNFO)-6)))
						TabNFO[TabNFO[,"nJ"]%in%nomJ,nom] = geneJun[nJ,"CountJun"]
					}
					
#					#	Ajoute le nouvel element a la carte
#					Map = rbind(Map, c(paste("Start_", nJ, sep=""), J_anorm[J,"StartJun"], matEV_RI[nrow(matEV_RI),"Type"]), c(paste("End_", nJ, sep=""), J_anorm[J,"EndJun"], matEV_RI[nrow(matEV_RI),"Type"]))
#					
#					#	Affichage
#					text(x=(J_anorm[J,"EndJun"]-J_anorm[J,"StartJun"])/2+J_anorm[J,"StartJun"], y=(Tr*10+rangJ*JunH)+4.7, label=nJ, cex=0.5, col=colRI)	#	Nom J
#					segments(x0=J_anorm[J,"StartJun"], y0=Tr*10+4+rangJ*JunH, x1=J_anorm[J,"EndJun"], col=colRI)	#	Segment horizontal
#					
#					#	Segments verticaux
#					if(!test_StartConnu)	segments(x0=J_anorm[J,"StartJun"], y0=Tr*10+4+rangJ*JunH, y1=0, col=colRI, lty="dotted")
#					if(!test_EndConnu)	segments(x0=J_anorm[J,"EndJun"], y0=Tr*10+4+rangJ*JunH, y1=0, col=colRI, lty="dotted")
				}
				write.table(resNFOjun, file=paste(pathResDataJunc, "/", nom, "_nfo_J_Gene.txt", sep=""), sep="\t", quote=FALSE, row.names=FALSE)
				
				if(nrow(matEV_RI)>0)
				{
					#	Calcul de score pour chaque Junction normale
					matEV_RI = cbind(matEV_RI, matrix(0, ncol=1, nrow=nrow(matEV_RI)))
					colnames(matEV_RI)[ncol(matEV_RI)] = "score"
					
					#	matEV_RI[matEV_RI[,"Junc_Normale"]=="","Junc_Normale"]="---"	#	remplace les vides par ---
					for(JN in names(table(matEV_RI[,"Junc_Normale"])))	#	JN = names(table(matEV_RI[,"J_N"]))[1]
					{
						tmp = matEV_RI[matEV_RI[,"Junc_Normale"]==JN,,drop=FALSE]
						
						if(as.numeric(tmp[1,"Junc_Normale_Count"])==0)
						{
							matEV_RI[matEV_RI[,"Junc_Normale"]==JN,"score"] = max(as.numeric(tmp[,"Junc_RI_Count"]))/1	#0.01
						}else{
							matEV_RI[matEV_RI[,"Junc_Normale"]==JN,"score"] = max(as.numeric(tmp[,"Junc_RI_Count"]))/as.numeric(tmp[1,"Junc_Normale_Count"])
						}
					}
					
					#	Nombre de reads junc / nombre max de junctions pour le gene
					maxNumbJ = 0
					minNumbJ = length(unlist(JuncT))
					for(V in 1:length(JuncT))
					{
						maxNumbJ = max(maxNumbJ, length(JuncT[[V]]))
						minNumbJ = min(minNumbJ, length(JuncT[[V]]))
					}
					
					MoyNcountParJunc = matrix(c(rep(sum(geneJun[,"CountJun"])/maxNumbJ, nrow(matEV_RI)),rep(sum(geneJun[,"CountJun"])/minNumbJ, nrow(matEV_RI))), ncol=2, byrow=FALSE)
					colnames(MoyNcountParJunc)=c("minMoyNcountParJunc", "maxMoyNcountParJunc")
					matEV_RI = cbind(matEV_RI, MoyNcountParJunc)
					write.table(matEV_RI, file=paste(pathResSample, "/", nom, "_RI.txt", sep=""), sep="\t", quote=FALSE, row.names=FALSE)
				}
			}
			
			###########################################################################################
#	SE
#			
			cols = c("Junc_SE", "ENSID", "Gene", "Chr", "Junc_SE_Start", "Junc_SE_End", "Junc_SE_Count", "NbreSkippedExons", "SkippedExonsID", "TranscritsID", "Junc_Normale", "Junc_Normale_Count")
			matEV_SE = matrix(0, nrow=0, ncol=length(cols))
			colnames(matEV_SE) = cols
			marge = limDecal
			
			tmpAlljunc = cbind(AllJunc, t(apply(AllJunc, 1, function(x) unlist(strsplit(x[2], "_")))))
			colnames(tmpAlljunc) = c("nomsJ", "AllJunc", "StartJun", "EndJun")
			
			tmpAlljunc = tmpAlljunc[tmpAlljunc[,"nomsJ"]%in%rownames(geneJun),,drop=FALSE]
			
			if(nrow(geneJun)>0)	for(J in 1:nrow(geneJun))
			{
				#	Liste des exons inclus dans la jonction
				SE_StartEnd = (exons[,"exon_chrom_start"]>=(geneJun[J,"StartJun"]-marge))&(exons[,"exon_chrom_end"]<=(geneJun[J,"EndJun"]+marge))
				
				if(sum(SE_StartEnd)>0)	#	Saut Exon
				{
					resSkippedExons = exons[SE_StartEnd,, drop=FALSE]
					tmpTranscrit = table(resSkippedExons[,"ensembl_transcript_id"])
					tmpSkippedExons = table(resSkippedExons[,"ensembl_exon_id"])
					nomSkippedExons = paste(names(tmpSkippedExons), collapse=",")
					nomTranscrits = paste(names(tmpTranscrit), collapse=",")
					
					if(min(tmpTranscrit)<max(tmpTranscrit))	#	si plusieurs comptages => Range
					{
						NSkippedExons = paste(min(tmpTranscrit), ";", max(tmpTranscrit), sep="")
					}else{
						NSkippedExons = min(tmpTranscrit)
					}
					
					if(sum(colnames(TabNFO)%in%nom)==0)	#	Si ech pas encore dans le tableau, Ajoute une colonne
					{
						TabNFO = cbind(TabNFO, rep(0, nrow(TabNFO)))
						colnames(TabNFO)[ncol(TabNFO)] = nom
					}
					
					coordsJ = rownames(geneJun)[J]
					nomJ = paste(coordsJ, "=>SE_", NSkippedExons, sep="")
					ColSE = "green"	#	SE anormal
					expl = paste("Saute ", NSkippedExons, " exon(s) => Transcrit Alternatif Connu", sep="")
					
					if(sum(AllJunc[,"nomsJ"]%in%coordsJ)<1)	#	SE inconnu, anormal
					{
						ColSE = "orange"
						expl = paste("Saute ", sum(SE_StartEnd), " exon(s)", sep="")
					}
					
					line = paste(nomJ, "\t", geneJun[J,"StartJun"], "\t", geneJun[J,"EndJun"], "\t", geneJun[J,"CountJun"], "\tSaute ", NSkippedExons, " exon(s) touche ", length(tmpTranscrit), " transcrit(s)", sep="")
					write(line, file="")
					write(line, file=nfoname, append=TRUE)
					EV_SE = exons[SE_StartEnd, c("ensembl_transcript_id", "ensembl_exon_id", "exon_chrom_start", "exon_chrom_end"), drop=FALSE]
					
					if(sum(TabNFO[,"nJ"]%in%nomJ)==0)	#	Junc inconnue
					{
						write(paste("Jonction inconnue (", nomJ, ") => ajoute une ligne de ", length(c(nomJ, geneJun[J,"StartJun"], geneJun[J,"EndJun"], expl, rep(0, ncol(TabNFO)-6))), " elements", sep=""), file="")
						
						TabNFO = rbind(TabNFO, c(gene, nomJ, Chr, geneJun[J,"StartJun"], geneJun[J,"EndJun"], expl, rep(0, ncol(TabNFO)-6)))
						TabNFO[TabNFO[,"nJ"]%in%nomJ, nom] = geneJun[coordsJ,"CountJun",drop=FALSE]
					}else{
						TabNFO[TabNFO[,"nJ"]%in%nomJ,nom] = geneJun[coordsJ,"CountJun",drop=FALSE]	#	geneJun[nJ,"CountJun"]
					}
					
					#	Liste les Junctions canoniques inclues dans le SE
					canonJunc = AllgeneJun[rownames(AllgeneJun)%in%rownames(TabAllJunc),,drop=FALSE]
					#	filtre la junction elle même
					canonJunc = canonJunc[!(canonJunc[,"Coords"]%in%paste(geneJun[J,"StartJun"], "_", geneJun[J,"EndJun"], sep="")),,drop=FALSE]
					
					test_JcanonInclue =
							((geneJun[J,"StartJun"]-marge)<=as.numeric(canonJunc[,"StartJun"]))&((geneJun[J,"EndJun"]+marge)>=as.numeric(canonJunc[,"EndJun"]))	#|	#	Junc canonique incluse dans le SE
							#((geneJun[J,"StartJun"]-marge)>=as.numeric(canonJunc[,"StartJun"]))&((geneJun[J,"StartJun"]-marge)<=as.numeric(canonJunc[,"EndJun"]))|	#	Start juncSE dans Canon Junc
							#((geneJun[J,"EndJun"]-marge)>=as.numeric(canonJunc[,"StartJun"]))&((geneJun[J,"EndJun"]-marge)<=as.numeric(canonJunc[,"EndJun"]))	#	End juncSE dans Canon Junc
					
					maxJuncNorm = 0
					if(sum(test_JcanonInclue)>0)
					{
						canonJuncInc = canonJunc[test_JcanonInclue,,drop=FALSE]

						if(sum(rownames(geneJun)%in%canonJuncInc[,"JunctionUniqueName"])>0)	#	Junc Canonique inclue avec comptage
						{
							tmpJuncNormCounts = geneJun[rownames(geneJun)%in%canonJuncInc[,"JunctionUniqueName"],,drop=FALSE]
							maxJuncNorm = max(as.numeric(tmpJuncNormCounts[,"CountJun"]))
							nomMaxJuncNorm = rownames(tmpJuncNormCounts)[as.numeric(tmpJuncNormCounts[,"CountJun"])==maxJuncNorm]
						}else{#	pas de comptages pour les junctions canoniques inclues dans le SE => 0
							maxJuncNorm = 0
							nomMaxJuncNorm = canonJuncInc[1,"JunctionUniqueName"]
						}
						if(length(nomMaxJuncNorm)>1)	nomMaxJuncNorm = paste(nomMaxJuncNorm, collapse=",")
						matEV_SE = rbind(matEV_SE, c(rownames(geneJun)[J], gene, nomGene, Chr, geneJun[J,c("StartJun", "EndJun", "CountJun"),drop=FALSE], NSkippedExons, nomSkippedExons, nomTranscrits, nomMaxJuncNorm, maxJuncNorm))
					}else{#	Pas de junction normale / canonique dans le SE
						matEV_SE = rbind(matEV_SE, c(rownames(geneJun)[J], gene, nomGene, Chr, geneJun[J,c("StartJun", "EndJun", "CountJun"),drop=FALSE], NSkippedExons, nomSkippedExons, nomTranscrits, "Pas_de_jonction_canonique_inclue_dans_le_SE", "0"))
					}
					
#					#	Affichage
#					text(x=(geneJun[coordsJ,"EndJun"]-geneJun[coordsJ,"StartJun"])/2+geneJun[coordsJ,"StartJun"], y=(Tr*10+J*JunH)+4.7, label=coordsJ, cex=0.5, col=ColSE)
#					segments(x0=geneJun[coordsJ,"StartJun"], y0=Tr*10+4+J*JunH, x1=geneJun[coordsJ,"EndJun"], col=ColSE)
#					segments(x0=geneJun[coordsJ,"StartJun"], y0=Tr*10+4+J*JunH, y1=0, col="orange", lty="dotted")
#					segments(x0=geneJun[coordsJ,"EndJun"], y0=Tr*10+4+J*JunH, y1=0, col="orange", lty="dotted")
					
					Map = rbind(Map, c(paste("Start_", coordsJ, sep=""), geneJun[coordsJ,"StartJun"], "SE"), c(paste("End_", coordsJ, sep=""), geneJun[coordsJ,"EndJun"], "SE"))
				}
			}
			Map = Map[order(as.numeric(Map[,"Pos"])),]
			write.table(Map, file=paste(pathResDataJunc, "/", nom, "_Map_J.txt", sep=""), sep="\t", quote=FALSE, row.names=FALSE)
			
			if(nrow(matEV_SE)>0)
			{
				#	Score
				matEV_SE = cbind(matEV_SE, matrix(0, ncol=1, nrow=nrow(matEV_SE)))
				colnames(matEV_SE)[ncol(matEV_SE)] = "score"
				for(SE in 1:nrow(matEV_SE))
				{
					if(!is.na(as.numeric(matEV_SE[SE,"Junc_Normale_Count"])))
					{
						if(max(as.numeric(matEV_SE[SE,"Junc_Normale_Count"]))==0)
						{
							matEV_SE[SE,"score"] = as.numeric(matEV_SE[SE,"Junc_SE_Count"])/1	#0.01
						}else{
							matEV_SE[SE,"score"] = as.numeric(matEV_SE[SE,"Junc_SE_Count"])/max(as.numeric(matEV_SE[SE,"Junc_Normale_Count"]))
						}
					}else{
						matEV_SE[SE,"score"] = 0
					}
				}
				#	Ajoute une colonne J canoniques
				matEV_SE = cbind(matEV_SE, matrix("---", ncol=1, nrow=nrow(matEV_SE)))
				colnames(matEV_SE)[ncol(matEV_SE)] = "Type"
				#	AllJunc = read.table(file=paste(pathEv, "/JonctionsNormales.txt", sep=""), sep="\t", header=TRUE)
				#	matEV_SE = tmpSE
				matEV_SE[matEV_SE[,"Junc_SE"]%in%AllJunc[,"nomsJ"],"Type"] = "SE_Canonique"
				
				#	Nombre de reads junc / nombre max de junctions pour le gene
				maxNumbJ = 0
				minNumbJ = length(unlist(JuncT))
				for(V in 1:length(JuncT))
				{
					maxNumbJ = max(maxNumbJ, length(JuncT[[V]]))
					minNumbJ = min(minNumbJ, length(JuncT[[V]]))
				}		
				MoyNcountParJunc = matrix(c(rep(sum(geneJun[,"CountJun"])/maxNumbJ, nrow(matEV_SE)),rep(sum(geneJun[,"CountJun"])/minNumbJ, nrow(matEV_SE))), ncol=2, byrow=FALSE)
				colnames(MoyNcountParJunc)=c("minMoyNcountParJunc", "maxMoyNcountParJunc")
				matEV_SE = cbind(matEV_SE, MoyNcountParJunc)
				tmpmatEV_SE = cbind(matEV_SE[,1], rep(Chr, nrow(matEV_SE)), matEV_SE[,2:ncol(matEV_SE), drop=FALSE])
				colnames(tmpmatEV_SE)[1:2] = c("Junc_SE", "Chr")
				write.table(tmpmatEV_SE, file=paste(pathResSample, "/", nom, "_SE.txt", sep=""), sep="\t", quote=FALSE, row.names=FALSE)
			}
#			dev.off()
		}

		if(nrow(TabNFO)>0)
		{
			ColI = grep("Interpretation", colnames(TabNFO))
			TotEchs = apply(TabNFO[,(ColI+1):ncol(TabNFO),drop=FALSE], 1, function(x) sum(as.numeric(x)>0))
			TabNFO = cbind(TabNFO[,1:ColI,drop=FALSE], apply(TabNFO[,(ColI+1):ncol(TabNFO),drop=FALSE], 1, function(x) sum(as.numeric(x)>0)), TabNFO[,(ColI+1):ncol(TabNFO),drop=FALSE])
			colnames(TabNFO)[(ColI+1)] = "TotEchs"
			colnames(TabNFO) = gsub(paste(gene, "_", Chr, "_", sep=""), "", colnames(TabNFO))
			
			#	Complete la matrice de comptage du NFO junc
			allSamples = gsub(".bam$", "", basename(bamsGene))
			manque = allSamples[!allSamples%in%colnames(TabNFO)[(ColI+2):ncol(TabNFO)]]
			if(length(manque)>0)
			{
				tmpMat = matrix(0, ncol=length(manque), nrow=nrow(TabNFO))
				colnames(tmpMat) = manque
				TabNFO = cbind(TabNFO, tmpMat)
			}
			
			TabNFO = TabNFO[order(as.numeric(TabNFO[,"TotEchs"]), decreasing=TRUE),,drop=FALSE]
			colnames(TabNFO) = gsub(paste(gene, "_", Chr, "_", sep=""), "", colnames(TabNFO))
			write.table(TabNFO, file=paste(pathEv, "tableNFOjunc.txt", sep=""), sep="\t", quote=FALSE, row.names=FALSE)
			
			colnames(countAllJunc) = gsub(paste(gene, "_", Chr, "_", sep=""), "", colnames(countAllJunc))
			write.table(countAllJunc, file=paste(pathEv, "AllCounts.txt", sep=""), sep="\t", quote=FALSE, row.names=FALSE)
			
			#	formatage excel table patients * evenements
			library(openxlsx)
			colTotEchs = grep("TotEchs", colnames(TabNFO))
			colSamples = c((colTotEchs+1):ncol(TabNFO))
			
			TempRes = as.data.frame(TabNFO)
			for(C in colSamples)	TempRes[,C] = as.numeric(as.character(TempRes[,C]))
			TempRes[,colTotEchs] = as.numeric(as.character(TempRes[,colTotEchs]))
			
			wb <- createWorkbook()
			addWorksheet(wb, "Patients_Evenements", gridLines = TRUE)
			boldHeader <- createStyle(textDecoration = 'bold') 
			writeData(wb, 1, TempRes, headerStyle = boldHeader)
			
			SpeStyle <- createStyle(fontColour = "black", bgFill = "red")
			NoCountStyle <- createStyle(fontColour = "black", bgFill = "grey")
			OneCountStyle <- createStyle(fontColour = "black", bgFill = "green")
			LowCountStyle <- createStyle(fontColour = "black", bgFill = "orange")
			HighCountStyle <- createStyle(fontColour = "black", bgFill = "red")
			
			conditionalFormatting(wb, "Patients_Evenements", cols=colTotEchs, rows=2:(nrow(TempRes)+1), rule="<=1", style = SpeStyle, gradient=TRUE)
			conditionalFormatting(wb, "Patients_Evenements", cols=colSamples, rows=2:(nrow(TempRes)+1), rule=">=20", style = HighCountStyle, gradient=TRUE)
			conditionalFormatting(wb, "Patients_Evenements", cols=colSamples, rows=2:(nrow(TempRes)+1), rule="<20", style = LowCountStyle, gradient=TRUE)
			conditionalFormatting(wb, "Patients_Evenements", cols=colSamples, rows=2:(nrow(TempRes)+1), rule="<=1", style = OneCountStyle, gradient=TRUE)
			conditionalFormatting(wb, "Patients_Evenements", cols=colSamples, rows=2:(nrow(TempRes)+1), rule="<1", style = NoCountStyle, gradient=TRUE)
			
			saveWorkbook(wb, file=paste(pathEv, "/tableNFOjunc.xlsx", sep=""), overwrite=TRUE)
			
			if(nrow(speJ)>0)
			{
				nomsJ = NULL
				for(E in 1:nrow(speJ))	if(sum(grepl(speJ[E,1], rownames(AllgeneJun)))>0){nomsJ = c(nomsJ, AllgeneJun[rownames(AllgeneJun)==speJ[E,1], "JunctionUniqueName"])}else{nomsJ = c(nomsJ, "")}
				rownames(speJ) = nomsJ
				write.table(speJ, file=paste(pathEv, "/JuncSPEcounts.txt", sep=""), sep="\t", quote=FALSE, col.names=c(paste("JunID\t", colnames(speJ)[1], sep=""), colnames(speJ)[2:length(colnames(speJ))]))
			}
			
			#	Regroupement, comparaison des evenements
			write(paste("\npathResSample=\"", pathResSample, "\"", sep=""), file="")
			write(paste("pathEv=\"", pathEv, "\"", sep=""), file="")
			write(paste("gene=\"", gene, "\"", sep=""), file="")
			write(paste("Chr=\"", Chr, "\"", sep=""), file="")
			
			AllRes = list.files(pathResSample, full.names=TRUE)
			
			#	Definition du nom du sample Pat, les autres res seront filtres
			SamplesTypes = as.matrix(read.table(pathtFileSamplesTypes, header=TRUE, sep="\t", quote=""))
			SamplesTypes[,"Sample"] = apply(SamplesTypes, 1, function(x) paste(x[1], "_", x[3], sep=""))
			SamplePat = SamplesTypes[grepl("Pat", SamplesTypes[,"Type"]), "Sample"]
			
			if(sum(grepl("_RI.txt", AllRes))>0)
			{
				All_RI = AllRes[grepl("_RI.txt", AllRes)]
				tmpCol = as.matrix(read.table(All_RI[1], header=TRUE, sep="\t", quote=""))
				tmpRI = matrix(0, ncol=(ncol(tmpCol)+1), nrow=0)
				for(Q in 1:length(All_RI))
				{
					tmp = as.matrix(read.table(All_RI[Q], header=TRUE, sep="\t", quote=""))
					tmp = cbind(tmp, rep(gsub("_RI.txt", "", basename(All_RI[Q])), nrow(tmp)))
					tmpRI = rbind(tmpRI, tmp)
					
				}
				#	table(tmpRI[,ncol(tmpRI)])
				colnames(tmpRI)[ncol(tmpRI)] = "Sample"
				#	filtre pour la comparaison
				#	tmpRI = tmpRI[grepl(SamplePat, tmpRI[,"Sample"]),, drop=FALSE]	#	Filtre les samples
				
				#	filtre les duplicats
				uniques = apply(tmpRI[,c("Chr", "Junc_RI_Start", "Junc_RI_End", "Sample"), drop=FALSE], 1, function(x) paste(x, sep="_", collapse="_"))
				tmpRI = tmpRI[!duplicated(uniques),,drop=FALSE]
				
				if(nrow(tmpRI)>0) write.table(tmpRI, file=paste(pathEv, "/resAll_RI.txt", sep=""), sep="\t", quote=FALSE, row.names=FALSE)
			}
			
			if(sum(grepl("_SE.txt", AllRes))>0)
			{
				All_SE = AllRes[grepl("_SE.txt", AllRes)]
				tmpCol = as.matrix(read.table(All_SE[1], header=TRUE, sep="\t", quote=""))
				tmpSE = matrix(0, ncol=(ncol(tmpCol)+1), nrow=0)
				for(Q in 1:length(All_SE))
				{
					tmp = as.matrix(read.table(All_SE[Q], header=TRUE, sep="\t", quote=""))
					tmp = cbind(tmp, rep(gsub("_SE.txt", "", basename(All_SE[Q])), nrow(tmp)))
					tmpSE = rbind(tmpSE, tmp)
				}
				colnames(tmpSE)[ncol(tmpSE)] = "Sample"
				#	tmpSE = tmpSE[grepl(SamplePat, tmpSE[,"Sample"]),, drop=FALSE]	#	Filtre les samples
	
				#	filtre les duplicats
				uniques = apply(tmpSE[,c("Chr", "Junc_SE_Start", "Junc_SE_End", "Sample"), drop=FALSE], 1, function(x) paste(x, sep="_", collapse="_"))
				tmpSE = tmpSE[!duplicated(uniques),,drop=FALSE]
	
				if(nrow(tmpSE)>0) write.table(tmpSE, file=paste(pathEv, "/resAll_SE.txt", sep=""), sep="\t", quote=FALSE, row.names=FALSE)
			}

			AllRes = list.files(pathEv, full.names=TRUE)	#	AllRes mis a jour apres filtrage

			########################################################################################################
			#	Filtrage et tests
			if(Comps)
			{
				AllCounts = as.matrix(read.table(paste(pathEv, "AllCounts.txt", sep=""), header=TRUE, check.names=FALSE, sep="\t", quote=""))
				#	N junc tot pour tous pour le score normalisé
				TotCountJuncs_All = apply(AllCounts[,(grep("Type", colnames(AllCounts))+1):ncol(AllCounts),drop=FALSE], 2, function(x) sum(as.numeric(x)))
				
				if(sum(grepl("_RI.txt", AllRes))>0)
				{
					All_RI = as.matrix(read.table(file=paste(pathEv, "/resAll_RI.txt", sep=""), header=TRUE, sep="\t", quote=""))
					All_RI[,"Sample"] = gsub(paste(gene, "_", Chr, "_", sep=""), "", All_RI[,"Sample"])
					
					#	Correction des scores => score * moyTotJuncCount / SampleTotJuncCount
					ScoreCorr = apply(All_RI, 1, function(x) as.numeric(x[grep("score", colnames(All_RI))])*(mean(as.numeric(TotCountJuncs_All))/as.numeric(TotCountJuncs_All[x[grep("Sample", colnames(All_RI))]])))
					All_RI = cbind(All_RI, apply(All_RI, 1, function(x) sum(as.numeric(AllCounts[,x[grep("Sample", colnames(All_RI))]]))))
					
					colnames(All_RI)[ncol(All_RI)] = "TotCountJuncs"
					All_RI = cbind(All_RI, ScoreCorr)
					colnames(All_RI)[ncol(All_RI)] = "ScoreCorrige"
					
					colSample = grep("Sample", colnames(All_RI))
					colENSID = grep("ENSID", colnames(All_RI))
					colChr = grep("Chr", colnames(All_RI))
					
					SamplesTypes = as.matrix(read.table(pathtFileSamplesTypes, header=TRUE, sep="\t", quote=""))
					SamplesTypes[,"Sample"] = apply(SamplesTypes, 1, function(x) paste(x[1], "_", x[3], sep=""))
					#	SampleNames = unlist(apply(SamplesTypes, 1, function(x) unique(All_RI[grepl(paste("_", x[1], sep=""), All_RI[,"Sample"]),"Sample"])))			
					SampleNames = unique(All_RI[,"Sample"])
					
					# renome les samples
					#for(W in 1:nrow(SamplesTypes)) SamplesTypes[W,"Sample"] = paste(All_RI[1,"ENSID"], "_", All_RI[1,"Chr"], "_", SamplesTypes[W,"Sample"], "_", SamplesTypes[W,"Proj"], sep="")
					
					Pat = SampleNames[SampleNames%in%SamplesTypes[grepl("Pat", SamplesTypes[,"Type"]), "Sample"]]
					Ctrl = SampleNames[SampleNames%in%SamplesTypes[grepl("Ctrl", SamplesTypes[,"Type"]), "Sample"]]

					Pat = Pat[Pat%in%unique(All_RI[,"Sample"])]	#	selectionne en fonction des res
					Ctrl = Ctrl[Ctrl%in%unique(All_RI[,"Sample"])]	#	selectionne en fonction des res
					All_Ctrl = All_RI[All_RI[,"Sample"]%in%Ctrl,,drop=FALSE]
					
					colJunc_RI_Start = grepl("Junc_RI_Start", colnames(All_Ctrl))
					colJunc_RI_End = grepl("Junc_RI_End", colnames(All_Ctrl))
					Ctrls = apply(All_Ctrl, 1, function(x) paste(as.numeric(x[colJunc_RI_Start]), "-", as.numeric(x[colJunc_RI_End]), sep=""))
					Ctrls = Ctrls[!grepl("0-0", Ctrls)]
					
#					resRI_f = matrix(0, ncol=(ncol(All_RI)+11), nrow=0)	#	Matrice de resultats des SE filtrés
#					colnames(resRI_f) = c(colnames(All_RI), "pvalFisher", "pval", "MoyScoreCtrls", "Nctrls", "ratioScores(P_vs_Ctrl)", "deltaScores(P_vs_Ctrl)", "pvalCorr", "MoyScoreCtrlsCorr", "ratioScores(P_vs_Ctrl)Corr", "deltaScores(P_vs_Ctrl)Corr", "Com_ou_Spe")
					
					if(length(Pat)>0)	for(P in 1:length(Pat))	#	Pat par pat
					{
						tmpPat = All_RI[All_RI[,"Sample"]==Pat[P],,drop=FALSE]
						colJunc_RI_Start = grep("Junc_RI_Start", colnames(tmpPat))
						colJunc_RI_End = grep("Junc_RI_End", colnames(tmpPat))
						
						#	Filtre les evenements avec read = 1
						tmpPat = tmpPat[as.numeric(tmpPat[,"Junc_RI_Count"])>=2,,drop=FALSE]
						if(nrow(tmpPat)>=1)	#	Si au moins un evenement avec read >1
						{
							#	Junctions normales=>evenements
							AllRIs = apply(tmpPat, 1, function(x) paste(as.numeric(x[colJunc_RI_Start]), "-", as.numeric(x[colJunc_RI_End]), sep=""))
							RIs = unique(AllRIs)
							RIs = RIs[!grepl("0-0", RIs)]
							
							for(E in 1:length(RIs))
							{
								if(sum(Ctrls==RIs[E])>=1)	#	Evenement retrouve chez les ctrls
								{
									evCtrl = All_Ctrl[Ctrls==RIs[E],, drop=FALSE]
									pval = try(t.test(x=as.numeric(evCtrl[,"score"]), mu=mean(as.numeric(tmpPat[AllRIs==RIs[E],"score"])))$p.value, silent=TRUE)
									if(class(pval) == "try-error")	pval = 1
									ratio = mean(as.numeric(tmpPat[AllRIs==RIs[E],"score"]))/mean(as.numeric(evCtrl[,"score"]))
									delta = abs(mean(as.numeric(tmpPat[AllRIs==RIs[E],"score"]))-mean(as.numeric(evCtrl[,"score"])))
									moyC = mean(as.numeric(evCtrl[,"score"]))
									
									pvalCorr = try(t.test(x=as.numeric(evCtrl[,"ScoreCorrige"]), mu=mean(as.numeric(tmpPat[AllRIs==RIs[E],"ScoreCorrige"])))$p.value, silent=TRUE)
									if(class(pvalCorr) == "try-error")	pvalCorr = 1
									ratioCorr = mean(as.numeric(tmpPat[AllRIs==RIs[E],"ScoreCorrige"]))/mean(as.numeric(evCtrl[,"ScoreCorrige"]))
									deltaCorr = abs(mean(as.numeric(tmpPat[AllRIs==RIs[E],"ScoreCorrige"]))-mean(as.numeric(evCtrl[,"ScoreCorrige"])))
									moyCcorr = mean(as.numeric(evCtrl[,"ScoreCorrige"]))
									
									#	test proportions Pat vs Ctrl
									
									RIpat = as.numeric(tmpPat[E,"Junc_RI_Count"])
									RIctrl = sum(sum(as.numeric(evCtrl[,"Junc_RI_Count"])))
									NormPat = as.numeric(tmpPat[E,"Junc_Normale_Count"])
									NormCtrl = sum(as.numeric(evCtrl[,"Junc_Normale_Count"]))
									
									#	http://www.sthda.com/english/wiki/two-proportions-z-test-in-r
									#res <- prop.test(x = c(RIpat, RIctrl), n = c(NormPat, NormCtrl))
									
									tab = matrix(c(NormPat, NormCtrl, RIpat, RIctrl), ncol=2)
									colnames(tab) = c("Norm", "RI")
									rownames(tab) = c("pat", "Ctrl")
									#resChi = chisq.test(tab, correct=TRUE)
									resFish = fisher.test(tab, conf.int=TRUE, conf.level=0.95)
									pvalFish = resFish$p.value
									
#									new = cbind(tmpPat[AllRIs==RIs[E],,drop=FALSE], matrix(rep(c(pvalFish, pval, moyC, nrow(evCtrl), ratio, delta, pvalCorr, moyCcorr, ratioCorr, deltaCorr, "ev_Com_C_P"), sum(AllRIs==RIs[E])), ncol=11, byrow=TRUE))
#									resRI_f = rbind(resRI_f, new)
								}else{					#	Evenement specifique du patient
									pval = try(t.test(x=rep(0, length(Ctrl)), mu=mean(as.numeric(tmpPat[AllRIs==RIs[E],"score"])))$p.value, silent=TRUE)
									if(class(pval) == "try-error")	pval = 1
									ratio = mean(as.numeric(tmpPat[AllRIs==RIs[E],"score"]))/0
									delta = mean(as.numeric(tmpPat[AllRIs==RIs[E],"score"]))
									moyC = 0
									
									pvalCorr = try(t.test(x=rep(0, length(Ctrl)), mu=mean(as.numeric(tmpPat[AllRIs==RIs[E],"ScoreCorrige"])))$p.value, silent=TRUE)
									if(class(pvalCorr) == "try-error")	pvalCorr = 1
									ratioCorr = mean(as.numeric(tmpPat[AllRIs==RIs[E],"ScoreCorrige"]))/0
									deltaCorr = mean(as.numeric(tmpPat[AllRIs==RIs[E],"ScoreCorrige"]))
									moyCcorr = 0
									
									pvalFish = "NA"
#									new = cbind(tmpPat[AllRIs==RIs[E],,drop=FALSE], matrix(rep(c(pvalFish, pval, moyC, length(Ctrl), ratio, delta, pvalCorr, moyCcorr, ratioCorr, deltaCorr, "ev_Spe_P"), sum(AllRIs==RIs[E])), ncol=11, byrow=TRUE))
#									resRI_f = rbind(resRI_f, new)
								}
							}
						}
					}
#					#	filtre score > 0.05
#					resRI_f = resRI_f[as.numeric(resRI_f[,"score"])>=0.05,,drop=FALSE]
#					
#					#	filtre les types d'evennement OK=> c("RIG_5p", "A5SS", "A3SS", "RI_aval", "RI_amont", "SE+RI_3p", "RI_5p+SE")
#					evOK = c("RIG_5p", "A5SS", "A3SS", "RI_aval", "RI_amont", "SE+RI_3p", "RI_5p+SE")
#					resRI_f = resRI_f[resRI_f[,"Type"]%in%evOK,,drop=FALSE]
#					
#					if(nrow(resRI_f)>0)	write.table(resRI_f, file=paste(pathEv, "/resRI_f.txt", sep=""), sep="\t", quote=FALSE, row.names=FALSE)
				}
				
				if(sum(grepl("_SE.txt", AllRes))>0)
				{
					All_SE = as.matrix(read.table(file=paste(pathEv, "/resAll_SE.txt", sep=""), header=TRUE, check.names=FALSE, sep="\t", quote=""))
					All_SE[,"Sample"] = gsub(paste(gene, "_", Chr, "_", sep=""), "", All_SE[,"Sample"])

					ScoreCorr = apply(All_SE, 1, function(x) as.numeric(x[grep("score", colnames(All_SE))])*(mean(as.numeric(TotCountJuncs_All))/as.numeric(TotCountJuncs_All[x[grep("Sample", colnames(All_SE))]])))
					All_SE = cbind(All_SE, apply(All_SE, 1, function(x) sum(as.numeric(AllCounts[,x[grep("Sample", colnames(All_SE))]]))))
					colnames(All_SE)[ncol(All_SE)] = "TotCountJuncs"
					All_SE = cbind(All_SE, ScoreCorr)
					colnames(All_SE)[ncol(All_SE)] = "ScoreCorrige"
					All_SE[is.na(All_SE[,"ScoreCorrige"]), "ScoreCorrige"] = 0	#	correction 2022/08/18 remplace les NA par un score nul
					
					SamplesTypes = as.matrix(read.table(pathtFileSamplesTypes, header=TRUE, sep="\t", quote=""))
					SamplesTypes[,"Sample"] = apply(SamplesTypes, 1, function(x) paste(x[1], "_", x[3], sep=""))
					SampleNames = unique(All_SE[,"Sample"])

					Pat = SampleNames[SampleNames%in%SamplesTypes[grepl("Pat", SamplesTypes[,"Type"]), "Sample"]]
					Ctrl = SampleNames[SampleNames%in%SamplesTypes[grepl("Ctrl", SamplesTypes[,"Type"]), "Sample"]]
					
					
					Pat = Pat[Pat%in%unique(All_SE[,"Sample"])]	#	selectionne en fonction des res
					Ctrl = Ctrl[Ctrl%in%unique(All_SE[,"Sample"])]	#	selectionne en fonction des res
					All_Ctrl = All_SE[All_SE[,"Sample"]%in%Ctrl,,drop=FALSE]
					
					colJunc_SE_Start = grepl("Junc_SE_Start", colnames(All_Ctrl))
					colJunc_SE_End = grepl("Junc_SE_End", colnames(All_Ctrl))
					Ctrls = apply(All_Ctrl, 1, function(x) paste(as.numeric(x[colJunc_SE_Start]), "-", as.numeric(x[colJunc_SE_End]), sep=""))
					
#					resSE_f = matrix(0, ncol=(ncol(All_SE)+11), nrow=0)	#	Matrice de resultats des SE filtrés
#					colnames(resSE_f) = c(colnames(All_SE),"pvalFisher",  "pval", "MoyScoreCtrls", "Nctrls", "ratioScores(P_vs_Ctrl)", "deltaScores(P_vs_Ctrl)", "pvalCorr", "MoyScoreCtrlsCorr", "ratioScores(P_vs_Ctrl)Corr", "deltaScores(P_vs_Ctrl)Corr", "Com_ou_Spe")
					
					if(length(Pat)>0)	for(P in 1:length(Pat))	#	Pat par pat
					{
						tmpPat = All_SE[All_SE[,"Sample"]==Pat[P],,drop=FALSE]
						tmpPat = tmpPat[as.numeric(tmpPat[,"Junc_SE_Count"])>=2,,drop=FALSE]	#	Filtre les evenements avec read = 1
						PcolJunc_SE_Start = grep("Junc_SE_Start", colnames(tmpPat))
						PcolJunc_SE_End = grep("Junc_SE_End", colnames(tmpPat))
						
						if(nrow(tmpPat)>=1)
						{
							AllSEs = apply(tmpPat, 1, function(x) paste(as.numeric(x[PcolJunc_SE_Start]), "-", as.numeric(x[PcolJunc_SE_End]), sep=""))
							SEs = unique(AllSEs)
							SEs = SEs[!grepl("0-0", SEs)]
							
							for(E in 1:length(SEs))	#	J310 E=45
							{
								if(sum(Ctrls==SEs[E])>=1)	#	Evenement retrouve chez les ctrls ################## comparer les noms J plutot que les coords
								{	
									evCtrl = All_Ctrl[Ctrls==SEs[E],, drop=FALSE]
									pval = try(t.test(x=as.numeric(evCtrl[,"score"]), mu=mean(as.numeric(tmpPat[AllSEs==SEs[E],"score"])))$p.value, silent=TRUE)
									if(class(pval) == "try-error")	pval = 1
									ratio = mean(as.numeric(tmpPat[AllSEs==SEs[E],"score"]))/mean(as.numeric(evCtrl[,"score"]))
									delta = abs(mean(as.numeric(tmpPat[AllSEs==SEs[E],"score"]))-mean(as.numeric(evCtrl[,"score"])))
									moyC = mean(as.numeric(evCtrl[,"score"]))
									
									pvalCorr = try(t.test(x=as.numeric(evCtrl[,"ScoreCorrige"]), mu=mean(as.numeric(tmpPat[AllSEs==SEs[E],"ScoreCorrige"])))$p.value, silent=TRUE)
									if(class(pvalCorr) == "try-error")	pvalCorr = 1
									ratioCorr = mean(as.numeric(tmpPat[AllSEs==SEs[E],"ScoreCorrige"]))/mean(as.numeric(evCtrl[,"ScoreCorrige"]))
									deltaCorr = abs(mean(as.numeric(tmpPat[AllSEs==SEs[E],"ScoreCorrige"]))-mean(as.numeric(evCtrl[,"ScoreCorrige"])))
									moyCcorr = mean(as.numeric(evCtrl[,"ScoreCorrige"]))
									
									####################################################################################
#	test proportions Pat vs Ctrl
									
									RIpat = as.numeric(tmpPat[E,"Junc_SE_Count"])
									RIctrl = sum(sum(as.numeric(evCtrl[,"Junc_SE_Count"])))
									NormPat = as.numeric(tmpPat[E,"Junc_Normale_Count"])
									NormCtrl = sum(as.numeric(evCtrl[,"Junc_Normale_Count"]))
									
									tab = matrix(c(NormPat, NormCtrl, RIpat, RIctrl), ncol=2)
									colnames(tab) = c("Norm", "SE")
									rownames(tab) = c("pat", "Ctrl")
									tab[is.na(tab)] = 0
									#resChi = chisq.test(tab, correct=TRUE)
									resFish = fisher.test(tab, conf.int=TRUE, conf.level=0.95)
									pvalFish = resFish$p.value
									
									####################################################################################
									
#									new = cbind(tmpPat[AllSEs==SEs[E],,drop=FALSE], matrix(rep(c(pvalFish, pval, moyC, nrow(evCtrl), ratio, delta, pvalCorr, moyCcorr, ratioCorr, deltaCorr, "ev_Com_C_P"), sum(AllSEs==SEs[E])), ncol=11, byrow=TRUE))
#									resSE_f = rbind(resSE_f, new)
								}else{					#	Evenement specifique du patient
									pval = try(t.test(x=rep(0, length(Ctrl)), mu=mean(as.numeric(tmpPat[AllSEs==SEs[E],"score"])))$p.value, silent=TRUE)
									if(class(pval) == "try-error")	pval = 1
									ratio = mean(as.numeric(tmpPat[AllSEs==SEs[E],"score"]))/0
									delta = mean(as.numeric(tmpPat[AllSEs==SEs[E],"score"]))
									moyC = 0
									
									pvalCorr = try(t.test(x=rep(0, length(Ctrl)), mu=mean(as.numeric(tmpPat[AllSEs==SEs[E],"ScoreCorrige"])))$p.value, silent=TRUE)
									if(class(pvalCorr) == "try-error")	pvalCorr = 1
									ratioCorr = mean(as.numeric(tmpPat[AllSEs==SEs[E],"ScoreCorrige"]))/0
									deltaCorr = mean(as.numeric(tmpPat[AllSEs==SEs[E],"ScoreCorrige"]))
									moyCcorr = 0
									pvalFish = "NA"
									
#									new = cbind(tmpPat[AllSEs==SEs[E],,drop=FALSE], matrix(rep(c(pvalFish, pval, moyC, length(Ctrl), ratio, delta, pvalCorr, moyCcorr, ratioCorr, deltaCorr, "ev_Spe_P"), sum(AllSEs==SEs[E])), ncol=11, byrow=TRUE))
#									resSE_f = rbind(resSE_f, new)
								}
							}
						}
					}
#					resSE_f = resSE_f[as.numeric(resSE_f[,"score"])>=0.05,,drop=FALSE]	#	filtre score > 0.05
#					resSE_f = resSE_f[!grepl("SE_Canonique", resSE_f[,"Type"]),,drop=FALSE]	#	 filtre les SE canoniques
#					
#					if(nrow(resSE_f)>0)	write.table(resSE_f, file=paste(pathEv, "/resSE_f.txt", sep=""), sep="\t", quote=FALSE, row.names=FALSE)
				}
			}else{
				write("!!! Fichier SamplesTypes.txt absent => pas de comparaisons  !!!", file="")
				write(paste(gene, "\tFichier SamplesTypes.txt absent => pas de comparaisons", sep=""), file=paste(basePath, "ListRes.txt", sep=""), append = TRUE)
			}
			pathURLzips = paste(basePath, "/zipRes/", sep="")
			if(!file.exists(pathURLzips))	dir.create(pathURLzips)
			saveWD = getwd()
			setwd(pathResGenes)
			system(paste("zip -r ", pathURLzips, basename(pathRes), ".zip ./", basename(pathRes), "/", basename(pathEv), sep=""))	#	Zip des res
			setwd(saveWD)
			
			write(paste(gene, "\tOK", sep=""), file=paste(basePath, "ListRes.txt", sep=""), append = TRUE)
		}else{
			write("!!! Aucun evennement detecte !!!", file="")
		}
	}else{
		write("!!! Pas de sample analysable !!!", file="")
		write(paste(gene, "\tNiet Pas de sample analysable", sep=""), file=paste(basePath, "ListRes.txt", sep=""), append = TRUE)
		unlink(rmdupPath, recursive=TRUE)
		unlink(pathPosGene, recursive=TRUE)
	}
}else{
	write("!!! Un seul exon !!!", file="")
	write(paste(gene, "\tNiet: Exon unique", sep=""), file=paste(basePath, "ListRes.txt", sep=""), append = TRUE)
}






