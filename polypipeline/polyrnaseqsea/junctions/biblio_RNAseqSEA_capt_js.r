
bibGCMaj<-function()  # test si les genes de la capture ENSg sont dans le gencode declare
{
  allGC = list.files("/Users/nicolas/Documents/RNAseqSEA/Refs/", recursive=FALSE, full.names=TRUE)
  allGC = allGC[!grepl("Junc", basename(allGC))]
  
  bibGC = matrix(0, ncol=4, nrow=0)
  for(N in 1:length(allGC))
  {
    write(paste("#\t", N, "/", length(allGC), " - ", gsub(".rds", "", basename(allGC[N])), sep=""), file="")
    tmp = readRDS(allGC[N])
    GCvers = unlist(strsplit(gsub(".rds", "", basename(allGC[N])), split ="_gencode"))
    bibGC = rbind(bibGC, cbind(tmp[,c("ensembl_gene_id", "external_gene_name")], rep(GCvers[1], nrow(tmp)), rep(GCvers[2], nrow(tmp))))
  }
  colnames(bibGC) = c("ensembl_gene_id", "external_gene_name", "Ref", "GCvers")
  saveRDS(bibGC, "/data-isilon/Cagnard/RNAseqSEA/Refs/bibGC.rds")
}

identifGC<-function(resPath)
{
  genes = as.matrix(read.table(paste(resPath, "/RNAseqSEA/ENSg.txt", sep=""), sep="\t"))
  bibGC = readRDS("/data-isilon/Cagnard/RNAseqSEA/Refs/bibGC.rds")
  
  genesOK = bibGC[bibGC[,"ensembl_gene_id"]%in%genes,,drop=FALSE]
  table(genesOK[,"GCvers"])
}

testENSgGC<-function(gcrds, resPath)  # teste les ID, met a jour si possible, sinon epure et fait une liste de rejet
{
  ENSGpath = list.files(paste(resPath, "/RNAseqSEA/", sep=""), pattern="^ENSg*.*.txt$", full.names=TRUE)
  genes = as.matrix(read.table(ENSGpath, sep="\t"))
  GCref = readRDS(gcrds)
  
  # test si ENS G ou T
  IDs_T2G = NULL
  IDs_T = NULL

  if(sum(GCref[,"ensembl_transcript_id"]%in%genes)>0) #  certains ENS sont des transcrits pas des genes !
  {
    IDs_T = genes[genes%in%GCref[,"ensembl_transcript_id"],,drop=FALSE] #  ID T ok pour la version GC
    write(paste("\t !!! ", length(IDs_T), " identifiants sont des ENS ID de transcrits !!!", sep=""), file="")
    write(paste("\t\t => Conversion en ENS genes ID", sep=""), file="")
    # conversion en genes
    IDs_T2G = unique(GCref[GCref[,"ensembl_transcript_id"]%in%IDs_T, "ensembl_gene_id"])
  }
  
  if(sum(GCref[,"ensembl_gene_id"]%in%genes)>0)
  {
    IDs_ok = genes[genes%in%GCref[,"ensembl_gene_id"],,drop=FALSE] #ID genes ok pour la version GC
  }
  
  ID_Not_ok = genes[!(genes%in%c(IDs_T, IDs_ok))]
  CorrIDs = NULL
  
  if((length(ID_Not_ok)>0))
  {
    write(paste("\t !!! ", length(ID_Not_ok), " identifiants absents de la référence Gencode ", esp, " v", gcvers, " !!!", sep=""), file="")
    write.table(ID_Not_ok, file=paste(resPath, "/RNAseqSEA/NotInRefGC_(", length(ID_Not_ok), ").txt", sep=""), sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
    file.rename(ENSGpath, paste(dirname(ENSGpath), "/_", basename(ENSGpath), sep=""))
  }
  
  if(length(ID_Not_ok)>0)
  {
    IDs_ok = unique(c(IDs_T2G, IDs_ok))
    
    if(exists(AnnotPath))
    {
     Annots = readRDS(AnnotPath)
     write(paste("\t\t => Tentative de recuperation des ENS genes ID par conversion des IDs en noms de genes", sep=""), file="")
     
     noms = unique(Annots[Annots[,"EnsID"]%in%ID_Not_ok,"Symbol"])
     CorrIDs = unique(unlist(apply(as.matrix(noms), 1, function(x) unique(GCref[GCref[,"external_gene_name"]%in%x, "ensembl_gene_id"]))))
     write(paste("\t\t => ", length(CorrIDs), " IDs recuperes", sep=""), file="")
    }
    IDs_ok = unique(c(IDs_ok, CorrIDs))
    write(paste("\t\t => ", length(IDs_ok), "/", length(genes), " IDs ok", sep=""), file="")
  }
  write.table(IDs_ok, file=paste(resPath, "/RNAseqSEA/ENSg.txt", sep=""), sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
}

ENSg<-function(idprojet, esp, biblioPath, resPath)
{
	source(biblioPath)
	filePath = paste(resPath, "/RNAseqSEA/", sep="")
	if(!file.exists(filePath))	dir.create(filePath)
	captDBB = system(paste(genelistPath, " -project=", idprojet, sep=""), intern=TRUE)
	tmp = do.call(rbind, apply(as.matrix(captDBB), 1, function(x) strsplit(x, "\t")))
	tmp = do.call(rbind, tmp)
	ENSGcaptDBB = unique(tmp[,2])

	ENSGcaptDBB = ENSGcaptDBB[ENSGcaptDBB%in%Annots[,"EnsID"]]
	
	write.table(ENSGcaptDBB, file=paste(resPath, "/RNAseqSEA/ENSg.txt", sep=""), row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
}

formatProjet<-function(resPath, align)
{
	formatResComp(resPath, align)
	formatResAll(resPath, align)
}

formatSamplesTypes_1vsAll<-function(idprojet, biblioPath, scriptPath, esp, gcvers, align, nCPUmax, limDecal, bamfiles, gcrds, resPath, sambambaPath, samtoolsPath, picardPath)
{
	filePath = paste(resPath, "/RNAseqSEA/", sep="")
	if(!file.exists(filePath))	dir.create(filePath)
	
	path_cmds_global = paste(filePath, "cmds_all.sh", sep="");
	unlink(path_cmds_global)
	
	if(length(bamfiles)>0)
	{
		cmds = matrix("", ncol=1, nrow=length(bamfiles))
		gabSamplesTypes = matrix(0, ncol=3, nrow=0)
		colnames(gabSamplesTypes) = c("Sample", "Type", "Proj")
		
		pathtFileENSg = list.files(paste(resPath, "/RNAseqSEA/", sep=""), full.names=TRUE)
		pathtFileENSg = pathtFileENSg[grepl("^ENSg.txt", basename(pathtFileENSg))]
		
		
		
		for(B in 1:length(bamfiles))
		{
			SamplesTypes = rbind(gabSamplesTypes, c(gsub(".bam", "", basename(bamfiles[B])), "Pat", idprojet))
			SamplesTypes = rbind(SamplesTypes, cbind(gsub(".bam", "", basename(bamfiles[!(c(1:length(bamfiles))%in%B)])), rep("Ctrl", (length(bamfiles)-1)), rep(idprojet, (length(bamfiles)-1))))
			write.table(SamplesTypes, file=paste(filePath, "SamplesTypes_", gsub(".bam", "", basename(bamfiles[B])), "_vs_all.txt", sep=""), row.names=FALSE, quote=FALSE, sep="\t")
			cmds[B] = paste("for i in `cat ", pathtFileENSg, "`;do echo Rscript ", scriptPath, " keepBams=TRUE gene=$i idprojet=", idprojet, " bamsProj=", idprojet, " titre=", gsub(".bam", "", basename(bamfiles[B])), "_vs_all esp=", esp, " gcvers=", gcvers, " align=", align, " nCPUmax=", nCPUmax, " limDecal=", limDecal, " biblioPath=", biblioPath, " gcrds=", gcrds, " resPath=", resPath, " sambambaPath=", sambambaPath, " samtoolsPath=", samtoolsPath, " picardPath=", picardPath, " >>" , path_cmds_global , ";done;", sep="")
		}
		write.table(cmds, file=paste(filePath, "/cmds.sh", sep=""), row.names=FALSE, col.names=FALSE, quote=FALSE, sep='\t')
		write(paste("\n\n\tLe formatage des comparaisons est ok, pour lancer les analyses:\n\t\tSe placer ici: ", resPath, "/\n\t\tlancer la commande: source \"./RNAseqSEA/cmds.sh\"", sep=""), file="")
	}else{
		write(paste("\n\n\t Les fichiers bams sont introuvables à cet endroit:\n\t\t ", resPath, "/align/", align, "/", sep=""), file="")
	}
}

formatSamplesTypes_1vsAll_bySampleGene<-function(idprojet, biblioPath, scriptPath, vers, RNAseqSEApath, esp, gcvers, align, nCPUmax, limDecal, bamfiles)
{
  filePath = paste(resPath, "/RNAseqSEA/", sep="")
  if(!file.exists(filePath))	dir.create(filePath)
	
	if(length(bamfiles)>0)
	{
		cmds = matrix("", ncol=1, nrow=0)
		gabSamplesTypes = matrix(0, ncol=3, nrow=0)
		colnames(gabSamplesTypes) = c("Sample", "Type", "Proj")
		
		pathtFileENSg = list.files(paste(resPath, "/RNAseqSEA/", sep=""), full.names=TRUE)
		pathtFileENSg = pathtFileENSg[grepl("^ENSg.txt", basename(pathtFileENSg))]
		
		listENSg = as.matrix(read.table(pathtFileENSg, sep="\t"))
		
		for(B in 1:length(bamfiles))
		{
			SamplesTypes = rbind(gabSamplesTypes, c(gsub(".bam", "", basename(bamfiles[B])), "Pat", idprojet))
			SamplesTypes = rbind(SamplesTypes, cbind(gsub(".bam", "", basename(bamfiles[!(c(1:length(bamfiles))%in%B)])), rep("Ctrl", (length(bamfiles)-1)), rep(idprojet, (length(bamfiles)-1))))
			write.table(SamplesTypes, file=paste(filePath, "SamplesTypes_", gsub(".bam", "", basename(bamfiles[B])), "_vs_all.txt", sep=""), row.names=FALSE, quote=FALSE, sep="\t")
			
			for(G in 1:nrow(listENSg))
			{
				cmds = rbind(cmds, matrix("", ncol=1, nrow=1))
				cmds[nrow(cmds), 1] = paste("/software/bin/run_cluster.pl -cpu=", nCPUmax, " -cmd=\"Rscript ", scriptPath, " keepBams=TRUE gene=", listENSg[G], " idprojet=", idprojet, " bamsProj=", idprojet, " titre=", gsub(".bam", "", basename(bamfiles[B])), "_vs_all esp=", esp, " align=", align, " nCPUmax=", nCPUmax, " limDecal=", limDecal, "biblioPath=", biblioPath, "\"", sep="")
			}
		}
		write.table(cmds, file=paste(filePath, "/cmds_bySampleGene.txt", sep=""), row.names=FALSE, col.names=FALSE, quote=FALSE, sep='\t')
		
		cmd = paste("mac2unix ", filePath, "cmds_bySampleGene.txt; dos2unix ", filePath, "cmds_bySampleGene.txt", sep="")
		system(paste(cmd, sep=""))
		
		write(paste("\n\n\tLe formatage des comparaisons est ok, pour lancer les analyses:\n\t\tSe placer ici: ", resPath, "/\n\t\tlancer la commande: source ", filePath, "cmds_bySampleGene.txt", sep=""), file="")
	}else{
		write(paste("\n\n\t Les fichiers bams sont introuvables à cet endroit:\n\t\t ", resPath, "/align/", align, "/", sep=""), file="")
	}
}

formatGTF<-function(gtfPath)
{
	library(rtracklayer)
	gtf <- rtracklayer::import(gtfPath)
	gtf_df = as.data.frame(gtf)
	gtf_exons  = gtf_df[gtf_df[,"type"] == "exon",]
	gtf_exons = gtf_exons[,c("seqnames", "strand", "gene_id", "transcript_id", "exon_id", "start", "end", "transcript_name", "level")]
	colnames(gtf_exons) = c("chromosome_name", "strand", "ensembl_gene_id", "ensembl_transcript_id", "ensembl_exon_id", "exon_chrom_start", "exon_chrom_end", "external_gene_name", "level")
	gtf_exons[,1] = gsub("chr", "", gtf_exons[,"chromosome_name"])
	gtf_exons[,"ensembl_gene_id"] = gsub("[.].*", "", gtf_exons[,"ensembl_gene_id"])
	gtf_exons[,"ensembl_transcript_id"] = gsub("[.].*", "", gtf_exons[,"ensembl_transcript_id"])
	gtf_exons[,"ensembl_exon_id"] = gsub("[.].*", "", gtf_exons[,"ensembl_exon_id"])

	gtf_exons[,"level"] = as.numeric(gtf_exons[,"level"])
	gtf_exons[,"exon_chrom_start"] = as.numeric(gtf_exons[,"exon_chrom_start"])
	gtf_exons[,"exon_chrom_end"] = as.numeric(gtf_exons[,"exon_chrom_end"])
	
	esp = gsub("/igv/.*", "", gsub("/data-isilon/public-data/", "", gtfPath))
	vers = gsub("gencode.", "", gsub(".gtf.gz", "", basename(gtfPath)))
	saveRDS(gtf_exons, paste("/data-isilon/Cagnard/RNAseqSEA/Refs/", esp, "_gencode", vers, ".rds", sep=""))
}

geneBams<-function(B)	#	Cree un bam pour les reads du gene
{
	write(paste("\t#\tSelection des reads du gene ", gene, sep=""), file="")
	
	bamProj = unlist(strsplit(dirname(bamfiles[B]), "/"))
	bamProj = bamProj[grepl("NGS", bamProj)]
	if(file.exists(paste(resPath, "./Echs.txt", sep="")))
	{
		newName = Echs[grepl(paste("^", gsub(".bam$", "", basename(bamfiles[B])), "$", sep=""), Echs[,1]),2]
		nomBam = paste(newName, "_", bamProj, ".bam", sep="")
	}else{
		nomBam = paste(gsub(".bam$", "", basename(bamfiles[B])), "_", bamProj, ".bam", sep="")
	}
	system(paste(samtoolsPath, " view -b ", bamfiles[B], " ", CHRbam, Chr, ":", Start-100, "-", End+100, " > ", pathBamsGene, "/", nomBam, sep=""))
	system(paste("java -jar ", picardPath, " BuildBamIndex -I ", pathBamsGene, "/", nomBam, sep=""))
}

rmdup<-function(I)
{
	write(paste("\t#\trmdup ", bamsGene[I], sep=""), file="")
	system(paste(sambambaPath, " markdup -r -t 4 ", bamsGene[I], " ", rmdupPath, substr(basename(bamsGene[I]), 0, (nchar(basename(bamsGene[I]))-4)), "_rmdup.bam", sep=""))
}

exJunBed<-function(J)	#	nom chromosome avec chr
{
	tmp = system(paste(sambambaPath, " view ",bamsRMdups[J], " ", CHRbam, Chr, ":", Start-100, "_", End+100, " | cut -f3,4,6  | awk '{if ($3 ~ /N/){printf(\"%s\",$1); start=$2; end=$2; split($3,a,\"[NIMDSH]\"); 
							split($3,b,\"[0-9]*\"); nb=length(b); for(i=2; i<=nb; i++){if(b[i] ~ /[MD]/){ end=end+a[i-1];} if(b[i] ~ /N/){end=end-1; printf(\"\t%s\t%s\",start, end); 
							start=end+a[i-1]+1; end=start;}} printf(\"\t%s\t%s\\n\", start,end-1);}}' > ", bedpath, gsub("_rmdup.bam", "", basename(bamsRMdups[J])), "_exons_junction_", Chr, ".mbed", sep=""), intern=TRUE)
}

JuncBam<-function(J)
{
	#	test du bed pour eviter les bams vides...
	if(file.size(paste(bedpath,gsub("_rmdup.bam", "", basename(bamsRMdups[J])), "_exons_junction_", Chr, ".mbed", sep="")) != 0L)
	{
		system(paste(sambambaPath, " view ",bamsRMdups[J], " | cut -f3,4,6  | awk '($3 ~ /N/)' > ", juncBamPath, gsub("_rmdup.bam", "", basename(bamsRMdups[J])), "_exons_junction_", Chr, ".bam", sep=""), intern=TRUE)
	}
}

filtrePairJunc<-function(x)
{
	x = x[!is.na(x)]
	x = as.numeric(x[3:(length(x)-1)])
	resPairs = matrix(x, ncol=2, byrow=TRUE)
	return(resPairs)
}

SelectJun<-function(S)	# S = beds[1]
{
	P = gsub("_exons_junction_.*", "", basename(S))	#	Patient
	write(paste("#\t Comptage des jonctions: ", P, sep=""), file="")
	
	tmp <- try(as.matrix(read.table(S, sep="\t", quote="", fill=TRUE, col.names = paste0("Col",seq_len(20)))))
	
	if(nrow(tmp)>1)
	{
		PairPos = cbind(apply(tmp, 1, function(x) filtrePairJunc(x)))
		if(is.list(PairPos[1]))
		{
			PairPos = do.call(rbind,PairPos)
		}else{
			PairPos = t(PairPos)
		}
		
		JunCounts = table(apply(PairPos, 1, function(x) paste(x, collapse= "_")))
		tmpC = t(apply(as.matrix(names(JunCounts)), 1, function(x) unlist(strsplit(x, "_"))))
		coords = cbind(as.numeric(tmpC[,1]), as.numeric(tmpC[,2]), as.numeric(JunCounts))
		coords = coords[order(coords[,1]),,drop=FALSE]
		colnames(coords) = c("StartJun", "EndJun", "CountJun")

		#	Selection des Junctions
		geneJun = coords[((coords[,"StartJun"]>=(Start-100))&(coords[,"StartJun"]<=(End+100)))|((coords[,"EndJun"]>=(Start-100))&(coords[,"EndJun"]<=(End+100))),,drop=FALSE]
		if(nrow(geneJun)>0)
		{
			saveRDS(geneJun, paste(pathPosGene, "/juncPairPos_", gene, "_", Chr, "_", P, ".rds", sep=""))
		}else{
			write(paste("#\t Aucun read pour ", P, sep=""), file="")
		}
	}else{
		write(paste("#\t Aucun read pour ", P, sep=""), file="")
	}
}

slogRNAseqSEA_captjs<-function(nfos)
{
  noms = c("idprojet", "esp", "align", "nodename", "login", "user", "effective_user", "wd")
  allnfos = matrix(0, ncol=length(noms), nrow=0)
  colnames(allnfos) = noms
  rdsPath = "/data-isilon/Cagnard/RNAseqSEA/RNAseqSEA_captjs_logs.rds"
  if(file.exists(rdsPath)) allnfos = readRDS(rdsPath)
  allnfos = rbind(allnfos, nfos)
  saveRDS(allnfos, rdsPath)
}

patRes<-function()
{	
	listProjs = list.files(projPath, pattern = paste("^RNAseqSEA_v", vers, "_", sep=""), full.names=TRUE)
	
	for(P in 1:length(listProjs))	#	Chaque projet
	{
		tmpwd = getwd()
		basePath = listProjs[P]
		write(paste("\n\t# projet ", basename(listProjs[P]), ": Formatage des resultst par patient", sep=""), file="")
		setwd(basePath)	#	Se place dans le projet
		SamplesResAll(basePath)
		setwd(tmpwd)
		#	liste des liens 
		fichs = list.files(basePath, full.names=TRUE)
		fichs = fichs[grepl(".zip", basename(fichs))]
		urls = paste("http://www.polyweb.fr/NGS/", idprojet, "/", esp, "/analisys/RNAseqSEA", if(!is.null(titre)){paste("_", titre, sep="")}, "/", basename(fichs), sep="")
		write.table(urls, file=paste(basePath, "/urls_", idprojet, if(!is.null(titre)){paste("_", titre, sep="")}, ".txt", sep=""), row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE)
	}
}

#	resAll
#	tableNFOjunc concatené et filtré de toute sles jonctions de tous les gènes avec annotations et table de comptages
#	Allres_SE, AllresRI, et filtrés pour chaque comparaison

formatResComp<-function(resPath, alignMethod)
{
  projPath = paste(resPath, "/analysis/", sep="")
  listProjs = list.files(projPath, pattern = paste("^RNAseqSEA_", sep=""), full.names=TRUE)
	
	for(P in 1:length(listProjs))	#	Chaque projet/comparaison
	{
		write(paste("# ", P, "/", length(listProjs), " ", basename(listProjs[P]), sep=""), file="")
		fichs = list.files(listProjs[P], full.names=TRUE, recursive=TRUE)
		fichs = fichs[grepl("SpliceRes", fichs)]
		fichs_NFOjunc = fichs[grepl("tableNFOjunc.txt$", fichs)]
		fichs_resAll_SE = fichs[grepl("resAll_SE.txt$", fichs)]
		fichs_resAll_RI = fichs[grepl("resAll_RI.txt$", fichs)]
#		fichs_resAll_SE_f = fichs[grepl("resSE_f.txt$", fichs)]
#		fichs_resAll_RI_f = fichs[grepl("resRI_f.txt$", fichs)]
		
		
		
		patName = basename(listProjs[P])
		patName = sub('RNAseqSEA_', '', patName)
		patName = sub('_vs_all', '', patName)
		
		resAllpath = paste(listProjs[P], "/AllRes/", sep="")
		if(!file.exists(resAllpath)) dir.create(resAllpath)
		pathJunctions = paste(resPath, "/junctions/", sep="")
		if(!file.exists(pathJunctions)) dir.create(pathJunctions)
		pathJunctionsRnaSeqSea = paste(pathJunctions, "/rnaseqsea/", sep="")
		if(!file.exists(pathJunctionsRnaSeqSea)) dir.create(pathJunctionsRnaSeqSea)
		pathJunctionsAlign = paste(pathJunctionsRnaSeqSea, "/", alignMethod, "/", sep="")
		if(!file.exists(pathJunctionsAlign)) dir.create(pathJunctionsAlign)

		if(length(fichs_NFOjunc)>0)
		{
			allNFOjunc = list()
			for(N in 1:length(fichs_NFOjunc))	allNFOjunc[[length(allNFOjunc)+1]] = as.matrix(read.table(fichs_NFOjunc[N], sep="\t", header=TRUE))
			allNFOjunc = do.call("rbind",allNFOjunc)
			write.table(allNFOjunc, file=paste(resAllpath, "allNFOjunc.txt", sep=""), sep="\t", quote=FALSE, row.names=FALSE)
		}

		if(length(fichs_resAll_SE)>0)
		{
			allResSE = list()
			for(N in 1:length(fichs_resAll_SE))	allResSE[[length(allResSE)+1]] = as.matrix(read.table(fichs_resAll_SE[N], sep="\t", header=TRUE))
			allResSE = do.call("rbind",allResSE)
			colnames(allResSE)[1] = paste("#", colnames(allResSE)[1], sep="")
			uniques = apply(allResSE[,c("Chr", "Junc_SE_Start", "Junc_SE_End", "Sample")], 1, function(x) paste(x, sep="_", collapse="_"))
			allResSE = allResSE[!duplicated(uniques),,drop=FALSE]
			write.table(allResSE, file=paste(resAllpath, "allResSE.txt", sep=""), sep="\t", quote=FALSE, row.names=FALSE)
			
			file1 <- paste(resAllpath, "allResSE.txt", sep="")
			file2 <- paste(pathJunctionsAlign, patName, "_SE.txt.gz", sep="")
			cmd1 <- paste(c('(grep "^#" ', file1, '; grep -v "^#" ', file1), collapse="")
			cmd2 <- paste(c('sort -t"`printf \'\t\'`" -k5,5 -k6,6n)'), collapse="")
			cmd3 <- paste(c('bgzip > ', file2,';'), collapse="")
			cmd_global <- paste(c( cmd1, cmd2, cmd3), collapse=" | ")
			system(cmd_global)
			cmd_tabix <- paste(c('tabix -s 5 -b 6 ', file2), collapse="")
			system(cmd_tabix)
		}

		if(length(fichs_resAll_RI)>0)
		{
			allResRI = list()
			for(N in 1:length(fichs_resAll_RI))	allResRI[[length(allResRI)+1]] = as.matrix(read.table(fichs_resAll_RI[N], sep="\t", header=TRUE))
			allResRI = do.call("rbind",allResRI)
			uniques = apply(allResRI[,c("Chr", "Junc_RI_Start", "Junc_RI_End", "Sample")], 1, function(x) paste(x, sep="_", collapse="_"))
			allResRI = allResRI[!duplicated(uniques),,drop=FALSE]
			colnames(allResRI)[1] = paste("#", colnames(allResRI)[1], sep="")
			write.table(allResRI, file=paste(resAllpath, "allResRI.txt", sep=""), sep="\t", quote=FALSE, row.names=FALSE)
			
			file_to_sort <- paste(resAllpath, "allResRI.txt", sep="")
			sort_bgzip_tabix_RI_file(file_to_sort);
			
			file1 <- paste(resAllpath, "allResRI.txt", sep="")
			file2 <- paste(pathJunctionsAlign, patName, "_RI.txt.gz", sep="")
			cmd1 <- paste(c('(grep "^#" ', file1, '; grep -v "^#" ', file1), collapse="")
			cmd2 <- paste(c('sort -t"`printf \'\t\'`" -k4,4 -k5,5n)'), collapse="")
			cmd3 <- paste(c('bgzip > ', file2,';'), collapse="")
			cmd_global <- paste(c( cmd1, cmd2, cmd3), collapse=" | ")
			system(cmd_global)
			cmd_tabix <- paste(c('tabix -s 4 -b 5 ', file2), collapse="")
			system(cmd_tabix)
		}

#		if(length(fichs_resAll_SE_f)>0)
#		{
#			allResSE_f = list()
#			for(N in 1:length(fichs_resAll_SE_f))	allResSE_f[[length(allResSE_f)+1]] = as.matrix(read.table(fichs_resAll_SE_f[N], sep="\t", header=TRUE))
#			allResSE_f = do.call("rbind",allResSE_f)
#			write.table(allResSE_f, file=paste(resAllpath, "allResSE_f.txt", sep=""), sep="\t", quote=FALSE, row.names=FALSE)
#		}
		
#		if(length(fichs_resAll_RI_f)>0)
#		{
#			allResRI_f = list()
#			for(N in 1:length(fichs_resAll_RI_f))	allResRI_f[[length(allResRI_f)+1]] = as.matrix(read.table(fichs_resAll_RI_f[N], sep="\t", header=TRUE))
#			allResRI_f = do.call("rbind",allResRI_f)
#			write.table(allResRI_f, file=paste(resAllpath, "allResRI_f.txt", sep=""), sep="\t", quote=FALSE, row.names=FALSE)
#		}

	}
}

sort_bgzip_tabix_SE_file<-function(file) {
	file_gz <- paste(file, ".gz", sep="")
	cmd1 <- paste(c('(grep "^#" ', file, '; grep -v "^#" ', file), collapse="")
	cmd2 <- paste(c('sort -t"`printf \'\t\'`" -k5,5 -k6,6n)'), collapse="")
	cmd3 <- paste(c('bgzip > ', file_gz,';'), collapse="")
	cmd_global <- paste(c( cmd1, cmd2, cmd3), collapse=" | ")
	system(cmd_global)
	cmd_tabix <- paste(c('tabix -s 5 -b 6 ', file_gz), collapse="")
	system(cmd_tabix)
}

sort_bgzip_tabix_RI_file<-function(file) {
	file_gz <- paste(file, ".gz", sep="")
	cmd1 <- paste(c('(grep "^#" ', file, '; grep -v "^#" ', file), collapse="")
	cmd2 <- paste(c('sort -t"`printf \'\t\'`" -k4,4 -k5,5n)'), collapse="")
	cmd3 <- paste(c('bgzip > ', file_gz,';'), collapse="")
	cmd_global <- paste(c( cmd1, cmd2, cmd3), collapse=" | ")
	system(cmd_global)
	cmd_tabix <- paste(c('tabix -s 4 -b 5 ', file_gz), collapse="")
	system(cmd_tabix)
}

formatResAll<-function(resPath, alignMethod)
{
	resAllpath = paste(resPath, "/junctions/rnaseqsea/", alignMethod, "/", sep="")
	fileout_RI = paste(c(resAllpath, 'allResRI.txt'), collapse="")
	fileout_SE = paste(c(resAllpath, 'allResSE.txt'), collapse="")
	fileout_RI_gz = paste(c(resAllpath, 'allResRI.txt.gz'), collapse="")
	fileout_SE_gz = paste(c(resAllpath, 'allResSE.txt.gz'), collapse="")
	cmd_zcat1 <- paste(c('zcat ', resAllpath, '*_RI.txt.gz | tail -n +2 >', fileout_RI), collapse="")
	system(cmd_zcat1)
	cmd_zcat2 <- paste(c('zcat ', resAllpath, '*_SE.txt.gz | tail -n +2 >', fileout_SE), collapse="")
	system(cmd_zcat2)
	sort_bgzip_tabix_RI_file(fileout_RI)
	sort_bgzip_tabix_SE_file(fileout_SE)
	pathAnalyse =  paste(resPath, "/analysis/", sep="")
	if(!file.exists(pathAnalyse)) dir.create(pathAnalyse)
	pathAnalyseAllRes =  paste(pathAnalyse, "/AllRes/", sep="")
	if(!file.exists(pathAnalyseAllRes)) dir.create(pathAnalyseAllRes)
	cmd_cp <- paste(c('cp ', fileout_RI_gz, '* ', fileout_SE_gz, '* ', pathAnalyseAllRes, '/.'), collapse="")
	system(cmd_cp)
	if(file.exists(fileout_RI)) { unlink(fileout_RI) }
	if(file.exists(fileout_SE)) { unlink(fileout_SE) }
}

formatResAll_OLD<-function(resPath)
{
  projPath = paste(resPath, "/analysis/", sep="")
        listProjs = list.files(projPath, pattern = paste("^RNAseqSEA_", sep=""), full.names=TRUE)
        resAllpath = paste(projPath, "/AllRes/", sep="")
        if(!file.exists(resAllpath)) dir.create(resAllpath)

        fichs = NULL
        for(P in 1:length(listProjs))   fichs = c(fichs, list.files(paste(listProjs[P], "/AllRes/", sep=""), full.names=TRUE, recursive=FALSE))

        fichs_NFOjunc = fichs[grepl("allNFOjunc.txt$", fichs)]
        fichs_resAll_SE = fichs[grepl("allResSE.txt$", fichs)]
        fichs_resAll_RI = fichs[grepl("allResRI.txt$", fichs)]
        fichs_resAll_SE_f = fichs[grepl("allResSE_f.txt$", fichs)]
        fichs_resAll_RI_f = fichs[grepl("allResRI_f.txt$", fichs)]

        if(length(fichs_NFOjunc)>0)
        {
                allNFOjunc = list()
                for(N in 1:length(fichs_NFOjunc))
                {
                        tmp = as.matrix(read.table(fichs_NFOjunc[N], sep="\t", header=TRUE))
                        tmp = cbind(tmp, rep(basename(dirname(dirname(listProjs[N]))), nrow(tmp)))
                        colnames(tmp)[ncol(tmp)] = "Comp"
                        allNFOjunc[[length(allNFOjunc)+1]] = tmp
                }
                tmp = do.call("rbind",lapply(allNFOjunc, ncol))

                if(length(table(tmp))==1)
                {
                        allNFOjunc = do.call("rbind",allNFOjunc)
                        indiceFiltre = apply(allNFOjunc, 1, function(x) paste(gsub(" ", "", x[3:5]), collapse="_", sep="_"))
                        allNFOjunc = allNFOjunc[duplicated(indiceFiltre),]
                        write.table(allNFOjunc, file=paste(resAllpath, "allNFOjunc.txt", sep=""), sep="\t", quote=FALSE, row.names=FALSE)
                        system(paste("bgzip -@ 40 ", resAllpath, "allNFOjunc.txt", sep=""))
                }else{
                        for(N in 1:length(fichs_NFOjunc))
                        {
                                tmpallNFOjunc = allNFOjunc[[N]]
                                indiceFiltre = apply(tmpallNFOjunc, 1, function(x) paste(gsub(" ", "", x[3:5]), collapse="_", sep="_"))
                                tmpallNFOjunc = tmpallNFOjunc[duplicated(indiceFiltre),]
                                write.table(tmpallNFOjunc, file=paste(resAllpath, "allNFOjunc_", gsub("/AllRes*.*", "", gsub(paste("*.*/RNAseqSEA_", sep=""), "", fichs_NFOjunc[N])), ".txt", sep=""), sep="\t", quote=FALSE, row.names=FALSE)
                        }
                }
        }
        
        if(length(fichs_resAll_RI)>0)
        {
                allResRI = list()
                for(N in 1:length(fichs_resAll_RI))
                {
                        tmp = as.matrix(read.table(fichs_resAll_RI[N], sep="\t", header=TRUE))
                        tmp = cbind(tmp, rep(basename(dirname(dirname(fichs_resAll_RI[N]))), nrow(tmp)))
                        colnames(tmp)[ncol(tmp)] = "Comp"
                        allResRI[[length(allResRI)+1]] = tmp
                }
                allResRI = do.call("rbind",allResRI)

                allResRI[,"Junc_RI_Start"] = as.numeric(allResRI[,"Junc_RI_Start"])
                allResRI[,"Junc_RI_End"] = as.numeric(allResRI[,"Junc_RI_End"])
                allResRI[,"Junc_RI_Count"] = as.numeric(allResRI[,"Junc_RI_Count"])
                allResRI[,"Junc_Normale_Count"] = as.numeric(allResRI[,"Junc_Normale_Count"])

                allResRI = allResRI[order(as.numeric(allResRI[,"Junc_RI_Start"])),]
                allResRI = allResRI[order(allResRI[,"Chr"]),]

                colnames(allResRI)[1] = paste("#", colnames(allResRI)[1], sep="")

                write.table(allResRI, file=paste(resAllpath, "allResRI.txt", sep=""), sep="\t", quote=FALSE, row.names=FALSE)
                system(paste("bgzip -@ 40 ", resAllpath, "allResRI.txt", sep=""))
                system(paste("tabix -S 1 -s ", grep("Chr", colnames(allResRI)), " -b ", grep("Junc_RI_Start", colnames(allResRI)), " -e ", grep("Junc_RI_End", colnames(allResRI)), " ", resAllpath, "allResRI.txt.gz", sep=""))
        }
}        
        