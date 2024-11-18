
#	/software/bin/run_cluster.pl -cpu=40 -cmd="Rscript /data-isilon/Cagnard/RNAseqSEA/RNAseqSEA_AllT/RNAseqSEA_AllT_formatRefs.r nCPU=20 esp=HG19 GCv=19 viaGTFgz=/data-isilon/public-data/HG19/igv/gencode.v19.annotation.gtf.gz"
#	/software/bin/run_cluster.pl -cpu=40 -cmd="Rscript /data-isilon/Cagnard/RNAseqSEA/RNAseqSEA_AllT/RNAseqSEA_AllT_formatRefs.r nCPU=20 esp=HG19 GCv=31 viaGTFgz=/data-isilon/public-data/HG19/igv/gencode.31.gtf.gz"
#	/software/bin/run_cluster.pl -cpu=40 -cmd="Rscript /data-isilon/Cagnard/RNAseqSEA/RNAseqSEA_AllT/RNAseqSEA_AllT_formatRefs.r nCPU=20 esp=HG19 GCv=34 viaGTFgz=/data-isilon/public-data/HG19/igv/gencode.34.gtf.gz"
#	/software/bin/run_cluster.pl -cpu=40 -cmd="Rscript /data-isilon/Cagnard/RNAseqSEA/RNAseqSEA_AllT/RNAseqSEA_AllT_formatRefs.r nCPU=20 esp=HG19 GCv=38 viaGTFgz=/data-isilon/public-data/HG19/igv/gencode.38.gtf.gz"
#	/software/bin/run_cluster.pl -cpu=40 -cmd="Rscript /data-isilon/Cagnard/RNAseqSEA/RNAseqSEA_AllT/RNAseqSEA_AllT_formatRefs.r nCPU=20 esp=HG19 GCv=39 viaGTFgz=/data-isilon/public-data/HG19/igv/gencode.39.gtf.gz"
#	/software/bin/run_cluster.pl -cpu=40 -cmd="Rscript /data-isilon/Cagnard/RNAseqSEA/RNAseqSEA_AllT/RNAseqSEA_AllT_formatRefs.r nCPU=20 esp=HG19 GCv=40 viaGTFgz=/data-isilon/public-data/HG19/igv/gencode.40.gtf.gz"
#	/software/bin/run_cluster.pl -cpu=40 -cmd="Rscript /data-isilon/Cagnard/RNAseqSEA/RNAseqSEA_AllT/RNAseqSEA_AllT_formatRefs.r nCPU=20 esp=HG19 GCv=41 viaGTFgz=/data-isilon/public-data/HG19/igv/gencode.41.gtf.gz"
#	/software/bin/run_cluster.pl -cpu=40 -cmd="Rscript /data-isilon/Cagnard/RNAseqSEA/RNAseqSEA_AllT/RNAseqSEA_AllT_formatRefs.r nCPU=20 esp=HG19 GCv=42 viaGTFgz=/data-isilon/public-data/HG19/igv/gencode.42.gtf.gz"
#	/software/bin/run_cluster.pl -cpu=40 -cmd="Rscript /data-isilon/Cagnard/RNAseqSEA/RNAseqSEA_AllT/RNAseqSEA_AllT_formatRefs.r nCPU=20 esp=HG19 GCv=43 viaGTFgz=/data-isilon/public-data/HG19/igv/gencode.43.gtf.gz"
#	/software/bin/run_cluster.pl -cpu=40 -cmd="Rscript /data-isilon/Cagnard/RNAseqSEA/RNAseqSEA_AllT/RNAseqSEA_AllT_formatRefs.r nCPU=20 esp=MM38 GCv=25 viaGTFgz=/data-isilon/public-data/MM38/igv/gencode.M25.gtf.gz"
#	/software/bin/run_cluster.pl -cpu=40 -cmd="Rscript /data-isilon/Cagnard/RNAseqSEA/RNAseqSEA_AllT/RNAseqSEA_AllT_formatRefs.r nCPU=20 esp=MM39 GCv=32 viaGTFgz=/data-isilon/public-data/MM39/igv/gencode.M32.gtf.gz"

#	/software/bin/run_cluster.pl -cpu=40 -cmd="Rscript /data-isilon/Cagnard/RNAseqSEA/RNAseqSEA_AllT/RNAseqSEA_AllT_formatRefs.r nCPU=40 esp=HG38 GCv=43 viaGTFgz=/data-isilon/Cagnard/RNAseqSEA/HG38_genecode43.gtf.gz"


#	esp = "HG19"; GCv = "43"; nCPU=20; viaGTFgz = "/data-isilon/public-data/HG19/igv/gencode.43.gtf.gz"

RNAseqSEApath = "/data-isilon/Cagnard/RNAseqSEA/"

args = commandArgs(trailingOnly=TRUE)
tmpArgs = t(rbind(apply(as.matrix(unlist(args)), 1, function(x) unlist(strsplit(x, "=")))))
if(sum(grepl("GCv", tmpArgs[,1]))>0)	GCv = tmpArgs[grep("GCv", tmpArgs[,1], ignore.case=TRUE),2]
if(sum(grepl("esp", tmpArgs[,1]))>0)	esp = tmpArgs[grep("esp", tmpArgs[,1], ignore.case=TRUE),2]
if(sum(grepl("viaGTFgz", tmpArgs[,1]))>0)	viaGTFgz = tmpArgs[grep("viaGTF", tmpArgs[,1], ignore.case=TRUE),2]
if(sum(grepl("nCPU", tmpArgs[,1]))>0)	nCPU = as.numeric(tmpArgs[grep("nCPU", tmpArgs[,1], ignore.case=TRUE),2])

system(command=paste("cp ", viaGTFgz, " ", RNAseqSEApath, "/Refs/", sep=""), intern=TRUE)
viaGTFz = paste(RNAseqSEApath, "/Refs/", basename(viaGTFgz), sep="")
system(command=paste("gunzip ", viaGTFz, sep=""), intern=TRUE)
viaGTF = gsub(".gz$", "", viaGTFz)
nomStruct = paste(esp, "_genecode", GCv, sep="")
library(rtracklayer)
tmpGTF = as.data.frame(readGFF(viaGTF))
cols = c("seqid", "strand", "gene_id", "transcript_id", "exon_id", "start", "end", "gene_name", "level")
if(sum(cols%in%colnames(tmpGTF))==length(cols))	refBED = tmpGTF[grepl("exon", tmpGTF[,"type"]),cols]
colnames(refBED) = c("chromosome_name", "strand", "ensembl_gene_id", "ensembl_transcript_id", "ensembl_exon_id", "exon_chrom_start", "exon_chrom_end", "external_gene_name", "level")
refBED[,"chromosome_name"] = gsub("chr", "", refBED[,"chromosome_name"])

if(sum(grepl("[_]", refBED[,"ensembl_gene_id"]))>0)
{
	refBED[,"ensembl_gene_id"] = gsub("_.*", "", refBED[,"ensembl_gene_id"])
	refBED[,"ensembl_transcript_id"] = gsub("_.*", "", refBED[,"ensembl_transcript_id"])
	refBED[,"ensembl_exon_id"] = gsub("_.*", "", refBED[,"ensembl_exon_id"])
}
if(sum(grepl("[.]", refBED[,"ensembl_gene_id"]))>0)
{
	refBED[,"ensembl_gene_id"] = gsub("[.].*", "", refBED[,"ensembl_gene_id"])
	refBED[,"ensembl_transcript_id"] = gsub("[.].*", "", refBED[,"ensembl_transcript_id"])
	refBED[,"ensembl_exon_id"] = gsub("[.].*", "", refBED[,"ensembl_exon_id"])
}
saveRDS(refBED, paste("/data-isilon/Cagnard/RNAseqSEA/Refs/", nomStruct, ".rds", sep=""))

#	II 
tmpAllJuncPath = paste(RNAseqSEApath, "/Refs/tmp_", nomStruct, sep="")
if(!file.exists(tmpAllJuncPath))	dir.create(tmpAllJuncPath)

junGene<-function(E)
{
	TabAllJunc = matrix(0, ncol=6, nrow=0)
	colnames(TabAllJunc) = c("chromosome_name", "ensembl_gene_id", "ensembl_transcript_id", "exons_junctions", "junction_chrom_start", "junction_chrom_end")
	
	Genes = partGenes[[E]]
	for(G in 1:length(Genes))
	{
		write(paste("\t#\tPart ", E, "-", G, "/", length(Genes), ": ", Genes[G], sep=""), file="")
		exons = refBED[refBED[,"ensembl_gene_id"]%in%Genes[G],,drop=FALSE]
		#	exons = exons[as.numeric(exons[,"level"])<=2,]	#	filtrage level Get level 1 & 2 annotation (manually annotated) only
		exons = exons[order(as.numeric(exons[,"exon_chrom_start"])),,drop=FALSE]
		exons = exons[order(exons[,"ensembl_transcript_id"]),,drop=FALSE]
		
		JuncT = list()	#	Liste des jonctions par transcrit
		for(Tid in unique(exons[,"ensembl_transcript_id"]))	#	Tid = unique(exons[,"ensembl_transcript_id"])[1]
		{
			Texons = exons[exons[,"ensembl_transcript_id"]==Tid,,drop=FALSE]
			if(nrow(Texons)>=2)
			{
				Texons = Texons[order(as.numeric(Texons[,"exon_chrom_start"])),,drop=FALSE]
				Tjunc = cbind(Texons[1:(nrow(Texons)-1),"exon_chrom_end"], Texons[2:(nrow(Texons)),"exon_chrom_start"])
				
				Juncs=NULL
				for(J in 1:nrow(Tjunc))	Juncs = c(Juncs, paste(J,(J+1), sep="-"))
				
				Tjunc = cbind(rep(Texons[1, "chromosome_name"], nrow(Tjunc)), rep(Texons[1, "ensembl_gene_id"], nrow(Tjunc)), rep(Texons[1, "ensembl_transcript_id"], nrow(Tjunc)), Juncs, Tjunc)
				TabAllJunc = rbind(TabAllJunc, Tjunc)
			}
		}
		saveRDS(TabAllJunc, file=paste(tmpAllJuncPath, "/tmpAllJunc_", E, ".rds", sep=""))
	}
}

ENSg = unique(refBED[,"ensembl_gene_id"])
chunk <- function(ENSg,nCPU) split(ENSg, factor(sort(rank(ENSg)%%nCPU)))
partGenes = chunk(ENSg,nCPU)

runif(1)
library(parallelMap)
parallelStart(mode = "multicore", cpus=nCPU, show.info=TRUE) 
f = function(E) junGene(E)
y = parallelMap(f, c(1:length(partGenes)))
parallelStop()

TabAllJunc = list()
AllrdsFiles = list.files(tmpAllJuncPath, full.names=TRUE)
for(N in AllrdsFiles)	TabAllJunc[[length(TabAllJunc)+1]] = readRDS(N)
TabAllJunc = do.call("rbind",TabAllJunc)
TabAllJunc = TabAllJunc[,c("chromosome_name", "junction_chrom_start", "junction_chrom_end", "ensembl_gene_id", "ensembl_transcript_id")]

TabAllJunc = TabAllJunc[order(as.numeric(TabAllJunc[,"junction_chrom_start"])),]
TabAllJunc = TabAllJunc[order(TabAllJunc[,"chromosome_name"]),]
write.table(TabAllJunc, file=paste(RNAseqSEApath, "/Refs/Junc_", nomStruct, ".bed", sep=""), sep="\t", quote=FALSE)
saveRDS(TabAllJunc, file=paste(RNAseqSEApath, "/Refs/Junc_", nomStruct, ".rds", sep=""))

unlink(tmpAllJuncPath, recursive=TRUE)

		