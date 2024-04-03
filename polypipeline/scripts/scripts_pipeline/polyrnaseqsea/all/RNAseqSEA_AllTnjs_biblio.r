###############################################################################################################
#	0- Sources

if(FALSE)
{
	library("rjson")
	allNFO <- fromJSON(file = "/data-isilon/Cagnard/RNAseqSEA/RNAseqSEA_AllT/NGS2019_2498.json")
	
	#	"outputs"   "analyse"   "project"   "softwares" "inputs"
	allNFO$outputs$ri	#	_RI.tab
	allNFO$analyse	#	date
	allNFO$project$gencode	#	"release" "rds"     "gtf"     "version"
	allNFO$softwares	#	"samtools" "sambamba" "bedtools"
	allNFO$inputs	#	bam path
	
	tabNFO = matrix(0, ncol=2, nrow=length(allNFO$outputs$ri))
	colnames(tabNFO) = c("out", "bam")
	
	tabNFO[,"out"] = dirname(unlist(allNFO$outputs$ri))
	tabNFO[,"bam"] = unlist(allNFO$inputs$bam)
	
	esp = allNFO$project$gencode$release
	GCv = allNFO$project$gencode$version
}

###############################################################################################################
#	I- Filtrage duplicats, selection des reads de junction et identification des bornes des reads

rmdup<-function(I)
{	
	write(paste("\t#\trmdup ", bams[I], sep=""), file="")
	#	system(paste("sambamba markdup -r -t 8 --tmpdir=",rmdupBamPath, " ", bams[I], " ", rmdupBamPath, basename(bams[I]), sep=""))
	
	system(paste("java -jar /software/bin/picard.jar MarkDuplicates --INPUT ", bams[I], " --OUTPUT ", rmdupBamPath, basename(bams[I]), " --METRICS_FILE ", gsub(".bam$", ".txt$", basename(bams[I])), sep=""))
}

splitBamSamTools<-function(B)	#	splitBamSamTools(B)
{
	chrs = system(paste("samtools view -@ 20 -H ", bams[B], " | perl -lne '/SN:(\\S+)/ and print $1'", sep=""), intern=TRUE)
	if(!is.null(ReName))
	{
		write(paste("\t#\t => Bams Rename", sep=""), file="")
		BamsReNames = as.matrix(read.table(ReName, sep="\t", header=TRUE))
		BamName = BamsReNames[BamsReNames[,1]%in%gsub(".bam$", "", basename(bams[B])),2]
	}else{
		BamName = gsub(".bam", "", basename(bams[B]))
	}
	
	BamPartsPath = paste(splitBamPath, "/", BamName, "/", sep="")
	if(!file.exists(BamPartsPath))	dir.create(BamPartsPath)
	
	for(C in 1:length(chrs))
	{
		write(paste("# Split de ", BamName, " [", B, "/", length(bams), "] chr:", chrs[C], sep=""), file="")
		system(paste("samtools view -@ 20 ", bams[B], " ", chrs[C], " -b > ", BamPartsPath, "/", chrs[C], ".bam", sep=""))
	}
}

###############################################################################################################
#	II Selection des reads de jonctions => bed

selectJunc<-function(B)
{
	nCPU = 4
	#avant = system(paste("samtools flagstat -@ ", nCPU, " ", bams[B], sep=""), intern=TRUE)
	#stats = rbind(stats, gsub(" .*", "", avant)[apply(as.matrix(nfo), 1, function(x) min(grep(x, avant)))])
	#	supprime les duplicates et selectionne les reads avec N
	write(paste("\t#\t", B, "/", length(bamsSplit), " - ", basename(dirname(bamsSplit[B])), "/", gsub(".bam", "", basename(bamsSplit[B])), sep=""), file="")
	system(paste("samtools view -@ ", nCPU, " -h -F 0x0400 ", bamsSplit[B], " | awk '$6 ~ /N/ || $1 ~ /^@/' | samtools view -@ ", nCPU, " -bS - > ", pathBamJuncs, "/", basename(dirname(bamsSplit[B])), "_", basename(bamsSplit[B]), sep=""))
	#apres = system(paste("samtools flagstat -@ ", nCPU, " ", pathBamJuncs, "/", basename(dirname(bamsSplit[B])), "_", basename(bamsSplit[B]), sep=""), intern=TRUE)
	#stats = rbind(stats, gsub(" .*", "", apres)[apply(as.matrix(nfo), 1, function(x) min(grep(x, apres)))])
	
	write(paste("\t#\t\t=> ", basename(dirname(bamsSplit[B])), "/", gsub(".bam", ".bed", basename(bamsSplit[B])), sep=""), file="")
	tmp = system(paste("samtools view -@ ", nCPU, " ", pathBamJuncs, "/", basename(dirname(bamsSplit[B])), "_", basename(bamsSplit[B]), " | cut -f3,4,6  | awk '{if ($3 ~ /N/){printf(\"%s\",$1); start=$2; end=$2; split($3,a,\"[NIMDSH]\"); 
							split($3,b,\"[0-9]*\"); nb=length(b); for(i=2; i<=nb; i++){if(b[i] ~ /[M]/){ end=end+a[i-1];} if(b[i] ~ /N/){end=end-1; printf(\"\t%s\t%s\",start, end); 
							start=end+a[i-1]+1; end=start;}} printf(\"\t%s\t%s\\n\", start,end-1);}}' > ", pathBed,basename(dirname(bamsSplit[B])), "_", gsub(".bam", "", basename(bamsSplit[B])), "_junctions.bed", sep=""), intern=TRUE)
	#	saveRDS(tmp, file=paste(pathBed,gsub(".bam", "", basename(bams[B])), "_junctions.rds", sep=""))
}

###############################################################################################################
#	III bed en RDS

bed2rds<-function(S)
{
	write(paste("\t#\t bed2rds: ", gsub(".bed", "", basename(listBed[S])), sep=""), file="")
	tmp <- try(read.table(listBed[S], sep="\t", quote="", fill=TRUE, col.names = paste0("Col",seq_len(20))))
	tmp = tmp[,2:ncol(tmp),drop=FALSE]
	saveRDS(tmp, file=gsub(".bed", ".rds", listBed[S]))
}

###############################################################################################################
#	IV calcul des positions des junctions et comptage

formatJunc<-function(x)
{
	posJ = c(x[!is.na(x)][c(2:(sum(!is.na(x))-1))])
	juncs = NULL
	for(J in 1:(length(posJ)-1))	juncs = c(juncs, posJ[J], posJ[J+1])
	resPairsPos = matrix(unlist(juncs), ncol=2, byrow=TRUE)
	return(resPairsPos)
	#	resPairsPos = apply(resPairsPos, 1, function(x) paste(x, sep="_", collapse="_"))
}

coordsJunc<-function(P)
{
	write(paste("\t#\t coordsJunc: ", P, "/", length(listBedRds), " - ", gsub(".rds", "", basename(listBedRds[P])), sep=""), file="")
	tmpBed = readRDS(listBedRds[P])
	tmpRes = apply(tmpBed, 1, function(x) formatJunc(x))
	if(is.list(tmpRes))
	{
		tmpRes = do.call("rbind",tmpRes)
	}else{
		tmpRes = t(tmpRes)
	}
	tmpRes = apply(tmpRes, 1, function(x) paste(x, collapse="_"))
	JunCounts = as.matrix(table(tmpRes))
	tmpC = t(apply(as.matrix(rownames(JunCounts)), 1, function(x) unlist(strsplit(x, "_"))))
	
	coords = cbind(as.numeric(tmpC[,1]), as.numeric(tmpC[,2]), JunCounts)
	coords = coords[order(as.numeric(coords[,2])),,drop=FALSE]
	coords = coords[order(coords[,1]),,drop=FALSE]
	colnames(coords) = c("StartJun", "EndJun", "CountJun")
	
	#	corrige le pb du ch dans nom de chrom
	chr = gsub(".*_", "", gsub("_junctions", "", gsub(".rds", "", basename(listBedRds[P]))))
	rownames(coords) = paste(chr, "_", rownames(coords), sep="")
	
	totRawJunc = nrow(coords)
	coords = coords[as.numeric(coords[,"CountJun"])>=limCountJun,,drop=FALSE]	#	limCountJun=2
	write(paste("Filtrage ",gsub(".rds", "", basename(listBedRds[P])),  " a limCountJun=", limCountJun, " :", totRawJunc, " => ", nrow(coords), "[", signif(100*((nrow(coords)-totRawJunc)/totRawJunc), digit=4), "%]", sep=""), file="")
	
	if(nrow(coords)>0) saveRDS(coords, paste(AllJuncPath, "/", basename(listBedRds[P]), sep=""))
}

###############################################################################################################
#	V Selection de la liste de genes
#	selection des ENSG ID des gènes qui ont des junctions

caracterizeJuncs<-function(C)	#	patient par patient
{
	nom = gsub(".rds", "", gsub("AllJunctionsCounts", "AllRes", basename(listCoordsBed[C])))
	write(paste("\t#\t", C, "/", length(listCoordsBed), "-", nom, sep=""), file="")

	bedtest = readRDS(listCoordsBed[C])
	
	Chr = gsub(addChr, "", gsub("_.*", "", rownames(bedtest)[1]))
	
	testINgene<-function(x)
	{
		bedRefTmp = bedRef[grepl(Chr, bedRef[,1]),]
		INgene = ((bedRefTmp[,2]<=as.numeric(x["StartJun"]))&
							(bedRefTmp[,3]>=as.numeric(x["StartJun"])))|
						((bedRefTmp[,2]<=as.numeric(x["EndJun"]))&
							(bedRefTmp[,3]>=as.numeric(x["EndJun"])))
		if(sum(INgene)>0)	return(as.character(bedRefTmp[INgene,4]))
		if(sum(INgene)==0)	return(NULL)
	}
	
	resIngene = apply(bedtest, 1, function(x) testINgene(x))
	genesList = unique(unlist(resIngene))
	saveRDS(genesList, paste(geneListsPath, gsub("AllJunctionsCounts", "ENSg", basename(listCoordsBed[C])), sep=""))
}

###############################################################################################################
#	VI	Caracterisation des jonctions

testanalysePatGene<-function(G)	#	gene="ENSG00000002746"; nom="cDNA_Fibro_DIE_MAG"
{
	nom = ToDoTab[G, 1]
	JunEch = readRDS(paste(AllJuncPath, "/", nom, "_junctions.rds", sep=""))
	JunEch = JunEch[!is.na(as.numeric(JunEch[,"EndJun"])),,drop=FALSE]
	Chr = geneRef[1,"chromosome_name"]

	if((nrow(JunEch)==0))
	{
		write(paste("\n\n#######################################################################", sep=""), file="")
		write(paste("# Analyse [", signif(100*(G/nrow(ToDoTab)), digit=4), "%] \t", G, "/", nrow(ToDoTab), "\t-\t", nom, "\t/\t", gene, sep=""), file="")
		write(paste("#\t!!!!!! nrow(JunEch)==0 !!!!!!", sep=""), file="")
		write(paste("#######################################################################\n", sep=""), file="")
	}
	#	junctions du gene
	if(sum(grepl(Chr, JunEch[,"Chr"]))==0)
	{
		write(paste("\n\n#######################################################################", sep=""), file="")
		write(paste("# Analyse [", signif(100*(G/nrow(ToDoTab)), digit=4), "%] \t", G, "/", nrow(ToDoTab), "\t-\t", nom, "\t/\t", gene, sep=""), file="")
		write(paste("#\t!!!!!! Pas de chr ", Chr, " dans le JunEch !!!!!!", sep=""), file="")
		write(paste("#######################################################################\n", sep=""), file="")
	}
	#geneJunEch = JunEch[JunEch[,"Chr"]==Chr,,drop=FALSE]
	#geneJunEch = geneJunEch[as.numeric(geneJunEch[,"EndJun"])>=Start,,drop=FALSE]
	#geneJunEch = geneJunEch[as.numeric(geneJunEch[,"StartJun"])<=End,,drop=FALSE]
}

#	http://www.polyweb.fr/NGS/NGS2023_6667/HG19_MT/RNAseqSEA/rmDupsBams/ARN1262_21MY1072_M.bam
analysePatGene<-function(G)	#	gene="ENSG00000002746"; nom="cDNA_Fibro_DIE_MAG"
{
	system(paste("clear", sep=""))
	gene = ToDoTab[G, 2]
	nom = ToDoTab[G, 1]
	
	IsResSE = IsResRI = FALSE
	
	pathResPat = paste(AllresPath, "/", nom, "/", sep="")
	if(!file.exists(pathResPat))	dir.create(pathResPat)
	
	TabNFO = matrix(0, ncol=6, nrow=0)	#	Taille des annots
	colnames(TabNFO) = c("ENSID", "nJ", "Chr", "Start", "End", "Interpretation")
	
	write(paste("\n\n################################################################", sep=""), file="")
	write(paste("# Analyse [", signif(100*(G/nrow(ToDoTab)), digit=4), "%] ", G, "/", nrow(ToDoTab), " - ", nom, " /", gene, sep=""), file="")
	write(paste("################################################################\n", sep=""), file="")
	
	JunEch = readRDS(paste(AllJuncPath, "/", nom, "_junctions.rds", sep=""))
	#	JunEch = JunEch[!is.na(as.numeric(JunEch[,"EndJun"])),,drop=FALSE]
	# nomGene = as.character(Annots[Annots[,"EnsID"]%in%gene,"Symbol"])
	nomGene = unique(as.character(refBED[refBED[,"ensembl_gene_id"]%in%gene,"external_gene_name"]))
	geneRef = refBED[refBED[,"ensembl_gene_id"]%in%gene,,drop=FALSE]
	geneJuncREF = juncBED[juncBED[,"ensembl_gene_id"]%in%gene,,drop=FALSE]
	geneJuncREF = cbind(as.numeric(geneJuncREF[,"junction_chrom_start"]), as.numeric(geneJuncREF[,"junction_chrom_end"]))	
	geneJuncREF = cbind(geneJuncREF, apply(geneJuncREF, 1, function(x) x[2]-x[1]))
	colnames(geneJuncREF) = c("StartJun", "EndJun", "JuncLength" )
	geneJuncREF = geneJuncREF[order(as.numeric(geneJuncREF[,"StartJun"])),,drop=FALSE]
	geneJuncREF = geneJuncREF[!duplicated(apply(geneJuncREF, 1, function(x) paste(gsub(" ", "", x[1]), "_", gsub(" ", "", x[2]), sep=""))),,drop=FALSE]
	Chr = geneRef[1,"chromosome_name"]
	rownames(geneJuncREF) = apply(geneJuncREF, 1, function(x) paste(Chr, "_", gsub(" ", "", x[1]), "_", gsub(" ", "", x[2]), sep="", collapse=""))
	
	if((nrow(geneJuncREF)>=1)&(nrow(JunEch)>=1))	#	Si au moins une jonction dans le gene
	{
		Chr = geneRef[1,"chromosome_name"]
		Start = min(geneRef[,c("exon_chrom_start", "exon_chrom_end")])
		End = max(geneRef[,c("exon_chrom_start", "exon_chrom_end")])
		
		########################################################################################
		#	Junctions possibles, junctions SPE de transcrit
		geneRef = geneRef[order(as.numeric(geneRef[,"exon_chrom_start"])),,drop=FALSE]
		geneRef = geneRef[order(geneRef[,"ensembl_transcript_id"]),,drop=FALSE]
		JuncT = list()	#	Liste des jonctions par transcrit
		for(Tid in unique(geneRef[,"ensembl_transcript_id"]))
		{
			Texons = geneRef[geneRef[,"ensembl_transcript_id"]==Tid,]
			Texons = Texons[order(as.numeric(Texons[,"exon_chrom_start"])),]
			Tjunc = cbind(Texons[1:(nrow(Texons)-1),"chromosome_name"], Texons[1:(nrow(Texons)-1),"exon_chrom_end"], Texons[2:(nrow(Texons)),"exon_chrom_start"])
			JuncT[[length(JuncT)+1]] = apply(Tjunc, 1, function(x) paste(gsub(" ", "", x[1]), "_", gsub(" ", "", x[2]), "_", gsub(" ", "", x[3]), sep=""))
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
		write.table(TabAllJunc, file=paste(pathResPat, "/TableAll_Junctions_Transcrits.txt", sep=""), sep="\t", quote=FALSE)
		
		#	Jonctions specifique de transcrit (parmis ceux pris en compte)
		JuncCount = apply(TabAllJunc, 1, function(x) sum(x))
		speJ = names(JuncCount)[JuncCount==1]
		speJ = cbind(speJ, apply(as.matrix(speJ), 1, function(x) names(JuncT)[grepl(x, JuncT)]))
		colnames(speJ) = c("JunctionSPE", "ensembl_transcript_id")
		
		#	junctions du gene
		if(sum(grepl(paste(addChr, Chr, "_", sep=""), rownames(JunEch)))==0)	write(paste("\t#\t!!!!!! Pas de chr ", Chr, " dans le JunEch !!!!!!", sep=""), file="")

		#	geneJunEch = JunEch[JunEch[,"Chr"]==Chr,,drop=FALSE]
		geneJunEch = JunEch[grepl(paste(addChr, Chr, "_", sep=""), rownames(JunEch)),,drop=FALSE]
		geneJunEch = geneJunEch[as.numeric(geneJunEch[,"EndJun"])>=Start,,drop=FALSE]
		geneJunEch = geneJunEch[as.numeric(geneJunEch[,"StartJun"])<=End,,drop=FALSE]
		
		if(nrow(geneJunEch)>0)
		{
			pathResGenePat = paste(pathResPat, "/", gene, "/", sep="")
			if(!file.exists(pathResGenePat))	dir.create(pathResGenePat)
			nfoname = paste(pathResGenePat, "/", nom, "_nfo_junc.txt", sep="")
			
			#	=> correction des pos
			geneJunEchCorr = geneJunEch
			for(C in 1:nrow(geneJunEch))
			{
				distStart = abs(as.numeric(geneJunEch[C,"StartJun"])-as.numeric(geneJuncREF[,"StartJun"]))
				distEnd = abs(as.numeric(geneJunEch[C,"EndJun"])-as.numeric(geneJuncREF[,"EndJun"]))

#	gestion du pb des corrections, recherche de la borne la plus proche dans l'intervalle de limDecal			
				if((min(distStart>0))&(min(distStart)<=limDecal))
				{
					if(length(unique((as.numeric(geneJuncREF[,"StartJun"])[distStart==(min(distStart))])))>1)
					{
						write(paste("C=", C, sep=""), file="")
						pb = c(gene, nom, paste("Correction bornes junction limDecal=", limDecal, paste(geneJunEch[C,], sep="_", collapse="_"), sep=""))
						write.table(pb, paste(ErrorPath, "/", nom, "_", gene, "_", paste(geneJunEch[C,], sep="_", collapse="_"), ".txt", sep=""), sep="\t", quote=FALSE)
					}
					geneJunEchCorr[C,"StartJun"] = min(as.numeric(geneJuncREF[,"StartJun"])[distStart==(min(distStart))])
				}
				if((min(distEnd>0))&(min(distEnd)<=limDecal)) geneJunEchCorr[C,"EndJun"] = min(as.numeric(geneJuncREF[,"EndJun"])[distEnd==(min(distEnd))])
				
				rownames(geneJunEchCorr) = apply(geneJunEchCorr, 1, function(x) paste(gsub(" ", "", x[1]), "_", gsub(" ", "", x[2]), "_", gsub(" ", "", x[3]), sep=""))
			}
			geneJunEch = geneJunEchCorr
			
			#geneJunEch = as.data.frame(geneJunEch)
			geneJunEch = cbind(as.numeric(geneJunEchCorr[,"StartJun"]), as.numeric(geneJunEchCorr[,"EndJun"]), as.numeric(geneJunEchCorr[,"CountJun"]))
			colnames(geneJunEch) = c("StartJun", "EndJun", "CountJun")
			rownames(geneJunEch) = apply(geneJunEch, 1, function(x) paste(Chr, "_", x[1], "_", x[2], sep=""))
			#for(DF in c("StartJun", "EndJun", "CountJun"))	geneJunEch[,DF] = as.numeric(as.character(geneJunEch[,DF]))
			
			cols = c("Junc_RI", "ENSID", "Gene", "Chr", "Junc_RI_Start", "Junc_RI_End", "Junc_RI_Count", "Junc_Normale", "Junc_Normale_Start", "Junc_Normale_End", "Junc_Normale_Count", "Type")
			matEV_RI = matrix(0, nrow=0, ncol=length(cols))
			colnames(matEV_RI) = cols
			
			#	Test si start et end des junctions sont connus
			test_RI = ((!(geneJunEch[,"StartJun"]%in%geneJuncREF[,"StartJun"])|!(geneJunEch[,"EndJun"]%in%geneJuncREF[,"EndJun"]))&(!(geneJunEch[,"StartJun"]%in%geneJunEch[,"EndJun"])|!(geneJunEch[,"EndJun"]%in%geneJunEch[,"StartJun"])))
			
			tmpExons = geneRef[order(as.numeric(geneRef[,"exon_chrom_start"])),]
			tmpExons = tmpExons[!duplicated(tmpExons[,"ensembl_exon_id"]),]
			Map = rbind(cbind(paste("Start_E", c(1:nrow(tmpExons)), sep="_"), tmpExons[,"exon_chrom_start"], tmpExons[,"ensembl_exon_id"]),
					cbind(paste("End_E", c(1:nrow(tmpExons)), sep="_"), tmpExons[,"exon_chrom_end"], tmpExons[,"ensembl_exon_id"]))
			colnames(Map) = c("N", "Pos", "EnsID_JuncType")
			
			if(sum(test_RI)>0)
			{
				J_anorm = geneJunEch[test_RI, c("StartJun", "EndJun", "CountJun"),drop=FALSE]	#	Junctions anormales				
				colresNFOjun = c("ENSID", "Junction", "Chr", colnames(J_anorm), "StartConnu", "EndConnu", "Start_in_exon", "End_in_exon", "Start_in_intron", "End_in_intron", "Start_in_3p",
						"End_in_3p", "Start_in_5p", "End_in_5p", "Start_in_exon1", "End_in_exon1", "Start_in_DerExon", "End_in_DerExon", "End_EQ_StartExon1", "Etart_EQ_EndDerExon",
						"EtartEnd_Meme_exon", "StartEnd_Meme_intron", "IE_contigus", "EI_contigus", "FinJrecouvreExon_SE", "DebutJrecouvreExon_SE")
				
				resNFOjun = matrix(0, ncol=length(colresNFOjun), nrow=0)
				colnames(resNFOjun) = colresNFOjun
				TabNFO = cbind(TabNFO, rep(0, nrow(TabNFO)))	#	Ajoute une colonne / ech
				colnames(TabNFO)[ncol(TabNFO)] = nom
				
				for(J in 1:nrow(J_anorm))
				{
					nJ = rownames(J_anorm)[J]
					
					#	Borne connue / inconnue
					test_StartConnu = J_anorm[J,"StartJun"]%in%geneRef[,"exon_chrom_end"]
					test_EndConnu = J_anorm[J,"EndJun"]%in%geneRef[,"exon_chrom_start"]
					
					test_EndJ_EndX = test_StartJ_StartX = FALSE
					if(!test_EndConnu)	test_EndJ_EndX = J_anorm[J,"EndJun"]%in%geneRef[,"exon_chrom_end"]	#	La junction se termine avec un exon sauté => SE
					if(!test_StartConnu)	test_StartJ_StartX = J_anorm[J,"StartJun"]%in%geneRef[,"exon_chrom_start"]	#	La junction debute avec un exon connu
					
					#	Jonction au moins partiellement dans un exon
					Jstart_in_exon = (J_anorm[J, "StartJun"]>geneRef[,"exon_chrom_start"])&(J_anorm[J, "StartJun"]<geneRef[,"exon_chrom_end"])
					test_Jstart_in_exon = sum(Jstart_in_exon)>0
					Jend_in_exon = (J_anorm[J, "EndJun"]>geneRef[,"exon_chrom_start"])&(J_anorm[J, "EndJun"]<geneRef[,"exon_chrom_end"])
					test_Jend_in_exon = sum(Jend_in_exon)>0
					
					#	Jonction au moins partiellement dans un intron
					test_Jstart_in_intron = !(test_Jstart_in_exon|test_StartConnu)&(J_anorm[J, "StartJun"]>min(geneRef[,"exon_chrom_start"]))&(J_anorm[J, "StartJun"]<max(geneRef[,"exon_chrom_end"]))
					test_Jend_in_intron = !(test_Jend_in_exon|test_EndConnu)&(J_anorm[J, "EndJun"]>min(geneRef[,"exon_chrom_start"]))&(J_anorm[J, "EndJun"]<max(geneRef[,"exon_chrom_end"]))
					
					#	Hors gene
					test_Jstart_in_3p = (J_anorm[J, "StartJun"]>max(geneRef[,"exon_chrom_end"]))
					test_Jend_in_3p = (J_anorm[J, "EndJun"]>max(geneRef[,"exon_chrom_end"]))
					
					test_Jend_in_5p = (J_anorm[J, "EndJun"]<min(geneRef[,"exon_chrom_start"]))
					test_Jstart_in_5p = (J_anorm[J, "StartJun"]<min(geneRef[,"exon_chrom_start"]))
					
					#	1 er exon
					test_Jstart_in_exon1 = (J_anorm[J, "StartJun"]>min(geneRef[,"exon_chrom_start"]))&(J_anorm[J, "StartJun"]<min(geneRef[,"exon_chrom_end"]))
					test_Jend_in_exon1 = (J_anorm[J, "EndJun"]>min(geneRef[,"exon_chrom_start"]))&(J_anorm[J, "EndJun"]<min(geneRef[,"exon_chrom_end"]))
					
					#	Der exon
					test_Jstart_in_Derexon = (J_anorm[J, "StartJun"]>max(geneRef[,"exon_chrom_start"]))&(J_anorm[J, "StartJun"]<max(geneRef[,"exon_chrom_end"]))
					test_Jend_in_Derexon = (J_anorm[J, "EndJun"]>max(geneRef[,"exon_chrom_start"]))&(J_anorm[J, "EndJun"]<max(geneRef[,"exon_chrom_end"]))
					
					test_Jend_EQ_Startexon1 = J_anorm[J, "EndJun"]==min(geneRef[,"exon_chrom_start"])	#	Fin Junc en debut exon1
					test_Jstart_EQ_EndDerexon = J_anorm[J, "StartJun"]== max(geneRef[,"exon_chrom_end"]) #	Debut Junc en fin dernier exon
					
					#	coords exon
					if(test_Jstart_in_exon)	# coords exon inclurant Jstart
					{
						Deb_Ex_Incl_Jstart = max(geneRef[geneRef[,"exon_chrom_start"]<J_anorm[J, "StartJun"], "exon_chrom_start"])
						Fin_Ex_Incl_Jstart = max(geneRef[geneRef[,"exon_chrom_end"]>J_anorm[J, "StartJun"], "exon_chrom_end"])
					}
					if(test_Jend_in_exon)	# coords exon inclurant Jend
					{
						Deb_Ex_Incl_Jend = max(geneRef[geneRef[,"exon_chrom_start"]<J_anorm[J, "EndJun"], "exon_chrom_start"])
						Fin_Ex_Incl_Jend = max(geneRef[geneRef[,"exon_chrom_end"]>J_anorm[J, "EndJun"], "exon_chrom_end"])
					}
					#	Start et End dans le meme exon
					test_JstartJend_Meme_exon = FALSE
					if(test_Jstart_in_exon&test_Jend_in_exon)	test_JstartJend_Meme_exon = (Deb_Ex_Incl_Jstart==Deb_Ex_Incl_Jend)&(Fin_Ex_Incl_Jstart==Deb_Ex_Incl_Jend)
					
					#	coords intron
					if(test_Jstart_in_intron)	# coords intron inclurant Jstart
					{
						Deb_Int_Incl_Jstart = max(geneRef[geneRef[,"exon_chrom_end"]<J_anorm[J, "StartJun"], "exon_chrom_end"])
						Fin_Int_Incl_Jstart = min(geneRef[geneRef[,"exon_chrom_start"]>J_anorm[J, "StartJun"], "exon_chrom_start"])
					}
					if(test_Jend_in_intron)	# coords intron inclurant Jend
					{
						Deb_Int_Incl_Jend = max(geneRef[geneRef[,"exon_chrom_end"]<J_anorm[J, "EndJun"], "exon_chrom_end"])
						Fin_Int_Incl_Jend = min(geneRef[geneRef[,"exon_chrom_start"]>J_anorm[J, "EndJun"], "exon_chrom_start"])
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
					canonGeneJun = geneJunEch[rownames(geneJunEch)%in%rownames(geneJuncREF),,drop=FALSE]
					
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
						expl = "Start in 3p ET end in 3p => RI en 3p"
						nomJ = paste(nJ, "=> RI_3p", sep="")
						line = paste(nomJ, "\t", gene, "\t", nomGene, "\t", "\t", geneJunEch[J,"StartJun"], "\t", geneJunEch[J,"EndJun"], "\t", geneJunEch[nJ,"CountJun"], "\t", expl, sep="")
						write(line, file=nfoname, append=TRUE)
						write(line, file="")
						matEV_RI = rbind(matEV_RI, cbind(newRI, as.matrix(rep("SE", nrow(newRI)))))
						
						if(sum(TabNFO[,"nJ"]%in%nomJ)==0)	TabNFO = rbind(TabNFO, c(gene, nomJ, Chr, geneJunEch[J,"StartJun"], geneJunEch[J,"EndJun"], expl, rep(0, ncol(TabNFO)-6)))
						TabNFO[TabNFO[,"nJ"]%in%nomJ,nom] = geneJunEch[nJ,"CountJun"]
					}
					
					if(test_Jstart_in_exon&test_Jend_in_3p)	#	Start in exon ET end in 3p
					{
						expl = "Start in exon ET end in 3p => raccourci exon en 3 et RI"
						nomJ = paste(nJ, "=> A3SS+RI", sep="")
						line = paste(nomJ, "\t", gene, "\t", nomGene, "\t", "\t", geneJunEch[J,"StartJun"], "\t", geneJunEch[J,"EndJun"], "\t", geneJunEch[nJ,"CountJun"], "\t", expl, sep="")
						write(line, file=nfoname, append=TRUE)
						write(line, file="")
						matEV_RI = rbind(matEV_RI, cbind(newRI, as.matrix(rep("SE", nrow(newRI)))))
						
						if(sum(TabNFO[,"nJ"]%in%nomJ)==0)	TabNFO = rbind(TabNFO, c(gene, nomJ, Chr, geneJunEch[J,"StartJun"], geneJunEch[J,"EndJun"], expl, rep(0, ncol(TabNFO)-6)))
						TabNFO[TabNFO[,"nJ"]%in%nomJ,nom] = geneJunEch[nJ,"CountJun"]
					}
					
					if(test_StartJ_StartX&test_Jend_in_exon1)	#	Start avec exon 1 ET end in exon 1
					{
						expl = "Start avec exon 1 et fin dans exon 1 => raccourci exon 1 en 5"
						nomJ = paste(nJ, "=> SE", sep="")
						line = paste(nomJ, "\t", gene, "\t", nomGene, "\t", "\t", geneJunEch[J,"StartJun"], "\t", geneJunEch[J,"EndJun"], "\t", geneJunEch[nJ,"CountJun"], "\t", expl, sep="")
						write(line, file=nfoname, append=TRUE)
						write(line, file="")
						matEV_RI = rbind(matEV_RI, cbind(newRI, as.matrix(rep("SE", nrow(newRI)))))
						
						if(sum(TabNFO[,"nJ"]%in%nomJ)==0)	TabNFO = rbind(TabNFO, c(gene, nomJ, Chr, geneJunEch[J,"StartJun"], geneJunEch[J,"EndJun"], expl, rep(0, ncol(TabNFO)-6)))
						TabNFO[TabNFO[,"nJ"]%in%nomJ,nom] = geneJunEch[nJ,"CountJun"]
					}
					
					if(test_StartJ_StartX&test_EndConnu)	#	(Start exon connu) ET (End connu)
					{
						expl = "Start avec exon et fin avec intron"
						nomJ = paste(nJ, "=> SE", sep="")
						line = paste(nomJ, "\t", gene, "\t", nomGene, "\t", "\t", geneJunEch[J,"StartJun"], "\t", geneJunEch[J,"EndJun"], "\t", geneJunEch[nJ,"CountJun"], "\t", expl, sep="")
						write(line, file=nfoname, append=TRUE)
						write(line, file="")
						matEV_RI = rbind(matEV_RI, cbind(newRI, as.matrix(rep("SE", nrow(newRI)))))
						
						if(sum(TabNFO[,"nJ"]%in%nomJ)==0)	TabNFO = rbind(TabNFO, c(gene, nomJ, Chr, geneJunEch[J,"StartJun"], geneJunEch[J,"EndJun"], expl, rep(0, ncol(TabNFO)-6)))
						TabNFO[TabNFO[,"nJ"]%in%nomJ,nom] = geneJunEch[nJ,"CountJun"]
					}
					
					#	A3SS+SE+RI_3p
					if(test_Jstart_in_exon1&test_Jend_in_3p)	#	(Start in premier exon) ET (End in 3p)
					{
						expl = "Start in premier exon et End in 3p"
						nomJ = paste(nJ, "=>A3SS+SE+RI_3p", sep="")
						line = paste(nomJ, "\t", gene, "\t", nomGene, "\t", "\t", geneJunEch[J,"StartJun"], "\t", geneJunEch[J,"EndJun"], "\t", geneJunEch[nJ,"CountJun"], "\t", expl, sep="")
						write(line, file=nfoname, append=TRUE)
						write(line, file="")
						matEV_RI = rbind(matEV_RI, cbind(newRI, as.matrix(rep("A3SS+SE+RI_3p", nrow(newRI)))))
						
						if(sum(TabNFO[,"nJ"]%in%nomJ)==0)	TabNFO = rbind(TabNFO, c(gene, nomJ, Chr, geneJunEch[J,"StartJun"], geneJunEch[J,"EndJun"], expl, rep(0, ncol(TabNFO)-6)))
						TabNFO[TabNFO[,"nJ"]%in%nomJ,nom] = geneJunEch[nJ,"CountJun"]
					}
					
					#	SE+RI_3p
					if(test_StartJ_StartX&test_Jend_in_3p)	#	(Start exon connu) ET (End in 3p)
					{
						expl = "Start exon connu et End in 3p"
						nomJ = paste(nJ, "=>SE+RI_3p", sep="")
						line = paste(nomJ, "\t", gene, "\t", nomGene, "\t", "\t", geneJunEch[J,"StartJun"], "\t", geneJunEch[J,"EndJun"], "\t", geneJunEch[nJ,"CountJun"], "\t", expl, sep="")
						write(line, file=nfoname, append=TRUE)
						write(line, file="")
						matEV_RI = rbind(matEV_RI, cbind(newRI, as.matrix(rep("RI_5p+SE", nrow(newRI)))))
						
						if(sum(TabNFO[,"nJ"]%in%nomJ)==0)	TabNFO = rbind(TabNFO, c(gene, nomJ, Chr, geneJunEch[J,"StartJun"], geneJunEch[J,"EndJun"], expl, rep(0, ncol(TabNFO)-6)))
						TabNFO[TabNFO[,"nJ"]%in%nomJ,nom] = geneJunEch[nJ,"CountJun"]
					}
					
					#	RI_5p+RI_3p
					if(test_Jstart_in_5p&test_Jend_in_3p)	#	(Start in 5p) ET (End in 5p)
					{
						expl = "Start en 5p et End in 3p"
						nomJ = paste(nJ, "=>RI_5p+RI_3p", sep="")
						line = paste(nomJ, "\t", gene, "\t", nomGene, "\t", "\t", geneJunEch[J,"StartJun"], "\t", geneJunEch[J,"EndJun"], "\t", geneJunEch[nJ,"CountJun"], "\t", expl, sep="")
						write(line, file=nfoname, append=TRUE)
						write(line, file="")
						matEV_RI = rbind(matEV_RI, cbind(newRI, as.matrix(rep("RI_5p+SE", nrow(newRI)))))
						
						if(sum(TabNFO[,"nJ"]%in%nomJ)==0)	TabNFO = rbind(TabNFO, c(gene, nomJ, Chr, geneJunEch[J,"StartJun"], geneJunEch[J,"EndJun"], expl, rep(0, ncol(TabNFO)-6)))
						TabNFO[TabNFO[,"nJ"]%in%nomJ,nom] = geneJunEch[nJ,"CountJun"]
					}
					
					#	RI_amont+SE
					if(test_Jstart_in_intron&test_EndJ_EndX)	#	(Start in intron) ET (End exon connu)
					{
						expl = "Start en 5p et end en fin d un exon connu"
						nomJ = paste(nJ, "=>SE+RI", sep="")
						line = paste(nomJ, "\t", gene, "\t", nomGene, "\t", "\t", geneJunEch[J,"StartJun"], "\t", geneJunEch[J,"EndJun"], "\t", geneJunEch[nJ,"CountJun"], "\t", expl, sep="")
						write(line, file=nfoname, append=TRUE)
						write(line, file="")
						matEV_RI = rbind(matEV_RI, cbind(newRI, as.matrix(rep("RI_5p+SE", nrow(newRI)))))
						
						if(sum(TabNFO[,"nJ"]%in%nomJ)==0)	TabNFO = rbind(TabNFO, c(gene, nomJ, Chr, geneJunEch[J,"StartJun"], geneJunEch[J,"EndJun"], expl, rep(0, ncol(TabNFO)-6)))
						TabNFO[TabNFO[,"nJ"]%in%nomJ,nom] = geneJunEch[nJ,"CountJun"]
					}
					
					#	RI_5p+SE
					if(test_Jstart_in_5p&test_EndJ_EndX)	#	(Start en 5p) ET (End exon connu)
					{
						expl = "Start en 5p et end en fin d un exon connu"
						nomJ = paste(nJ, "=>SE+RI", sep="")
						line = paste(nomJ, "\t", gene, "\t", nomGene, "\t", "\t", geneJunEch[J,"StartJun"], "\t", geneJunEch[J,"EndJun"], "\t", geneJunEch[nJ,"CountJun"], "\t", expl, sep="")
						write(line, file=nfoname, append=TRUE)
						write(line, file="")
						matEV_RI = rbind(matEV_RI, cbind(newRI, as.matrix(rep("RI_5p+SE", nrow(newRI)))))
						
						if(sum(TabNFO[,"nJ"]%in%nomJ)==0)	TabNFO = rbind(TabNFO, c(gene, nomJ, Chr, geneJunEch[J,"StartJun"], geneJunEch[J,"EndJun"], expl, rep(0, ncol(TabNFO)-6)))
						TabNFO[TabNFO[,"nJ"]%in%nomJ,nom] = geneJunEch[nJ,"CountJun"]
					}
					
					#	SE+RI
					if(test_StartJ_StartX&test_Jend_in_intron)	#	(Start exon connu) ET (End in intron)
					{
						expl = "Start au debut d un exon et Retention intergenique"
						nomJ = paste(nJ, "=>SE+RI", sep="")
						line = paste(nomJ, "\t", gene, "\t", nomGene, "\t", "\t", geneJunEch[J,"StartJun"], "\t", geneJunEch[J,"EndJun"], "\t", geneJunEch[nJ,"CountJun"], "\t", expl, sep="")
						write(line, file=nfoname, append=TRUE)
						write(line, file="")
						matEV_RI = rbind(matEV_RI, cbind(newRI, as.matrix(rep("RI_5p", nrow(newRI)))))
						
						if(sum(TabNFO[,"nJ"]%in%nomJ)==0)	TabNFO = rbind(TabNFO, c(gene, nomJ, Chr, geneJunEch[J,"StartJun"], geneJunEch[J,"EndJun"], expl, rep(0, ncol(TabNFO)-6)))
						TabNFO[TabNFO[,"nJ"]%in%nomJ,nom] = geneJunEch[nJ,"CountJun"]
					}
					
					#	RI_5p
					if(test_Jstart_in_5p&test_Jend_in_5p&!test_EndConnu&!test_StartConnu)	#	(Start in 5p) ET (End in 5p)
					{
						expl = "Retention intergenique en 5p"
						nomJ = paste(nJ, "=>RI_5p", sep="")
						line = paste(nomJ, "\t", gene, "\t", nomGene, "\t", "\t", geneJunEch[J,"StartJun"], "\t", geneJunEch[J,"EndJun"], "\t", geneJunEch[nJ,"CountJun"], "\t", expl, sep="")
						write(line, file=nfoname, append=TRUE)
						write(line, file="")
						matEV_RI = rbind(matEV_RI, cbind(newRI, as.matrix(rep("RI_5p", nrow(newRI)))))
						
						if(sum(TabNFO[,"nJ"]%in%nomJ)==0)	TabNFO = rbind(TabNFO, c(gene, nomJ, Chr, geneJunEch[J,"StartJun"], geneJunEch[J,"EndJun"], expl, rep(0, ncol(TabNFO)-6)))
						TabNFO[TabNFO[,"nJ"]%in%nomJ,nom] = geneJunEch[nJ,"CountJun"]
					}
					
					#	RIG_5p
					if(test_Jstart_in_5p&test_EndConnu&test_Jend_EQ_Startexon1)	#	(Start in 5p) ET (End connu) ET (End = start 1er exon)
					{
						expl = "Retention intergenique en 5p"
						nomJ = paste(nJ, "=>RIG_5p", sep="")
						line = paste(nomJ, "\t", gene, "\t", nomGene, "\t", "\t", geneJunEch[J,"StartJun"], "\t", geneJunEch[J,"EndJun"], "\t", geneJunEch[nJ,"CountJun"], "\t", expl, sep="")
						write(line, file=nfoname, append=TRUE)
						write(line, file="")
						matEV_RI = rbind(matEV_RI, cbind(newRI, as.matrix(rep("RIG_5p", nrow(newRI)))))
						
						if(sum(TabNFO[,"nJ"]%in%nomJ)==0)	TabNFO = rbind(TabNFO, c(gene, nomJ, Chr, geneJunEch[J,"StartJun"], geneJunEch[J,"EndJun"], expl, rep(0, ncol(TabNFO)-6)))
						TabNFO[TabNFO[,"nJ"]%in%nomJ,nom] = geneJunEch[nJ,"CountJun"]
					}
					
					#	EL
					if(test_Jstart_in_exon&test_Jend_in_exon&test_JstartJend_Meme_exon)	#	(Start in exon) ET (End in exon) ET (Start et End dans le meme exon)
					{
						expl = "Exon Loss"
						nomJ = paste(nJ, "=>EL", sep="")
						line = paste(nomJ, "\t", gene, "\t", nomGene, "\t", ChrgeneJunEch[J,"StartJun"], "\t", geneJunEch[J,"EndJun"], "\t", geneJunEch[nJ,"CountJun"], "\t", expl, sep="")
						write(line, file=nfoname, append=TRUE)
						write(line, file="")
						matEV_RI = rbind(matEV_RI, cbind(newRI, as.matrix(rep("EL", nrow(newRI)))))
						
						if(sum(TabNFO[,"nJ"]%in%nomJ)==0)	TabNFO = rbind(TabNFO, c(gene, nomJ, Chr, geneJunEch[J,"StartJun"], geneJunEch[J,"EndJun"], expl, rep(0, ncol(TabNFO)-6)))
						TabNFO[TabNFO[,"nJ"]%in%nomJ,nom] = geneJunEch[nJ,"CountJun"]
					}
					
					#	A5SS+A3SS
					if(test_Jstart_in_exon&test_Jend_in_exon&(!test_JstartJend_Meme_exon))	#	(Start in exon) ET (End in exon) ET (Start et end dans 2 exons differents)
					{
						expl = "Exon1 raccourci en aval + Exon2 raccourci en amont"
						nomJ = paste(nJ, "=>A5SS+A3SS", sep="")
						line = paste(nomJ, "\t", gene, "\t", nomGene, "\t", "\t", geneJunEch[J,"StartJun"], "\t", geneJunEch[J,"EndJun"], "\t", geneJunEch[nJ,"CountJun"], "\t", expl, sep="")
						write(line, file=nfoname, append=TRUE)
						write(line, file="")
						matEV_RI = rbind(matEV_RI, cbind(newRI, as.matrix(rep("A5SS+A3SS", nrow(newRI)))))
						
						if(sum(TabNFO[,"nJ"]%in%nomJ)==0)	TabNFO = rbind(TabNFO, c(gene, nomJ, Chr, geneJunEch[J,"StartJun"], geneJunEch[J,"EndJun"], expl, rep(0, ncol(TabNFO)-6)))
						TabNFO[TabNFO[,"nJ"]%in%nomJ,nom] = geneJunEch[nJ,"CountJun"]
					}
					
					#	A5SS
					if(test_StartConnu&test_Jend_in_exon)	#	(start connu) ET (End in exon)
					{
						expl = "Exon2 raccourci en amont"
						nomJ = paste(nJ, "=>A5SS", sep="")
						line = paste(nomJ, "\t", gene, "\t", nomGene, "\t", "\t", geneJunEch[J,"StartJun"], "\t", geneJunEch[J,"EndJun"], "\t", geneJunEch[nJ,"CountJun"], "\t", expl, sep="")
						write(line, file=nfoname, append=TRUE)
						write(line, file="")
						matEV_RI = rbind(matEV_RI, cbind(newRI, as.matrix(rep("A5SS", nrow(newRI)))))
						
						if(sum(TabNFO[,"nJ"]%in%nomJ)==0)	TabNFO = rbind(TabNFO, c(gene, nomJ, Chr, geneJunEch[J,"StartJun"], geneJunEch[J,"EndJun"], expl, rep(0, ncol(TabNFO)-6)))
						TabNFO[TabNFO[,"nJ"]%in%nomJ,nom] = geneJunEch[nJ,"CountJun"]
					}
					
					#	A3SS
					if(test_Jstart_in_exon&test_EndConnu)	#	(start in exon) ET (End connu)
					{
						expl = "Exon1 raccourci en aval"
						nomJ = paste(nJ, "=>A3SS", sep="")
						line = paste(nomJ, "\t", gene, "\t", nomGene, "\t", "\t", geneJunEch[J,"StartJun"], "\t", geneJunEch[J,"EndJun"], "\t", geneJunEch[nJ,"CountJun"], "\t", expl, sep="")
						write(line, file=nfoname, append=TRUE)
						write(line, file="")
						matEV_RI = rbind(matEV_RI, cbind(newRI, as.matrix(rep("A3SS", nrow(newRI)))))
						
						if(sum(TabNFO[,"nJ"]%in%nomJ)==0)	TabNFO = rbind(TabNFO, c(gene, nomJ, Chr, geneJunEch[J,"StartJun"], geneJunEch[J,"EndJun"], expl, rep(0, ncol(TabNFO)-6)))
						TabNFO[TabNFO[,"nJ"]%in%nomJ,nom] = geneJunEch[nJ,"CountJun"]
					}
					
					#	RI_aval
					if(test_Jstart_in_intron&test_EndConnu)	#	(Start in intron) ET (End connu)
					{
						expl = "Jonction en aval de RI"
						nomJ = paste(nJ, "=>RI_aval", sep="")
						line = paste(nomJ, "\t", gene, "\t", nomGene, "\t", "\t", geneJunEch[J,"StartJun"], "\t", geneJunEch[J,"EndJun"], "\t", geneJunEch[nJ,"CountJun"], "\t", expl, sep="")
						write(line, file=nfoname, append=TRUE)
						write(line, file="")
						matEV_RI = rbind(matEV_RI, cbind(newRI, as.matrix(rep("RI_aval", nrow(newRI)))))
						
						if(sum(TabNFO[,"nJ"]%in%nomJ)==0)	TabNFO = rbind(TabNFO, c(gene, nomJ, Chr, geneJunEch[J,"StartJun"], geneJunEch[J,"EndJun"], expl, rep(0, ncol(TabNFO)-6)))
						TabNFO[TabNFO[,"nJ"]%in%nomJ,nom] = geneJunEch[nJ,"CountJun"]
					}
					
					#	RI_amont
					if(test_StartConnu&test_Jend_in_intron)	#	(Start connu) ET (End in intron)
					{
						expl = "Jonction en amont de RI"
						nomJ = paste(nJ, "=>RI_amont", sep="")
						line = paste(nomJ, "\t", gene, "\t", nomGene, "\t", "\t", geneJunEch[J,"StartJun"], "\t", geneJunEch[J,"EndJun"], "\t", geneJunEch[nJ,"CountJun"], "\t", expl, sep="")
						write(line, file=nfoname, append=TRUE)
						write(line, file="")
						matEV_RI = rbind(matEV_RI, cbind(newRI, as.matrix(rep("RI_amont", nrow(newRI)))))
						
						if(sum(TabNFO[,"nJ"]%in%nomJ)==0)	TabNFO = rbind(TabNFO, c(gene, nomJ, Chr, geneJunEch[J,"StartJun"], geneJunEch[J,"EndJun"], expl, rep(0, ncol(TabNFO)-6)))
						TabNFO[TabNFO[,"nJ"]%in%nomJ,nom] = geneJunEch[nJ,"CountJun"]
					}
					
					#	RIs
					if(test_Jstart_in_intron&test_Jend_in_intron&test_JstartJend_Meme_intron)	#	(Start in intron) ET (End in intron) ET (Start et End dans le meme intron)
					{
						expl = "Junction entre RI"
						nomJ = paste(nJ, "=>RIs", sep="")
						line = paste(nomJ, "\t", gene, "\t", nomGene, "\t", "\t", geneJunEch[J,"StartJun"], "\t", geneJunEch[J,"EndJun"], "\t", geneJunEch[nJ,"CountJun"], "\t", expl, sep="")
						write(line, file=nfoname, append=TRUE)
						write(line, file="")
						matEV_RI = rbind(matEV_RI, cbind(newRI, as.matrix(rep("RIs", nrow(newRI)))))
						
						if(sum(TabNFO[,"nJ"]%in%nomJ)==0)	TabNFO = rbind(TabNFO, c(gene, nomJ, Chr, geneJunEch[J,"StartJun"], geneJunEch[J,"EndJun"], expl, rep(0, ncol(TabNFO)-6)))
						TabNFO[TabNFO[,"nJ"]%in%nomJ,nom] = geneJunEch[nJ,"CountJun"]
					}
					
					#	SE+RIs
					if(test_Jstart_in_intron&test_Jend_in_intron&(!test_JstartJend_Meme_intron))	#	(Start in intron) ET (End in intron) ET (Start et End dans 2 introns differents) 
					{
						expl = "Junction introns / RIs avec saut exon1"
						nomJ = paste(nJ, "=>SE+RIs", sep="")
						line = paste(nomJ, "\t", gene, "\t", nomGene, "\t", "\t", geneJunEch[J,"StartJun"], "\t", geneJunEch[J,"EndJun"], "\t", geneJunEch[nJ,"CountJun"], "\t", expl, sep="")
						write(line, file=nfoname, append=TRUE)
						write(line, file="")
						matEV_RI = rbind(matEV_RI, cbind(newRI, as.matrix(rep("SE+RIs", nrow(newRI)))))
						
						if(sum(TabNFO[,"nJ"]%in%nomJ)==0)	TabNFO = rbind(TabNFO, c(gene, nomJ, Chr, geneJunEch[J,"StartJun"], geneJunEch[J,"EndJun"], expl, rep(0, ncol(TabNFO)-6)))
						TabNFO[TabNFO[,"nJ"]%in%nomJ,nom] = geneJunEch[nJ,"CountJun"]
					}
					
					#	A5SS+RI_amont
					if(test_Jstart_in_intron&test_Jend_in_exon&IE_contig)	#	(start intron) ET (End in exon) ET (IE contigus)
					{
						expl = "Jonction entre RI et exon raccourci en amont + start dans exon suivant end dans intron"
						nomJ = paste(nJ, "=>A5SS+RI_amont", sep="")
						line = paste(nomJ, "\t", gene, "\t", nomGene, "\t", "\t", geneJunEch[J,"StartJun"], "\t", geneJunEch[J,"EndJun"], "\t", geneJunEch[nJ,"CountJun"], "\t", expl, sep="")
						write(line, file=nfoname, append=TRUE)
						write(line, file="")
						matEV_RI = rbind(matEV_RI, cbind(newRI, as.matrix(rep("A5SS+RI_amont", nrow(newRI)))))
						
						if(sum(TabNFO[,"nJ"]%in%nomJ)==0)	TabNFO = rbind(TabNFO, c(gene, nomJ, Chr, geneJunEch[J,"StartJun"], geneJunEch[J,"EndJun"], expl, rep(0, ncol(TabNFO)-6)))
						TabNFO[TabNFO[,"nJ"]%in%nomJ,nom] = geneJunEch[nJ,"CountJun"]
					}
					
					#	A5SS+SE+RI_amont
					if(test_Jstart_in_intron&test_Jend_in_exon&(!IE_contig))	#	(start intron) ET (End in exon) ET (IE NON contigus)
					{
						expl = "Jonction entre RI et exon raccourci en amont avec start dans exon, end dans intron NON contigus"
						nomJ = paste(nJ, "=>A5SS+SE+RI_amont", sep="")
						line = paste(nomJ, "\t", gene, "\t", nomGene, "\t", "\t", geneJunEch[J,"StartJun"], "\t", geneJunEch[J,"EndJun"], "\t", geneJunEch[nJ,"CountJun"], "\t", expl, sep="")
						write(line, file=nfoname, append=TRUE)
						write(line, file="")
						matEV_RI = rbind(matEV_RI, cbind(newRI, as.matrix(rep("A5SS+SE+RI_amont", nrow(newRI)))))
						
						if(sum(TabNFO[,"nJ"]%in%nomJ)==0)	TabNFO = rbind(TabNFO, c(gene, nomJ, Chr, geneJunEch[J,"StartJun"], geneJunEch[J,"EndJun"], expl, rep(0, ncol(TabNFO)-6)))
						TabNFO[TabNFO[,"nJ"]%in%nomJ,nom] = geneJunEch[nJ,"CountJun"]
					}
					
					#	RI_5p+A5SS
					if(test_Jstart_in_5p&test_Jend_in_exon1)	#	(start dans 5p) ET (End dans 1er exon)
					{
						expl = "Junction entre 5p et Exon raccourci en amont"
						nomJ = paste(nJ, "=>RI_5p+A5SS", sep="")
						line = paste(nomJ, "\t", gene, "\t", nomGene, "\t", "\t", geneJunEch[J,"StartJun"], "\t", geneJunEch[J,"EndJun"], "\t", geneJunEch[nJ,"CountJun"], "\t", expl, sep="")
						write(line, file=nfoname, append=TRUE)
						write(line, file="")
						matEV_RI = rbind(matEV_RI, cbind(newRI, as.matrix(rep("RI_5p+A5SS", nrow(newRI)))))
						
						if(sum(TabNFO[,"nJ"]%in%nomJ)==0)	TabNFO = rbind(TabNFO, c(gene, nomJ, Chr, geneJunEch[J,"StartJun"], geneJunEch[J,"EndJun"], expl, rep(0, ncol(TabNFO)-6)))
						TabNFO[TabNFO[,"nJ"]%in%nomJ,nom] = geneJunEch[nJ,"CountJun"]
					}
					
					#	RI_5p+SE+A5SS
					if(test_Jstart_in_5p&test_Jend_in_exon&(!test_Jend_in_exon1))	#	(start dans 5p) ET (End dans exon) ET (End pas in exon1)
					{
						expl = "Junction entre 5p, saut exon et Exon raccourci en amont"
						nomJ = paste(nJ, "=>RI_5p+SE+A5SS", sep="")
						line = paste(nomJ, "\t", gene, "\t", nomGene, "\t", "\t", geneJunEch[J,"StartJun"], "\t", geneJunEch[J,"EndJun"], "\t", geneJunEch[nJ,"CountJun"], "\t", expl, sep="")
						write(line, file=nfoname, append=TRUE)
						write(line, file="")
						matEV_RI = rbind(matEV_RI, cbind(newRI, as.matrix(rep("RI_5p+SE+A5SS", nrow(newRI)))))
						
						if(sum(TabNFO[,"nJ"]%in%nomJ)==0)	TabNFO = rbind(TabNFO, c(gene, nomJ, Chr, geneJunEch[J,"StartJun"], geneJunEch[J,"EndJun"], expl, rep(0, ncol(TabNFO)-6)))
						TabNFO[TabNFO[,"nJ"]%in%nomJ,nom] = geneJunEch[nJ,"CountJun"]
					}
					
					#	A3SS+RI_aval
					if(test_Jstart_in_exon&test_Jend_in_intron&EI_contig)	#	(Start in exon) ET (End in intron) ET (EI contigus)
					{
						expl = "Exon raccourci en aval + RI + start dans exon suivant end dans intron"
						nomJ = paste(nJ, "=>A3SS+RI_aval", sep="")
						line = paste(nomJ, "\t", gene, "\t", nomGene, "\t", "\t", geneJunEch[J,"StartJun"], "\t", geneJunEch[J,"EndJun"], "\t", geneJunEch[nJ,"CountJun"], "\t", expl, sep="")
						write(line, file=nfoname, append=TRUE)
						write(line, file="")
						matEV_RI = rbind(matEV_RI, cbind(newRI, as.matrix(rep("A3SS+RI_aval", nrow(newRI)))))
						
						if(sum(TabNFO[,"nJ"]%in%nomJ)==0)	TabNFO = rbind(TabNFO, c(gene, nomJ, Chr, geneJunEch[J,"StartJun"], geneJunEch[J,"EndJun"], expl, rep(0, ncol(TabNFO)-6)))
						TabNFO[TabNFO[,"nJ"]%in%nomJ,nom] = geneJunEch[nJ,"CountJun"]
					}
					
					#	A3SS+SE+RI_aval
					if(test_Jstart_in_exon&test_Jend_in_intron&(!EI_contig))	#	(Start in exon) ET (End in intron) ET (EI NON contigus)
					{
						expl = "Exon raccourci en aval + RI + start dans exon, end dans intron NON contigus"
						nomJ = paste(nJ, "=>A3SS+SE+RI_aval", sep="")
						line = paste(nomJ, "\t", gene, "\t", nomGene, "\t", "\t", geneJunEch[J,"StartJun"], "\t", geneJunEch[J,"EndJun"], "\t", geneJunEch[nJ,"CountJun"], "\t", expl, sep="")
						write(line, file=nfoname, append=TRUE)
						write(line, file="")
						matEV_RI = rbind(matEV_RI, cbind(newRI, as.matrix(rep("A3SS+SE+RI_aval", nrow(newRI)))))
						
						if(sum(TabNFO[,"nJ"]%in%nomJ)==0)	TabNFO = rbind(TabNFO, c(gene, nomJ, Chr, geneJunEch[J,"StartJun"], geneJunEch[J,"EndJun"], expl, rep(0, ncol(TabNFO)-6)))
						TabNFO[TabNFO[,"nJ"]%in%nomJ,nom] = geneJunEch[nJ,"CountJun"]
					}
					
					#	A3SS+RI_3p
					if(test_Jstart_in_Derexon&test_Jend_in_3p)	#	(Start in Dernier exon) ET (End in 3p)
					{
						expl = "Dernier Exon raccourci en aval + fin de junction en 5p"
						nomJ = paste(nJ, "=>A3SS+RI_3p", sep="")
						line = paste(nomJ, "\t", gene, "\t", nomGene, "\t", "\t", geneJunEch[J,"StartJun"], "\t", geneJunEch[J,"EndJun"], "\t", geneJunEch[nJ,"CountJun"], "\t", expl, sep="")
						write(line, file=nfoname, append=TRUE)
						write(line, file="")
						matEV_RI = rbind(matEV_RI, cbind(newRI, as.matrix(rep("A3SS+RI_3p", nrow(newRI)))))
						
						if(sum(TabNFO[,"nJ"]%in%nomJ)==0)	TabNFO = rbind(TabNFO, c(gene, nomJ, Chr, geneJunEch[J,"StartJun"], geneJunEch[J,"EndJun"], expl, rep(0, ncol(TabNFO)-6)))
						TabNFO[TabNFO[,"nJ"]%in%nomJ,nom] = geneJunEch[nJ,"CountJun"]
					}
					
					#	SE+RI_3p
					if(test_StartConnu&test_Jend_in_3p)	#	(Start connu) ET (End in 3p)
					{
						expl = "Jonction prolongee en 3p avec saut exon"
						nomJ = paste(nJ, "=>SE+RI_3p", sep="")
						line = paste(nomJ, "\t", gene, "\t", nomGene, "\t", "\t", geneJunEch[J,"StartJun"], "\t", geneJunEch[J,"EndJun"], "\t", geneJunEch[nJ,"CountJun"], "\t", expl, sep="")
						write(line, file=nfoname, append=TRUE)
						write(line, file="")
						matEV_RI = rbind(matEV_RI, cbind(newRI, as.matrix(rep("SE+RI_3p", nrow(newRI)))))
						
						if(sum(TabNFO[,"nJ"]%in%nomJ)==0)	TabNFO = rbind(TabNFO, c(gene, nomJ, Chr, geneJunEch[J,"StartJun"], geneJunEch[J,"EndJun"], expl, rep(0, ncol(TabNFO)-6)))
						TabNFO[TabNFO[,"nJ"]%in%nomJ,nom] = geneJunEch[nJ,"CountJun"]
					}
					
					#	RI_5p+SE
					if(test_Jstart_in_5p&test_EndConnu)	#	(Start in 5p) ET (End connu)
					{
						expl = "Jonction prolongee en 5p avec saut exon"
						nomJ = paste(nJ, "=>RI_5p+SE", sep="")
						line = paste(nomJ, "\t", gene, "\t", nomGene, "\t", "\t", geneJunEch[J,"StartJun"], "\t", geneJunEch[J,"EndJun"], "\t", geneJunEch[nJ,"CountJun"], "\t", expl, sep="")
						write(line, file=nfoname, append=TRUE)
						write(line, file="")
						matEV_RI = rbind(matEV_RI, cbind(newRI, as.matrix(rep("RI_5p+SE", nrow(newRI)))))
						
						if(sum(TabNFO[,"nJ"]%in%nomJ)==0)	TabNFO = rbind(TabNFO, c(gene, nomJ, Chr, geneJunEch[J,"StartJun"], geneJunEch[J,"EndJun"], expl, rep(0, ncol(TabNFO)-6)))
						TabNFO[TabNFO[,"nJ"]%in%nomJ,nom] = geneJunEch[nJ,"CountJun"]
					}
					
					#	RI_5p+SE+RI
					if(test_Jstart_in_5p&test_Jend_in_intron)	#	(Start in 5p) ET (End in intron)
					{
						expl = "Jonction prolongée en 5p avec saut d’exon et RI"
						nomJ = paste(nJ, "=>RI_5p+SE+RI", sep="")
						line = paste(nomJ, "\t", gene, "\t", nomGene, "\t", "\t", geneJunEch[J,"StartJun"], "\t", geneJunEch[J,"EndJun"], "\t", geneJunEch[nJ,"CountJun"], "\t", expl, sep="")
						write(line, file=nfoname, append=TRUE)
						write(line, file="")
						matEV_RI = rbind(matEV_RI, cbind(newRI, as.matrix(rep("RI_5p+SE+RI", nrow(newRI)))))
						
						if(sum(TabNFO[,"nJ"]%in%nomJ)==0)	TabNFO = rbind(TabNFO, c(gene, nomJ, Chr, geneJunEch[J,"StartJun"], geneJunEch[J,"EndJun"], expl, rep(0, ncol(TabNFO)-6)))
						TabNFO[TabNFO[,"nJ"]%in%nomJ,nom] = geneJunEch[nJ,"CountJun"]
					}
					
					#	RIG_3p
					if(test_Jstart_EQ_EndDerexon&test_Jend_in_3p)	#	(Start = End dernier exon)) ET (End in 3p)
					{
						expl = "Retention intergenique 3p"
						nomJ = paste(nJ, "=>RIG_3p", sep="")
						line = paste(nomJ, "\t", gene, "\t", nomGene, "\t", "\t", geneJunEch[J,"StartJun"], "\t", geneJunEch[J,"EndJun"], "\t", geneJunEch[nJ,"CountJun"], "\t", expl, sep="")
						write(line, file=nfoname, append=TRUE)
						write(line, file="")
						matEV_RI = rbind(matEV_RI, cbind(newRI, as.matrix(rep("RIG_3p", nrow(newRI)))))
						
						if(sum(TabNFO[,"nJ"]%in%nomJ)==0)	TabNFO = rbind(TabNFO, c(gene, nomJ, Chr, geneJunEch[J,"StartJun"], geneJunEch[J,"EndJun"], expl, rep(0, ncol(TabNFO)-6)))
						TabNFO[TabNFO[,"nJ"]%in%nomJ,nom] = geneJunEch[nJ,"CountJun"]
					}
					
					#	RI_aval+SE+RIG_3p
					if(test_Jstart_in_intron&test_Jend_in_3p)	#	(Start in intron) ET (End in 3p)
					{
						expl = "Retention intron en amont, saut exon et Retention intergenique 3p"
						nomJ = paste(nJ, "=>RI_amont+SE+RIG_3p", sep="")
						line = paste(nomJ, "\t", gene, "\t", nomGene, "\t", "\t", geneJunEch[J,"StartJun"], "\t", geneJunEch[J,"EndJun"], "\t", geneJunEch[nJ,"CountJun"], "\t", expl, sep="")
						write(line, file=nfoname, append=TRUE)
						write(line, file="")
						matEV_RI = rbind(matEV_RI, cbind(newRI, as.matrix(rep("RI_aval+SE+RIG_3p", nrow(newRI)))))
						
						if(sum(TabNFO[,"nJ"]%in%nomJ)==0)	TabNFO = rbind(TabNFO, c(gene, nomJ, Chr, geneJunEch[J,"StartJun"], geneJunEch[J,"EndJun"], expl, rep(0, ncol(TabNFO)-6)))
						TabNFO[TabNFO[,"nJ"]%in%nomJ,nom] = geneJunEch[nJ,"CountJun"]
					}
					
					#	Ajoute le nouvel element a la carte
					Map = rbind(Map, c(paste("Start_", nJ, sep=""), J_anorm[J,"StartJun"], matEV_RI[nrow(matEV_RI),"Type"]), c(paste("End_", nJ, sep=""), J_anorm[J,"EndJun"], matEV_RI[nrow(matEV_RI),"Type"]))
				}
				write.table(resNFOjun, file=paste(pathResGenePat, "/", nom, "_nfo_J_Gene.txt", sep=""), sep="\t", quote=FALSE, row.names=FALSE)
				
				if(nrow(matEV_RI)>0)
				{
					IsResRI = TRUE
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
					
					MoyNcountParJunc = matrix(c(rep(sum(as.numeric(geneJunEch[,"CountJun"]))/maxNumbJ, nrow(matEV_RI)),rep(sum(as.numeric(geneJunEch[,"CountJun"]))/minNumbJ, nrow(matEV_RI))), ncol=2, byrow=FALSE)
					colnames(MoyNcountParJunc)=c("minMoyNcountParJunc", "maxMoyNcountParJunc")
					matEV_RI = cbind(matEV_RI, MoyNcountParJunc)
					matEV_RI = cbind(matEV_RI, rep(nom, nrow(matEV_RI)))
					colnames(matEV_RI)[ncol(matEV_RI)] = "Sample"
					write.table(matEV_RI, file=paste(pathResGenePat, "/", nom, "_RI.txt", sep=""), sep="\t", quote=FALSE, row.names=FALSE)
					saveRDS(matEV_RI, file=paste(pathResRI, "/", nom, "/", nom, "_", gene, "_RI.rds", sep=""))
				}
			}
			
			###########################################################################################
			#	SE
			#			
			cols = c("Junc_SE", "ENSID", "Gene", "Chr", "Junc_SE_Start", "Junc_SE_End", "Junc_SE_Count", "NbreSkippedExons", "SkippedExonsID", "TranscritsID", "Junc_Normale", "Junc_Normale_Count")
			matEV_SE = matrix(0, nrow=0, ncol=length(cols))
			colnames(matEV_SE) = cols
			
			if(nrow(geneJunEch)>0)	for(J in 1:nrow(geneJunEch))
			{
				#	Liste des exons inclus dans la jonction
				SE_StartEnd = (geneRef[,"exon_chrom_start"]>=(as.numeric(geneJunEch[J,"StartJun"])-limDecal))&(geneRef[,"exon_chrom_end"]<=(as.numeric(geneJunEch[J,"EndJun"])+limDecal))
				
				if(sum(SE_StartEnd)>0)	#	Saut Exon
				{
					resSkippedExons = geneRef[SE_StartEnd,, drop=FALSE]
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
					
					coordsJ = rownames(geneJunEch)[J]
					nomJ = paste(coordsJ, "=>SE_", NSkippedExons, sep="")
					expl = paste("Saute ", NSkippedExons, " exon(s) => Transcrit Alternatif Connu", sep="")
					
					if(sum(AllJunc%in%coordsJ)<1)	expl = paste("Saute ", sum(SE_StartEnd), " exon(s)", sep="") #	SE inconnu, anormal
					
					line = paste(nomJ, "\t", geneJunEch[J,"StartJun"], "\t", geneJunEch[J,"EndJun"], "\t", geneJunEch[J,"CountJun"], "\tSaute ", NSkippedExons, " exon(s) touche ", length(tmpTranscrit), " transcrit(s)", sep="")
					write(line, file="")
					write(line, file=nfoname, append=TRUE)
					EV_SE = geneRef[SE_StartEnd, c("ensembl_transcript_id", "ensembl_exon_id", "exon_chrom_start", "exon_chrom_end"), drop=FALSE]
					
					if(sum(TabNFO[,"nJ"]%in%nomJ)==0)	#	Junc inconnue
					{
						write(paste("Jonction inconnue (", nomJ, ") => ajoute une ligne de ", length(c(nomJ, geneJunEch[J,"StartJun"], geneJunEch[J,"EndJun"], expl, rep(0, ncol(TabNFO)-6))), " elements", sep=""), file="")
						
						TabNFO = rbind(TabNFO, c(gene, nomJ, Chr, geneJunEch[J,"StartJun"], geneJunEch[J,"EndJun"], expl, rep(0, ncol(TabNFO)-6)))
						TabNFO[TabNFO[,"nJ"]%in%nomJ, nom] = geneJunEch[coordsJ,"CountJun",drop=FALSE]
					}else{
						TabNFO[TabNFO[,"nJ"]%in%nomJ,nom] = geneJunEch[coordsJ,"CountJun",drop=FALSE]	#	geneJun[nJ,"CountJun"]
					}
					
					#	Liste les Junctions canoniques inclues dans le SE
					canonGeneJun = geneJunEch[rownames(geneJunEch)%in%rownames(geneJuncREF),,drop=FALSE]
					test_JcanonInclue = ((as.numeric(geneJunEch[J,"StartJun"])-limDecal)<=as.numeric(canonGeneJun[,"StartJun"]))&((as.numeric(geneJunEch[J,"EndJun"])+limDecal)>=as.numeric(canonGeneJun[,"EndJun"]))	#|	#	Junc canonique incluse dans le SE
					
					maxJuncNorm = 0
					if(sum(test_JcanonInclue)>0)
					{
						canonJuncInc = canonGeneJun[test_JcanonInclue,,drop=FALSE]
						
						if(sum(rownames(geneJunEch)%in%rownames(canonJuncInc))>0)	#	Junc Canonique inclue avec comptage
						{
							tmpJuncNormCounts = geneJunEch[rownames(geneJunEch)%in%rownames(canonJuncInc),,drop=FALSE]
							maxJuncNorm = max(as.numeric(tmpJuncNormCounts[,"CountJun"]))
							nomMaxJuncNorm = rownames(tmpJuncNormCounts)[as.numeric(tmpJuncNormCounts[,"CountJun"])==maxJuncNorm]
						}else{#	pas de comptages pour les junctions canoniques inclues dans le SE => 0
							maxJuncNorm = 0
							nomMaxJuncNorm = rownames(canonJuncInc)[1]
						}
						if(length(nomMaxJuncNorm)>1)	nomMaxJuncNorm = paste(nomMaxJuncNorm, collapse=",")
						matEV_SE = rbind(matEV_SE, c(rownames(geneJunEch)[J], gene, nomGene, Chr, geneJunEch[J,c("StartJun", "EndJun", "CountJun"),drop=FALSE], NSkippedExons, nomSkippedExons, nomTranscrits, nomMaxJuncNorm, maxJuncNorm))
					}else{#	Pas de junction normale / canonique dans le SE
						matEV_SE = rbind(matEV_SE, c(rownames(geneJunEch)[J], gene, nomGene, Chr, geneJunEch[J,c("StartJun", "EndJun", "CountJun"),drop=FALSE], NSkippedExons, nomSkippedExons, nomTranscrits, "Pas_de_jonction_canonique_inclue_dans_le_SE", "0"))
					}
					
					Map = rbind(Map, c(paste("Start_", coordsJ, sep=""), geneJunEch[coordsJ,"StartJun"], "SE"), c(paste("End_", coordsJ, sep=""), geneJunEch[coordsJ,"EndJun"], "SE"))
				}
			}
			Map = Map[order(as.numeric(Map[,"Pos"])),]
			write.table(Map, file=paste(pathResGenePat, "/", nom, "_Map_J.txt", sep=""), sep="\t", quote=FALSE, row.names=FALSE)
			
			if(nrow(matEV_SE)>0)
			{
				IsResSE = TRUE
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
				matEV_SE[matEV_SE[,"Junc_SE"]%in%AllJunc,"Type"] = "SE_Canonique"
				
				#	Nombre de reads junc / nombre max de junctions pour le gene
				maxNumbJ = 0
				minNumbJ = length(unlist(JuncT))
				for(V in 1:length(JuncT))
				{
					maxNumbJ = max(maxNumbJ, length(JuncT[[V]]))
					minNumbJ = min(minNumbJ, length(JuncT[[V]]))
				}		
				MoyNcountParJunc = matrix(c(rep(sum(as.numeric(geneJunEch[,"CountJun"]))/maxNumbJ, nrow(matEV_SE)),rep(sum(as.numeric(geneJunEch[,"CountJun"]))/minNumbJ, nrow(matEV_SE))), ncol=2, byrow=FALSE)
				colnames(MoyNcountParJunc)=c("minMoyNcountParJunc", "maxMoyNcountParJunc")
				matEV_SE = cbind(matEV_SE, MoyNcountParJunc)
				tmpmatEV_SE = cbind(matEV_SE[,1], matEV_SE[,2:ncol(matEV_SE), drop=FALSE])
				colnames(tmpmatEV_SE)[1] = "Junc_SE"
				
				tmpmatEV_SE = cbind(tmpmatEV_SE, rep(nom, nrow(tmpmatEV_SE)))
				colnames(tmpmatEV_SE)[ncol(tmpmatEV_SE)] = "Sample"
				write.table(tmpmatEV_SE, file=paste(pathResGenePat, "/", nom, "_SE.txt", sep=""), sep="\t", quote=FALSE, row.names=FALSE)
				saveRDS(tmpmatEV_SE, file=paste(pathResSE, "/", nom, "/", nom, "_", gene, "_SE.rds", sep=""))
			}
			
			if((!IsResSE)|(!IsResRI))
			{
				#	Pas de res, suppression du dossier gène ?
				unlink(pathResGenePat, recursive=TRUE)
			}else{
				write.table(TabNFO, file=paste(pathResGenePat, "tableNFOjunc.txt", sep=""), sep="\t", quote=FALSE, row.names=FALSE)
				#	saveRDS(TabNFO, file=paste(pathResGenePat, "tableNFOjunc.rds", sep="")
			}
		}else{
			write(paste("\n################################################################################################################################", sep=""), file="")
			write(paste("################################################################################################################################", sep=""), file="")
			write(paste("#\t\t", nom, " / ", gene, "\t\t\t\t#", sep=""), file="")
			write(paste("################################################################################################################################", sep=""), file="")
			write(paste("################################################################################################################################\n", sep=""), file="")
		}
	}
}

###############################################################################################################
#	VII	Condenser les résultats
condenseResPat<-function(R)
{
	AllResRDSfiles = list.files(paste(pathRes, "/Allres", tabResCondens[R,2], "/", tabResCondens[R,1], "/", sep=""), full.names=TRUE)
	
	if(length(AllResRDSfiles)>0)
	{
		allRes = list()
		for(N in 1:length(AllResRDSfiles)){allRes[[length(allRes)+1]] = readRDS(AllResRDSfiles[N]); write(paste(tabResCondens[R,2], "_", tabResCondens[R,1], " (", N, "/", length(AllResRDSfiles), ")", sep=""), file="")}
		allRes = do.call("rbind",allRes)
		allRes = allRes[order(allRes[,min(grep("_Start", colnames(allRes)))]),,drop=FALSE]
		allRes = allRes[order(allRes[,grep("^Chr$", colnames(allRes))]),,drop=FALSE]
		
		#evs = apply(allRes, 1, function(x) paste(x[1], x[length(x)], sep="_", collapse="_"))
		#evs[duplicated(evs)]
		#allRes = allRes[!duplicated(evs),,drop=FALSE]
		
		#write.table(allRes, file=paste(AllresPath, "allRes_", tabResCondens[R,2], "_", tabResCondens[R,1], ".txt", sep=""), sep="\t", quote=FALSE, row.names=FALSE)
		saveRDS(allRes, file=paste(AllresPath, "allRes_", tabResCondens[R,2], "_", tabResCondens[R,1], ".rds", sep=""))
		write.table(allRes, file=paste(pathAllResAnalysis, "allRes_", tabResCondens[R,2], "_", tabResCondens[R,1], ".txt", sep=""), sep="\t", quote=FALSE, row.names=FALSE)
	}else{
		write(paste("\n/Allres", tabResCondens[R,2], "/", tabResCondens[R,1], "/ est vide !", sep=""), file="")
	}
}

###############################################################################################################
#	Stats

Stats<-function(idprojet, esp, nCPU, align)
{
	pathRes = paste("/data-isilon/sequencing/ngs/", idprojet, "/", esp, "/RNAseqSEA/", sep="")
	pathBamJuncs = paste(pathRes, "/BamJuncs/", sep="")
	
	resStatsPath = paste(pathRes, "/Stats/", sep="")
	if(!file.exists(resStatsPath))	dir.create(resStatsPath)
	
#	nreads / ech
	bams = list.files(paste("/data-isilon/sequencing/ngs/", idprojet, "/", esp, "/align/", align, "/", sep=""), full.names=TRUE, pattern="\\.bam$")
	nMappedReads = NULL
	for(B in 1:length(bams))
	{
		write(paste("#\t", B, "/", length(bams), "\tcomptage des Mapped reads de ", basename(bams[B]), sep=""), file="")
		nMappedReads = c(nMappedReads, system(paste("samtools view -@ ", nCPU, " -c -F 260 ", bams[B], sep=""), intern=TRUE))
	}
	if(!is.null(ReName))	names(nMappedReads) = BamsReNames[BamsReNames[,1]%in%gsub(".bam$", "", basename(bams)),2]
	if(is.null(ReName))	names(nMappedReads) = gsub(".bam", "", basename(bams))
	saveRDS(as.matrix(nMappedReads), paste(resStatsPath, "/nMappedReads_ech.rds", sep=""))
	
#	n ReadsJunc / ech / chr 
	bamsJuncs = list.files(pathBamJuncs, full.names=TRUE, pattern="\\.bam$")
	
	nJuncMappedReads = NULL
	for(B in 1:length(bamsJuncs))
	{
		write(paste("#\t", B, "/", length(bamsJuncs), "\tcomptage des Mapped reads de ", basename(bamsJuncs[B]), sep=""), file="")
		nJuncMappedReads = c(nJuncMappedReads, system(paste("samtools view -@ ", nCPU, " -c -F 260 ", bamsJuncs[B], sep=""), intern=TRUE))
	}
	names(nJuncMappedReads) = gsub(".bam$", "", basename(bamsJuncs))
	saveRDS(as.matrix(nJuncMappedReads), paste(resStatsPath, "/nJuncMappedReads_ech_chr.rds", sep=""))
	
	nJuncMappedReadsEch = NULL
	for(S in 1:length(samples))	nJuncMappedReadsEch = c(nJuncMappedReadsEch, sum(as.numeric(nJuncMappedReads[grepl(paste("^", samples[S], "_", sep=""), names(nJuncMappedReads))])))
	
	names(nJuncMappedReadsEch) = samples
	saveRDS(as.matrix(nJuncMappedReadsEch), paste(resStatsPath, "/nJuncMappedReads_ech.rds", sep=""))
	
#	RI SE
	AllresPath = paste(pathRes, "/AllRes/", sep="")
	AllRes_RI = list.files(AllresPath, full.names=TRUE, pattern="\\.rds$")
	
	RISE = matrix(0, ncol=4, nrow=length(samples))
	rownames(RISE) = samples
	colnames(RISE) = c("nSE", "nSE_Genes", "nRI", "nRI_Genes")
	for(S in 1:length(samples))
	{
		tmpSE = readRDS(paste(AllresPath, "allRes_SE_", samples[S], ".rds", sep=""))
		RISE[S, 1] = nrow(tmpSE)
		RISE[S, 2] = length(unique(tmpSE[,"ENSID"]))
		
		tmpRI = readRDS(paste(AllresPath, "allRes_RI_", samples[S], ".rds", sep=""))
		RISE[S, 3] = nrow(tmpRI)
		RISE[S, 4] = length(unique(tmpRI[,"ENSID"]))
	}
	saveRDS(RISE, paste(resStatsPath, "/RISE_ech.rds", sep=""))
	
	stats = cbind(nMappedReads, nJuncMappedReadsEch, RISE)
	saveRDS(stats, paste(resStatsPath, "/stats_ech.rds", sep=""))
	
	write.table(stats, file=paste(resStatsPath, "/RISE_ech.txt", sep=""), sep="\t", quote=FALSE)
}

slogRNAseqSEA_AllTnjs<-function(nfos)
{
  noms = c("idprojet", "esp", "align", "nodename", "login", "user", "effective_user", "wd")
  allnfos = matrix(0, ncol=length(noms), nrow=0)
  colnames(allnfos) = noms
  rdsPath = "/data-isilon/Cagnard/RNAseqSEA/RNAseqSEA_AllTnjs_logs.rds"
  if(file.exists(rdsPath)) allnfos = readRDS(rdsPath)
  allnfos = rbind(allnfos, nfos)
  saveRDS(allnfos, rdsPath)
}

GenVer<-function(idprojet)
{
	cfg = cfgParse()
	library(RMySQL)
	con <- dbConnect(MySQL(), user="polyweb", password=cfg$polyprod["pw",], dbname="PolyprojectNGS", host="10.200.27.108")	
	res<-dbSendQuery(con, paste('SELECT `releases`.`name` FROM `PolyprojectNGS`.`releases`, `PolyprojectNGS`.`project_release`, `PolyprojectNGS`.`projects` WHERE `releases`.`release_id`= `project_release`.`release_id`AND `projects`.`name` = "', idprojet, '" AND `project_release`.`project_id`=`projects`.`project_id`;', sep=""))
	on.exit(dbDisconnect(con))
	type = fetch(res, n=-1)
	return(type)
}

cfgParse<-function()
{
	if(file.exists("/software/polyweb/poly-src/GenBo/lib/obj-nodb/genbo.cfg"))
	{
		con <- file("/software/polyweb/poly-src/GenBo/lib/obj-nodb/genbo.cfg",open="r")
	}else{
		viaScripts = "/data-isilon/Cagnard/MagicMorgan/Methode_RNAseq_DevL/"
		con <- file(paste(viaScripts, "/genbo.cfg", sep=""),open="r")
	}
	line <- readLines(con)
	close(con)
	
	line = line[grepl('\\[', line)|grepl('\\:', line)]	#	Selectionne les lignes utiles
	Sections = c(grep('\\[', line), length(line)+1)
	
	cfg = list()
	noms = c()
	for(N in 1:(length(Sections)-1))
	{
		section = line[Sections[N]:(Sections[N+1]-1)]
		nom = gsub("\\[", "", section[1])
		nom = gsub("\\]", "", nom)
		noms = c(noms, nom)
		
		section = t(rbind(apply(as.matrix(section[2:length(section), drop=FALSE]), 1, function(x)	as.matrix(c(unlist(strsplit(x, '\\:')), "")[1:2]))))
		rownames(section) = section[,1]
		section = section[,2, drop=FALSE]
		colnames(section) = nom
		
		cfg[[length(cfg)+1]] = section
	}
	names(cfg) = noms
	
	return(cfg)
}

###############################################################################################################
#	Formatage des res

formatResAll<-function(pathAllResAnalysis)
{
  write(paste("#\tFormatage des resultats", sep=""), file="")
  fichs = list.files(pathAllResAnalysis, full.names=TRUE, recursive=FALSE)
	
	fichs_resAll_SE = fichs[grepl("allResSE.rds$", fichs)]
	fichs_resAll_RI = fichs[grepl("allResRI.rds$", fichs)]

	if(length(fichs_resAll_SE)>0)
	{
	  write(paste("#\t\t Compression et indexation du resAll SE", sep=""), file="")
	  allResSE = readRDS(fichs_resAll_SE)
		allResSE[,"Junc_SE_Start"] = format(as.numeric(allResSE[,"Junc_SE_Start"]), scientific=FALSE)
		allResSE[,"Junc_SE_End"] = format(as.numeric(allResSE[,"Junc_SE_End"]), scientific=FALSE)
		allResSE[,"Junc_SE_Count"] = as.numeric(allResSE[,"Junc_SE_Count"])
		allResSE[,"Junc_Normale_Count"] = as.numeric(allResSE[,"Junc_Normale_Count"])

		allResSE = allResSE[order(as.numeric(allResSE[,"Junc_SE_Start"])),,drop=FALSE]
		allResSE = allResSE[order(allResSE[,"Chr"]),,drop=FALSE]
		
		#colnames(allResSE)[1] = paste("#", colnames(allResSE)[1], sep="")
		write.table(allResSE, file=paste(pathAllResAnalysis, "allResSE.txt", sep=""), sep="\t", quote=FALSE, row.names=FALSE)
		#system(paste("bgzip ", pathAllResAnalysis, "allResSE.txt", sep=""))
		#system(paste("tabix -S 1 -s ", grep("Chr", colnames(allResSE)), " -b ", grep("Junc_SE_Start", colnames(allResSE)), " -e ", grep("Junc_SE_End", colnames(allResSE)), " ", pathAllResAnalysis, "allResSE.txt.gz", sep=""))
		#unlink(fichs_resAll_SE)
	}
	
	if(length(fichs_resAll_RI)>0)
	{
	  write(paste("#\t\t Compression et indexation du resAll RI", sep=""), file="")
	  allResRI = readRDS(fichs_resAll_RI)
		allResRI[,"Junc_RI_Start"] = format(as.numeric(allResRI[,"Junc_RI_Start"]), scientific=FALSE)
		allResRI[,"Junc_RI_End"] = format(as.numeric(allResRI[,"Junc_RI_End"]), scientific=FALSE)
		allResRI[,"Junc_RI_Count"] = as.numeric(allResRI[,"Junc_RI_Count"])
		allResRI[,"Junc_Normale_Count"] = as.numeric(allResRI[,"Junc_Normale_Count"])
		
		allResRI = allResRI[order(as.numeric(allResRI[,"Junc_RI_Start"])),]
		allResRI = allResRI[order(allResRI[,"Chr"]),]
		
		#colnames(allResRI)[1] = paste("#", colnames(allResRI)[1], sep="")
		write.table(allResRI, file=paste(pathAllResAnalysis, "allResRI.txt", sep=""), sep="\t", quote=FALSE, row.names=FALSE)
		#system(paste("bgzip ", pathAllResAnalysis, "allResRI.txt", sep=""))
		#system(paste("tabix -S 1 -s ", grep("Chr", colnames(allResRI)), " -b ", grep("Junc_RI_Start", colnames(allResRI)), " -e ", grep("Junc_RI_End", colnames(allResRI)), " ", pathAllResAnalysis, "allResRI.txt.gz", sep=""))
		#unlink(fichs_resAll_RI)
	}
	
#	system(paste("bgzip ", pathAllResAnalysis, "allResRI.txt", sep=""))
#	system(paste("bgzip ", pathAllResAnalysis, "allResSE.txt", sep=""))
	
	write(paste("#\t\t Compression et suppression res RI / SE par sample", sep=""), file="")
	resPat = list.files(pathAllResAnalysis, full.names=TRUE, recursive=FALSE, pattern="\\.txt$")
	resPatPath = paste(pathAllResAnalysis, "/resPat/", sep="")
	dir.create(resPatPath)
	setwd(pathAllResAnalysis)
	for(L in 1:length(resPat))
	{
		system(paste("tar -zcvf ", resPatPath, "/", basename(resPat[L]), ".gz ./", basename(resPat[L]), sep=""))
		unlink(resPat[L], recursive=TRUE)
	}
	
	write(paste("#\t\t Suppression des split bams", sep=""), file="")
	unlink(splitBamPath, recursive=TRUE)
	write(paste("#\t\t Suppression des Junc", sep=""), file="")
	unlink(AllJuncPath, recursive=TRUE)
}




