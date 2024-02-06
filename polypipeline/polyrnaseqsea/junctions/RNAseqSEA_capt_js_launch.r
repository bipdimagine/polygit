
# /software/bin/run_cluster.pl -cpu=40 -cmd="Rscript  /data-isilon/Cagnard/RNAseqSEA/RNAseqSEA_capt_js/RNAseqSEA_capt_js_launch.r idprojet=NGS2023_7143"
# /software/polyweb/poly-disk/poly-src/polygit/polymorphism-cgi/rnaseq/create_config_splices_analyse_file.pl -project=NGS2023_6989 force=FALSE

# idprojet="NGS2021_3899"

 
args = commandArgs(trailingOnly=TRUE)
tmpArgs = t(rbind(apply(as.matrix(unlist(args)), 1, function(x) unlist(strsplit(x, "=")))))
if(sum(grepl("^idprojet$", tmpArgs[,1]))>0)	idprojet = tmpArgs[grep("^idprojet$", tmpArgs[,1], ignore.case=TRUE),2]

nCPUmax = 4
limDecal=1
CorrAnnot=TRUE
# if(file.exists())

write(paste("\t#\t1- Cree le fichier de config js", sep=""), file="")
#configPath = system(paste("/software/polyweb/poly-disk/poly-src/polygit/polymorphism-cgi/rnaseq/create_config_splices_analyse_file.pl -force=1 -project=", idprojet, sep=""), intern=TRUE)
configPath = system(paste("/data-isilon/bipd-src/mbras/git_repository/polymorphism-cgi/rnaseq/create_config_splices_analyse_file.pl -force=1 -project=", idprojet, sep=""), intern=TRUE)

if(file.exists(configPath))
{
  write(paste("\t#\t2- Parse le fichier de config et recupere les infos", sep=""), file="")
  library("rjson")
  allNFO <- fromJSON(file = configPath)
  
  scriptPath = allNFO$rnaseqseapaths$script
  biblioPath = allNFO$rnaseqseapaths$biblio
  
  # AnnotPath = "/Users/nicolas/Documents/MagicMorgan/Methode_RNAseq_DevL/Annots/AllEns.Rds"
  AnnotPath = "/data-isilon/Cagnard/MagicMorgan/Methode_RNAseq_DevL/Annots/AllEns.Rds"
  genelistPath = "/software/polyweb/poly-disk/poly-src/polygit/polyscripts/polyrnaseq/genes_list.pl"
  runcluster = "/software/bin/run_cluster.pl"
  
  Annots = readRDS(AnnotPath)
  source(biblioPath)
  
  esp = allNFO$project$gencode$release
  gcvers = allNFO$project$gencode$version
  resPath = dirname(allNFO$analyse$path)
  align = basename(dirname(unlist(allNFO$inputs$bam)[1]))
  gcrds = allNFO$project$gencode$rds
  gtf = allNFO$project$gencode$gtf
  juncrds = allNFO$project$gencode$junctions_canoniques_rds
  bamfiles = unlist(allNFO$inputs$bam)
  bedtoolsPath = allNFO$softwares$bedtools
  sambambaPath = allNFO$softwares$sambamba
  samtoolsPath = allNFO$softwares$samtools
  picardPath = allNFO$softwares$picard
  capture = allNFO$analyse$capture
  
  write(paste("\t#\t3- Recupere la liste de capture", sep=""), file="")
  try(ENSg(idprojet, esp, biblioPath, resPath), silent=TRUE)
  
  if(sum(grepl("^ENSg*.*.txt", list.files(paste(resPath, "/RNAseqSEA/", sep=""))))==1)
  {
    write(paste("\t#\t\tListe de capture ok", sep=""), file="")
    write(paste("\t#\t4- Valide la liste de capture", sep=""), file="")
    testENSgGC(gcrds, resPath)
    
    write(paste("\t#\t4- Cree les fichier samplefiles des comparaisons", sep=""), file="")
    formatSamplesTypes_1vsAll(idprojet, biblioPath, scriptPath, esp, gcvers, align, nCPUmax=4, limDecal=1, bamfiles, gcrds, resPath, sambambaPath, samtoolsPath, picardPath)
    
    write(paste("\t#\t5- Exectute les commandes du cmd.txt", sep=""), file="")
	cmd = paste("sh ", resPath, "/RNAseqSEA/cmds.sh", sep="");
    system(cmd, wait = TRUE);
	cmd = paste("cat ", resPath, "/RNAseqSEA/cmds_all.sh | ", runcluster," -limit=300 -cpu=", nCPUmax, " ", sep="");
    system(cmd, wait = TRUE);
    
    write(paste("\t#\t6- Formate les resultats", sep=""), file="")
    formatProjet(resPath, align)
    
    write(paste("\t#\t7- Supprime les fichiers temporaires", sep=""), file="")
    tmp = list.files(paste(resPath, "/analysis/", sep=""), pattern = "^RNAseqSEA_", recursive = FALSE, full.names=TRUE)
    runif(1)
    library(parallelMap)
    parallelStart(mode = "multicore", cpus=nCPUmax, show.info=TRUE) 
    f = function(resGenes) unlink(paste(resGenes, "/resGenes/", sep=""), recursive = TRUE)
    y = parallelMap(f, tmp)
    parallelStop()
  }else{
    write.table(paste("La capture n'est pas renseignee dans la base: Placer la liste des ID Ensembl dans le fichier ", resPath, "/RNAseqSEA/ENSg.txt", sep=""), file=paste(resPath, "/RNAseqSEA/ReadMe.txt", sep=""), sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)
    write(paste("Placer la liste des ID Ensembl dans le fichier: ", resPath, "/RNAseqSEA/ENSg.txt", sep=""), file="")
  }
}else{
  write(paste("Fichier de config js introuvable: ", configPath, sep=""), file="")
}
