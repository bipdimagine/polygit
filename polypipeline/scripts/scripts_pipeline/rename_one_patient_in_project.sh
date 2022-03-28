#!/bin/bash


####exemple de ligne de commande pour l'execution : /bip-d/perl/GenBo/script/ngs_exome/pipeline/rename_one_patient_in_project.sh NGS2014_0606 abc ABC

#source /bip-d/soft/etc/analysis.cfg


project=$1
patient_actual=$2
patient_new=$3


perl /bip-d/perl/GenBo/script/ngs_exome/pipeline/rename_one_patient_in_project.pl -project=${project} -patient_actual=${patient_actual} -patient_new=${patient_new}


###check if files has been correctly moved

DIR_PROJECT=/data-isilon/sequencing/ngs

cd $DIR_PROJECT/${project}/HG19
check=$(find -name ${patient_actual}* |wc -l)

if [ $check != 0 ] ;
	then
	echo "ERROR : $check files not replaced : " 
	echo $(find -name ${patient_actual}*)
else
	echo "OK : All files replaced."
fi









