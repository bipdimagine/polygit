#!/bin/bash
echo $1
add_calling_method.sh -project=$1 -method=canvas,manta,wisecondor
../../../bds_cache.pl -project=$1 -yes=1 -force=1 > $1.log

