#!/bin/zsh +x
neventPerFile=500
for e in 10;do
#for e in 3 4 5 6;do
#for e in 3 4 5 6 8 10 15;do
#for e in 3;do
#for e in 8;do
#for e in 6 8;do
#for e in 8 10 15;do
#for e in 2;do
#    for i in {0..0};do
    for i in {0..49999};do
#    for i in {50000..99999};do
#    for i in {0..9999};do
	f=$((${i}/$neventPerFile))
	q=$(($i - ${f}*${neventPerFile}))
	if [ $q -eq 0 ];
	then
	    echo $i
	    #	    qsub -q all -e ${LEAFDIR}/jobs/fitLE_e${e}_start${i}.err -o ${LEAFDIR}/jobs/fitLE_e${e}_start${i}.log ${LEAFDIR}/jobs/fitLE_e${e}_start${i}.sh
	    #mv /disk01/usr5/bquilain/LEAF_hkhybridmpmt10pc14374100Hz_4.2kHzbl_10MeV_tmp${i}.root /disk01/usr5/bquilain/LEAF_hkhybridmpmt10pc14374100Hz_4.2kHzbl_10MeV_tmpdir_${i}.root 
	    qsub -q all -e ${LEAFDIR}/jobs/fitLE_e${e}_start${i}.err -o ${LEAFDIR}/jobs/fitLE_e${e}_start${i}.log ${LEAFDIR}/jobs/fitLE_e${e}_start${i}.sh 
	fi
    done
done
