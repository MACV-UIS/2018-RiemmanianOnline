
BASE=`pwd`"/"
LIB="/home/macv/DATOS/gustavo/libsvm-3.22/tools/"

SALIDA=$BASE"res-experimento03.txt"
TEMPORAL=$LIB"temporal/"

DIR="/home/macv/DATOS/fabio/2018-RiemmanianOnline/c++_Code/RecursiveMesCov_KTH/resultados/03_comp_Theta/"

cp $DIR"Gtest_concMes_IFs_24_mean_var_selected.txt" $TEMPORAL"test.txt"
cp $DIR"Gtrain_concMes_IFs_25_mean_var_selected.txt" $TEMPORAL"train.txt"

cd $LIB
export OMP_NUM_THREADS=11
/usr/bin/python easy.py $TEMPORAL"train.txt" $TEMPORAL"test.txt" >> $SALIDA
rm test.txt.scale
rm train.txt.model
rm train.txt.range
rm train.txt.scale
rm train.txt.scale.out
rm train.txt.scale.png


cd $BASE

echo ""
echo "LISTO!"
echo ""
