#!/bin/bash


make -j distFunc=L2


./range -l ../data/selectivity100/5blobs/70k50Dpoints_mnist_l2.txt ../data/selectivity100/5blobs/mnist50_sample_workload_1k_l2.txt > out_l.txt
./range -s ../data/selectivity100/5blobs/70k50Dpoints_mnist_l2.txt ../data/selectivity100/5blobs/mnist50_sample_workload_1k_l2.txt > out_s.txt
./range -m ../data/selectivity100/5blobs/70k50Dpoints_mnist_l2.txt ../data/selectivity100/5blobs/mnist50_sample_workload_1k_l2.txt > out_m.txt
./range -m -v 128 ../data/selectivity100/5blobs/70k50Dpoints_mnist_l2.txt ../data/selectivity100/5blobs/mnist50_sample_workload_1k_l2.txt > out_m128.txt
./range -c -v 128 ../data/selectivity100/5blobs/70k50Dpoints_mnist_l2.txt ../data/selectivity100/5blobs/mnist50_sample_workload_1k_l2.txt > out_c128.txt

python3 ../AdaptiveIndexNewVersion_Clean/compareFiles.py out_l.txt out_s.txt
python3 ../AdaptiveIndexNewVersion_Clean/compareFiles.py out_l.txt out_m.txt
python3 ../AdaptiveIndexNewVersion_Clean/compareFiles.py out_l.txt out_m128.txt
python3 ../AdaptiveIndexNewVersion_Clean/compareFiles.py out_l.txt out_c128.txt


./knn -l -n 5 ../data/selectivity100/5blobs/70k50Dpoints_mnist_l2.txt ../data/selectivity100/5blobs/mnist50_sample_workload_1k_l2.txt > knn_l.txt 2> results_l.txt
./knn -k -n 5 ../data/selectivity100/5blobs/70k50Dpoints_mnist_l2.txt ../data/selectivity100/5blobs/mnist50_sample_workload_1k_l2.txt > knn_k.txt 2> results_k.txt
./knn -o -n 5 ../data/selectivity100/5blobs/70k50Dpoints_mnist_l2.txt ../data/selectivity100/5blobs/mnist50_sample_workload_1k_l2.txt > knn_o.txt 2> results_o.txt
./knn -k -v 128 -n 5 ../data/selectivity100/5blobs/70k50Dpoints_mnist_l2.txt ../data/selectivity100/5blobs/mnist50_sample_workload_1k_l2.txt > knn_k128.txt 2> results_k128.txt


python3 ../AdaptiveIndexNewVersion_Clean/checkResultsKNN.py results_l.txt results_k.txt
python3 ../AdaptiveIndexNewVersion_Clean/checkResultsKNN.py results_l.txt results_o.txt
python3 ../AdaptiveIndexNewVersion_Clean/checkResultsKNN.py results_l.txt results_k128.txt


make -j distFunc=ED


./range -l -n 5 ../data/words/wordsNewNoSpaces.txt ../data/words/queries/2letterwords2nn.txt > out_lstring.txt
./range -s -n 5 ../data/words/wordsNewNoSpaces.txt ../data/words/queries/2letterwords2nn.txt > out_sstring.txt
./range -m -n 5 ../data/words/wordsNewNoSpaces.txt ../data/words/queries/2letterwords2nn.txt > out_mstring.txt
./range -m -v 128 -n 5 ../data/words/wordsNewNoSpaces.txt ../data/words/queries/2letterwords2nn.txt > out_m128string.txt
./range -c -v 128 -n 5 ../data/words/wordsNewNoSpaces.txt ../data/words/queries/2letterwords2nn.txt > out_c128string.txt


python3 ../AdaptiveIndexNewVersion_Clean/compareFiles.py out_lstring.txt out_sstring.txt
python3 ../AdaptiveIndexNewVersion_Clean/compareFiles.py out_lstring.txt out_mstring.txt
python3 ../AdaptiveIndexNewVersion_Clean/compareFiles.py out_lstring.txt out_m128string.txt
python3 ../AdaptiveIndexNewVersion_Clean/compareFiles.py out_lstring.txt out_c128string.txt


./knn -l -n 5 ../data/words/wordsNewNoSpaces.txt ../data/words/queries/2letterwords2nn.txt > knn_lstring.txt 2> results_lstring.txt
./knn -k -n 5 ../data/words/wordsNewNoSpaces.txt ../data/words/queries/2letterwords2nn.txt > knn_kstring.txt 2> results_kstring.txt
./knn -k -v 128 -n 5 ../data/words/wordsNewNoSpaces.txt ../data/words/queries/2letterwords2nn.txt > knn_k128string.txt 2> results_k128string.txt

python3 ../AdaptiveIndexNewVersion_Clean/checkResultsKNNStrings.py results_lstring.txt resulresults_kstringts_k.txt
python3 ../AdaptiveIndexNewVersion_Clean/checkResultsKNNStrings.py results_lstring.txt results_k128string.txt