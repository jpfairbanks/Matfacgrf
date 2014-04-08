edges="../data/sipscan"
batchfile=$edges.batch.mtx
initfile=$edges.initial.mtx
p="0.2"
echo Splitting $edges.mtx into $initfile + $batchfile randomly with probability $p of being in the batch
python randedges.py -p $p $edges.mtx > $initfile 2> $batchfile
