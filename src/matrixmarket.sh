PREFIX="$HOME/.julia/v0.3/Matfacgrf/src"
mmpy="$PREFIX/matrixmarket.py"
tail -n +2 |tr "," " "|  cut -d" " -f2,3 |\
    sort -n | uniq -c | tr --squeeze " " |\
    awk '{print $2,$3,$1}' | python $mmpy
