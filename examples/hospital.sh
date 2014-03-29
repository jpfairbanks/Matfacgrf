julia=~/julia/julia
$julia examples/visualize.jl -f data/hospital_edges.csv -k 6 --batchsize 3000 --max_vertices 75
printf "Running correlation study\n"
$julia -e 'reload("src/correlation.jl"); test();'
