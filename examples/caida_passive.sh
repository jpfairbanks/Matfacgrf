#caida_passive.sh: processes the files produced by ../src/pcapextract.sh into timestamped edge list
# the output files have a header row that looks like
#        time,src_id,dst_ip
#
# and data rows which look like
#        2.3232e-5,3,102
# the time is relative to the first packet arrival time.
#example file name equinix-chicago.dirA.20130529-125710.UTC.anon
filename=$1
PROJHOME="$HOME/Projects/Matfacgrf/examples/"
renamevertices="$PROJHOME/caida.py"
awkscript="$PROJHOME/caida_passive.awk"
cat  $filename.tsv  | awk -v FS="\t" -f $awkscript |\
    python $renamevertices > $filename.edges
    
#> $filename.csv
