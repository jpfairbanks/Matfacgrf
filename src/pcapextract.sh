#pcapextract.sh
#Author: James Fairbanks
#Date: 2014-04-08
#Extracts a pcap file into a csv edgelist with
# time_relative	ip.src	ip.dst
# call with the first argument a pcap file
#infile=equinix-chicago.dirA.20130529-125710.UTC.anon
infile=$1
params="-e frame.time_relative -e ip.src -e ip.dst"
tshark -T fields -r $infile.pcap $params > $infile.tsv
