#pcapextract.sh
#Author: James Fairbanks
#Date: 2014-04-08
#Extracts a pcap.gz file into a csv edgelist with
# time_epoch	ip.src	ip.dst
# call with the first argument a pcap file
# example: pcapextract.sh equinix-chicago.dirA.20130529-125710.UTC.anon
# will read equinix-chicago.dirA.20130529-125710.UTC.anon.pcap.gz and
# produce equinix-chicago.dirA.20130529-125710.UTC.anon.tsv
infile=$1
pcapsuffix="pcap"
params="-e frame.time_epoch -e ip.src -e ip.dst"
tshark -T fields -r $infile.$pcapsuffix $params > $infile.tsv
