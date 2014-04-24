BEGIN{OFS=",";print "time","src_id","dst_ip";}
$0 !~ /,/ {print $1, $2, $3}
/,/{errcount+=1; print errcount, $0> /dev/stderr }
