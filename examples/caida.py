#caida.py: this program converts pcap data from CAIDA into a sparse matrix COO format
# it can also write the transformed_pcap data to stdout
#Author: James Fairbanks
#Date: 2014-04-04 16:54

import scipy as sp
import scipy.io as spio
import pandas as pd
import numpy as np
from sys import stdin
from sys import stdout

verbose = False
#fp = "/home/users/jfairbanks/remotes/shared/DATA/CAIDA/sipscan/sipscan.csv"
fp = stdin
chunksize = None
matrix_market = False
transformed_pcap = True

#Load all the data at once.
df = pd.read_csv(fp, header=0, chunksize=chunksize)

#compute a unique integer label for each vertex.
hash = dict()
i = 0
for src,dst in zip(df.src_id, df.dst_ip):
	#print(src,dst)
	if src not in hash:
		hash[src] = i
		i += 1
	if dst not in hash:
		hash[dst] = i
		i += 1
if verbose:
	print(hash)

#convert ip addresses to vertex labels
I = np.zeros(len(df.src_id))
J = np.zeros(len(df.dst_ip))
data = np.zeros(len(df.src_id))
edges = zip(df.src_id, df.dst_ip)

if matrix_market:
	for k, e in enumerate(edges):
		#print(k, e)
		src = e[0]
		dst = e[1]
		I[k] = hash[src]
		J[k] = hash[dst]
		data[k] += 1
	#write the matrix market format
	A = sp.sparse.coo_matrix((data, (I, J)))
	if verbose:
		print(A)
	spio.mmwrite("sipscan", A)


#write out original table with new vertex labels.
if transformed_pcap:
	df.src_id = df.src_id.map(hash)
	df.dst_ip = df.dst_ip.map(hash)
	df.to_csv(stdout, index=False)
