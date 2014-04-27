"""matrixmarket.py: Converts an edgelist file to matrix market format"""
from __future__ import print_function
import sys
import scipy as sp
import scipy.io as spio
import pandas as pd

if __name__ == '__main__':
    df = pd.read_csv(sys.stdin, header=None, sep=' ')
    I = df[0]
    J = df[1]
    data = df[2]
    m = max(I.max(), J.max()) +1
    print(m, file=sys.stderr)
    A = sp.sparse.coo_matrix((data, (I, J)), shape=(m, m))
    B = (A + A.T)/2
    #write the matrix market format
    spio.mmwrite(sys.stdout, B)
