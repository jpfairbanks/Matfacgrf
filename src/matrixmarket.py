"""matrixmarket.py: Converts a csv timestamped edgelist file to matrix market format"""

import sys
import scipy as sp
import scipy.io as spio
import pandas as pd

if __name__ == '__main__':
    df = pd.read_csv(sys.stdin, header=None,)
    I = df[1]
    J = df[2]
    data = [1 for i in range(len(I))]
    #write the matrix market format
    A = sp.sparse.coo_matrix((data, (I, J)))
    spio.mmwrite(sys.stdout, A)
