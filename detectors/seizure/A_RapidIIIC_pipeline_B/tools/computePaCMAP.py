import scipy.io as sio
import sys
import pacmap
import numpy as np

def pyumap(fn_in, fn_out):
  random_seed = 1986
  mat_contents = sio.loadmat(fn_in)
  X = mat_contents["X"]

  embedding = pacmap.PaCMAP(n_components=2, n_neighbors=None, MN_ratio=0.5, FP_ratio=2.0) 
  # fit the data (The index of transformed data corresponds to the index of the original data)
  X_pacmap = embedding.fit_transform(X, init="pca")
  sio.savemat(fn_out+".mat", {"X_pacmap":X_pacmap})

if __name__ == '__main__':
    x = str(sys.argv[1])
    y = str(sys.argv[2])
    pyumap(x, y) 
