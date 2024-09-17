import umap
import scipy.io as sio
import sys

def pyumap(fn, fout):
  random_seed = 1986
  mat_contents = sio.loadmat(fn+".mat")

  X = mat_contents["X"]
  Y = mat_contents["Y_model"]
  y = mat_contents["y_model"]
  idx_epoch = mat_contents["idx_epoch"]


  X_umap = umap.UMAP(n_neighbors = 15, min_dist = 0.1, metric='euclidean').fit_transform(X)
  sio.savemat(fout+".mat", {"X_umap":X_umap, "Y":Y, "y":y, "idx_epoch":idx_epoch})

if __name__ == '__main__':
    x = str(sys.argv[1])
    y = str(sys.argv[2])

    pyumap(x, y) 

######################################################
# import umap
# from time import time
# import scipy.io as sio

# random_seed = 1986

# fn = "step3_output"
# print(fn)

# mat_contents = sio.loadmat(fn+".mat")
# Xn = mat_contents["Xn"]
# XdBn = mat_contents["XdBn"]
# y = mat_contents["y_"]
# z = mat_contents["z_"]

# t3 = time()
# XdBn_umap = umap.UMAP(n_neighbors = 15, min_dist = 0.1, metric='euclidean').fit_transform(XdBn)
# t4 = time()
# print("%s: %.2g sec" % ("dBnorm UMAP", t4 - t3))

# sio.savemat(fn+'_embed_', {"XdBn_umap":XdBn_umap, "y":y, "z":z})