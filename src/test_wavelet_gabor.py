from sm_make_wavelet_bank import compute_wavelet_gabor
import numpy as np
from scipy.io import savemat,loadmat

def gen_and_save_raw_signal(fname="raw_signal.mat"):
    # create random signal with std = 1
    signal = np.random.normal(0,1,1000)

    # serialize to .m file so that matlab script can read it
    savemat(fname,{"signal":signal})
    print(f"Generated and saved raw signal to {fname}")
    return 

def load_ml_cellmat(fname="raw_signal.mat",field="signal") -> np.ndarray:
    # Also use this method with fname="wt_ml.m", field="wt"
    matfile_dict = loadmat(fname)
    return np.asarray(matfile_dict[field]) # assume what we want is stored in "signal" field

def read_wt_matlab():
    # Returns output of matlab's awt_freq transform

    return # dummy


def compare_signals(wt_py:np.ndarray, wt_ml:np.ndarray): # -> bool:
    import matplotlib.pyplot as plt
    # Display them alongside one another for Visual comparison
    plt.subplots(2,1,figsize=(12,7))
    plt.suptitle()
    plt.subplot(2,1,1)
    plt.title("Python")
    plt.plot(wt_py,)

    plt.subplot(2,1,2)
    plt.title("Matlab")

    plt.show(block=True)

    # Numerically test if the signals correspond pointwise (to 5 sig figs)
    py_rounded = wt_py.round(5)
    ml_rounded = wt_ml.round(5)
    print("The following should be all true.")
    for pyrow,mlrow in zip(py_rounded,ml_rounded):
        print(f"All same: {(pyrow==mlrow).all()}")


if __name__=="__main__":
    gen_and_save_raw_signal("raw_signal.mat")
    print("generated and saved raw signal")

    # # Get raw signal
    signal = np.squeeze(load_ml_cellmat(fname="raw_signal.mat",field="signal"))
    # fs = 16
    # freqs = [0.5,1.0,2.0,4.0]
    # # Compute wavelt transform with python module
    # wt_py = compute_gabor_wavelet(signal,fs,freqs)
    # 
    # # Get computed wavelet transform from matlab awt_freqs
    # wt_ml = load_ml_cellmat(fname="wt_ml.m",fieldname="wt")








