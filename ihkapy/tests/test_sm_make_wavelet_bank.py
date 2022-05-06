from ihkapy.sm_make_wavelet_bank import *

#[test]
def test_compute_wavelet_gabor(plot=False):
    signal = np.random.normal(0,1000,np.power(2,16))
    fs = 16.0 # Hz
    freqs = np.power(10,np.linspace(np.log10(0.5),np.log10(8),5))
    # logger.debug(f"Test compute_wavelet_gabor\nfreqs = {freqs}")
    xi = 5
    wt = compute_wavelet_gabor(signal,fs,freqs,xi)
    print(f"\n\twt.shape={wt.shape}\n\tInput signal shape={signal.shape}\n\tLen freqs={len(freqs)}")
    assert len(signal),len(freqs) == wt.shape
    if plot==True:
        logging.getLogger('matplotlib').setLevel(logging.WARNING)
        import matplotlib.pyplot as plt
        # plt.figure(figsize=(12,8))
        plt.subplots(2,1,figsize=(12,7))
        plt.suptitle(f"Wavelet Transforms of Gaussian Random Noise (sigma=1000)\nSample Frequency = {fs}")

        plt.subplot(2,1,1)
        plt.title("REAL")
        labels = [f"freq = {i:.2f}" for i in freqs]
        plt.plot(np.real(wt[100:200,:]),".",alpha=0.5,label=labels)
        plt.plot(np.real(wt[100:200,:]),"--",color="k",linewidth=0.5,alpha=0.1)
        plt.legend()

        plt.subplot(2,1,2)
        plt.title("IMAGINARY")
        labels = [f"freq = {i:.2f}" for i in freqs]
        plt.plot(np.imag(wt[100:200,:]),".",alpha=0.5,label=labels)
        plt.plot(np.imag(wt[100:200,:]),"--",color="k",linewidth=0.5,alpha=0.1)
        plt.legend()

        plt.show(block=True)
test_compute_wavelet_gabor(plot=True) # set plot=True to display plots in test 
print("TEST PASSED: compute_wavelet_gabor()") 
# strange logging bug, see https://stackoverflow.com/questions/72127312/python-logger-setlevel-bug 



#[test]
"""Warning this test is not pure, it relies on data being in certain files."""
def test_make_wavelet_bank():
    # TODO: implement this test
    return
test_make_wavelet_bank()
# uncomment below once implemented
# logger.info("TEST PASSED: make_wavelet_bank()") 



