
config = load_toml()    # load the toml config file
RAW_DATA_PATH = config.fileio.RAW_DATA_PATH # somthing like this
WAVELET_BINARIES_PATH = config.fileio.WAVELET_BINARIES_PATH


def get_wavelet_bank(
        path_to_edf_file,
        n_chan,
        now,
        list_,
        all,
        the,
        options,
        its,
        unfortunately,
        the,
        least,
        of_evils): # get power per channel
    return

def wavelet_convolve(signal, fs, freq):
    return 


def main():
    # Find all raw.edf and raw.txt pairs inside RAW_DATA_PATH
    # raw_data_files = [(patient_1_edf_path,patient_1_txt_path),(p2_edf,p2_txt),...]

    # For each file, get the power on each channel
    kwargs = {"n_chan":xyz,"abc":config.params.data.etc} # varargin
    for edf,path in raw_data_files:
        get_wavelet_bank(edf,**kwargs)






