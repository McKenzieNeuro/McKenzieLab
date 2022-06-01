# IHKA python
IntraHippocampal Kainic Acid

Python modules for ihka data analysis, based off of MatLab code from [McKenzieLab's repository](https://github.com/McKenzieNeuro/McKenzieLab/tree/main/IHKA).

## Install & Setup 
To access the methods and their modules, clone this repo with `git clone https://github.com/dcxSt/IHKApy`, enter it's root directory `cd IHKApy`, then pip install with `pip install .`

## Dependencies
- numpy
- scipy
  - scipy.io (only used in tests to read matlab files)
  - scipy.stats (for zscore-ing)
- pyedf (TODO:check to see if this one works okay)
- tqmd (for progress bar)
- re (regex library for tests and asserts, make sure verything formated well)

## TODO
*This list is ordered from highest to lowest priority.*
- [x] Implement Options.toml config file as user input interface
- [x] Implement wavelet transform module, for turning raw edf into 
  - [x] Bullet proof with tests
  - [x] Test compare output with matlab scripts
    - [ ] also, implement this test in the test suit (not a priority, but should be done eventually)
  - [x] Implement and test `make_all` method
- [x] Implement file I/O for reading from binary files
  - [x] Implement reading multiple segments at once
  - [x] Bullet proof with tests
  - [x] Test compare output with matlab scripts
- [x] Figure out how to connect to Azure blob storage via the API
- [ ] Figure out how to put this pipeline onto the Azure cloud
  - [x] create VM and ssh into it
  - [ ] DEBUG: read and write to blob storage unit via python azure api fromt he VM
  - [ ] Mout storage from VM (see [azure doc](https://docs.google.com/document/d/1lXst8D3eh3-yyND3NJNg4Wx16jm18QAJQrwL0WdpC-Q/edit?usp=sharing) for links to microsoft documentation on how to do this)
  - [ ] access and manipulate data in the blob storage from the azure api dexterously (i.e. run IHKApy code over the connection)
- [ ] Implement featurizing module
  - [ ] Tim selection algorithm
- [ ] Implement train model module
- [ ] Generate docs
- [ ] Make sure dependencies are well listed
- [ ] Draw dependency graph
- [ ] Detailed sketch of pipeline 
- [ ] Sort out dependencies

- [ ] Little todos, formatting tasks
  - [ ] rename data ops FREQ_NUM etc. instead of what, add "FREQ" to SPACING param too
- [ ] Change `fileio` the dictionary to `fio_ops` because it's ambiguous and can eaily get confused with the module by the same name. Change this everywhere in the codebase including the `Options.toml`. 
- [ ] Make sure precision isn't hard coded anywhere, get it from ops "int16"
- [ ] Do a consistency check of variable names

*Questions / Concerns*
- [ ] The raw data that gets serialized is not rescaled before it is cast into int16, in the case of the mouse data this causes non-negligible quantization because the RMS is around 12 (11.861 for the first channel), most of samples are in the single and low double digits. They could probably benefit from being scaled by a similar scale factor as the power and phase convolutions. 

*Discussion*
Although Pynapple's abstractions are very nice, they're not exactly what we want. If we wanted to make use of the interval restrict functionality, we would have to load large chunks of data (13GB) into RAM just to sample a small fraction of it. Instead, we implement methods to read small specific windows of binary files at a time. We could load those windows into a pynapple Tsd object, but I think this will not be necessary; it would be like using a sledgehammer to insert a pin into a corkboard. The only place something like the Interval Sets object might be handy is when defining those windows which to sample, but here I think the same can be accomplished with either an array or a pandas data frame. 

There is an intuitive but rigorous naming convention for binary files: each file  must match `${basename}_ch_\d\d\d.dat` where the three digits are the zero-padded index of the channel, and basename must match the `.edf` data file and the `.txt` metadata (seizure timestamps) file. 


## Technical notes and comments
Currently we sample wavelet frequencies from a log space spanning 0.5 Hz to 200 Hz; our data's sampling frequency is 2000 Hz. 

Under the hood, scipy's coherence function (and mscohere) computes the mean of the square of the cross spectral density divided by the psd of each signal, and averages over as many 256-sample windows that it can fit into our 5-second long segment (=10000 samples). The smallest frequency it samples is only half a wavelength, the second smallest is a full wavelength. These low frequencies (the lowest of which is about 7Hz) are already much greater than our low wavelet frequencies.  


## Conventions
- The `edf` raw data file and it's `txt` metadata file must have the same basename and both be located in the `RAW_DATA_PATH` directory. 

## Expo

Displayed below is a 100-datapoint chunk of the output of `sm_make_wavelet_bank.compute_wavelet_gabor(signal,fs,freqs,xi)` where `signal` is a sample of random gaussian noise with the standard deviation set to 1000 i.e. `signal=np.random.normal(0,1000,500)`, `fs=16.0`, `freqs=[0.5,1,2,4,8]`, and `xi=5`. 

![sample from sm_make_wavelet_bank.make_wavelet_bank](https://raw.githubusercontent.com/dcxSt/IHKApy/main/img/wavelet_transforms_demo.png)

### Refs
- Article on IHKA [https://www.sciencedirect.com/science/article/abs/pii/S001448862030323X](https://www.sciencedirect.com/science/article/abs/pii/S001448862030323X) 
- Writing good pythonic python: [PEP8](https://pep8.org/#break-before-or-after-binary-operator)
- [Guide to python packaging (the docs)](https://python-packaging.readthedocs.io/en/latest/dependencies.html)




