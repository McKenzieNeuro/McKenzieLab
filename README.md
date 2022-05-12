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
  - [ ] Bullet proof with tests
  - [x] Test compare output with matlab scripts
    - [ ] also, implement this test in the test suit (not a priority, but should be done eventually)
- [x] Implement file I/O for reading from binary files
  - [x] Bullet proof with tests
  - [x] Test compare output with matlab scripts
- [x] Figure out how to connect to Azure blob storage via the API
- [ ] Figure out how to put this pipeline onto the Azure cloud
  - [x] create VM and ssh into it
  - [ ] DEBUG: read and write to blob storage unit via python azure api fromt he VM
  - [ ] Mout storage from VM (see [azure doc](https://docs.google.com/document/d/1lXst8D3eh3-yyND3NJNg4Wx16jm18QAJQrwL0WdpC-Q/edit?usp=sharing) for links to microsoft documentation on how to do this)
  - [ ] access and manipulate data in the blob storage from the azure api dexterously (i.e. run IHKApy code over the connection)
- [ ] Implement featurizing module
- [ ] Implement train model module
- [ ] Generate docs
- [ ] Make sure dependencies are well listed
- [ ] Draw dependency graph
- [ ] Detailed sketch of pipeline 
- [ ] Sort out dependencies

*Questions / Concerns*
- [ ] The raw data that gets serialized is not rescaled before it is cast into int16, in the case of the mouse data this causes non-negligible quantization because the RMS is around 12 (11.861 for the first channel), most of samples are in the single and low double digits. They could probably benefit from being scaled by a similar scale factor as the power and phase convolutions. 

## Expo

Displayed below is a 100-datapoint chunk of the output of `sm_make_wavelet_bank.compute_wavelet_gabor(signal,fs,freqs,xi)` where `signal` is a sample of random gaussian noise with the standard deviation set to 1000 i.e. `signal=np.random.normal(0,1000,500)`, `fs=16.0`, `freqs=[0.5,1,2,4,8]`, and `xi=5`. 

![sample from sm_make_wavelet_bank.make_wavelet_bank](https://raw.githubusercontent.com/dcxSt/IHKApy/main/img/wavelet_transforms_demo.png)

### Refs
- Article on IHKA [https://www.sciencedirect.com/science/article/abs/pii/S001448862030323X](https://www.sciencedirect.com/science/article/abs/pii/S001448862030323X) 
- Writing good pythonic python: [PEP8](https://pep8.org/#break-before-or-after-binary-operator)
- [Guide to python packaging (the docs)](https://python-packaging.readthedocs.io/en/latest/dependencies.html)




