# IHKA python
IntraHippocampal Kainic Acid

Python modules for ihka data analysis, based off of MatLab code from [McKenzieLab's repository](https://github.com/McKenzieNeuro/McKenzieLab/tree/main/IHKA).

## Dependencies
- numpy
- pyedf (check to see if this one works okay)
- 
- 

## TODO
*This list is ordered from highest to lowest priority.*
- [x] Implement Options.toml config file as user input interface
- [x] Implement wavelet transform module, for turning raw edf into 
  - [ ] Bullet proof with tests
  - [ ] Test compare output with matlab scripts
- [x] Implement file I/O for reading from binary files
  - [x] Bullet proof with tests
  - [ ] Test compare output with matlab scripts
- [x] Figure out how to connect to Azure blob storage via the API
- [ ] Figure out how to put this pipeline onto the Azure cloud
  - [x] create VM and ssh into it
  - [ ] DEBUG: read and write to blob storage unit via python azure api fromt he VM
  - [ ] access and manipulate data in the blob storage from the azure api dexterously (i.e. run IHKApy code over the connection)
- [ ] Implement featurizing module
- [ ] Implement train model module
- [ ] Generate docs
- [ ] Make sure dependencies are well listed
- [ ] Draw dependency graph
- [ ] Detailed sketch of pipeline 
- [ ] Sort out dependencies

## Expo

Displayed below is a 100-datapoint chunk of the output of `sm_make_wavelet_bank.compute_wavelet_gabor(signal,fs,freqs,xi)` where `signal` is a sample of random gaussian noise with the standard deviation set to 1000 i.e. `signal=np.random.normal(0,1000,500)`, `fs=16.0`, `freqs=[0.5,1,2,4,8]`, and `xi=5`. 

![sample from sm_make_wavelet_bank.make_wavelet_bank](https://raw.githubusercontent.com/dcxSt/IHKApy/main/img/wavelet_transforms_demo.png)

### Refs
- Article on IHKA [https://www.sciencedirect.com/science/article/abs/pii/S001448862030323X](https://www.sciencedirect.com/science/article/abs/pii/S001448862030323X) 
- Writing good pythonic python: [PEP8](https://pep8.org/#break-before-or-after-binary-operator)

