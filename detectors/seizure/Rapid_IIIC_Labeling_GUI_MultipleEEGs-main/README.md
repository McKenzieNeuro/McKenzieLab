# Rapid_IIIC_Labeling_GUI_MultipleEEGs
GUI for rapid labeling of segments from multiple EEGs, and instructions for preparing data for these labeling tasks

This repository describes a GUI developed by Dr. Jin Jing ("JJ"), PhD that enables rapid labeling of IIIC patterns. By following these instructions you can set up the GUI to annotate your own data. Note that this GUI allows labeling of isolated samples from multiple EEGs in one task. This is different from the GUI that allows labeling of a single EEG exhaustively. 

Files and data (~8GB) can be downloaded [here](https://www.dropbox.com/scl/fo/9u457f1usu6zkjr45s6vy/h?dl=0&rlkey=8z5vhcv2x0u2yi9qzsot3xqom).

**Requirements:** MATLAB, EEGLAB (https://sccn.ucsd.edu/eeglab/index.php), and Python (Anaconda)

**Input data:** Raw EDF files inside .\Data\EDF\; scalp monopolar/C2 EEG that contains a full set of 19 channels + 1 EKG (optional), with channels named as follows:

Fp1 F3 C3 P3 F7 T3 T5 O1 Fz Cz Pz Fp2 F4 C4 P4 F8 T4 T6 O2 (EKG)  

<img width="193" alt="image" src="https://user-images.githubusercontent.com/10371730/217501280-143a4953-4c6b-4f2e-b4b5-18ecb03e2db1.png">

**Step1:** Read EDF to MAT using EEGLAB toolbox for MATLAB. Run script step1_readEDF2MAT.m, which converts EDF format to MAT format in .\Data\MAT\ that contains the following variables:

- _data_: EEG array
- _channels_: list of channel names in data
- _Fs_: the sampling rate of data
- _startTime_: the start time vector of data
<img width="247" alt="image" src="https://user-images.githubusercontent.com/10371730/217501324-572925b8-504f-4ad3-b580-9921e4768144.png">

**Step2:** Preprocess MAT to select/rearrange channels, resample to 200Hz , and denoise with 0.5 40Hz band-pass and 5Hz band-stop centered at the power-line frequency (US: 60Hz, UK: 50Hz). Output files are saved in .\Data\processed\
<img width="247" alt="image" src="https://user-images.githubusercontent.com/10371730/217501379-2c7a14e8-ef35-460b-977e-7630f0b421de.png">

**Step3:** run SPaRCNet (Python backend)

==Configure Python==

Install anaconda3 and open a terminal

`$ conda create -n iiic python=3.6`  
`$ activate iiic`  
`$ conda install -c conda-forge hdf5storage`  
`$ pip install mne`  
`$ pip install torch==1.5.0+cpu torchvision==0.6.0+cpu -f https://download.pytorch.org/whl/torch_stable.html`  

Run MATLAB wrapper step3_runSPaRCNat.m

- CSV score table will be export to .\Data\iiic\
- Each row is the probabilities for 6 classes: Other, Seizure, LPD, GPD, LRDA, and GRDA
- Starting from the 1^st^ 10sec EEG segment and moving at 2sec step in time. Eg. row #1: scores for 0 10sec, row #2: scores for 2 12sec, ...

<img width="340" alt="image" src="https://user-images.githubusercontent.com/10371730/217501490-c9d58c3c-3d37-422f-ad01-ed00bcfd7225.png">

**Step4:** Run step4_readCSV.m to read CSV to MAT to make sure every 2sec segment has scores. The output files are saved in .\Data\iiic\model_prediction\.

**Step5:** Run step5_computeSpectrograms to get regional average spectrograms in .\Data\Spectrograms\, which contains the following variables:

- _Sdata_: 4 regional average spectrograms
- _stimes_: time coordinates
- _sfreqs_: frequency coordinates
- _params_: spectrogram parameters

**Step6:** Run step6_segementEEG.m to divide EEG into stationary periods via change point (CP) detection. The output look-up-table is exported to .\Data\CPDs\ which contains the following variables:

-	isCPcenters: 0 (not CP center) or 1 (is CP center)
-	isCPs: 0 (not CP) pr 1 (is CP)
-	lut_cpd: index of each CP center (column #1) and its range [start (column #2),  end (column #3)]; lut = lookup table
<img width="272" alt="image" src="https://user-images.githubusercontent.com/10371730/217501727-fc8cdb60-188e-4e10-ba63-12abc9a93b81.png">

**Step7:** 
Run step7_parseCPcenters.m to parse data of CP centers for labelling GUI in .\Data\CP_centers\ each contains the following variables:
-	SEG: 14sec EEG 
-	Sdata: 10min spectrograms
-	fileKey: EEG file token
-	idx_CPcenter: the index of CP center (2sec unit) in cEEG
-	idx_CPrange: the range of CP segment (2sec unit) in cEEG
-	scores: model predicted probabilities, in order of Other, Seizure, LPD, GPD, LRDA, GRDA
-	sfregs: frequency coordinates of spectrograms
<img width="470" alt="image" src="https://user-images.githubusercontent.com/10371730/217500593-34a296c3-b575-4cab-ac9f-ab81d63f8193.png">

**Step8:**
Run `step8_getLUT.m` to get the global look-up-table (LUT) for labelling GUI:
-	Column #1: EEG name
-	Column #2: sample index in 2sec segment (CP center)
-	Column #3: CP range 
-	Column #4: model prediction
-	Column #5: model probabilities for Other, Seizure, LPD, GPD, LRDA, GRDA

<img width="470" alt="image" src="https://user-images.githubusercontent.com/10371730/217500532-57e6e6eb-5f3b-4023-8e96-de9743545ab1.png">

**Step9:**
Run `step9_getBoW.m` to get the bag of word (BoW) model using spectrograms with 500 words using K-means clustering method and compute the normalized distribution of words for each sample as BoW feature. This is further used in labelling GUI to do similarity search in chi -square distance.

**Step10:**
-	Configure PaCMAP (Python library:  https://github.com/YingfanWang/PaCMAP). 
-	Run step10_getPaCMap.m to get the global embedding.

**Step11:**
Run labeling GUI step11_IIICGUI_mPatients.m   

input rater initials to store scores  
<img width="135" alt="image" src="https://user-images.githubusercontent.com/10371730/217502524-d3c6aea2-ae40-4348-b115-eef124a5c4c2.png">

click Start to continue   
<img width="470" alt="image" src="https://user-images.githubusercontent.com/10371730/217502639-33f1a7aa-c8c4-4192-8b1c-222060313aa7.png">

Label 30 to 50 samples selected by GUI per interaction (samples at the class boundaries).  
<img width="467" alt="image" src="https://user-images.githubusercontent.com/10371730/217505443-c2b3304f-2e3c-45f7-bc6e-81002e2f38e7.png">

Press Yes to update labels (spreading in PaCMAP by nearest labeled points) and enter next iteration.  
Press Done button to seal and export the labels.  
<img width="470" alt="image" src="https://user-images.githubusercontent.com/10371730/217521077-be6e0511-3bd3-4f18-8123-662e8aa5f0df.png">
