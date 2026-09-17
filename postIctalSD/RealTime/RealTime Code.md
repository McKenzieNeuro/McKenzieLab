# Closed-Loop Neural Recording & Stimulation System

A MATLAB-based real-time electrophysiology platform for closed-loop brain stimulation experiments. The system records neural signals via Intan hardware, computes theta/delta band power ratios, runs a machine learning seizure-prediction model, and triggers stimulation based on configurable criteria — all in parallel across multiple CPU cores.

---

## Table of Contents

- [Overview](#overview)
- [Requirements](#requirements)
- [System Architecture](#system-architecture)
- [Configuration](#configuration)
- [Experiment Programs](#experiment-programs)
- [Stimulation Trigger Modes](#stimulation-trigger-modes)
- [Data Outputs](#data-outputs)
- [Helper Functions](#helper-functions)
- [Getting Started](#getting-started)

---

## Overview

This script manages multi-subject, multi-day neural recording and stimulation sessions. It coordinates several parallel processes:

- Acquiring raw neural data over TCP from an Intan RHS or RHD recording system
- Computing real-time theta (5–12 Hz) and delta (1–4 Hz) band amplitudes via Hilbert transform
- Running a trained random forest classifier to predict seizure onset
- Sending stimulation triggers to Intan based on configurable closed-loop criteria
- Recording synchronized video from multiple webcams
- Saving all data streams to disk concurrently

---

## Requirements

- MATLAB with the following toolboxes:
  - Parallel Computing Toolbox
  - Image Processing Toolbox
  - Computer Vision Toolbox (for pose/frame labeling)
  - Data Acquisition Toolbox (for Arduino)
  - Signal Processing Toolbox
- Intan RHS and/or RHD recording system with TCP enabled
- Arduino Uno (for LED/timestamp sync pulses)
- Logitech HD Pro Webcam C920 (or compatible)
- A trained `.mat` model file (RUSBoost random forest) for seizure prediction

---

## System Architecture

The system uses MATLAB's `spmd` (Single Program, Multiple Data) parallel computing to distribute work across CPU cores. Each core has a dedicated role:

| Core | Role | Description |
|---|---|---|
| `intanIO` | **Intan I/O** | Reads TCP data, controls stimulation, orchestrates all other cores |
| `thetaCalc` | **Theta Calculator** | Bandpass filters (5–12 Hz) and computes mean Hilbert amplitude |
| `deltaCalc` | **Delta Calculator** | Bandpass filters (1–4 Hz) and computes mean Hilbert amplitude |
| `predictModel(1..N)` | **Seizure Predictor** | One core per active subject; extracts features and runs the classifier |
| `dataSaver` | **Data Writer** | Writes all binary `.dat` files to disk |
| `camReader` | **Camera Reader** | Captures frames from all webcams |
| `camSaver` | **Video Writer** | Synchronizes and saves frames as `.mp4` files |
| `frameLabeler` | **Pose Estimator** | Optionally applies a CNN to label animal pose in video |
| `setup` | **Setup** | Configuration and outer-loop coordination |
| `DBGplotData` | **Debug Plotter** | Optional real-time data visualization |

---

## Configuration

All experiment parameters are set at the top of the script. Key sections:

### Subjects & Boxes

```matlab
stimParam(1).config(1).subject = '';       % Box 1 subject ID (empty = inactive)
...
stimParam(1).config(4).subject = 'PTP6.0'; % Box 4 subject ID
```

Boxes with empty subject fields are automatically excluded from recording and stimulation.  Third dimension corresponds to *n* subject.

### Recording Start Time

```matlab
StartTime = [8, 00]; % 24-hour format [hr, min]
nextDay = 0;         % Set to 1 to delay start until tomorrow
```

### Intan Hardware

```matlab
useRHS = 0;          % Enable Intan RHS system
useRHD = 1;          % Enable Intan RHD system
recIntan = 1;        % 1 = record to disk, 0 = run only
intanIP_RHD = '127.0.0.1';
intanPort1_RHD = 5000; % Command port
intanPort2_RHD = 5001; % Data output port
```

### Analysis & Downsampling

```matlab
sF   = 20000; % Raw sample rate (Hz)
dSF  = 1250;  % Downsampled rate (Hz)
```

### Prediction Model

```matlab
usePrediction = 1;
modelPath{1,1} = 'R:\path\to\model.mat'; % One path per active box
```

The model `.mat` file must contain a `rusTree` field (RUSBoost classifier) and an `ops` struct with fields `freqs`, `durFeat`, `features`, `ch_subj`.

### Debug Mode

To enable verbose logging for specific cores:

```matlab
debug_cores = [predictModel, dataSaver, frameLabeler, DBGplotData];
```

---

## Experiment Programs

Set `program` to select the session structure:

| Value | Program | Description |
|---|---|---|
| `1` | **Baseline** | Baseline recording with no stimulation |
| `2` | **Open Loop** | Baseline → open-loop stimulation → baseline |
| `3` | **Kindling** | Baseline → induction stimulation → baseline |
| `4` | **Closed Loop** | Intelligent prediction-triggered stimulation |
| `5` | **Induction Peak** | Open-loop S1 followed by peak-theta-triggered S2 |
| `6` | **Induction Decay** | Open-loop S1 followed by decay-theta-triggered S2 |

Each program defines a `schedule` (block types per day) and `duration` (seconds per block). The schedule iterates across recording days automatically.

---

## Stimulation Trigger Modes

Each block label encodes which stimulator and trigger type to use:

| Label | Meaning |
|---|---|
| `BL` | Baseline — no stimulation |
| `S1_O` | Stim 1, open loop (fixed ISI) |
| `S2_ptC` | Stim 2, triggered at peak theta/delta ratio |
| `S2_dtC` | Stim 2, triggered at decay of theta/delta ratio |
| `S1_tdC` | Stim 1, triggered when theta/delta falls below threshold |
| `S1_ipC` | Stim 1, intelligent prediction closed loop |
| `S1X_O` | Sham version of S1 open loop (no current delivered) |
| `SX` | Sham for this box only |

### Closed-Loop Trigger Parameters

```matlab
threshold_tdC = 1;  % T/D ratio threshold for tdC mode
threshold_ptC = 2;  % T/D ratio threshold for ptC mode
sustain_ptC   = 4;  % Seconds T/D must stay above threshold (ptC)
threshold_dtC = 1;  % T/D threshold for decay detection
sustain_dtC   = 4;  % Seconds T/D must stay below threshold after peak (dtC)
stimClass     = 2;  % Model output class that triggers stimulation (ipC)
numConsec     = 5*dSF; % Consecutive samples of stimClass needed (ipC)
```

---

## Data Outputs

All files are saved to a timestamped session directory under the selected datastore path.  All variables which are dynamically named (e.g. to save sampling rate in the file name) are denoted as {variable}.

| File | Format | Description |
|---|---|---|
| `amplifier_{numChan}Ch_{dsF}KHz_int16.dat` | int16 binary | Downsampled amplifier data |
| `timestamps_2Ch_{dsF}KHz_uint8.dat` | uint8 binary | Sample timestamps |
| `ToD_2Ch_{dsF}_uint8.dat` | uint8 binary | Theta/delta ratio (×25 scaled) |
| `Labels_2Ch_{dsF}_uint8.dat` | uint8 binary | Per-sample classifier labels |
| `Stim_{numChan}Ch_{dsF}KHz_uint8.dat` | uint8 binary | Stimulation event markers |
| `digitalin_1Ch_{dsF}KHz_.dat` | uint16 binary | Downsampled digital input |
| `stim_{numChan}Ch_{dsF}KHz_.dat` | int16 binary | Processed stimulation timestamps |
| `*.mp4` | MPEG-4 video | Synchronized webcam recordings (10 fps default) |
| `recInfo.mat` | MATLAB struct | Session metadata, block timing, stim parameters |
| `config.mat` | MATLAB struct | Full configuration snapshot for reproducibility |

---

## Helper Functions

| Function | Description |
|---|---|
| `setStim(TCP, enable, config)` | Configures Intan stim parameters over TCP |
| `byte2doubleConv(...)` | Decodes raw TCP byte stream into amplifier data and timestamps |
| `convLPFiltFix(...)` | FFT-based windowed sinc low-pass filter with edge padding |
| `plotFrame(input)` | Real-time video display callback |
| `plotData(input)` | Real-time theta/delta ratio plot callback |
| `plotLabel(input)` | Real-time classifier label plot callback |
| `makeBox(input)` | Creates stop dialog box to interrupt recording |
| `getTime()` | Returns current time in seconds since midnight |
| `returnSubfolder(dirIn)` | Lists subdirectories of a given path |
| `RemapPoint(...)` | Remaps image coordinates between resolutions |
| `symsepchar(StrIn, Sym)` | Splits a string on a delimiter character |
| `sm_wavelet(x, Fs, freqlist)` | Morlet continuous wavelet transform |

---

## Getting Started

1. Set subject IDs and active boxes in `stimParam(1).config`.
2. Set `program` to the desired experiment type.
3. Set `StartTime` to the desired recording start.
4. Verify Intan TCP ports and IP addresses match the recording software.
5. Set `modelPath` entries to point to your trained `.mat` classifier file.
6. Run the script. You will be prompted to:
   - Label webcams by cage number
   - Select active channels for each port
   - Choose a datastore directory
   - Optionally copy data to a backup directory
7. The system will wait until `StartTime`, then begin recording automatically.

To stop a session early, close the **Stop Closed Loop** dialog that appears on screen.