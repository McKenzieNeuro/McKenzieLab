# -*- coding: utf-8 -*-
"""
The University of New Mexico Hospital EEG Corpus Builder.

@author: Nick Aase (aase@cs.unm.edu)

Requirements
------------
pandas
numpy
PyMuPDF
pyedflib
mne
mne_bids

Notes
-----
22/06/30 (NEJA): Initial build.
22/07/11 (NEJA): Augmented REs used. Began to play with NLP software (nltk)
22/07/13 (NEJA): Made a class for documents that integrates anonymization.
                 Need to do some refactoring, though.
22/07/15 (NEJA): Implemented BIDS saving and labeling. Patients are anonymized
                 (for the most part), and dates are offset (but kept uniform
                 from patient to patient so as to improve longitudinal use).
                 More refactoring performed, and more to come.
22/07/22 (NEJA): Moved Clinical_Notes class to unmhcb_io.py as well as static
                 I/O methods from the Corpus class.
22/07/26 (NEJA): Fixed conversion from NKT to EDF.
22/07/28 (NEJA): Adjusted birthdates and recording dates for EDF output files.
                 Now BDs are rolled back to Jan. 1st of the year patients were
                 born, and recording dates are relative to that.
22/07/29 (NEJA): Changed glob to walk in sysScan, which eliminates the need for
                 manual recursion.
22/08/03 (NEJA): Substantial tweaks now that more than one patient's data
                 are being processed.
"""
from path import Path
import pandas as pd
import numpy as np
import filecmp
import os
import sys
import time
from dateutil.relativedelta import relativedelta
import mne
import mne_bids
import warnings
import logging
from grabchart import grabchart, RELATIVE_PATH
from unmhcb_io import get_PNT_metadata, get_MRN_from_PNT, copyNKT, \
    PNT_FIELDS, Clinical_Notes, convert_NKT2EDF


# Constants
BIDS_DB_PATH = Path("../test")
SHARED_EEG_DRIVE = Path("//EEG-MAS-GW/EEG/")
STAGING_DIR = Path('nk_staging_data')
CSV_FILENAME = 'PNT_metadata.csv'
SUB_PAD = '05'
SES_PAD = '02'


class CorpusBuilder:
    """The main datastructure for building the EEG Corpus."""

    def sysScan(self, fullPath: str, depth=None):
        """
        Scan a file system for EEG data corresponding to a list of MRNs.

        Parameters
        ----------
        fullPath : str
            Top-most directory to begin searching.
        depth : int | None, optional
            If *not* None, specify a recursion depth. The default is None.

        Returns
        -------
        None.

        """
        curDir = Path(fullPath)
        # NEJA Maybe put something to avoid re-scanning directories...
        for file in curDir.walk('*.PNT'):  # Just consider PNT files
            if os.stat(file).st_size == 0:
                continue  # Skip blank files
            mrn = get_MRN_from_PNT(file)
            if mrn in self.mrns:
                logging.info(f"Found {mrn} at {file}")
                repeat = False
                for saved_file in self.df.Path:
                    if filecmp.cmp(file, saved_file, shallow=False):
                        # Avoid writing the same data multiple times
                        repeat = True
                        warnings.warn(f"Duplicate files found: "
                                      f"{file} & {saved_file})!!")
                        break
                if not repeat:
                    metadata = [file] + get_PNT_metadata(file)
                    session = {i: j for i, j in zip(self.df.columns, metadata)}
                    self.df = self.df.append(session, ignore_index=True)
        """ NEJA handled by the walk method
        # Recurse to the appropriate depth
        depth = None if depth is None else depth - 1
        if depth is None or depth > 0:
            for subDir in curDir.dirs():
                self.sysScan(subDir, depth)
        """

    def generateRE(self):
        """
        Create a list of regular expressions to use in anonymizing data.

        Returns
        -------
        A list of lookahead regular expressions based upon PNT file data.
        """
        regexNames = []
        nameList = list(PNT_FIELDS.keys())[2:]
        for idx in nameList:
            names = self.df[idx][self.df[idx].notna()].unique()
            for name in names:
                splitName = name.split(',')
                regexNames += [f'(?=.*{name}).+' for name in splitName]
        return regexNames

    def addPatient(self, patient):
        """Add a patient object to the corpus."""
        self.patients.append(patient)

    def __init__(self, mrnPath='mrn.txt'):
        with open(mrnPath, 'r') as f:
            self.mrns = [mrn.strip() for mrn in f]
        assert len(self.mrns) > 0

        # Build up DataFrame of original recorded runs
        try:
            self.df = pd.read_csv(CSV_FILENAME)
        except FileNotFoundError:
            warnings.warn("No previous data appear to have been recorded.\n"
                          "Defaulting to a blank DataFrame.")
            colns = ['Path', 'LocalPath', 'Notes'] + list(PNT_FIELDS)
            self.df = pd.DataFrame(columns=colns)
        self.patients = []

    def __iter__(self):
        """Iterate over patients."""
        self.pntIdx = 0
        return self

    def __next__(self):
        """Get the next patient."""
        if self.pntIdx < len(self.patients):
            retVal = self.patients[self.pntIdx]
            self.pntIdx += 1
            return retVal
        else:
            raise StopIteration

    def __len__(self):
        """Return number of unique patients."""
        return len(self.patients)


class Patient:
    """Each patient in the Corpus has a unique object representing them."""

    def anonymizePDF(self, pdf, name):
        """Redact data from the clinical notes."""
        if isinstance(pdf, Clinical_Notes):
            pdf.redact()
            # Clear metadata from file
            pdf.set_metadata({i: '' for i in pdf.metadata.keys()})
            pdf.save(name)
            logging.info(f"{pdf} anonymized as {name}.")

    def correlate(self):
        """Correlate which clinical notes go with which EEG recordings."""
        # NEJA what should the consequences be for not finding a pdf?
        for pdf in self.notes:
            start = pdf.sessionDates[0]
            end = pdf.sessionDates[1]
            validEntries = (self.df.Time > start) & (self.df.Time < end)
            recordIdx = self.df.index[validEntries].tolist()
            for idx in recordIdx:
                self.df.loc[idx, 'Notes'] = pdf
        # Sort by increasing dates for easier saving in BIDS format
        self.df.sort_values('Time', inplace=True, ignore_index=True)

    def __init__(self, mrn, recording_info, notes):
        self.mrn = mrn
        self.df = recording_info.reset_index()
        self.notes = []
        for note in notes:
            self.notes.append(Clinical_Notes(mrn, note))
        self.anonID = 0  # Anonymous UID
        birthday = self.df.DOB[0]
        self.birthday = [int(i) for i in birthday.split('/')]
        self.nlpNames = pd.Series(dtype=str)

    def __repr__(self):
        """Return Patient w/ MRN."""
        return f"Patient object (MRN: {self.mrn})"


if __name__ == '__main__':
    corpus = CorpusBuilder()
    # NEJA TEMP grabchart(corpus.mrns)  # Pull reports from PowerChart
    corpus.sysScan(SHARED_EEG_DRIVE)  # Scan for PNT files

    # Potentially reduce the MRNs to those that were found during the scan?
    corpus.mrns = corpus.df.MRN.unique().tolist()
    corpus.df.to_csv(CSV_FILENAME, index=False)

    # Keep sensitive names from PNT files for all patients
    Clinical_Notes.sensitiveNames = corpus.generateRE()
    startTime = time.perf_counter()
    for mrn in corpus.mrns:
        notes_path = RELATIVE_PATH / str(mrn)
        note_names = notes_path.abspath().glob('*.pdf')
        recording_data = corpus.df.loc[corpus.df.MRN == mrn]
        patient = Patient(mrn, recording_data, note_names)
        patient.correlate()
        corpus.addPatient(patient)
        activeStage = STAGING_DIR / str(mrn)
        if not activeStage.exists():
            Path.makedirs(activeStage)
        startTime = time.perf_counter()
        for idx, nkrecording in enumerate(patient.df.Path):
            nkrecording = Path(nkrecording)
            recordingName = nkrecording.basename()[:-4]
            nkPath = nkrecording.dirname()
            # NEJA TEMP copyNKT(nkPath, activeStage, recordingName)
            patient.df.loc[idx, 'LocalPath'] = activeStage / recordingName
        stopTime = time.perf_counter() - startTime
        print(f"Saved NK files locally in {stopTime}sec.")
    # Construct the filenames and data structures for the BIDS format
    # NEJA NOTE: Consider parallel processing by patient?
    for i, patient in enumerate(corpus, 1):  # BIDS indices start at 1
        patient.anonID = i
        subject = f"{patient.anonID:{SUB_PAD}}"
        for j, pdf in enumerate(patient.df.Notes.unique(), 1):
            session = f"{j:{SES_PAD}}"
            patientSession = patient.df[patient.df.Notes == pdf].reset_index()
            for run, entry in patientSession.iterrows():
                localPath = entry.LocalPath
                nkFName = localPath / '.'.join([localPath.basename(), 'EEG'])
                nkData = mne.io.read_raw_nihon(nkFName, verbose=True)
                if run == 0:  # Modify the record date to be Jan 1st
                    delta = relativedelta(months=patient.birthday[1],
                                          days=patient.birthday[2])
                newDate = nkData.info['meas_date'] - delta
                # NEJA some recordings overlap (and are redundant). Why?
                data = convert_NKT2EDF(nkData, newDate, patient.birthday[0])
                bidsPath = mne_bids.BIDSPath(subject=subject,
                                             session=session,
                                             task='rest',  # NEJA see below
                                             run=run+1,
                                             root=BIDS_DB_PATH)
                # BIDSPath requires a task-<label> note as per spec sheet
                # From BIDS spec:
                #  Note that the TaskName field does not have to be a
                #  ”behavioral task” that subjects perform, but can reflect
                #  some information about the conditions present when the data
                #  was acquired (for example, "rest", "sleep", or "seizure")
                mne_bids.write_raw_bids(data, bidsPath)

                print('---\n\n---')

            pdfName = f"{bidsPath.basename[:16]}_clinical_notes.pdf"
            patient.anonymizePDF(pdf, bidsPath.directory / pdfName)
    runtime = time.perf_counter() - startTime
    print(f"Wrote BIDS data in {runtime} seconds.")

"""
TODO Remove short clips if they are indeed redundant.
TODO Determine why the generated EDF file is being read twice!
     Happens in mne_bids.write_raw_bids...
TODO Better redaction for pdfs (e.g., currently missing a random name)
TODO BIDS annotation anonymization
TODO Where to add parallelization?
DONE Determine if NKT files that overlap are redundant.
     May be clips to look at specific events.
DONE Adjust data so month, day, and year are relative to patient BD
DONE Adjust Patient BD: YY/01/01 (where YY is their *actual* year of birth)
DONE Adjust anonymized data so year is relative to patient BD
DONE Save clinical notes as plaintext...?
    PowerChart doesn't really allow for this (can be html, rtx, etc. but isn't
    consistent). Plaintext could be dumped by pyMuPDF, but headers/footers
    will make it ugly.
DONE Save redacted PDFs along with BIDS datasets
DONE Recording should be 12hrs? Appears to be 6
DONE Define proper start times for contiguous files (currently always Jan 1st)
    All edf files created are assigned a date relative to the start of the
    session.

Interesting sidenote: it is encouraged to keep the orig format in a directory
    called ./sourcedata
"""
