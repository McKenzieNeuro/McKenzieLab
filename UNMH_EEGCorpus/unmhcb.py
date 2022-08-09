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
22/07/22 (NEJA): Moved ClinicalNotes class to unmhcb_io.py as well as static
                 I/O methods from the Corpus class.
22/07/26 (NEJA): Fixed conversion from NKT to EDF.
22/07/28 (NEJA): Adjusted birthdates and recording dates for EDF output files.
                 Now BDs are rolled back to Jan. 1st of the year patients were
                 born, and recording dates are relative to that.
22/07/29 (NEJA): Changed glob to walk in sysScan, which eliminates the need for
                 manual recursion.
22/08/03 (NEJA): Substantial tweaks now that more than one patient's data
                 are being processed. Corpus's addPatient now handles the
                 transfer of the NK data.
"""
from path import Path
import pandas as pd
import filecmp
import os
from dateutil.relativedelta import relativedelta
import mne
import mne_bids
import random
import warnings
import logging
from grabchart import grabchart, NOTES_PATH
from unmhcb_io import get_PNT_metadata, get_MRN_from_PNT, copyNKT, \
    PNT_FIELDS, ClinicalNotes, convert_NKT2EDF

random
# Constants
BIDS_DB_PATH = Path("EEG_Corpus")
SHARED_EEG_DRIVE = Path("//EEG-MAS-GW/EEG/")
STAGING_DIR = Path('nk_staging_data')
CSV_FILENAME = 'PNT_metadata.csv'
SUB_PAD = '05'
SES_PAD = '02'

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class CorpusBuilder:
    """The main datastructure for building the EEG Corpus."""

    def sysScan(self, fullPath: str):
        """
        Scan a file system for EEG data corresponding to a list of MRNs.

        Parameters
        ----------
        fullPath : str
            Top-most directory to begin searching.

        Returns
        -------
        None.

        """
        logger.info(f"Scanning {fullPath} for patients' Nihon Koden files")
        curDir = Path(fullPath)
        # NEJA NOTE: Maybe put something to avoid re-scanning directories...?
        for file in curDir.walk('*.PNT'):  # Just consider PNT files
            if os.stat(file).st_size == 0:
                continue  # Skip blank files
            mrn = get_MRN_from_PNT(file)
            if mrn in self.mrns:
                logger.info(f"Found {mrn} at {file}")
                repeat = False
                for saved_file in self.df.Path:
                    if filecmp.cmp(file, saved_file, shallow=False):
                        # Avoid writing the same data multiple times
                        repeat = True
                        warnings.warn(f"Duplicate files found: "
                                      f"{file} & {saved_file})!!")
                        break
                if not repeat:
                    metadata = get_PNT_metadata(file) + [file]
                    session = {i: j for i, j in zip(self.df.columns, metadata)}
                    self.df = self.df.append(session, ignore_index=True)

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

    def addPatient(self, mrn, overwrite=False):
        """Add a patient object to the corpus."""
        notes_path = NOTES_PATH / str(mrn)
        note_names = notes_path.abspath().glob('*.pdf')
        recording_data = self.df.loc[self.df.MRN == mrn]
        patient = Patient(mrn, recording_data, note_names)
        patient.anonID = len(self.patients) + 1  # BIDS indices start at 1
        patient.correlate()
        activeStage = STAGING_DIR / str(mrn)
        if not activeStage.exists():
            Path.makedirs(activeStage)
        for idx, nkrecording in enumerate(patient.df.Path):
            nkrecording = Path(nkrecording)
            recordingName = nkrecording.basename()[:-4]
            # Copy the recordings off of the current DB to lighten the load.
            # (This may ultimately be unnecessary.)
            nkPath = nkrecording.dirname()
            copyNKT(nkPath, activeStage, recordingName, overwrite)
            logger.info(f"Copied NK dataset {idx+1} of {len(patient.df.Path)}")
            patient.df.loc[idx, 'LocalPath'] = activeStage / recordingName
        self.patients.append(patient)

    def writeBIDS(self):
        """
        Construct the filenames and data structures in the BIDS format.

        Will read each patient's data, perform anonymization, then finally
        and add to the BIDS-formatted DB.

        Note: BIDS indices start at 1

        Returns
        -------
        None.

        """
        # NEJA NOTE: Consider PP here
        for patient in self.patients:
            subject = f"{patient.anonID:{SUB_PAD}}"
            for j, pdf in enumerate(patient.df.Notes.unique(), 1):
                session = f"{j:{SES_PAD}}"
                if pdf is not None:
                    mask = patient.df.Notes == pdf
                else:
                    mask = patient.df.Notes.isnull()
                patientSession = patient.df[mask].reset_index(drop=True)
                for run, entry in patientSession.iterrows():
                    localPath = entry.LocalPath
                    nkFName = localPath / '.'.join([
                                                    localPath.basename(),
                                                    'EEG'
                                                    ])
                    nkData = mne.io.read_raw_nihon(nkFName, verbose=True)
                    if run == 0:  # Modify the record date to be Jan 1st
                        delta = relativedelta(months=patient.birthday[1],
                                              days=patient.birthday[2])
                    newDate = nkData.info['meas_date'] - delta
                    bidsPath = mne_bids.BIDSPath(subject=subject,
                                                 session=session,
                                                 task='rest',  # NEJA see below
                                                 run=run+1,
                                                 root=BIDS_DB_PATH)
                    # BIDSPath requires a task-<label> note as per spec sheet
                    # From BIDS spec:
                    #  Note that the TaskName field does not have to be a
                    #  ”behavioral task” that subjects perform, but can reflect
                    #  some information about the conditions present when
                    #  the data was acquired (for example, "rest", "sleep",
                    #  or "seizure")
                    data = convert_NKT2EDF(nkData,
                                           newDate,
                                           patient.birthday[0])
                    mne_bids.write_raw_bids(data, bidsPath)

                print('---\n\n---')
                if pdf is not None:
                    pdfName = f"{bidsPath.basename[:16]}_ClinicalNotes.pdf"
                    patient.anonymizePDF(pdf, bidsPath.directory / pdfName)

    def __init__(self, mrnPath='mrn.txt'):
        with open(mrnPath, 'r') as f:
            self.mrns = [mrn.strip() for mrn in f]
        assert len(self.mrns) > 0

        # Build up DataFrame of original recorded runs
        try:
            self.df = pd.read_csv(CSV_FILENAME)
        except FileNotFoundError:
            warnings.warn("No previous data appear to have been recorded. "
                          "Defaulting to a blank DataFrame.")
            colns = list(PNT_FIELDS) + ['Path', 'LocalPath', 'Notes']
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

    def __repr__(self):
        """Includes info about how many patients and recordings."""
        return f"CorpusBuilder object (containing {len(self.patients)}"\
            f" patients and {self.df.shape[0]} remote recordings)"


class Patient:
    """Each patient in the Corpus has a unique object representing them."""

    def anonymizePDF(self, pdf, name):
        """Redact data from the clinical notes."""
        if isinstance(pdf, ClinicalNotes):
            pdf.redact()
            # Clear metadata from file
            pdf.set_metadata({i: '' for i in pdf.metadata.keys()})
            pdf.save(name)
            logger.info(f"{pdf} anonymized as {name}.")

    def correlate(self):
        """Correlate which clinical notes go with which EEG recordings."""
        for pdf in self.notes:
            start = pdf.sessionDates[0]
            end = pdf.sessionDates[1]
            validEntries = (self.df.Time >= start) & (self.df.Time <= end)
            recordIdx = self.df.index[validEntries].tolist()
            for idx in recordIdx:
                self.df.loc[idx, 'Notes'] = pdf
        # Sort by increasing dates for easier saving in BIDS format
        self.df.sort_values('Time', inplace=True, ignore_index=True)
        for idx, row in self.df.iterrows():
            if pd.isna(row.Notes):
                warnings.warn(f"Patient {self.mrn} has no known notes "
                              f"associated with recording {row.Path} "
                              "Consider checking PowerChart again.")
                self.df.loc[idx, 'Notes'] = None

    def __init__(self, mrn, recording_info, notes):
        self.mrn = mrn
        self.df = recording_info.reset_index(drop=True)
        self.notes = []
        for note in notes:
            clinicalDoc = ClinicalNotes(mrn, note)
            if len(clinicalDoc.sessionDates) == 2:
                self.notes.append(clinicalDoc)
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
    # NEJA TEMP corpus.sysScan(SHARED_EEG_DRIVE)  # Scan for PNT files

    # Potentially reduce the MRNs to those that were found during the scan?
    corpus.mrns = corpus.df.MRN.unique().tolist()
    corpus.df.to_csv(CSV_FILENAME, index=False)

    # Keep sensitive names from PNT files for all patients
    ClinicalNotes.sensitiveNames = corpus.generateRE()
    for mrn in corpus.mrns:
        if mrn != 5397382:
            corpus.addPatient(mrn)

    corpus.writeBIDS()

"""
TODO Randomize BIDs subject ID
TODO Remove short clips if they are indeed redundant.
TODO Better redaction for pdfs (e.g., currently missing a random name)
TODO Annotation anonymization?
TODO Where to add parallelization? (Should be procs due to I/O load)
DONE Get rid of double-indexed DF in patient data
DONE Inhibit double writes / reads of EDF data
     Don't think it's going to happen with mne_bids the way it is.
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
