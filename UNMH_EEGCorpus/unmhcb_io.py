# -*- coding: utf-8 -*-
"""
Handle I/O functionality for the UNMH Corpus Builder.

@author: Nick Aase (naase@cs.unm.edu)

Notes
-----
22/07/22 (NEJA): Initial build
22/07/26 (NEJA): Added NKT to EDF file conversion.
22/07/28 (NEJA): Fixed conversion to anonymize the start date of a session,
                 but maintain time consistancy between runs. (e.g., a 6hr
                 run that starts at 10PM on Aug. 23rd will now show as running
                 from 10PM Jan. 1st to 4AM Jan. 2nd)
"""
import os
import re
from sys import stderr
from datetime import datetime, date
import mne
import numpy as np
import pyedflib
from path import Path
try:
    import fitz  # Used for redacting PDFs
except ImportError:
    raise ImportError('Make sure PyMuPDF is installed.')


PNT_FIELDS = {'Time': (0x40, 14),
              'MRN': (0x600, 16),
              'Name': (0x62C, 24),
              'DOB': (0x65F, 13),
              'RefName': (0x682, 24),
              'Physician': (0x6A5, 24),
              'Operator': (0xB00, 24)}
REGEXES = [r"Result date\s*:\s*.*",
           r"Performed by\s*:\s*.*",
           r"Verified by\s*:\s*.*",
           r"Encounter info\s*:\s*.*",
           r"Printed by:.*",
           r"Printed on:.*",
           r"FIN\s*#*\s*:\s*\d*",
           r"MRN[\s#]*\s*:\s*\d*",
           r"MR\s*#\s*:\s*\d*",
           r"MRN Patient.*[\n\s]+.*",  # Mind Research Network num.
           r"Physician[\s\(\)s]*:.*",
           r"Technologist[\s\(\)s]*:.*\n?.*",
           r"(?=.*M\.?D\.?)[^T].*",
           r"(?=.*Ph\.?D\.?).*",
           r"\d{1,2}/\d{1,2}/\d{2,4}"]
BLACK = (0, 0, 0)
WHITE = (255, 255, 255)


def _readField(f, fieldInfo: tuple) -> str:
    """Read a single field of a PNT file."""
    pos, fieldLen = fieldInfo
    f.seek(pos)
    try:
        field = np.fromfile(f, f'|S{fieldLen}', 1).astype(f'U{fieldLen}')[0]
    except Exception as e:
        print('NEJA DEBUG could only get this error to happen once...',
              file=stderr)
        raise e
    field = field.replace('\x00', '')
    return field.replace(' ', '')


def get_PNT_metadata(file: Path or str, fields=PNT_FIELDS) -> [str]:
    """
    Get metadata from a PNT file.

    Parameters
    ----------
    file : Path or str
        Path to PNT file.
    fields : {str:(int, int)}, optional
        One or more key/value pairs. Key is the name of the field
        (as defined here), and the value is a tuple corresponding to the
        number of bytes into the PNT file the data is and the length of
        the field. The default is PNT_FIELDS.

        Valid keys for fields are 'Time', 'MRN', 'Name', 'DOB', 'RefName',
                                  'Physician', and 'Operator'

    Returns
    -------
    [str]
        Each relevant field's data specified by the fields parameter.
    """
    data = []
    with open(file, 'r') as f:
        for field in fields:
            data.append(_readField(f, PNT_FIELDS[field]))
    return data


def get_MRN_from_PNT(file: Path):
    """Return only the MRN from a valid PNT file."""
    return get_PNT_metadata(file, fields=['MRN'])[0]


def copyNKT(sourceDir, destDir, recName):
    """Copy all NKT data for a run from the EEG DB."""
    if not isinstance(sourceDir, Path):
        sourceDir = Path(sourceDir)
    recordingData = sourceDir.glob(f'*{recName}*')
    destDir = Path(destDir) / recName
    if not destDir.exists():
        os.mkdir(destDir)
    for item in recordingData:
        if item.isdir():
            if item.basename().endswith('VOR'):  # Skip video files
                continue
            else:
                try:
                    item.copytree(destDir / item.basename(),
                                  copy_function=Path.copy)
                except FileExistsError as e:  # NEJA is an error the right choice?
                    raise e("Cannot recreate subdirectories!")
        else:
            item.copy(destDir)


def convert_NKT2EDF(nkData, sessionStart=None, birthYear=1947):
    """
    Save a Nihon Koden-formatted file as EDF.

    After a Nihon Koden-formatted recording has been read by MNE it can be
    passed into this function. It will generate a new temporary EDF file
    with an anonymized name, birthdate and startdate, then return an MNE
    object that is representative of the new file.

    Care should be given, as a new file with potentially sensitive PII is
    created locally!

    Parameters
    ----------
    nkData : mne.io.nihon.nihon.RawNihon
        Loaded Raw MNE data from a NK-formatted file.
    sessionStart : datetime.datetime (optional)
        Offset session start time (helpful in anonymizing). Defaults to None.
    birthYear : int (optional)
        Patient's year of birth. Defaults to 1947. (Because 47.)

    Returns
    -------
    mne.io.edf.edf.RawEDF
        Raw MNE data pulled from a newly created EDF.

    """
    if not isinstance(nkData, mne.io.nihon.nihon.RawNihon):
        raise TypeError('NKT file not in valid MNE form!')
    if not nkData.preload:
        nkData.load_data()
    if sessionStart is None:
        sessionStart = nkData.info['meas_date']
    info = nkData.describe(data_frame=True)
    f = pyedflib.EdfWriter('temp.edf', nkData.info['nchan'])
    chType = nkData.get_channel_types()
    for ch_num, ch_data in info.iterrows():
        f.setLabel(ch_num, f"{chType[ch_num]} {ch_data['name']}")
        f.setPhysicalDimension(ch_num, ch_data['unit'])
        f.setSamplefrequency(ch_num, nkData.info['sfreq'])
        if ch_data['min'] < ch_data['max']:
            f.setPhysicalMinimum(ch_num, ch_data['min'])
            f.setPhysicalMaximum(ch_num, ch_data['max'])
    del info

    f.setStartdatetime(sessionStart)
    f.setPatientName('X')
    f.setBirthdate(date(birthYear, 1, 1))  # All patients were born on Jan 1st
    f.writeSamples(nkData.get_data())

    for annotation in nkData.annotations:
        onset = annotation['onset']
        duration = annotation['duration']
        description = annotation['description']
        f.writeAnnotation(onset, duration, description)

    f.close()
    return mne.io.read_raw_edf('temp.edf', infer_types=True)


class Clinical_Notes(fitz.Document):
    """Subclass of fitz.Document that implements anonymizing."""

    sensitiveNames = ''

    def getSensitiveTxt(self, text: str) -> [str]:
        """
        Find sensitive data in a string.

        Parameters
        ----------
        text : str
            Text to search. Typically one page of a PDF report in this context.

        Returns
        -------
        [str]
            A list of all the substrings in text that match:
                The REGEXES global list,
                The patient's MRN,
                Any names identified in the preamble of the PDF,
                Any names identified in the PNT files.
        """
        # First, try to get names from the preamble
        preambleRE = [r"Performed by\s*:\s*.*",
                      r"Verified by\s*:\s*.*"]
        for pattern in preambleRE:
            search = re.search(pattern, text)
            if search:
                reLength = len(pattern[:-8])
                names = search.group()[reLength:]
                search = re.search(' on', names)
                if search:
                    names = names[:search.start()].strip().split(',')
                    self.regexes += [f'(?=.*{name.strip()}).+'
                                     for name in names]

        instances = []
        # Find all instances specified in the regex list
        for pattern in self.regexes:
            instances += re.findall(pattern, text, re.IGNORECASE)

        # Find all sensitive names pulled from the PNT files
        for name in Clinical_Notes.sensitiveNames:
            instances += re.findall(name, text, re.IGNORECASE)

        # Find instances following the Action list segment
        search = re.search(r"Completed Action List", text, re.IGNORECASE)
        if search:
            startIdx = search.end()
            instances += re.findall(r'^\*.*\n??.*',
                                    text[startIdx:],
                                    re.MULTILINE)

        # Clean up whitespace for better accuracy
        instances = [i.strip() for i in instances]
        return instances

    def redact(self):
        """Redact sensitive data from the file page by page."""
        for page in self:
            # wrap_Contents is needed for fixing alignment issues with rect
            # boxes in some cases where there are alignment issues
            page.wrap_contents()
            pageText = page.get_text()
            sensitive = self.getSensitiveTxt(pageText)

            for data in sensitive:
                areas = page.search_for(data)

                # drawing outline over sensitive datas
                for area in areas:
                    page.add_redact_annot(area, fill=BLACK)

            page.apply_redactions()

    def __init__(self, mrn, path):
        self.mrn = mrn
        self.regexes = REGEXES + [str(self.mrn)]
        self.path = str(path)
        super().__init__(self.path)

        # Get date interval of encounter. All related runs should fall
        # between the interval and are considered part of one session.
        self.sessionDates = []
        page = self[0]  # Patient object is an array of pages
        pageText = page.get_text()
        search = re.search(r"(Encounter).*\n.*", pageText, re.IGNORECASE)
        if search:
            date_re = REGEXES[-1]
            self.sessionDates = re.findall(date_re, search.group())
            assert len(self.sessionDates) == 2
            for i, sesDate in enumerate(self.sessionDates):
                try:
                    dt = datetime.strptime(sesDate, '%m/%d/%Y')
                    self.sessionDates[i] = 1000000
                    self.sessionDates[i] *= int(dt.strftime('%Y%m%d'))
                except ValueError:  # Year was not 4 digits in PDF
                    dt = datetime.strptime(sesDate, '%m/%d/%y')
                    self.sessionDates[i] = 20000000000000
                    self.sessionDates[i] += int(dt.strftime('%y%m%d'))
