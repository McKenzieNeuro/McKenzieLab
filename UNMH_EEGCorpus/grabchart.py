# -*- coding: utf-8 -*-
"""
A script to select Neuroclinical data from a list of patients in PowerChart.

(The script name hails from a play on words of LabChart, though...)

@author: Nick Aase (aase@cs.unm.edu)

Requirements
------------
pyautogui
OpenCV (Not 100% necessary, but REALLY helps.)

Notes
-----
    22/06/09 (NEJA): Added optional confidence in image recognition, as it
        seemed to fail frequently. Began implementation of appending functions.
    22/06/10 (NEJA): Modularizing everything, making the navigator more
        navigatory, corrected for some edge cases where windows weren't
        selectable or files already existed.
    22/06/13: Saved off some snapshots of screens where the script fails.
        They could be integrated into the code to make the functionality
        more robust.
    22/06/30 (NEJA): Consolidated the helper functions that used to be in
        another file.
    22/07/21 (NEJA): Added further precision to finding the Neurodiagnostics
        folder.
"""
from sys import stderr
import os
import time
from path import Path
from navigator import Navigator
import warnings
import logging
# Check dependencies
try:
    import pyautogui
except ImportError:
    raise ImportError("Make sure pyautogui is installed!")

psVersion = [int(i) for i in pyautogui.pyscreeze.__version__.split('.')]
if psVersion[1] <= 1 and psVersion[2] < 28:
    # Pixel detection is borked in old versions
    print("Old pyscreeze version detected!"
          "Try running pip install -U pyscreeze",
          file=stderr)
    raise ImportError("pyscreeze version must be at least 0.1.28!")


# Constants
MIN_RECORD_DATE = '08/01/2007'
F_MIN_RECORD_DATE = 1185951600.0  # Float representation of 8/1/2007
PDF_PATH = "\\".join(["\\\\Client", Path.getcwd(), "clinical_notes"])
RELATIVE_PATH = Path('clinical_notes')
WINDOW_TITLES = {'main': 'PowerChart Organizer for',
                 'relationship': 'Assign a Relationship',
                 'print': 'Print',
                 'save': 'Save Print Output As',
                 'error': 'Confirm Save As'}


def findPatient(mrn):
    """Bring up patient search."""
    time.sleep(1)
    almostTopright = (nav.win.topright[0]-1, nav.win.topright[1]+1)
    pyautogui.moveTo(almostTopright)
    pyautogui.moveRel(-100, 110)
    nav.click()
    pyautogui.typewrite(mrn)
    pyautogui.moveRel(65, 0)
    nav.click()
    time.sleep(1)
    nav.moveToWin("Patient Search")
    if nav.win:
        pyautogui.moveRel(0, 50)
    else:
        nav.moveToWin(WINDOW_TITLES['main'])
        pyautogui.moveTo(nav.win.center)
        pyautogui.moveRel(0, 150)
    nav.click(doubleClick=True, leadDelay=1, tailDelay=3)


def loadClinicalNotes(mrn):
    """Select Clinical Notes."""
    if mrn not in nav.win.title:  # Check to see if a new window is open
        time.sleep(4)
        nav.moveToWin(f'{mrn} Opened by')
        nav.maxWin()
    pyautogui.moveTo(213, 37)  # Chart dropdown (spacing is not relative)
    nav.click(leadDelay=2)
    pyautogui.moveRel(0, 270)
    nav.click(leadDelay=2)   # NEJA may need to be much longer
    #  Establish new directory for patient
    try:
        os.mkdir(f"{RELATIVE_PATH}\\{mrn}")
    except FileExistsError:
        warnings.warn(
            f"Patient {mrn}'s data appears to have been saved already.")

    logging.info('Adjust search dates')
    pyautogui.moveRel(100, -25)
    nav.click(button='right', leadDelay=2)
    pyautogui.moveRel(5, 5)
    nav.click(leadDelay=2)
    pyautogui.moveTo(nav.win.center)
    pyautogui.moveRel(100, -20)
    nav.click(doubleClick=True)
    latest_date = MIN_RECORD_DATE
    patientPath = f'{RELATIVE_PATH}/{mrn}/'
    files = os.listdir(patientPath)
    if len(files) > 0:  # Patient's data has been seen before
        floatTime = os.path.getmtime(patientPath + files[-1])
        mostRecent = max(floatTime, F_MIN_RECORD_DATE)
        timeStruct = time.strptime(time.ctime(mostRecent))
        latest_date = f"{timeStruct[1]:02}/{timeStruct[2]:02}/{timeStruct[0]}"
    pyautogui.typewrite(latest_date)
    pyautogui.hold('enter')
    pyautogui.moveRel(-50, 123)
    nav.click(tailDelay=2)


def selectNeurodiagnostics(mrn):
    """Select Neurodiagnostics folder image."""
    success = True
    fieldLocation = nav.moveToImg('img/neuro_folder.png')
    counter = 0
    while fieldLocation is None and counter < 5:  # Check a max of 5 fields
        # Check another image if the first check fails
        time.sleep(1.5)
        logging.info("Trying another screenshot for navigation")
        fieldLocation = nav.moveToImg('img/radio_buttons.png')
        if fieldLocation:
            break
        pyautogui.moveTo(nav.win.left, nav.win.centery)
        pyautogui.moveRel(250, 0)
        nav.click()
        pyautogui.press('n')
        fieldLocation = nav.moveToImg('img/neuro_selected0.png')
        if not fieldLocation:
            fieldLocation = nav.moveToImg('img/neuro_selected1.png')
        if not fieldLocation:
            fieldLocation = nav.moveToImg('img/neuro_selected2.png')
        counter += 1

    if fieldLocation is None:  # Skip, but note, any unprocessed MRNs
        warnings.warn("Cannot find Neurodiagnostics folder via screenshot."
                      "Are there any records for this time interval?")
        success = False
    else:
        #  Bring up Print dialog
        nav.click(button='right', leadDelay=2, tailDelay=2)
        pyautogui.moveRel(20, 40)
        nav.click(tailDelay=5)
    return success


def printToPDF(mrn):
    """Choose print to PDF."""
    idx = 0
    notePath = f'{RELATIVE_PATH}{mrn}/'
    if os.path.exists(notePath):  # Patient's data has been seen before
        files = os.listdir(notePath)
        idx = len(files)
    nav.moveToWin(WINDOW_TITLES['print'], ignoreMissing=True)
    while nav.win is not None:
        validDialog = False
        if nav.win is not None:
            validDialog = True
            pyautogui.moveRel(300, 0)
            fieldLocation = None
            for _ in range(5):
                time.sleep(2)
                pyautogui.press('m')
                fieldLocation = nav.moveToImg("img/print_to_pdf0.png",
                                              confidence=0.7)
                if fieldLocation is None:
                    logging.info(
                        "Initial image movement failed, trying another.")
                    fieldLocation = nav.moveToImg("img/print_to_pdf2.png",
                                                  confidence=0.7)

                if fieldLocation is None:
                    logging.info(
                        "Secondary image movement failed. Trying another")
                    fieldLocation = nav.moveToImg("img/print_to_pdf3.png",
                                                  confidence=0.7)

                if fieldLocation:
                    break

            if not fieldLocation:  # Try just hitting enter
                warnings.warn("Microsoft PDF print option is not confirmed...")
                pyautogui.press('enter')
            else:
                pyautogui.moveTo(fieldLocation)
                time.sleep(1)
                nav.click(doubleClick=True)

        while True:
            nav.win = None
            nav.moveToWin(WINDOW_TITLES['save'], ignoreMissing=True)
            if nav.win is None and validDialog:
                time.sleep(1)  # The save dialog can take a while to appear
            else:
                break

        # Save documents
        if validDialog:
            filename = f"{PDF_PATH}\\{mrn}\\{mrn}_{idx}"
            logging.info(f"Writing {filename}.pdf")
            pyautogui.typewrite(filename)
            pyautogui.moveTo(nav.win.bottomright)
            pyautogui.moveRel(-50, -50)
            nav.click(tailDelay=3)
            # Perhaps we'll want to overwrite files, but throw error for now
            nav.moveToWin(WINDOW_TITLES['error'], ignoreMissing=True)
            if nav.win is not None:  # Odd way to check for overwriting, but...
                raise FileExistsError("Trying to overwrite an existing file!")
            else:
                idx += 1
                nav.moveToWin(WINDOW_TITLES['print'], ignoreMissing=True)
    return idx


nav = Navigator()


def grabchart(mrns: list):
    """Perform the primary processing."""
    global nav
    try:
        os.mkdir(f"{RELATIVE_PATH}")
    except FileExistsError:
        # This isn't a bad thing
        logging.info("Path to notes exists.")

    for i, mrn in enumerate(mrns):
        logging.info(f"Recording notes for patient {mrn} "
                     f"({i+1} of {len(mrns)})")
        try:
            nav.moveToWin(WINDOW_TITLES['main'])
        except NameError:
            raise RuntimeError("PowerChart can't be found! Is it running?")
        mrn = str(mrn)
        nav.win.activate()
        nav.maxWin()
        findPatient(mrn)

        #  Choose research role if the dialog appears (it frequently doesn't)
        nav.moveToWin(WINDOW_TITLES['relationship'], ignoreMissing=True)
        if nav.win is not None:
            pyautogui.moveRel(0, -20)
            nav.click(leadDelay=1)
            pyautogui.moveRel(75, 150)
            nav.click()

        nav.moveToWin(WINDOW_TITLES['main'])
        loadClinicalNotes(mrn)
        if not selectNeurodiagnostics(mrn):
            with open('unsaved_mrns.txt', 'a+', encoding='utf-8') as f:
                f.write(f"{mrn}\n")
            continue
        numSaved = printToPDF(mrn)
        logging.info(f"Successfully saved {numSaved} notes for patient {mrn}")
        
    logging.info("Process Complete")


"""
Issues:
    Opened 22/06/10:
        If main window is obscured, it won't be brought to the front
    Closed 22/06/10:
        Added a check to see if window is minimized; seems to have taken
        care of any oddities.
    Opened 22/06/10:
        Directory for patient created before an open window is confirmed...
    Closed 22/06/13:
        Flipped the ordering while modularizing the code.
    Opened 22/07/21:
        Since this is finicky, add a feature that allows for the user to
        pick up where they left off?
"""
