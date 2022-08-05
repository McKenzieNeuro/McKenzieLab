"""
A navigation class to aid in getting around GUIs.

This is particularly relevent for PowerChart.

@author: Nick Aase (aase@cs.unm.edu)

Notes
-----
    22/05/09 (NEJA): Approximate date of creation
    22/07/18 (NEJA): Trying to make code more flexible on different screens.
    22/08/02 (NEJA): Removed sleep wrapper--was not helpful.
"""
import pyautogui
import sys
from time import sleep
import pygetwindow


class Navigator:
    """Class for performing more coordinated activities with pyautogui."""

    @staticmethod
    def click(button='left', doubleClick=False, leadDelay=0, tailDelay=0):
        """
        Add optional delay before and/or after each click.

        Parameters
        ----------
        button : str, optional
            Specifies left or right click. The default is 'left'.
        doubleClick : bool, optional
            Perform a double-click. The default is False.
        leadDelay : int, optional
            Time to delay (in sec) before clicking. The default is 0.
        tailDelay : int, optional
            Time to delay (in sec) after clicking. The default is 0.

        Raises
        ------
        SyntaxError
            If button parameter isn't 'left' or 'right'.

        Returns
        -------
        None.

        """
        sleep(leadDelay)
        if doubleClick:
            pyautogui.doubleClick()
        else:
            if button.lower() == 'left':
                pyautogui.click()
            elif button.lower() == 'right':
                pyautogui.rightClick()
            else:
                raise SyntaxError(f"'{button}' undefined in this scope.")
        sleep(tailDelay)

    def moveToImg(self,
                  img: str,
                  dx=0,
                  dy=0,
                  confidence=0.999) -> pyautogui.Point:
        """
        Find a location based off an image and possibly move from there.

        Parameters
        ----------
        img : str
            Filename of an image to find and click on the screen.
        dx : int, optional
            Move relative to the center of the image. The default is 0.
        dy : int, optional
            Move relative to the center of the image. The default is 0.
        confidence: float, optional
            Confidence in image matching. The default is 0.999 as per OpenCV.

        Returns
        -------
        Point
            Coordinates of cursor on screen if successful.
        """
        try:
            self.coords = pyautogui.locateCenterOnScreen(
                img, confidence=confidence)
        except NotImplementedError:
            self.coords = pyautogui.locateCenterOnScreen(img)
        pyautogui.moveTo(self.coords)
        pyautogui.moveRel(dx, dy)
        return self.coords

    def moveToWin(self, name: str, dx=0, dy=0, ignoreMissing=False):
        """
        Find a location from the substring of the window name.

        Window manipulation relies on pygetwindow, which doesn't yet
        support linux as of June 2022 (not sure about mac).

        Parameters
        ----------
        name : str
            The name of a window to find and click on.
        dx : int, optional
            Move relative to the center of the image. The default is 0.
        dy : int, optional
            Move relative to the center of the image. The default is 0.
        ignoreMissing : bool, optional
            Suppress error if a window isn't found. The default is False.

        Raises
        ------
        NameError
            When a window with the specified name isn't present.
        AssertionError
            If the navigator is not being used on a windows machine.

        Returns
        -------
        Window object
            First window found corresponding to the name argument.
        """
        assert sys.platform == 'win32'
        self.win = None
        for window in pyautogui.getAllWindows():
            if name in window.title:
                self.win = window
                self.coords = self.win.center
                if self.win.isMinimized:
                    self.win.restore()
                    sleep(0.5)
                pyautogui.moveTo(self.coords)
                pyautogui.moveRel(dx, dy)
                try:
                    self.win.activate()
                except pygetwindow.PyGetWindowException:
                    self.click()
                break

        if not ignoreMissing and self.win is None:
            raise NameError(f'No "{name}" window present!')
        else:
            return self.win

    def maxWin(self):
        """
        Maximize the current window.

        Pyautogui's maximize() method is not always functional, so this
        circumvents the problem by maximizing using the upper-right icons.
        """
        if pyautogui.size()[0] > self.win.topright[0]:
            # Window is not centered and/or maximized
            pyautogui.moveTo(self.win.topright)
            pyautogui.moveRel(-70, 20)
            self.click()

        elif pyautogui.size()[0] < self.win.topright[0]:
            # Window is off the screen
            pyautogui.moveTo(self.win.topleft)
            pyautogui.moveRel(70, 20)
            self.click(button='right')
            pyautogui.moveRel(10, 100)
            self.click()

    def closeWin(self):
        """Close the current window via mouse movement."""
        pyautogui.moveTo(self.win.topright)
        pyautogui.moveRel(-25, 20)
        self.click()

    def __init__(self):
        self.win = None
        self.coords = pyautogui.Point(None, None)
