import pyvisa
import numpy as np

class TLDC2200:
    def __init__(self):
        # Opens a resource manager
        self.rm = pyvisa.highlevel.ResourceManager()
        # Opens the connection to the TLDC2200 LED controller. The variable instr is the handle for the device.
        self.instr = self.rm.open_resource('USB0::0x1313::0x80C8::M00804408::INSTR')
        
        self.instr.write("*CLS") # clears any instructions and/or errors from previous run
        self.instr.write("SOURCE1:MODE CC") # sets the device to constant current mode
        
        self._current: float = 10
        
    def generate_current_list(self, min: float, max: float, num: int, spacing_type: int) -> np.ndarray:
        """
        Simple function that returns an evenly spaced list of the light intensities in Amps.
        ## Inputs
        - min: The minimum light intensity in mA.
        - max: The maximum light intensity in mA.
        - num: The number of light intensity points.
        - spacing_type: The type of spacing between light intensity points.
            - 0: linear spacing
            - 1: logarithmic spacing
        ### Note
        This function calls the function `numpy.geomspace()` requires numpy version >= 1.16.0.
        """

        if spacing_type==0:
            return np.linspace(min/1000, max/1000, num)
        else:
            return np.geomspace(min/1000, max/1000, num)

    def close(self):
        """
        Close the handle to the device and the resource manager.
        
        The device will remain in REMOTE mode until the USB is disconnected or LOCAL button is pressed.
        When in REMOTE mode, most buttons on the console will be locked.
        """
        
        self.instr.close()
        self.rm.close()

    def on(self):
        self.instr.write("OUTPUT1:STATE ON") # type: ignore

    def off(self):
        self.instr.write("OUTPUT1:STATE OFF") # type: ignore

    @property
    def current(self) -> float:
        return float(self.instr.query("SENSE3:CURRENT:DATA?")) # type: ignore
    
    @current.setter
    def current(self, value: float) -> None:
        self._current = value
        self.instr.write(f'SOURCE1:CCURENT:CURRENT {value}') # type: ignore
        