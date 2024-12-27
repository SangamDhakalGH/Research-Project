import pyvisa
import time

class SignalGenerator:
    """ This class is used to control the signal generator. """

    def __init__(self) -> None:
        # establishes connection to the instrument
        self.rm    = pyvisa.highlevel.ResourceManager()
        self.instr = self.rm.open_resource('USB0::5710::4002::TW00050162::INSTR')
        
        # clears previous instructions
        self.instr.write('*CLR') # type: ignore

        # set the signal generator to default conditions
        self._waveform     = "PULSE"
        self._pulse_width  = 160e-6
        self._low_voltage  = 0
        self._high_voltage = 1.7
        self._edge_time    = 5e-9
        self._frequency    = 1000

    def open(self) -> None:
        """
        Apply the default conditions to the signal generator.
        """
        
        self.waveform     = "PULSE"
        self.pulse_width  = 160e-6 # seconds
        self.low_voltage  = 0
        self.high_voltage = 1.7
        self.edge_time    = 5e-9 # seconds
        self.frequency    = 1000 # Hz

    def close(self) -> None:
        self.instr.close()
        self.rm.close()
    
    @property
    def waveform(self) -> str:
        self._waveform = self.instr.query('FUNCTION?') # type: ignore
        time.sleep(0.1)
        return self._waveform
    
    @waveform.setter
    def waveform(self, value: str) -> None:
        self.instr.write(f'FUNCTION {value}') # type: ignore
        time.sleep(0.1)
        self._waveform = value

    @property
    def pulse_width(self) -> float:
        self._pulse_width = float(self.instr.query('PULSE:WIDTH?')) # type: ignore
        time.sleep(0.1)
        return self._pulse_width
    
    @pulse_width.setter
    def pulse_width(self, value: float) -> None:
        self.instr.write(f'PULSE:WIDTH {value}') #type: ignore
        time.sleep(0.1)
        self._pulse_width = value

    @property
    def edge_time(self) -> float:
        self._edge_time = float(self.instr.query('PULSE:TRANSITION?')) # type: ignore
        time.sleep(0.1)
        return self._edge_time

    @edge_time.setter
    def edge_time(self, value: float) -> None:
        self.instr.write(f'PULSE:TRANSITION {value}') # type: ignore
        time.sleep(0.1)
        self._edge_time = value

    @property
    def frequency(self) -> float:
        self._frequency = float(self.instr.query('FREQUENCY?')) # type: ignore
        time.sleep(0.1)
        return self._frequency

    @frequency.setter
    def frequency(self, value: int) -> None:
        self.instr.write(f'FREQUENCY {value}') # type: ignore
        time.sleep(0.1)
        self._frequency = value

    @property
    def high_voltage(self) -> float:
        self._high_voltage = float(self.instr.query('VOLTAGE:HIGH?'))
        time.sleep(0.1)
        return self._high_voltage

    @high_voltage.setter
    def high_voltage(self, value: float) -> None:
        self.instr.write(f'VOLTAGE:HIGH {value}') # type: ignore
        time.sleep(0.1)
        self._high_voltage = value

    @property
    def low_voltage(self) -> float:
        self._low_voltage = float(self.instr.query('VOLTAGE:LOW?'))
        time.sleep(0.1)
        return self._low_voltage

    @low_voltage.setter
    def low_voltage(self, value: float) -> None:
        self.instr.write(f'VOLTAGE:LOW {value}') # type: ignore
        time.sleep(0.1)
        self._low_voltage = value

    def output_on(self) -> None:
        self.instr.write('OUTPUT ON') # type: ignore
        time.sleep(0.1)

    def output_off(self):
        self.instr.write('OUTPUT OFF') # type: ignore
        time.sleep(0.1)