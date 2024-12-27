import ctypes
import time

class Motor():
    """
    Class for controlling the motor.
    """
    def __init__(self, serial_num:str):
        self.motor = ctypes.cdll.LoadLibrary("C:/Program Files/Thorlabs/Kinesis/Thorlabs.MotionControl.KCube.DCServo.dll")
        
        self.serial_num = ctypes.c_char_p(serial_num)
        if self.motor.TLI_BuildDeviceList() == 0:
            self.motor.CC_Open(self.serial_num)
            self.motor.CC_StartPolling(self.serial_num, ctypes.c_int(200))
        else: raise Exception()
        
        self._position: float = 0.0
        # Set up the device to convert real units to device units.
        # These numbers are from the manual.
        STEPS_PER_REV = ctypes.c_double(512.00)  # for the mts50-z8
        gbox_ratio    = ctypes.c_double( 67.49)  # gearbox ratio
        pitch         = ctypes.c_double(  1.00)
        self.motor.CC_SetMotorParamsExt(self.serial_num, STEPS_PER_REV, gbox_ratio, pitch)

        # Home the device if not already homed
        if self.position == 0:
            # the position will be exactly 0 when the device is first turned on.
            # if it has already been homed, it will not home again unless is it at its home position.
            self.motor.CC_Home(self.serial_num)

    def close(self) -> None:
        """ Close the connection to the KCube. """
        self.motor.CC_Close(self.serial_num)
        
        # this is for when testing the motor when it is not connected.
        #motor.TLI_UninitializeSimulations()
    
    @property
    def position(self) -> float:
        """ Returns the position of the motor as a float value, in millimeters. """
        
        # get position in device units
        self.motor.CC_RequestPosition(self.serial_num)
        time.sleep(0.2)
        dev_pos = ctypes.c_int(self.motor.CC_GetPosition(self.serial_num))
        
        # convert device units to real units
        real_pos = ctypes.c_double()
        self.motor.CC_GetRealValueFromDeviceUnit(self.serial_num, dev_pos, ctypes.byref(real_pos), 0)

        # correct _position to match real units
        self._position = real_pos.value
    
        return self._position
    
    @position.setter
    def position(self,value: float) -> None:
        """
        Increments the stage to a new position.
        ## Parameters
        - value (float): The new position in millimeters.
        """
        
        # find new desired new position in real units
        new_position_mm = ctypes.c_double(value) 
        
        # find new direst position in device units
        new_position_device_units = ctypes.c_int()
        self.motor.CC_GetDeviceUnitFromRealValue(self.serial_num, new_position_mm, ctypes.byref(new_position_device_units), 0)

        # move to new position (time.sleep(0.2) is to give time for the command to transfer)
        self.motor.CC_SetMoveAbsolutePosition(self.serial_num, new_position_device_units)
        time.sleep(0.2)
        self.motor.CC_MoveAbsolute(self.serial_num)
        time.sleep(0.2)
        
        # adjust _position to match the new position in real units
        self._position = new_position_mm.value

