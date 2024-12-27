import numpy as np
from ctypes import *
import time
from picosdk.ps5000a import ps5000a as ps
from picosdk.functions import adc2mV, assert_pico_ok, mV2adc
from picosdk.constants import PICO_STATUS

class picoscope:
    def __init__(self):
        # enum definitions:
        #channels
        self.PS5000a_CHANNEL_A=0
        self.PS5000a_CHANNEL_B=1
        self.PS5000a_CHANNEL_C=2
        self.PS5000a_CHANNEL_D=3
        #self.PS5000a_EXTERNAL=4
        self.PS5000a_MAX_CHANNELS=4   # self.PS5000_EXTERNAL
        self.PS5000a_TRIGGER_AUX=3
        self.PS5000a_MAX_TRIGGER_SOURCES=2
        
        #voltage ranges
        self.PS5000a_10MV=0
        self.PS5000a_20MV=1
        self.PS5000a_50MV=2
        self.PS5000a_100MV=3
        self.PS5000a_200MV=4
        self.PS5000a_500MV=5
        self.PS5000a_1V=6
        self.PS5000a_2V=7
        self.PS5000a_5V=8
        self.PS5000a_10V=9
        self.PS5000a_20V=10
        self.range_fctr=np.array(([1e-2, 2e-2, 5e-2, 0.1, 0.2, 0.5, 1, 2, 5, 10, 20]),dtype=float)
        
        # trigger direction
        self.PS5000a_TRIG_ABOVE=0
        self.PS5000a_TRIG_BELOW=1
        self.PS5000a_TRIG_RISING=2
        self.PS5000a_TRIG_FALLING=3
        
        # ac/dc coupling
        self.PS5000a_COUPLING_DC=1
        self.PS5000a_COUPLING_AC=0
        
        self.tic=time.time()
        
        # default behavior
        self.fail = False # 
        self.sel_rapid_block = False             # is RapidBlockMode desired?
        self.nseg = 100                    # only for rapid block mode: number of segments
        self.sel_chon = [1, 1]                   # only chA [1,0] or both channels [1,1] on
        self.sel_v = 8                            # selected voltage range
        self.sel_dt = 1e-6                        # sample interval in secs
        self.sel_nsamp = 6000                      # number of samples
        self.sel_oversample = 10                   # number of oversamples
        self.sel_threshold = 1000                  # threshold in mV
        self.sel_autotrig_ms = 200
        self.sel_direction = self.PS5000a_TRIG_RISING
        self.sel_trig_pos_rel = 0.1    # relative position of trigger in collected data
        self.sel_trig_ch = self.PS5000a_CHANNEL_A              
        self.sel_seg_index=0        # the index ofthe memory segment to be used
        self.DATA_READY = WINFUNCTYPE(c_void_p, c_short, c_short, POINTER(c_void_p))
        self.data_ready = self.DATA_READY(self.py_data_ready)
        self.flag_newdata = False
                
    def py_data_ready(self,a,b,c):
        #self.tic=time()
        #print"Data ready has been called"
        #print a 
        #print b 
        #print c 
        self.set_data_buffer()
        if not self.fail:
            self.get_values()
        else:
            print('get values abandoned due to previous failure')
        #telaps=time()-self.tic
        #print "elapsed time for data collection: " + str(telaps) + " s."
        
    def stop_sampling(self):
        #print "stopping the sampling process"
        hh=ps.ps5000aStop(self.chandle)
        #if hh==0:
        #    "Sampling successfully stopped"
        #else:
        #    "A problem occurred while closing board. Error code: " + str(hh)

    def open_device(self):
        self.status = {}
        self.resolution = ps.PS5000A_DEVICE_RESOLUTION["PS5000A_DR_12BIT"]
        self.chandle = c_int16()
        
        self.status["openunit"] = ps.ps5000aOpenUnit(byref(self.chandle),None,self.resolution)

        try:
            assert_pico_ok(self.status["openunit"])
        except:
            powerStatus = self.status["openunit"]
            print(powerStatus)
            if powerStatus == 286:
                self.status["changePowerSource"] = ps.ps5000aChangePowerSource(self.chandle, powerStatus)
            elif powerStatus == 282:
                self.status["changePowerSource"] = ps.ps5000aChangePowerSource(self.chandle, powerStatus)
            else:
                raise

            assert_pico_ok(self.status["changePowerSource"])

    def flash_led(self):
        print("trying to flash the led")
        hh=ps.ps5000aFlashLed(self.chandle,5)
        print(hh)
        
    def close_device(self):
        hh=ps.ps5000aCloseUnit(self.chandle)
        if hh==0:
            print("board successfully closed")
        else:
            print("A problem occurred while closing board. Error code: " + str(hh))
            
    def set_channels(self):
        self.status["setChA"] = ps.ps5000aSetChannel(self.chandle,
                                                    ps.PS5000A_CHANNEL['PS5000A_CHANNEL_A'],
                                                    1, # enabled
                                                    ps.PS5000A_COUPLING['PS5000A_DC'],
                                                    ps.PS5000A_RANGE['PS5000A_2V'],
                                                    0) # analogue offset
        print(self.status["setChA"])

        self.status["setChD"] = ps.ps5000aSetChannel(self.chandle,
                                                    ps.PS5000A_CHANNEL['PS5000A_CHANNEL_D'],
                                                    1, # enabled
                                                    ps.PS5000A_COUPLING['PS5000A_DC'],
                                                    ps.PS5000A_RANGE['PS5000A_2V'],
                                                    0) # analogue offset
        print(self.status["setChD"])
                 
    def get_timebase(self):
        self.sel_timebase=int(self.sel_dt*1e9-1)
        print(self.sel_timebase)
        self.sel_nsamp_pretrig=int(np.floor(self.sel_nsamp*self.sel_trig_pos_rel)) # number of samples before trigger RRRRRRREMOVE * NSEG!!!!!
        self.sel_nsamp_posttrig=self.sel_nsamp-self.sel_nsamp_pretrig       # number of samples after trigger

        time_interval_ns=(c_long)()
        max_samples=(c_long)()
        hh=ps.ps5000aGetTimebase(self.chandle,
                                 self.sel_timebase,
                                 self.sel_nsamp,
                                 byref(time_interval_ns),
                                 #self.sel_oversample,
                                 byref(max_samples),
                                 self.sel_seg_index)
        if hh==0:
            print("Timebase obtained. Interval: " + str(time_interval_ns) + " ns; No of samples: " + str(max_samples))
        else:
            print("Error while requesting time base. Error code: " + str(hh))
            
    def set_trigger(self):
        trig_delay=0
        
        self.maxADC = c_int16()
        self.status["maximumValue"] = ps.ps5000aMaximumValue(self.chandle, byref(self.maxADC))
        assert_pico_ok(self.status["maximumValue"])

        #thresh_adc =int(self.sel_threshold / 5000 * 32767) # threshold in counts 
        thresh_adc = int(mV2adc(self.sel_threshold,ps.PS5000A_RANGE["PS5000A_10V"], self.maxADC))
        #print(f'trig ch: {self.sel_trig_ch}')
        print(f'threshold in cts: {thresh_adc}')
        #print(f'direction: {self.sel_direction}')
        #print(f'autotrig: {self.sel_autotrig_ms}')
        
        hh=ps.ps5000aSetSimpleTrigger(self.chandle,
                                      1,
                                      ps.PS5000A_CHANNEL["PS5000A_CHANNEL_A"],
                                      thresh_adc,
                                      self.sel_direction,
                                      trig_delay,
                                      self.sel_autotrig_ms)
        if hh==0:
            print("Simple trigger set." )
        else:
            print("Error while setting simple trigger. Error code: " + str(hh))
            
    def set_memory_segments(self):
        self.n_max_samples=(c_long)() # function will return the maximum number of samples here
        hh=ps.ps5000aMemorySegments(self.chandle, self.nseg, byref(self.n_max_samples))
        if hh==0:
            print("Number of memory segments set to " + str(self.nseg) + '.' )
        else:
            print("Error while trying to set memory segments. Error code: " + str(hh))
                        
    def set_no_of_captures(self):
        hh=ps.ps5000aSetNoOfCaptures(self.chandle, self.nseg)
        if hh==0:
            print("Number of captures set to " + str(self.nseg) + '.' )
        else:
            print("Error while trying to set number of captures. Error code: " + str(hh))
            
    def run_block(self):
        self.time_busy_ms=(c_long)()
        block_ready_fct=POINTER(c_int)()
        pparm=(c_void_p)()
        
        #hh=ps.ps5000RunBlock(self.chandle, self.sel_nsamp_pretrig, self.sel_nsamp_posttrig,self.sel_timebase,self.sel_oversample,
        #                         byref(self.time_busy_ms),self.sel_seg_index,block_ready_fct,byref(pparm))
        hh=ps.ps5000aRunBlock(self.chandle,
                              self.sel_nsamp_pretrig,
                              self.sel_nsamp_posttrig,
                              self.sel_timebase,
                              #self.sel_oversample,
                              byref(self.time_busy_ms),
                              self.sel_seg_index,
                              self.data_ready,
                              byref(pparm))
        #if hh==0:
        #    print "Block mode acquisition started. Board will be busy for at least " + str(self.time_busy_ms) + " ms." 
        #else:
        #    print "Error while trying to start block acquisition. Error code: " + str(hh)

    def set_data_buffer(self):
        if self.sel_rapid_block:
            #self.bufA = [] # trying to use a list of ctype arrays because self.bufA[ii] always points to the same address
            #self.bufA = [(c_short * self.sel_nsamp)() for r in range(self.sel_nsamp)] 
            self.bufA = (c_int16* self.sel_nsamp * self.nseg)()
            for ii in range(self.nseg):
                hh1 = ps.ps5000aSetDataBufferBulk(self.chandle,0,byref(self.bufA[ii]),self.sel_nsamp,ii)
                if hh1 == 0:
                    pass
                    #print 'CHA: buffer for segment ' + str(ii) + 'ok' 
                else:
                    print('CHA: error while allocating buffer for segment ' + str(ii) + '. Error code =' + str(hh1))
                    self.fail = True
        else:
            self.bufA = (c_int16*self.sel_nsamp)()
            hh1 = ps.ps5000aSetDataBuffer(self.chandle, # handle
                                          ps.PS5000A_CHANNEL['PS5000A_CHANNEL_A'], # channel
                                          byref(self.bufA), # buffer pointer
                                          self.sel_nsamp, # buffer length
                                          0,
                                          ps.PS5000A_RATIO_MODE['PS5000A_RATIO_MODE_NONE'],)
       
        if self.sel_chon[1]:
            if self.sel_rapid_block:
                self.bufB = (c_int16 * self.sel_nsamp * self.nseg)()
                for ii in range(self.nseg):
                    hh2 = ps.ps5000aSetDataBufferBulk(self.chandle,1,byref(self.bufB[ii]),self.sel_nsamp,ii)
                    if hh1 == 0:
                        pass
                        #print 'CHB: buffer for segment ' + str(ii) + 'ok' 
                    else:
                        print('CHB: error while allocating buffer for segment ' + str(ii) + '. Error code =' + str(hh1))
                        self.fail = True

            else:
                self.bufB = (c_int16*self.sel_nsamp)()
                self.bufB_min = (c_int16*self.sel_nsamp)()
                hh2 = ps.ps5000aSetDataBuffers(self.chandle, # handle
                                               ps.PS5000A_CHANNEL['PS5000A_CHANNEL_D'], # channel
                                               byref(self.bufB), # buffer pointer
                                               byref(self.bufB_min),
                                               self.sel_nsamp, # buffer length
                                               0,
                                               0)
           
        #if hh1==0:
        #    print "Buffer of Channel A allocated successfully."
        #else:
        #    print "Error while allocating buffer of channel A. Error code: " + str(hh1)
        #if hh2==0:
        #    print "Buffer of Channel A allocated successfully."
        #else:
        #    print "Error while allocating buffer of channel A. Error code: " + str(hh1)

    def block_callback(handle, statusCallback, param):
        global wasCalledBack, ready
        wasCalledBack = True
        if statusCallback != PICO_STATUS['PICO_CANCELLED']:
            ready = True
    
    def get_values(self):
        fctr=self.range_fctr[self.sel_v]/32768.
        #print(f'fctr: {fctr}')
        nsamp=(c_ulong)(self.sel_nsamp)
        #print nsamp
        
        self.maxADC = c_int16()
        self.status["maximumValue"] = ps.ps5000aMaximumValue(self.chandle, byref(self.maxADC))
        assert_pico_ok(self.status["maximumValue"])
        
        
        # Convert the python function into a C function pointer.
        cFuncPtr = ps.BlockReadyType(self.block_callback)

        if self.sel_rapid_block:
            ofl=(c_short * self.nseg)()
            hh=ps.ps5000aGetValuesBulk(self.chandle,byref(nsamp),0,self.nseg-1,byref(ofl))
        else:
            ofl=(c_short)()
            hh=ps.ps5000aGetValues(self.chandle,0,byref(nsamp),0,0,0,byref(ofl))
            
        #if self.sel_chon[1]:
        #    hh=ps.ps5000aGetValues(self.chandle,1,byref(nsamp),0,0,0,byref(ofl)) # NOT NEEDED???
        
        if hh==0:
            #print "Buffer filled with " + str(nsamp) + " samples." 
            #self.ADCA=np.array(self.bufA,dtype=float)*fctr
            self.ADCA = np.array(adc2mV(self.bufA, ps.PS5000A_RANGE["PS5000A_10V"] , self.maxADC))/1000
            if self.sel_chon[1]:
                self.ADCB = np.array(adc2mV(self.bufB, ps.PS5000A_RANGE["PS5000A_10V"], self.maxADC))/1000
                
            self.flag_newdata = True
            #print self.ADCA.shape
        else:
            print("Error while trying to transfer to buffer. Error code: " + str(hh))
            self.fail = True
        
    def set_sel_threshold(self, val) -> None:
        self.sel_threshold = val

class new_picoscope:
    def __init__(self):
        # channels
        self.channelA = ps.PS5000A_CHANNEL["PS5000A_CHANNEL_A"]
        self.channelB = ps.PS5000A_CHANNEL["PS5000A_CHANNEL_B"]
        self.channelC = ps.PS5000A_CHANNEL["PS5000A_CHANNEL_C"]
        self.channelD = ps.PS5000A_CHANNEL["PS5000A_CHANNEL_D"]
        
        self.trigger_channel = self.channelA # defaults to channel A
        self.data_channel    = self.channelD # defaults to channel D
        
        self.coupling_type = ps.PS5000A_COUPLING["PS5000A_DC"] # not sure what this does
        
        # device resolution modes, lower 'bit-ness' allows for faster time intervals
        self.resolution16 = ps.PS5000A_DEVICE_RESOLUTION["PS5000A_DR_16BIT"]
        self.resolution15 = ps.PS5000A_DEVICE_RESOLUTION["PS5000A_DR_15BIT"]
        self.resolution14 = ps.PS5000A_DEVICE_RESOLUTION["PS5000A_DR_14BIT"]
        self.resolution12 = ps.PS5000A_DEVICE_RESOLUTION["PS5000A_DR_12BIT"]
        self.resolution8  = ps.PS5000A_DEVICE_RESOLUTION["PS5000A_DR_8BIT" ]    
        self.resolution   = 12
        
        self.timeInterval_ns = (c_int32)() # nanoseconds
        self.sel_nsamp       = 100000
        self.sel_nseg        = 2
        self.nseg            = 100
        self.sel_threshold   = 2000 # mV
        self.chandle = c_int16()
        self.status  = {}

        self.sel_timebase = 3
        self.sel_dt        = 1e-8
        self.sel_direction = 2 # rising
        self.trigger_delay = 0
        self.range_fctr = np.array(([1e-2, 2e-2, 5e-2, 0.1, 0.2, 0.5, 1, 2, 5, 10, 20]),dtype=float)
        self.sel_v = 9  # defaults to `10V`
        self.offset = 0 # analog offset
        self.sel_trig_pos_rel = 0.05
        self.sel_rapid_block = False
        self.sel_chon = [1,1]
        self.flag_newdata = False
        self.sel_seg_index = 0
        self.DATA_READY = WINFUNCTYPE(c_void_p, c_short, c_short, POINTER(c_void_p))
        self.data_ready = self.DATA_READY(self.py_data_ready)
        self.fail = False
        self.sel_autotrig_ms = 200
        
    def py_data_ready(self,a,b,c):
        try:    self.set_data_buffer()
        except: print("Error while trying to transfer to buffer" )
        
        try:
            #print(self.fail)
            if not self.fail: self.get_values()
            else:             print('get values abandoned due to previous failure')
        except:               print("Get values failed" )

    def open_device(self):
        """ Opens the picoscope device. """
        
        match self.resolution: # chooses the right device resolution mode
            case  8: resolution = self.resolution8
            case 12: resolution = self.resolution12
            case 14: resolution = self.resolution14
            case 15: resolution = self.resolution15
            case 16: resolution = self.resolution16
        
        self.status["openunit"] = ps.ps5000aOpenUnit(byref(self.chandle), None, resolution)
        assert_pico_ok(self.status["openunit"])

        # this section checks that the picoscope is connected to the power source. Not sure if it is needed but it was in the examples from the manufacturer.
        try: assert_pico_ok(self.status["openunit"])
        except: # PicoNotOkError:
            powerStatus = self.status["openunit"]
            if   powerStatus == 286: self.status["changePowerSource"] = ps.ps5000aChangePowerSource(self.chandle, powerStatus)
            elif powerStatus == 282: self.status["changePowerSource"] = ps.ps5000aChangePowerSource(self.chandle, powerStatus)
            else: raise
            assert_pico_ok(self.status["changePowerSource"])

    def close_device(self):
        """  Closes the picoscope. """
        self.status["closed"] = ps.ps5000aCloseUnit(self.chandle)
        assert_pico_ok(self.status["closed"])

    def enable_channel(self, channel, voltage_range, offset, name):
        print(voltage_range, offset)
        self.status[f'Channel {name}'] = ps.ps5000aSetChannel(self.chandle, channel, 1, self.coupling_type, voltage_range, offset)
        
    def set_channels(self):
        """
        Sets up the channels to be used. The others are disabled.
        The trigger channel is hard coded to have a range of 10V, but the data channel is chosen by the user.
        If the data channel is chosen to have a range of 1V, then the signal is offset by -0.5V.
        This allows an actual voltage range of -0.5 to 1.5V.
        """
        list_of_channels = [self.channelA, self.channelB, self.channelC, self.channelD]

        for channel in list_of_channels:
            match channel:
                case self.trigger_channel:
                    # enable the trigger channel. It is forced to have a range of 10V
                    self.enable_channel(channel, 9, 0, "Trigger")
                case self.data_channel:
                    # enable the data channel
                    # if the chosen voltage range is 1V, then the signal is offset by -0.5V, otherwise the offset is 0.
                    if self.sel_v == 6: self.enable_channel(channel, 6,         -0.5, "Data")              
                    else:               self.enable_channel(channel, self.sel_v, 0,   "Data")
                case _:
                    # disable the channel
                    self.status[f'Disabled {channel}'] = ps.ps5000aSetChannel(self.chandle, channel, 0, self.coupling_type, self.sel_v, 0)

    def set_data_buffer(self):
        if self.sel_rapid_block:
            self.bufA = (c_int16* self.sel_nsamp * self.nseg)()
            for ii in range(self.nseg):
                self.status["BufferA"] = ps.ps5000aSetDataBufferBulk(self.chandle, 0, byref(self.bufA[ii]), self.sel_nsamp, ii)
                if self.status["BufferA"] == 0: continue
                else:
                    assert_pico_ok(self.status["BufferA"])
                    self.fail = True
        else:
            self.bufA = (c_int16* self.sel_nsamp)()
            self.status["BufferA"] = ps.ps5000aSetDataBuffer(self.chandle, self.trigger_channel, byref(self.bufA), self.sel_nsamp, 0, ps.PS5000A_RATIO_MODE['PS5000A_RATIO_MODE_NONE'])
        assert_pico_ok(self.status["BufferA"])
        if self.sel_chon[1]:
            if self.sel_rapid_block:
                self.bufB = (c_int16 * self.sel_nsamp * self.nseg)()
                for ii in range(self.nseg):
                    self.status["BufferD"] = ps.ps5000aSetDataBufferBulk(self.chandle, 1, byref(self.bufB[ii]), self.sel_nsamp, ii)
                    if self.status["BufferD"] == 0: continue
                    else:
                        assert_pico_ok(self.status["BufferD"])
                        self.fail = True
            else:
                self.bufB = (c_int16*self.sel_nsamp)()
                self.bufB_min = (c_int16*self.sel_nsamp)()
                self.status["BufferD"] = ps.ps5000aSetDataBuffers(self.chandle, self.data_channel, byref(self.bufB), byref(self.bufB_min), self.sel_nsamp, 0, 0)
            assert_pico_ok(self.status["BufferD"])
    
    def set_memory_segments(self):
        self.n_max_samples=(c_long)() # function will return the maximum number of samples here
        self.status["memorySegments"]=ps.ps5000aMemorySegments(self.chandle, self.nseg, byref(self.n_max_samples))
        assert_pico_ok(self.status["memorySegments"])
                        
    def set_no_of_captures(self):
        self.status["numberCaptures"] = ps.ps5000aSetNoOfCaptures(self.chandle, self.nseg)
        assert_pico_ok(self.status["numberCaptures"])
    
    def get_timebase(self):
        """ Gets the timebase. """
        print("Getting timebase...")
        
        self.returnedMaxSamples = (c_int32)()
        
        self.preTriggerSamples  = int(np.floor(self.sel_nsamp * self.sel_trig_pos_rel))
        self.postTriggerSamples = self.sel_nsamp - self.preTriggerSamples
        print(f"{self.sel_timebase = }")
        self.status["getTimebase"] = ps.ps5000aGetTimebase2(self.chandle, self.sel_timebase, self.sel_nsamp, byref(self.timeInterval_ns), byref(self.returnedMaxSamples), 0)
        assert_pico_ok(self.status["getTimebase"])
        print(f'Timebase obtained. Interval: {self.timeInterval_ns.value} ns; No of samples: {self.returnedMaxSamples.value}')
    
    def set_trigger(self):
        self.maxADC = c_int16()
        self.status["maximumValue"] = ps.ps5000aMaximumValue(self.chandle, byref(self.maxADC))
        assert_pico_ok(self.status["maximumValue"])

        thresh_adc = int(mV2adc(self.sel_threshold,9, self.maxADC)) # the 9 is for a 10V range
        print(f'trigger delay {int(self.trigger_delay)}')
        self.status["simpleTrigger"] = ps.ps5000aSetSimpleTrigger(self.chandle, 1, self.trigger_channel, thresh_adc, self.sel_direction, int(self.trigger_delay), self.sel_autotrig_ms)
        assert_pico_ok(self.status["simpleTrigger"])

    def block_callback(self, handle, statusCallback): # idk if this is needed
        wasCalledBack = True
        if statusCallback != PICO_STATUS['PICO_CANCELLED']:
            self.ready = True

    def run_block(self):
        pparm=(c_void_p)()
        self.time_busy_ms = (c_int32)()
        # Convert the python function into a C function pointer.
        cFuncPtr = ps.BlockReadyType(self.block_callback) # might not be needed

        self.status["runBlock"] = ps.ps5000aRunBlock(self.chandle, self.preTriggerSamples, self.postTriggerSamples, self.sel_timebase,
                                                     byref(self.time_busy_ms), self.sel_seg_index, self.data_ready, byref(pparm))
        
        assert_pico_ok(self.status["runBlock"])
        assert_pico_ok(self.status["getTimebase"])
    
    def get_values(self):
        self.nsamp = (c_ulong)(self.sel_nsamp)
        try:
            self.maxADC = c_int16()
            self.status["maximumValue"] = ps.ps5000aMaximumValue(self.chandle, byref(self.maxADC))
            assert_pico_ok(self.status["maximumValue"])
        except:
            print("failure trying to get the maximum value")

        if self.sel_rapid_block:
            ofl=(c_short * self.nseg)()
            self.status["getValuesBulk"] = ps.ps5000aGetValuesBulk(self.chandle, byref(self.nsamp), 0, self.nseg-1, byref(ofl))
            print(self.status)
            assert_pico_ok(self.status["getValuesBulk"])
        else:
            ofl=(c_short)()
            self.status["getValues"] = ps.ps5000aGetValues(self.chandle, 0, byref(self.nsamp), 0, 0, 0, byref(ofl))
            assert_pico_ok(self.status["getValues"])
    
        if self.status["getValues"]==0 or self.status["getValuesBulk"]==0:
            #print(self.sel_v)
            self.ADCA = np.array(adc2mV(self.bufA,          9, self.maxADC))/1000                # convert to V from mV
            self.ADCB = np.array(adc2mV(self.bufB, self.sel_v, self.maxADC))/1000 - self.offset # in V

            self.flag_newdata = True
        else:
            assert_pico_ok(self.status["getValues"])
            self.fail = True
    
    def set_sel_threshold(self, val) -> None:
        self.sel_threshold = val