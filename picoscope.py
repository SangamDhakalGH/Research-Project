import numpy as np
#import serial
from ctypes import *
import time
import os

class picoscope:
    def __init__(self):
        # enum definitions:
        #channels
        self.PS4000_CHANNEL_A=0
        self.PS4000_CHANNEL_B=1
        self.PS4000_EXTERNAL=4
        self.PS4000_MAX_CHANNELS=self.PS4000_EXTERNAL
        self.PS4000_TRIGGER_AUX=3
        self.PS4000_MAX_TRIGGER_SOURCES=2
        
        #voltage ranges
        self.PS4000_10MV=0
        self.PS4000_20MV=1
        self.PS4000_50MV=2
        self.PS4000_100MV=3
        self.PS4000_200MV=4
        self.PS4000_500MV=5
        self.PS4000_1V=6
        self.PS4000_2V=7
        self.PS4000_5V=8
        self.PS4000_10V=9
        self.PS4000_20V=10
        self.range_fctr=np.array(([1e-2, 2e-2, 5e-2, 0.1, 0.2, 0.5, 1, 2, 5, 10, 20]),dtype=float)
        
        # trigger direction
        self.PS4000_TRIG_ABOVE=0
        self.PS4000_TRIG_BELOW=1
        self.PS4000_TRIG_RISING=2
        self.PS4000_TRIG_FALLING=3
        
        # ac/dc coupling
        self.PS4000_COUPLING_DC=1
        self.PS4000_COUPLING_AC=0
        
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
        self.sel_threshold = 0                    # threshold in adc counts
        self.sel_autotrig_ms = 200
        self.sel_direction = self.PS4000_TRIG_RISING
        self.sel_trig_pos_rel = 0.1    # relative position of trigger in collected data
        self.sel_trig_ch = self.PS4000_CHANNEL_A              
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
        hh=self.L.ps4000Stop(self.unit_hd)
        #if hh==0:
        #    "Sampling successfully stopped"
        #else:
        #    "A problem occurred while closing board. Error code: " + str(hh)

    def open_device(self):
        #os.chdir(r"C:\Program Files\Pico Technology\SDK\lib")
        #self.L = cdll.LoadLibrary("ps4000.dll")
        self.L = cdll.LoadLibrary("C:\Program Files\Pico Technology\SDK\lib\ps4000.dll")
        #self.L=windll.PS4000
        print(self.L)
        self.unit_hd=(c_ushort)()# the handle to talk to the device
        hh=self.L.ps4000OpenUnit(byref(self.unit_hd))
        if hh==0:
            print("Picoscope successfully opened. Handle:")
            print(self.unit_hd)
        else:
            print("board could not be opened, Error code: " + str(hh))
        
    def flash_led(self):
        print("trying to flash the led")
        hh=self.L.ps4000FlashLed(self.unit_hd,5)
        print(hh)
        
    def close_device(self):
        hh=self.L.ps4000CloseUnit(self.unit_hd)
        if hh==0:
            print("board successfully closed")
        else:
            print("A problem occurred while closing board. Error code: " + str(hh))
            
    def set_channels(self):
        for jj in range(2):
            hh=self.L.ps4000SetChannel(self.unit_hd,jj,self.sel_chon[jj],self.PS4000_COUPLING_DC,self.sel_v)
            if hh==0:
                print("Channel " + str(jj+1) + " set successfully.")
            else:
                print("Error while setting channel " + str(jj+1) + ". Error code: " + str(hh))
                
    def get_timebase(self):
        self.sel_timebase=int(self.sel_dt*1e7-1) # according to timebase rule on pg17 in programmers guide
        self.sel_nsamp_pretrig=int(np.floor(self.sel_nsamp*self.sel_trig_pos_rel)) # number of samples before trigger RRRRRRREMOVE * NSEG!!!!!
        self.sel_nsamp_posttrig=self.sel_nsamp-self.sel_nsamp_pretrig       # number of samples after trigger

        time_interval_ns=(c_long)()
        max_samples=(c_long)()
        hh=self.L.ps4000GetTimebase(self.unit_hd,self.sel_timebase,self.sel_nsamp,byref(time_interval_ns),
                                    self.sel_oversample,byref(max_samples),self.sel_seg_index)
        if hh==0:
            print("Timebase obtained. Interval: " + str(time_interval_ns) + " ns; No of samples: " + str(max_samples))
        else:
            print("Error while requesting time base. Error code: " + str(hh))
            
    def set_trigger(self):
        trig_delay=0
        
        thresh_adc =int(self.sel_threshold / 5000. * 32767) # threshold in counts 
        
        print('trig ch:')
        print(self.sel_trig_ch)
        print('threshold in cts')
        print(thresh_adc)
        print('direction')
        print(self.sel_direction)
        print('autotrig')
        print(self.sel_autotrig_ms)
        
        hh=self.L.ps4000SetSimpleTrigger(self.unit_hd,1,self.sel_trig_ch,
                                         thresh_adc,self.sel_direction,trig_delay,self.sel_autotrig_ms)
        if hh==0:
            print("Simple trigger set." )
        else:
            print("Error while setting simple trigger. Error code: " + str(hh))
            
    def set_memory_segments(self):
        self.n_max_samples=(c_long)() # function will return the maximum number of samples here
        hh=self.L.ps4000MemorySegments(self.unit_hd, self.nseg, byref(self.n_max_samples))
        if hh==0:
            print("Number of memory segments set to " + str(self.nseg) + '.' )
        else:
            print("Error while trying to set memory segments. Error code: " + str(hh))
                        
    def set_no_of_captures(self):
        hh=self.L.ps4000SetNoOfCaptures(self.unit_hd, self.nseg)
        if hh==0:
            print("Number of captures set to " + str(self.nseg) + '.' )
        else:
            print("Error while trying to set number of captures. Error code: " + str(hh))
            
    def run_block(self):
        self.time_busy_ms=(c_long)()
        block_ready_fct=POINTER(c_int)()
        pparm=(c_void_p)()
        
        #hh=self.L.ps4000RunBlock(self.unit_hd, self.sel_nsamp_pretrig, self.sel_nsamp_posttrig,self.sel_timebase,self.sel_oversample,
        #                         byref(self.time_busy_ms),self.sel_seg_index,block_ready_fct,byref(pparm))
        hh=self.L.ps4000RunBlock(self.unit_hd, self.sel_nsamp_pretrig, self.sel_nsamp_posttrig,self.sel_timebase,self.sel_oversample,
                                 byref(self.time_busy_ms),self.sel_seg_index,self.data_ready,byref(pparm))
        #if hh==0:
        #    print "Block mode acquisition started. Board will be busy for at least " + str(self.time_busy_ms) + " ms." 
        #else:
        #    print "Error while trying to start block acquisition. Error code: " + str(hh)

    def set_data_buffer(self):
        if self.sel_rapid_block:
            #self.bufA = [] # trying to use a list of ctype arrays because self.bufA[ii] always points to the same address
            #self.bufA = [(c_short * self.sel_nsamp)() for r in range(self.sel_nsamp)] 
            self.bufA = (c_short * self.sel_nsamp * self.nseg)()
            for ii in range(self.nseg):
                hh1 = self.L.ps4000SetDataBufferBulk(self.unit_hd,0,byref(self.bufA[ii]),self.sel_nsamp,ii)
                if hh1 == 0:
                    pass
                    #print 'CHA: buffer for segment ' + str(ii) + 'ok' 
                else:
                    print('CHA: error while allocating buffer for segment ' + str(ii) + '. Error code =' + str(hh1))
                    self.fail = True
        else:
            self.bufA = (c_short*self.sel_nsamp)()
            hh1 = self.L.ps4000SetDataBuffer(self.unit_hd,0,byref(self.bufA),self.sel_nsamp)
       
        if self.sel_chon[1]:
            if self.sel_rapid_block:
                self.bufB = (c_short * self.sel_nsamp * self.nseg)()
                for ii in range(self.nseg):
                    hh2 = self.L.ps4000SetDataBufferBulk(self.unit_hd,1,byref(self.bufB[ii]),self.sel_nsamp,ii)
                    if hh1 == 0:
                        pass
                        #print 'CHB: buffer for segment ' + str(ii) + 'ok' 
                    else:
                        print('CHB: error while allocating buffer for segment ' + str(ii) + '. Error code =' + str(hh1))
                        self.fail = True

            else:
                self.bufB = (c_short*self.sel_nsamp)()
                hh2 = self.L.ps4000SetDataBuffer(self.unit_hd,1,byref(self.bufB),self.sel_nsamp)
           
        #if hh1==0:
        #    print "Buffer of Channel A allocated successfully."
        #else:
        #    print "Error while allocating buffer of channel A. Error code: " + str(hh1)
        #if hh2==0:
        #    print "Buffer of Channel A allocated successfully."
        #else:
        #    print "Error while allocating buffer of channel A. Error code: " + str(hh1)

    def get_values(self):
        fctr=self.range_fctr[self.sel_v]/32768.
        #print fctr
        nsamp=(c_ulong)(self.sel_nsamp)
        #print nsamp
        
        if self.sel_rapid_block:
            ofl=(c_short * self.nseg)()
            hh=self.L.ps4000GetValuesBulk(self.unit_hd,byref(nsamp),0,self.nseg-1,byref(ofl))
        else:
            ofl=(c_short)()
            hh=self.L.ps4000GetValues(self.unit_hd,0,byref(nsamp),0,0,0,byref(ofl))
            
        #if self.sel_chon[1]:
        #    hh=self.L.ps4000GetValues(self.unit_hd,1,byref(nsamp),0,0,0,byref(ofl)) # NOT NEEDED???
        
        if hh==0:
            #print "Buffer filled with " + str(nsamp) + " samples." 
            self.ADCA=np.array(self.bufA,dtype=float)*fctr
            if self.sel_chon[1]:
                self.ADCB=np.array(self.bufB,dtype=float)*fctr
            self.flag_newdata = True
            #print self.ADCA.shape
        else:
            print("Error while trying to transfer to buffer. Error code: " + str(hh))
            self.fail = True
        
    def set_sel_threshold(self, val) -> None:
        self.sel_threshold = val