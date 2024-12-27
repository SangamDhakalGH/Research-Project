import sys
import os
import time
import locale
import PyQt6.QtGui     as qtg
import PyQt6.QtWidgets as qtw
import matplotlib      as mp
mp.use('QtAgg')
import numpy  as np
import pandas as pd
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from picoscope_5000   import new_picoscope as picoscope
from MatplotlibCanvas import MatplotlibCanvas
from OutputBox        import OutputBox
from MotorControl     import Motor
from TLDC2200_pyvisa  import TLDC2200
from Berkeley645      import SignalGenerator
from femtofunc        import save_matrix_milano_convention #, binnen

mp.rcParams['font.size'] = 8
locale.setlocale(locale.LC_NUMERIC, 'C')

class muspp_ui(qtw.QMainWindow):    
    def __init__(self) -> None:
        self.setupUi()

    def setupUi(self) -> None:
        ''' Setting up the interface. '''
        qtw.QMainWindow.__init__(self)

        layout      = qtw.QGridLayout()
        layout.setSpacing(1)
        layout_psh  = qtw.QVBoxLayout()
        layout_cmb  = qtw.QVBoxLayout()
        self.widget = qtw.QWidget()
        font_norm   = qtg.QFont()
        font_norm.setPointSize(9)

        font_bold = qtg.QFont()
        font_bold.setBold(True)
        font_bold.setPointSize(9)
        self.widget.setFont(font_norm)

        ##### Define GUI Elements on the Left Side #####
        # Minimum for LED current
        self.min_LED_current_label = qtw.QLabel(   text='Min Bias Light Current (mA)', maximumWidth=200)
        self.min_LED_current       = qtw.QLineEdit(text='10',                          maximumWidth=100)
        # Maximum for LED current
        self.max_LED_current_label = qtw.QLabel(   text='Max Bias Light Current (mA)', maximumWidth=200)
        self.max_LED_current       = qtw.QLineEdit(text='10000',                       maximumWidth=100)
        # Number of points to set LED at
        self.steps_LED_current_label = qtw.QLabel(   text='Number of Bias Values', maximumWidth=150)
        self.steps_LED_current       = qtw.QLineEdit(text='15',                    maximumWidth=100)
        # Drop-down Menu for Lin or Log evenly spaced light intensities
        self.label_LED_spacing = qtw.QLabel(   text='Bias Intensity Spacing', maximumWidth=150)
        self.combo_LED_spacing = qtw.QComboBox(maximumWidth=100)
        self.combo_LED_spacing.addItems(["Linear","Logarithmic"])
        # Voltage Combo Box
        self.label_V = qtw.QLabel(   text='Y Scale', maximumWidth=100)
        self.combo_V = qtw.QComboBox(maximumWidth=125)
        # options for voltage combo box
        voltage_options = ['10 mV',  '20 mV', '50 mV', '100 mV', '200 mV',
                           '500 mV', '1 V (offset of -0.5V)', '2 V', '5 V', '10 V', '20 V']
        self.combo_V.addItems(voltage_options)
        self.combo_V.setCurrentIndex(9)  # sets default to '10 V'
        # Time Interval Input Box
        self.label_dt = qtw.QLabel(text='Sampling Interval', maximumWidth=100)
        self.combo_dt = qtw.QComboBox(maximumWidth=125, toolTip='This must be chosen before opening the picoscope.')
        self.time_interval_options = [ '2 ns | 8 bit' ,  '4 ns | 12 bit' , '8 ns | 12 bit',   '8 ns | 15 bit',  '16 ns | 15 bit',  '32 ns | 16 bit',  '48 ns | 16 bit',
                                      '64 ns | 16 bit', '80 ns | 16 bit', '96 ns | 16 bit', '112 ns | 16 bit', '128 ns | 16 bit', '144 ns | 16 bit']
        self.combo_dt.addItems(self.time_interval_options)
        self.combo_dt.setCurrentIndex(3) # sets default to '4 ns'
        # Frequency Input Box
        self.label_freq = qtw.QLabel(text='Frequency (Hz)', maximumWidth=100)
        self.line_freq  = qtw.QLineEdit(text='1000',        maximumWidth=100)
        # Number of Segments Input Box
        self.label_nseg = qtw.QLabel(text='Number of Segments', maximumWidth=100)
        self.line_nseg  = qtw.QLineEdit(text='300',             maximumWidth=100)
        # Delta(V)/V ratio box
        self.label_desired_ratio = qtw.QLabel(text='\u0394V/V ratio', maximumWidth=200)
        self.line_desired_ratio  = qtw.QLineEdit(text='0.05',         maximumWidth=100)
        # Noise or Time Combo Box
        self.label_Time_Noise = qtw.QLabel(text='Ends when?', maximumWidth=100)
        self.combo_Time_Noise = qtw.QComboBox(maximumWidth=130)
        self.combo_Time_Noise.addItems(["Time","Noise","Whichever Comes First"])
        # Time Input Box
        self.label_time = qtw.QLabel(text='Time per Measurement (s)', maximumWidth=200)
        self.line_time  = qtw.QLineEdit(text='30',                    maximumWidth=100)
        # Target Noise Input Box
        self.label_target_noise = qtw.QLabel(text='Target Noise', maximumWidth=100)
        self.line_target_noise  = qtw.QLineEdit(text='2e-3',      maximumWidth=100)
        # Data Channel Combo Box
        self.label_data_channel = qtw.QLabel(text='Data Channel', maximumWidth=100)
        self.combo_data_channel = qtw.QComboBox(maximumWidth=100, toolTip='Channel where the voltage you want to measure would come in. (Default is Channel D)')
        channels = ['Channel A', 'Channel B', 'Channel C', 'Channel D', 'EXT']
        self.combo_data_channel.addItems(channels)
        self.combo_data_channel.setCurrentIndex(3)  # sets default to `Channel D`
        # Trigger Channel Combo Box, the trigger direction is hard coded to be `rising`
        self.label_trig_channel    = qtw.QLabel(text='Trigger Channel', maximumWidth=100)
        self.combo_trigger_channel = qtw.QComboBox(maximumWidth=100, toolTip='Channel where the laser trigger would come in. (Default is Channel A)')
        self.combo_trigger_channel.addItems(channels)
        self.combo_trigger_channel.setCurrentIndex(0)  # sets default to `Channel A`
        # Trigger Threshold Input Box
        self.label_trigger_threshold = qtw.QLabel(text='Trig Thresh.(mV)', maximumWidth=100)
        self.line_trigger_threshold  = qtw.QLineEdit(text='3000',          maximumWidth=100)

        ##### Place GUI Elements on the Left Side #####
        # adds the list of widgets utilizing a loop (was originally adding elements explicitely which caused lots of repitition)
        widget_list = [self.min_LED_current_label,   self.min_LED_current,         self.max_LED_current_label,   self.max_LED_current,
                       self.steps_LED_current_label, self.steps_LED_current,       self.label_LED_spacing,       self.combo_LED_spacing,       self.label_V,
                       self.combo_V,  self.label_dt, self.combo_dt,                self.label_freq,               self.line_freq,              self.label_nseg,
                       self.line_nseg,               self.label_desired_ratio,     self.line_desired_ratio,      self.label_Time_Noise,        self.combo_Time_Noise,
                       self.label_time,              self.line_time,               self.label_target_noise,      self.line_target_noise,       self.label_data_channel,
                       self.combo_data_channel,      self.label_trig_channel,      self.combo_trigger_channel,   self.label_trigger_threshold, self.line_trigger_threshold]
        for i in widget_list: layout_cmb.addWidget(i)

        ##### Define and Place Graphical Elements #####
        layout_graph1 = qtw.QVBoxLayout()
        layout_graph2 = qtw.QVBoxLayout()
        
        self.left_graph   = MatplotlibCanvas(self,height=4)
        self.left_toolbar = NavigationToolbar(self.left_graph, self)
        layout_graph1.addWidget(self.left_toolbar)
        layout_graph1.addWidget(self.left_graph)
        
        self.right_graph   = MatplotlibCanvas(self)
        self.right_toolbar = NavigationToolbar(self.right_graph, self)
        self.right_graph.axes.set_ylim(-0.2,1.5)
        self.right_graph.axes.set_ylabel("Voltage (V)")
        layout_graph2.addWidget(self.right_toolbar)
        layout_graph2.addWidget(self.right_graph)

        ##### Define a Output Box #####
        self.output_box = OutputBox()
        layout_graph1.addWidget(self.output_box)

        ##### Define GUI Elements on the Right Side #####
        self.psh_open_devices  = qtw.QPushButton(text='Open Devices',           maximumWidth=150, toolTip='Opens the connects to the Picoscope, LED, and motor.')
        self.psh_voc_scan      = qtw.QPushButton(text='Scan Voc Range',         maximumWidth=150, toolTip='Check the Voc across the bias light range in order to determine the range that is worth measuring.')
        self.psh_find_trig_lev = qtw.QPushButton(text='Find Trigger Level',     maximumWidth=150)
        self.psh_close_devices = qtw.QPushButton(text='Close Devices',          maximumWidth=150, toolTip='Closes the connects to the Picoscope, LED, and motor.')
        self.psh_pp            = qtw.QPushButton(text='Pump-Probe 1',           maximumWidth=150, toolTip='TPV/TPC w/o the ability to control the delta(V)/V ratio.')
        self.psh_pp_filter     = qtw.QPushButton(text='Pump-Probe 1 w/ Filter', maximumWidth=150, toolTip='TPV w/ the ability to control the delta(V)/V ratio.')
        self.psh_pp_filter_TPC = qtw.QPushButton(text='TPC w/ Position from file', maximumWidth=150)
        self.psh_pp_filter_TPC_check = qtw.QCheckBox(text='Retain Filter Pos for TPC?',           toolTip='Saves the positions of the filter in millimeters to have identical conditions for TPC measurments.')
        self.psh_pp2        = qtw.QPushButton(text='Pump-Probe 2', maximumWidth=150)
        self.psh_pp3        = qtw.QPushButton(text='Pump-Probe 3', maximumWidth=150)
        self.check_pumponly = qtw.QCheckBox(  text='Pump only?')
        self.check_TATPC    = qtw.QCheckBox(  text='TA/TPC?',     toolTip='Connect PD to CHA for TA, connect Solar Cell to CHB over a 50 Hz resistor for TPC')
        self.check_chref    = qtw.QCheckBox(  text='Reference',   toolTip='Connect signal PD to CHA for TA, connect reference PD to CHB.')
        self.psh_pp_fr      = qtw.QPushButton(text='PP Free Run', maximumWidth=150)
        self.psh_stop       = qtw.QPushButton(text='Continue',    maximumWidth=150, toolTip='Stop sampling and continue to next step')
        self.psh_abort      = qtw.QPushButton(text='Abort',       maximumWidth=150, toolTip='Abort current measurement')
        
        #### Place GUI Elements on Right Side #####
        ### again, this was converted to a loop by Steven ###
        widget_list2 = [self.psh_open_devices, self.psh_voc_scan,  self.psh_find_trig_lev, self.psh_close_devices,
                        self.psh_pp,           self.psh_pp_filter, self.psh_pp_filter_TPC_check,self.psh_pp_filter_TPC,
                        self.psh_pp2,          self.psh_pp3,       self.check_pumponly,    self.check_TATPC,
                        self.check_chref,      self.psh_pp_fr,     self.psh_stop,          self.psh_abort]
        for i in widget_list2: layout_psh.addWidget(i)

        # place layouts
        layout.addLayout(layout_cmb,    0, 0, 1, 1)
        layout.addLayout(layout_graph1, 0, 1, 1, 1)
        layout.addLayout(layout_graph2, 0, 2, 1, 1)
        layout.addLayout(layout_psh,    0, 3, 1, 1)
        
        self.widget.setLayout(layout)
        self.setCentralWidget(self.widget)
        self.setWindowTitle("TPV/TPC Measurements")

        self.pfad_open     = os.getcwd()
        self.stop_request  = False
        self.abort_request = False
        self.mode = 'freerun'
        # if true, pump-probe will have two parts and require user to block the probe for the second part
        self.measure_pumponly = True
        self.ncore = 20  # nr of iteration for core loop before plotting is done
        # noise which is considered acceptable so sampling loop will go to next step
        self.target_noise = 2e-6
        self.fig2_yrange  = 1e-3
        self.fig1_yrange  = 1e-3
        self.fig2_ycenter = 0.0
        self.fig1_ycenter = 0.0
        self.plot_nbin    = 500
        
        # connect the appropriate functions to each button
        self.psh_open_devices .clicked.connect(self.open_devices)
        self.psh_voc_scan     .clicked.connect(self.Voc_Scan)
        self.psh_close_devices.clicked.connect(self.close_devices)
        self.psh_find_trig_lev.clicked.connect(self.find_trigger_level)
        self.psh_pp           .clicked.connect(lambda mode='measure': self.pump_probe_one(mode))
        self.psh_pp_filter    .clicked.connect(lambda mode='measure': self.pump_probe_one_filter(mode))
        self.psh_pp_filter_TPC.clicked.connect(lambda mode='measure': self.pump_probe_one_filter_TPC(mode))
        self.psh_pp2          .clicked.connect(lambda mode='measure': self.pump_probe_two(mode))
        self.psh_pp3          .clicked.connect(lambda mode='measure': self.pump_probe_three(mode))
        self.psh_pp_fr        .clicked.connect(lambda mode='freerun': self.pump_probe_freerun(mode)) # add lambda fct to indicate that this is free-run
        self.psh_stop         .clicked.connect(self.stop)
        self.psh_abort        .clicked.connect(self.abort)
        self.line_trigger_threshold.returnPressed.connect(self.refresh_trig)

    def find_trigger_level(self) -> None:
        ''' 
        Finds the voltage level at which to place the trigger for the picoscope.
        Can the trigger level just be fixed to something like 3 volts? I will need to test this.
        '''

        self.output_box.write_to_output_box('Setting Trigger Level...')

        self.setup_experiment()
        self.ps.trigger_delay = 0
        self.ps.set_trigger()
        self.signal_generator.output_on()

        A, _ = self.pump_probe_core(20)

        self.left_graph.draw_plot(self.t,A) # this is for testing purposes

        # set the value for the new trigger
        new_trigger = int(1.75*1000) #int((max(A)-1)*1000) # A is in volts, so multiply by 1000 to get mV

        self.ps.sel_threshold = new_trigger
        self.ps.set_trigger() # set the trigger to the proper amount
        self.signal_generator.output_off()
        self.output_box.write_to_output_box(f"setting trigger level to {new_trigger} mV.")
        self.line_trigger_threshold.setText(str(new_trigger))
        del self.ps.ADCA # clears this so it can be used in the future without issue

    def refresh_trig(self) -> None:
        self.output_box.write_to_output_box('calling refresh trig')
        self.ps.sel_threshold = int(str(self.line_trigger_threshold.text()))
        self.ps.set_trigger()

    def read_user_settings(self) -> None:
        self.output_box.write_to_output_box('Reading User Settings...')
        
        # Voltage Range
        self.ps.sel_v = self.combo_V.currentIndex()
        # Sampling Interval (the amount of time between each point, in seconds)
        self.ps.sel_dt     = int(str(self.time_interval_options[self.combo_dt.currentIndex()]).split(' ')[0])*1e-9
        self.ps.resolution = int(str(self.time_interval_options[self.combo_dt.currentIndex()]).split(' ')[3])
        match self.ps.sel_dt*1e9:
            # case {time interval in ns}: self.sel_timebase = {corresponding timebase from page 22 of Pico500a manual}
            case 2:   self.ps.sel_timebase = 1 
            case 4:   self.ps.sel_timebase = 2
            case 8:   self.ps.sel_timebase = 3
            case 16:  self.ps.sel_timebase = 4
            case 32:  self.ps.sel_timebase = 5
            case 48:  self.ps.sel_timebase = 6
            case 64:  self.ps.sel_timebase = 7
            case 80:  self.ps.sel_timebase = 8
            case 96:  self.ps.sel_timebase = 9
            case 112: self.ps.sel_timebase = 10
            case 128: self.ps.sel_timebase = 11
            case 144: self.ps.sel_timebase = 12
        # Number of Samples and Frequency
        self.frequency    = int(str(self.line_freq.text()))
        self.signal_generator.frequency = self.frequency
        self.ps.sel_nsamp = int(1/(self.ps.sel_dt * self.frequency))
        # Number of Segments
        self.ps.nseg      = int(str(self.line_nseg.text()))
        # Target Noise
        self.target_noise = float(str(self.line_target_noise.text()))
        # Pump Only? true/false checkbox
        self.measure_pumponly = self.check_pumponly.isChecked()
        # TA/TPC checkbox
        self.flag_TATPC = self.check_TATPC.isChecked()
        # Data Channel
        self.ps.data_channel    = self.combo_data_channel.currentIndex()
        # Trigger Channel
        self.ps.trigger_channel = self.combo_trigger_channel.currentIndex()        
        # Trigger Threshold Value
        self.ps.sel_threshold   = int(str(self.line_trigger_threshold.text()))
        # Minimum and Maximum LED current
        self.min_current = float(str(self.min_LED_current.text()))
        self.max_current = float(str(self.max_LED_current.text()))
        # Number of Steps for LED
        self.steps_LED   = int(str(self.steps_LED_current.text()))
        # Choice for time/noise to end measurement
        self.time_noise  = self.combo_Time_Noise.currentIndex()
        # Time Limit
        self.time_limit_value = int(str(self.line_time.text()))
        # Checking Ratio time limit
        self.checking_ratio_time_limit = 5
        # Desired delta(V)/V ratio
        self.voltage_ratio_value = float(str(self.line_desired_ratio.text()))
        # Space LED points linearly or log-wise
        self.LED_spacing = self.combo_LED_spacing.currentIndex()
        # Checkbox for doing TPV and TPC in one go
        self.retain_filter_positions = self.psh_pp_filter_TPC_check.isChecked()

        self.output_box.write_to_output_box("User Settings Read Successfully")

    def open_devices(self) -> None:
        ''' Function for Open Devices button. '''    
        self.output_box.write_to_output_box("Trying to Open Picoscope...")
        try:
            self.ps = picoscope()
            self.ps.resolution = int(str(self.time_interval_options[self.combo_dt.currentIndex()]).split(' ')[3]) # determines what resolution to open the picoscope at.
            self.ps.open_device()
        except:
            self.output_box.write_to_output_box("Connection to the Picoscope was UNSUCCESSFUL")  
        else:
            self.output_box.write_to_output_box("Connection to the Picoscope was SUCCESSFUL")

        self.output_box.write_to_output_box("Trying to open the LED connection...")
        try:
            self.LED = TLDC2200()
        except:
            self.output_box.write_to_output_box("Connection to the LED was UNSUCCESSFUL")
        else:
            self.output_box.write_to_output_box("Connection to the LED was SUCCESSFUL")

        self.output_box.write_to_output_box("Trying to open connection to the motor...")
        try:
            self.motor = Motor(b"27601355")
        except:
            self.output_box.write_to_output_box("Connection to the motor was UNSUCCESSFUL")
        else:
            self.output_box.write_to_output_box("Connection to the motor was SUCCESSFUL. Please allow the motor to home itself.")
        
        self.output_box.write_to_output_box("Trying to open connection to the signal generator...")
        try:
            self.signal_generator = SignalGenerator()
            self.signal_generator.open()
        except:
            self.output_box.write_to_output_box("Connection to the signal generator was UNSUCCESSFUL")
        else:
            self.output_box.write_to_output_box("Connection to the signal generator was SUCCESSFUL")

    def stop(self) -> None:
        ''' Function for Stop button. '''
        self.stop_request = True

    def abort(self) -> None:
        ''' Function for Abort button. '''
        #self.stop_request = True
        self.abort_request = True

    def close_devices(self) -> None:
        ''' Function for Close Devices button '''
        self.ps.close_device()
        self.LED.close()
        self.motor.close()
        self.signal_generator.close()
        self.output_box.write_to_output_box("All device connections closed.")

    def setup_experiment(self) -> None:
        self.stop_request = False
        self.abort_request = False
        self.read_user_settings()
        #self.prepare_plots()
        self.total_time = self.ps.sel_dt * self.ps.sel_nsamp 
        self.left_graph.clear_lines()
        self.right_graph.clear_lines()
        self.t = np.linspace(self.ps.sel_dt, self.total_time, self.ps.sel_nsamp) - (self.ps.sel_trig_pos_rel * self.total_time)

        self.ps.set_channels()
        self.ps.get_timebase()
        self.ps.set_trigger()
        if self.ps.sel_rapid_block:
            self.ps.set_memory_segments()
            self.ps.set_no_of_captures()
        if self.check_chref.isChecked():
            self.ps.sel_chon[1] = 1
        else:
            self.ps.sel_chon[1] = 0

    def Voc_Scan(self) -> None:
        """
        This function sweeps the whole range of the bias light in order to gather the Voc (or Jsc).
        The user then determines the appropriate bias light range to use for `pump_probe`
        """
        
        # control the number of samples
        self.read_user_settings()
        self.ps.sel_nsamp = 1000
        #self.prepare_plots()
        self.total_time = self.ps.sel_dt * self.ps.sel_nsamp 
        self.ps.set_channels()
        self.ps.get_timebase()
        if self.check_chref.isChecked():
            self.ps.sel_chon[1] = 1
        else:
            self.ps.sel_chon[1] = 0
            
        # to save the resulting table
        file_name = qtw.QFileDialog.getSaveFileName(self, 'Save File', 'D:\\', '*.csv')
        
        file_name = str(file_name[0])
        self.output_box.write_to_output_box(file_name)
        
        # if file name is empty, return nothing
        if file_name == '':
            return
        
        # Define the values for the LED
        self.LED_current_list = self.LED.generate_current_list(self.min_current, self.max_current, self.steps_LED, self.LED_spacing)  # in Amps
        
        # list for the Voc values
        Vocs = [np.nan] * len(self.LED_current_list)

        self.right_graph.wipe()
        self.right_graph.axes.set_title("$V_{OC}$ Values")
        self.right_graph.axes.set_xscale('log')
        self.right_graph.axes.set_ylabel('$V_{OC}$')
        self.right_graph.axes.set_xlabel('Bias Light Current (A)')
        self.right_graph.axes.set_xlim(self.min_current/2000,self.max_current/500) 

        for id, LED_power in enumerate(np.flip(self.LED_current_list)):
            self.left_graph.axes.set_title('$V_{OC}$ for LED set to ' + f'{LED_power:.4f}')
            
            # sets the LED current
            self.LED.current = LED_power
            # turn on LED and wait a short time
            self.LED.on()
            self.output_box.write_to_output_box("LED on")
            time.sleep(0.5)
            # output applied current
            self.output_box.write_to_output_box(f"Applied LED current (Amps): {self.LED.current:.4f}") # type: ignore
            
            # starting time
            start_time = time.time()
            elapsed_time = 0
            
            # find the proper voltage range
            self.find_voltage_range()

            # initial reading
            _, CB = self.pump_probe_core(20)
            self.left_graph.wipe()
            self.left_graph.axes.set_ylim(np.min(CB)-0.05, np.max(CB)+0.05)
            self.left_graph.axes.set_ylabel('$V_{OC}$ (V)')
            self.left_graph.axes.set_xlabel('Sample')
            left_line = self.left_graph.draw_plot(y=CB,color='tab:blue')
            
            Voc = self.find_voc_or_jsc(left_line,start_time,elapsed_time,CB)
            
            self.output_box.write_to_output_box(f'The Voc for {LED_power:.4f} Amps is {Voc:.4f} Volts.')
            Vocs[-(id+1)] = Voc
            
            # turn off LED and wait a short time
            self.LED.off()
            self.output_box.write_to_output_box("LED off")

            # update the graph to show a scatter of all Voc values on the right
            self.right_graph.axes.set_ylim(0, np.nanmax(np.array(Vocs))+0.1)
            self.right_graph.clear_lines()
            self.right_graph.draw_scatter(self.LED_current_list, Vocs, color="tab:blue")
            qtw.QApplication.processEvents()
            if self.stop_request:
                return None
            if self.abort_request:
                return None
            
            self.wait_for_cooling(LED_power,15)
            
        self.output_box.write_to_output_box("All Voc checks complete")

        # a dialog window pops up prompting user to place 50 ohm resistor in line if they want to also measure Jsc
        Jsc_dialog_box = qtw.QMessageBox(self)
        Jsc_dialog_box.setWindowTitle("TPV to TPC Switch")
        Jsc_dialog_box.setIcon(qtw.QMessageBox.Icon.Question)
        Jsc_dialog_box.setStandardButtons(qtw.QMessageBox.StandardButton.Yes | qtw.QMessageBox.StandardButton.No)
        Jsc_dialog_box.setText("If you would like to measure Jsc as well, place the 50 ohm resistor inline and press 'Yes'. \n Otherwise, press 'No'.")
        
        button_clicked = Jsc_dialog_box.exec()

        if button_clicked == qtw.QMessageBox.StandardButton.Yes:
            #self.setup_experiment()
            
            # list for the Voc (or Jsc) values
            Jscs = [np.nan] * len(self.LED_current_list)
            self.right_graph.wipe()
            self.right_graph.axes.set_title("$J_{SC}$ Values")
            self.right_graph.axes.set_xscale('log')
            self.right_graph.axes.set_yscale('log')
            self.right_graph.axes.set_ylabel('$J_{SC}$')
            self.right_graph.axes.set_xlabel('Bias Light Current (A)')
            self.right_graph.axes.set_xlim(self.min_current/2000,self.max_current/500) 

            for id, LED_power in enumerate(np.flip(self.LED_current_list)):
                self.left_graph.axes.set_title('$V_{OC}$ for LED set to ' + f'{LED_power:.4f}')
                
                # sets the LED current
                self.LED.current = LED_power
                # turn on LED and wait a short time
                self.LED.on() 
                self.output_box.write_to_output_box("LED on")
                time.sleep(0.5)
                # output applied current
                self.output_box.write_to_output_box(f"Applied LED current (Amps): {self.LED.current:.4f}") # type: ignore
                
                # starting time
                start_time = time.time()
                elapsed_time = 0
                
                # find the proper voltage range
                self.find_voltage_range()

                # initial reading
                _, CB = self.pump_probe_core(20)
                self.left_graph.wipe()
                self.left_graph.axes.set_ylim(np.min(CB)-0.01, np.max(CB)+0.01)
                self.left_graph.axes.set_ylabel('$J_{SC}$ (V)')
                self.left_graph.axes.set_xlabel('Sample')
                left_line = self.left_graph.draw_plot(y=CB,color='tab:blue')
                
                Jsc = self.find_voc_or_jsc(left_line,start_time,elapsed_time, CB)
                
                self.output_box.write_to_output_box(f'The Jsc for {LED_power:.4f} Amps is {Jsc:.4f} Volts.')
                Jscs[-(id+1)] = Jsc
                
                # turn off LED and wait a short time
                self.LED.off() 
                self.output_box.write_to_output_box("LED off")

                # update the graph to show a scatter of all Voc values on the right
                self.right_graph.axes.set_ylim(np.nanmin(np.array(Jscs)/10),np.nanmax(np.array(Jscs))*10)
                self.right_graph.clear_lines()
                self.right_graph.draw_scatter(self.LED_current_list, Jscs,color="tab:blue")
                qtw.QApplication.processEvents()
                if self.stop_request:
                    return None
                
                self.wait_for_cooling(LED_power,15)
                
            self.output_box.write_to_output_box("All Jsc checks complete")
            # save just the Voc data
            data = pd.DataFrame({'Bias (A)': self.LED_current_list, 'Voc (V)': Vocs, 'Jsc (V)': Jscs})
            data.to_csv(file_name,index=False)
        else:
            # save just the Voc data
            data = pd.DataFrame({'Bias': self.LED_current_list, 'Voc': Vocs})
            data.to_csv(file_name,index=False)

    def wait_for_cooling(self, LED_power: float, wait_time: float) -> None:
        if LED_power != self.LED_current_list[0]:
            # if not the final light intensity, sleep
            self.wait(wait_time)

    def wait(self, wait_time: float) -> None:
        self.output_box.write_to_output_box(f"Waiting for the sample to cool for {wait_time} seconds.")
        delay_begin_time = time.time()
        delay_end_time = delay_begin_time + wait_time
        while time.time() < delay_end_time:
            qtw.QApplication.processEvents() # this allows the system to remain "usable while waiting"
            if self.abort_request:
                return
            continue
        self.output_box.write_to_output_box("Waiting complete.")

    def find_voc_or_jsc(self, left_line: list, start_time: float, elapsed_time: float, CB: np.ndarray) -> float:
        ii=1
        
        qtw.QApplication.processEvents()
        while (elapsed_time <= 5):
            D = self.pump_probe_core(20)

            CB = (ii*CB + D[1])/(ii+1) # running average
            left_line[0].set_ydata(CB)
            #self.left_graph.axes.set_ylim(np.min(CB)-0.05, np.max(CB)+0.05)
            self.left_graph.draw()
            qtw.QApplication.processEvents()
            ii += 1
            # checks elapsed time. If it is greater than 5, the loop stops.
            elapsed_time = time.time() - start_time
            
        return np.mean(CB) # calculate the Voc or Jsc for that light intensity

    def pump_probe_three(self, mode='measure') -> None:
        '''wrap around pump probe 
        so pump probe can be written such that it can be called with different time windows
        to be stitched together
        '''

        # to save the resulting table
        file_name = qtw.QFileDialog.getSaveFileName(self, 'Save File', 'D:\\', '*.csv')
        
        file_name = str(file_name[0])
        self.output_box.write_to_output_box(file_name)
        
        # if file name is empty, return nothing
        if file_name == '':
            return
        
        self.mode = mode

        # 1. SHORT SCAN
        self.setup_experiment()
        inz = int(self.ps.sel_nsamp * self.ps.sel_trig_pos_rel * 0.8)  # up to here, average is calculated for baseline
        D = self.pump_probe(inz, mode='measure')  # call the main loop
        if self.abort_request: return

        DTA1 = D[0]
        DTD1 = D[1]
        nza1 = D[2]
        nzb1 = D[3]
        t1 = self.t * 1.0

        self.line_dt.setText(str(self.ps.sel_dt * 10.)) # this won't work right bc of the changes to how the time interval is chosen.
        self.setup_experiment()
        inz = int(self.ps.sel_nsamp * self.ps.sel_trig_pos_rel * 0.8)  # up to here, average is calculated for baseline
        D = self.pump_probe(inz, mode='measure')  # call the main loop
        if self.abort_request: return

        DTA2 = D[0]
        DTD2 = D[1]
        nza2 = D[2]
        nzb2 = D[3]
        t2 = self.t * 1.0

        self.line_dt.setText(str(self.ps.sel_dt * 10.)) # this won't work right bc of the changes to how the time interval is chosen.
        self.setup_experiment()
        inz = int(self.ps.sel_nsamp * self.ps.sel_trig_pos_rel * 0.8)  # up to here, average is calculated for baseline
        D = self.pump_probe(inz, mode='measure')  # call the main loop
        if self.abort_request: return

        DTA3 = D[0]
        DTD3 = D[1]
        nza3 = D[2]
        nzb3 = D[3]
        t3 = self.t * 1.0

        if self.measure_pumponly:
            # THE SAMPLE LOOP WITH ONLY PUMP
            msc = qtw.QMessageBox(None)
            msc.setText('Please block the probe beam(s)')
            msc.setStandardButtons(qtw.QMessageBox.StandardButton.Ok | qtw.QMessageBox.StandardButton.Cancel)
            msc.setDefaultButton(qtw.QMessageBox.StandardButton.Ok)
            r = msc.exec()
            msc.accept()
            if r == qtw.QMessageBox.StandardButton.Cancel:
                self.abort_request = True
            if r == qtw.QMessageBox.StandardButton.Ok:
                self.abort_request = False

            self.stop_request = False  # re-arm the stop request

            self.line_dt.setText(str(self.ps.sel_dt * 0.01)) # this won't work right bc of the changes to how the time interval is chosen.
            self.setup_experiment()
            D = self.pump_only(inz, nza1, nzb1, DTD1, mode='measure')  # call the main loop
            if self.abort_request: return

            DTAx1 = D[0]
            DTDx1 = D[1]
            Apu1 = D[2]
            Bpu1 = D[3]

            self.line_dt.setText(str(self.ps.sel_dt * 10.0)) # this won't work right bc of the changes to how the time interval is chosen.
            self.setup_experiment()
            D = self.pump_only(inz, nza2, nzb2, DTD2, mode='measure')  # call the main loop
            if self.abort_request: return

            DTAx2 = D[0]
            DTDx2 = D[1]
            Apu2 = D[2]
            Bpu2 = D[3]

            self.line_dt.setText(str(self.ps.sel_dt * 10.0)) # this won't work right bc of the changes to how the time interval is chosen.
            self.setup_experiment()
            D = self.pump_only(inz, nza3, nzb3, DTD3, mode='measure')  # call the main loop
            if self.abort_request: return

            DTAx3 = D[0]
            DTDx3 = D[1]
            Apu3 = D[2]
            Bpu3 = D[3]

            self.line_dt.setText(str(self.ps.sel_dt * 0.01)) # this won't work right bc of the changes to how the time interval is chosen.
        else:
            DTAx1 = DTA1
            DTDx1 = DTD1
            DTAx2 = DTA2
            DTDx2 = DTD2
            DTAx3 = DTA3
            DTDx3 = DTD3
            Apu1 = Apu2 = Apu3 = Bpu1 = Bpu2 = Bpu3 = DTAx1 * 0.0

        clip2 = (t2 > t1[-1])
        clip3 = (t3 > t2[-1])

        self.output_box.write_to_output_box(clip2)
        self.output_box.write_to_output_box(clip3)

        t2part = t2[clip2]
        t3part = t3[clip3]

        tall = np.hstack((t1, t2part, t3part))
        DTAx = np.hstack((DTAx1, DTAx2[clip2], DTAx3[clip3]))
        DTDx = np.hstack((DTDx1, DTDx2[clip2], DTDx3[clip3]))
        DTA0x = np.hstack((DTA1, DTA2[clip2], DTA3[clip3]))
        DTD0x = np.hstack((DTD1, DTD2[clip2], DTD3[clip3]))
        Apux = np.hstack((Apu1, Apu2[clip2], Apu3[clip3]))
        Bpux = np.hstack((Bpu1, Bpu2[clip2], Bpu3[clip3]))

        raus = DTAx.reshape(-1, 1)
        if self.ps.sel_chon[1]:
            raus = DTD0x.reshape(-1, 1)
            raus = np.hstack((raus, (Apux-Bpux).reshape(-1, 1), DTDx.reshape(-1, 1)))
        else:
            raus = DTA0x.reshape(-1, 1)
            raus = np.hstack((raus, Apux.reshape(-1, 1), DTAx.reshape(-1, 1)))
        self.output_box.write_to_output_box(raus.shape)
        self.output_box.write_to_output_box(tall.shape)
        save_matrix_milano_convention(range(raus.shape[1]), tall, raus, file_name, suff='')

        self.output_box.write_to_output_box(' measurement finished.')

    def pump_probe_two(self, mode='measure') -> None:
        '''wrap around pump probe 
        so pump probe can be written such that it can be called with different time windows
        to be stitched together
        '''

        # to save the resulting table
        file_name = qtw.QFileDialog.getSaveFileName(self, 'Save File', 'D:\\', '*.csv')
        
        file_name = str(file_name[0])
        self.output_box.write_to_output_box(file_name)
        
        # if file name is empty, return nothing
        if file_name == '':
            return
        
        self.mode = mode

        # 1. SHORT SCAN
        self.setup_experiment()
        inz = int(self.ps.sel_nsamp * self.ps.sel_trig_pos_rel * 0.8)  # up to here, average is calculated for baseline
        D = self.pump_probe(inz, mode='measure')  # call the main loop
        if self.abort_request: return

        DTA1 = D[0]
        DTD1 = D[1]
        nza1 = D[2]
        nzb1 = D[3]
        t1 = self.t * 1.0

        self.line_dt.setText(str(self.ps.sel_dt * 10.)) # this won't work right bc of the changes to how the time interval is chosen.
        self.line_nsamp.setText(str(self.ps.sel_nsamp * 10))
        self.setup_experiment()
        inz = int(self.ps.sel_nsamp * self.ps.sel_trig_pos_rel * 0.8)  # up to here, average is calculated for baseline
        D = self.pump_probe(inz, mode='measure')  # call the main loop
        if self.abort_request: return

        DTA2 = D[0]
        DTD2 = D[1]
        nza2 = D[2]
        nzb2 = D[3]
        t2 = self.t * 1.0

        self.line_dt.setText(str(self.ps.sel_dt * 0.1)) # this won't work right bc of the changes to how the time interval is chosen.
        self.line_nsamp.setText(str(int(self.ps.sel_nsamp * 0.1)))

        if self.measure_pumponly:
            # THE SAMPLE LOOP WITH ONLY PUMP
            msc = qtw.QMessageBox(None)
            msc.setText('Please block the probe beam(s)')
            msc.setStandardButtons(qtw.QMessageBox.StandardButton.Ok | qtw.QMessageBox.StandardButton.Cancel)
            msc.setDefaultButton(qtw.QMessageBox.StandardButton.Ok)
            r = msc.exec()
            msc.accept()
            if r == qtw.QMessageBox.StandardButton.Cancel:
                self.abort_request = True
            if r == qtw.QMessageBox.StandardButton.Ok:
                self.abort_request = False

            self.setup_experiment()
            inz = int(self.ps.sel_nsamp * self.ps.sel_trig_pos_rel * 0.8)  # up to here, average is calculated for baseline
            D = self.pump_only(inz, nza1, nzb1, DTD1, mode='measure')  # call the main loop
            if self.abort_request: return

            DTAx1 = D[0]
            DTDx1 = D[1]
            Apu1 = D[2]
            Bpu1 = D[3]

            self.line_dt.setText(str(self.ps.sel_dt * 10.0)) # this won't work right bc of the changes to how the time interval is chosen.
            self.line_nsamp.setText(str(self.ps.sel_nsamp * 10))
            self.setup_experiment()
            inz = int(self.ps.sel_nsamp * self.ps.sel_trig_pos_rel * 0.8)  # up to here, average is calculated for baseline
            D = self.pump_only(inz, nza2, nzb2, DTD2, mode='measure')  # call the main loop
            if self.abort_request: return

            DTAx2 = D[0]
            DTDx2 = D[1]
            Apu2 = D[2]
            Bpu2 = D[3]

            self.line_dt.setText(str(self.ps.sel_dt * 0.1)) # this won't work right bc of the changes to how the time interval is chosen.
            self.line_nsamp.setText(str(int(self.ps.sel_nsamp * 0.1)))

        else:
            DTAx1 = DTA1
            DTDx1 = DTD1
            DTAx2 = DTA2
            DTDx2 = DTD2

            Apu1 = Bpu1 = DTAx1 * 0.0
            Apu2 = Bpu2 = DTAx2 * 0.0

        clip2 = (t2 > t1[-1])

        self.output_box.write_to_output_box(clip2)

        t2part = t2[clip2]

        tall = np.hstack((t1, t2part))
        DTAx = np.hstack((DTAx1, DTAx2[clip2]))
        DTDx = np.hstack((DTDx1, DTDx2[clip2]))
        DTA0x = np.hstack((DTA1, DTA2[clip2]))
        DTD0x = np.hstack((DTD1, DTD2[clip2]))
        Apux = np.hstack((Apu1, Apu2[clip2]))
        Bpux = np.hstack((Bpu1, Bpu2[clip2]))

        raus = DTAx.reshape(-1, 1)
        if self.ps.sel_chon[1]:
            raus = DTD0x.reshape(-1, 1)
            raus = np.hstack((raus, (Apux-Bpux).reshape(-1, 1), DTDx.reshape(-1, 1)))
        else:
            raus = DTA0x.reshape(-1, 1)
            raus = np.hstack((raus, Apux.reshape(-1, 1), DTAx.reshape(-1, 1)))
        self.output_box.write_to_output_box(raus.shape)
        self.output_box.write_to_output_box(tall.shape)
        save_matrix_milano_convention(range(raus.shape[1]), tall, raus, file_name, suff='')

        self.output_box.write_to_output_box(' measurement finished.')

    def pump_probe_freerun(self, mode='freerun') -> None:
        self.mode = mode
        self.setup_experiment()
        inz = int(self.ps.sel_nsamp * self.ps.sel_trig_pos_rel * 0.8)  # up to here, average is calculated for baseline
        # run a block to obtain array size
        D = self.pump_probe(inz, mode='freerun')  # call the main loop

    def turn_on_devices(self,LED_power) -> None:
        self.LED.current = LED_power
        self.LED.on()
        self.signal_generator.output_on()
        self.output_box.write_to_output_box("LED and Laser On")
        time.sleep(0.5)
        self.output_box.write_to_output_box(f"Applied LED current (Amps): {self.LED.current:.4f}") # type: ignore

        self.find_voltage_range()

    def turn_off_devices(self) -> None:
        self.LED.off() 
        self.signal_generator.output_off()
        self.output_box.write_to_output_box("LED and Laser Off")

    def find_pulse_frequency(self, time, CB: np.ndarray) -> bool:
        """
        This function checks the slope of the data channel prior to the pulse.
        If the slope is negative and large in magnitude, decrease the signal generator frequency.
        Otherwise, proceed with the measurement.
        
        # Inputs
        - time: The time axis
        - CB: The data channel
        
        # Outputs
        - True if the frequency is fine. The measurement will proceed.
        - False if the frequency is not fine. The frequency will be decreased and the check will be repeated.
        """
        
        # find the values of the data channel prior to the peak
        index = 9*CB.argmax()//10
        CB_subset   = CB[:index]
        time_subset = time[:index]
        voc = np.mean(CB_subset)
        # linear regression of the subsets
        slope, intercept = np.polyfit(time_subset, CB_subset, 1)
        
        # if slope is positive enough, proceed with the measurment, else decrease the frequency and repeat
        if slope/voc > -1.5:
            self.output_box.write_to_output_box(f"Slope of the data channel prior to the pulse is {slope:.4f}. Proceeding with the measurement.")
            return True
        else:
            self.output_box.write_to_output_box(f"Slope of the data channel prior to the pulse is {slope:.4f}. Decreasing the pulse frequency.")
            # adjust the necessary values to change the frequency and related variables.
            self.signal_generator.frequency = int(self.signal_generator.frequency * 0.75) # change the frequency
            self.ps.sel_nsamp = int(1/(self.ps.sel_dt * self.signal_generator.frequency)) # update the number of samples
            self.line_freq.setText(str(self.signal_generator.frequency))                  # change the value in the text box to reflect the new value
            self.ps.set_channels()
            self.ps.get_timebase()
            self.ps.set_trigger()
            self.output_box.write_to_output_box(f"New pulse frequency: {self.signal_generator.frequency} Hz")
            return False

    def pump_probe_one(self, mode='measure') -> None:
        ''' 
        Used for TPV/TPC measurements.
        This function iterates through a list of a bias light intensities as defined by `self.LED.generate_current_list(min,max,num,spacing)`.
        It saves the resulting numpy arrays to `.dat` files.
        The files are named according to `file_name` and a suffix added for the bias light current.
        The first column is the time in seconds and the second column is the data in volts.
        The other two columns can be ignored.
        '''

        # time axis label for the graph
        self.right_graph.axes.set_xlabel("Time (s)")
        self.right_graph.refresh()
        qtw.QApplication.processEvents()

        self.doing_ratio_checking = False

        file_name = self.save_dialog()
        if file_name == '': return

        self.mode = mode
        self.ps.offset = 0
        self.setup_experiment()
        # inz = number of samples * relative position of trigger in collected data * 0.8. Not Sure the purpose of this
        inz = int(self.ps.sel_nsamp * self.ps.sel_trig_pos_rel * 0.8)  # up to here, average is calculated for baseline
        # run a block to obtain array size

        # Define the values for the LED
        self.LED_current_list = self.LED.generate_current_list(self.min_current, self.max_current, self.steps_LED, self.LED_spacing)  # in Amps

        # Set the frequency for the laser
        self.signal_generator.frequency = self.frequency

        # grab current time to see how long the measurements take
        start_time = time.time()

        # this for loop was added by Steven
        for LED_power in np.flip(self.LED_current_list):
            self.turn_on_devices(LED_power)
            
            # calls the main loop to take measurement
            D = self.pump_probe(inz, mode)

            if self.stop_request: continue # moves to the next iteration of the for loop
            if self.abort_request:
                self.turn_off_devices()
                break # stops the for loop, function then ends

            # split apart data from measurement
            DTA = D[0]
            DTD = D[1]
            nza = D[2]
            nzb = D[3]

            if self.measure_pumponly: # this isn't evaluated as true under normal operation.
                # THE SAMPLE LOOP WITH ONLY PUMP
                msc = qtw.QMessageBox(None)
                msc.setText('Please block the probe beam(s)')
                msc.setStandardButtons(qtw.QMessageBox.StandardButton.Ok | qtw.QMessageBox.StandardButton.Cancel)
                msc.setDefaultButton(qtw.QMessageBox.StandardButton.Ok)
                r = msc.exec()
                msc.accept()
                if r == qtw.QMessageBox.StandardButton.Cancel:
                    self.abort_request = True
                if r == qtw.QMessageBox.StandardButton.Ok:
                    self.abort_request = False

                self.stop_request = False  # re-arm the stop request
                nza0 = nza * 1.0
                nzb0 = nzb * 1.0

                # call the main loop
                D = self.pump_only(inz, nza0, nzb0, DTD, mode='measure')
                if self.abort_request: return

                DTAx = D[0]
                DTDx = D[1]
                Apu  = D[2]
                Bpu  = D[3]
            else:
                DTAx = DTA
                DTDx = DTD
                Apu = Bpu = DTA * 0.0

            raus = DTAx.reshape(-1, 1)
            if self.ps.sel_chon[1]:
                raus = DTD.reshape(-1, 1)
                raus = np.hstack((raus, (Apu-Bpu).reshape(-1, 1), DTDx.reshape(-1, 1)))
            else:
                raus = DTA.reshape(-1, 1)
                raus = np.hstack((raus, Apu.reshape(-1, 1), DTAx.reshape(-1, 1)))

            self.turn_off_devices()
            
            self.save_data(DTDx, self.tA, file_name, LED_power)
            
            self.wait_for_cooling(LED_power,30)
            
        self.total_runtime(start_time)

    def pump_probe_one_filter(self, mode='measure') -> None:
        """
        Operation is similar to the `pump_probe_one` function.
        This function does TPV measurements over a series of bias light intensities while maintaining a ratio of delta(V)/V by incrementing the motor controlled with `MotorControl.py`.
        If desired, the function prompts the user to change the setup to do TPC measurements immediately after TPV.
        """

        filter_positions = []
        self.early_return = False
        
        file_name = self.save_dialog()
        if file_name == '': return

        self.mode = mode
        self.ps.offset = 0
        self.setup_experiment()
        # inz = number of samples * relative position of trigger in collected data * 0.8. Not Sure the purpose of this
        inz = int(self.ps.sel_nsamp * self.ps.sel_trig_pos_rel * 0.8)  # up to here, average is calculated for baseline
        # run a block to obtain array size

        # Define the values for the LED
        self.LED_current_list = self.LED.generate_current_list(self.min_current, self.max_current, self.steps_LED, self.LED_spacing)  # in Amps

        # Set the frequency for the laser
        self.signal_generator.frequency = self.frequency

        # grab current time to see how long the measurements take
        start_time = time.time()
        
        # this for loop was added by Steven
        for LED_power in np.flip(self.LED_current_list):
            self.doing_ratio_checking = True
            final_pass = False

            self.turn_on_devices(LED_power)
            
            # while loop to check delta(V)/V ratio
            ratio_satisfied = False # initialize break condition for when ratio is met
            while True:
                # adjusts the voltage range of the picoscope
                self.find_voltage_range()
                
                # calls the main loop to take measurement
                D = self.pump_probe(inz, mode)

                if self.stop_request: break # stops the while loop
                if self.abort_request: break # stops the while loop

                # split apart data from measurement
                DTA = D[0]
                DTD = D[1]
                nza = D[2]
                nzb = D[3]

                if self.measure_pumponly:
                    # THE SAMPLE LOOP WITH ONLY PUMP
                    msc = qtw.QMessageBox(None)
                    msc.setText('Please block the probe beam(s)')
                    msc.setStandardButtons(qtw.QMessageBox.StandardButton.Ok | qtw.QMessageBox.StandardButton.Cancel)
                    msc.setDefaultButton(qtw.QMessageBox.StandardButton.Ok)
                    r = msc.exec()
                    msc.accept()
                    if r == qtw.QMessageBox.StandardButton.Cancel:
                        self.abort_request = True
                    if r == qtw.QMessageBox.StandardButton.Ok:
                        self.abort_request = False

                    self.stop_request = False  # re-arm the stop request
                    nza0 = nza * 1.0
                    nzb0 = nzb * 1.0

                    # call the main loop
                    D = self.pump_only(inz, nza0, nzb0, DTD, mode='measure')
                    if self.abort_request: return

                    DTAx = D[0]
                    DTDx = D[1]
                    Apu  = D[2]
                    Bpu  = D[3]
                else:
                    DTAx = DTA
                    DTDx = DTD
                    Apu = Bpu = DTA * 0.0

                raus = DTAx.reshape(-1, 1)
                if self.ps.sel_chon[1]:
                    raus = DTD.reshape(-1, 1)
                    raus = np.hstack((raus, (Apu-Bpu).reshape(-1, 1), DTDx.reshape(-1, 1)))
                else:
                    raus = DTA.reshape(-1, 1)
                    raus = np.hstack((raus, Apu.reshape(-1, 1), DTAx.reshape(-1, 1)))
                self.output_box.write_to_output_box(raus.shape)

                if final_pass:
                    break
                self.turn_off_devices()
                ratio_satisfied = self.filter_laser_power(D[1])
                self.turn_on_devices(LED_power)
                if ratio_satisfied:
                    self.doing_ratio_checking = False
                    final_pass = True
            
            self.turn_off_devices()
            
            # grab the position of the filter
            current_position = self.motor.position
            self.output_box.write_to_output_box(f'Current Position: {current_position:.3f} mm')
            filter_positions.append(current_position)

            self.save_data(DTDx, self.tA, file_name, LED_power)
            
            if self.abort_request:
                break # breaks the for loop

            if self.early_return == False:
                # move motor small amount in preperation for next light intensity
                current_position = self.motor.position
                self.output_box.write_to_output_box(f'Current Position: {current_position:.3f} mm')
                self.motor.position = current_position - 1.0
            else: self.early_return = False # re-arms the boolean variable and moves on. will only be reached when early_return is true
            
            self.wait_for_cooling(LED_power,30)

        self.total_runtime(start_time)

        if self.retain_filter_positions == True:
            np.savetxt('filter_positions.csv', filter_positions, delimiter = ',')
            self.retain_filter_positions = False # changes it to false to prevent infinite looping

            # a dialog window pops up prompting user to place 50 ohm resistor in line
            TPC_dialog_box = qtw.QMessageBox(self)
            TPC_dialog_box.setWindowTitle("TPV to TPC Switch")
            TPC_dialog_box.setIcon(qtw.QMessageBox.Icon.Information)
            TPC_dialog_box.setText("Please insert the 50 ohm resistor inline to the picoscope.\nThen hit the button.")
            
            button_clicked = TPC_dialog_box.exec()

            if button_clicked == qtw.QMessageBox.StandardButton.Ok:
                self.pump_probe_one_filter_TPC(mode="measure")

    def filter_laser_power(self,voltage_values: np.ndarray) -> bool:
        """
        This function checks the `delta(V)/V` ratio against the desired ratio.
        
        ## Paramters:
        1. self
        2. voltage_values: A 1-D numpy array of the voltage measurements.
        
        ## Returns:
        A boolean value which indicates whether the ratio condition is satifisied.
        Once it is satisfied, `True` is returned, the while loop stops. and the measurement is saved.
        
        ## Process:
        - Find the baseline voltage and maximum voltage.
        - Define an accpetable range surrounding the deisred voltage ratio.
        - Adjust the position of the filter a small amount trying to get closer to the desired voltage ratio.
        """
        
        # find the location and amplitude of the peak voltage value
        index_peak  = np.argmax(voltage_values)
        max_voltage = np.amax(voltage_values)
        
        # take half of the data before that peak and take the mean
        baseline_voltage: float = np.mean(voltage_values[20:index_peak//2]) # the first 20 are ignored to eliminate any boundary effects from `np.convolve` (not yet in place but could be in the future)
        # find the ratio between change in voltage from baseline to peak.
        actual_ratio    : float = (max_voltage - baseline_voltage) / baseline_voltage
        # if ratio is within 20% (for example) of desired ratio then move on to save the data
        ratio_leniency = 0.2
        minimum_acceptable_ratio: float = (1 - ratio_leniency) * self.voltage_ratio_value
        maximum_acceptable_ratio: float = (1 + ratio_leniency) * self.voltage_ratio_value
        
        self.output_box.write_to_output_box(f'Baseline Voltage:        {baseline_voltage:.3f} V')
        self.output_box.write_to_output_box(f'Max Voltage:             {max_voltage:.3f} V')
        self.output_box.write_to_output_box(f'Acceptable Ratio Range: ({minimum_acceptable_ratio:.2%}, {maximum_acceptable_ratio:.2%})')
        self.output_box.write_to_output_box(f'Actual Ratio:            {actual_ratio:.2%}')
        
        if (minimum_acceptable_ratio <= actual_ratio) and (actual_ratio <= maximum_acceptable_ratio):
            # if voltage ratio is in the "acceptable region" then end the while loop and move on to saving.
            self.output_box.write_to_output_box(f'Actual Ratio: {actual_ratio:.2%} is within {ratio_leniency:2.0%} of the desired ratio of {self.voltage_ratio_value:.2%}.')
            return True # ends the while loop
            
        # get the current position of the motor
        current_position = self.motor.position
        self.output_box.write_to_output_box(f'Current Position: {current_position:.3f} mm')
        
        difference_in_ratio = abs(self.voltage_ratio_value - actual_ratio)/self.voltage_ratio_value
        
        if (difference_in_ratio >= 0.1):
            increment_mm = 1.0
            wait = 1.00
        if (difference_in_ratio >= 0.2):
            increment_mm = 2.0 
            wait = 2.00
        if (difference_in_ratio >= 0.3):
            increment_mm = 3.0
            wait = 3.00
        if (difference_in_ratio >= 0.4):
            increment_mm = 4.0
            wait = 4.00
        if (difference_in_ratio >= 0.5):
            increment_mm = 5.0
            wait = 5.00 
        
        if (actual_ratio < self.voltage_ratio_value):
            # if the current ratio is less than the desired ratio
            if (current_position > 53):
                # saftey condition to stop when the motor is at its maximum extent
                self.output_box.write_to_output_box("Desired ratio is unachievable. Motor is at its extreme position.")
                self.early_return = True
                return True
            
            self.motor.position = current_position + increment_mm
            self.output_box.write_to_output_box(f'The stage has attempted to move to {current_position+increment_mm:.3f} mm.')
        
        if (actual_ratio > self.voltage_ratio_value):
            # if the current ratio is larger than the desired ratio
            if (current_position < 0.05):
                # saftey condition to stop when the motor is at its minimum extent
                self.output_box.write_to_output_box("Desired ratio is unachievable. Motor is at its extreme position.")
                self.early_return = True
                return True

            self.motor.position = current_position - increment_mm
            self.output_box.write_to_output_box(f'The stage has attempted to move to {current_position-increment_mm:.3f} mm.')
    
        self.wait(wait)
        return False

    def pump_probe_one_filter_TPC(self, mode="measure") -> None:
        """
        This function is called at the end of `pump_probe_one_filter` when `retain_filter_positions` is checked as `True`.
        """
        
        self.doing_ratio_checking = False

        # import list of filter positions and reverses it
        filter_positions_mm: list[float] = np.loadtxt("filter_positions.csv").tolist()

        file_name = self.save_dialog()
        if file_name == '': return
        
        self.mode = mode
        self.ps.offset = 0
        self.setup_experiment()
        # inz = number of samples * relative position of trigger in collected data * 0.8. Not Sure the purpose of this
        inz = int(self.ps.sel_nsamp * self.ps.sel_trig_pos_rel * 0.8)  # up to here, average is calculated for baseline
        # run a block to obtain array size

        # Define the values for the LED
        self.LED_current_list = self.LED.generate_current_list(self.min_current, self.max_current, self.steps_LED, self.LED_spacing)  # in Amps

        # grab current time to see how long the measurements take
        start_time = time.time()

        # this for loop was added by Steven
        for id, LED_power in enumerate(np.flip(self.LED_current_list)):
            
            # move motor to correct position
            current_position = self.motor.position
            self.output_box.write_to_output_box(f'Current Position: {current_position:.3f} mm')
            if type(filter_positions_mm)==float:
                desired_position = filter_positions_mm
            else:
                desired_position = filter_positions_mm[id]

            self.motor.position = desired_position
            self.output_box.write_to_output_box(f'The stage has attempted to move to {desired_position:.3f} mm.')
            self.wait(30)

            self.turn_on_devices(LED_power)
            
            # calls the main loop to take measurement
            D = self.pump_probe(inz, mode)

            if self.stop_request: continue
            if self.abort_request: break

            # split apart data from measurement
            DTA = D[0]
            DTD = D[1]
            nza = D[2]
            nzb = D[3]

            if self.measure_pumponly:
                # THE SAMPLE LOOP WITH ONLY PUMP
                msc = qtw.QMessageBox(None)
                msc.setText('Please block the probe beam(s)')
                msc.setStandardButtons(qtw.QMessageBox.StandardButton.Ok | qtw.QMessageBox.StandardButton.Cancel)
                msc.setDefaultButton(qtw.QMessageBox.StandardButton.Ok)
                r = msc.exec()
                msc.accept()
                if r == qtw.QMessageBox.StandardButton.Cancel:
                    self.abort_request = True
                if r == qtw.QMessageBox.StandardButton.Ok:
                    self.abort_request = False

                self.stop_request = False  # re-arm the stop request
                nza0 = nza * 1.0
                nzb0 = nzb * 1.0

                # call the main loop
                D = self.pump_only(inz, nza0, nzb0, DTD, mode='measure')
                if self.abort_request: return

                DTAx = D[0]
                DTDx = D[1]
                Apu = D[2]
                Bpu = D[3]
            else:
                DTAx = DTA
                DTDx = DTD
                Apu = Bpu = DTA * 0.0

            raus = DTAx.reshape(-1, 1)
            if self.ps.sel_chon[1]:
                raus = DTD.reshape(-1, 1)
                raus = np.hstack((raus, (Apu-Bpu).reshape(-1, 1), DTDx.reshape(-1, 1)))
            else:
                raus = DTA.reshape(-1, 1)
                raus = np.hstack((raus, Apu.reshape(-1, 1), DTAx.reshape(-1, 1)))

            self.turn_off_devices()
            
            self.save_data(DTDx, self.tA, file_name, LED_power)
            
        self.total_runtime(start_time)
 
    def pump_probe(self, inz: int, mode='measure') -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """
        Calls `pump_probe_core` and takes a running average of the data.
        Returns the data after averaging together calls to `pump_probe_core` until an end condition is met (time or noise).
        Also plots the raw data and trigger channel on the left graph, and the averaged data and smoothed data on the right graph.
        
        # Inputs:
        - inz: Designates the portion of the data used for calculating the noise level.
        - mode: the mode of operation, either 'measure' or 'freerun'
        
        # Returns:
        A tuple containing the trigger channel data, the averaged data, 2 variables whose purpose is related to the noise level, and the smoothed data
        """

        # THE SAMPLE LOOP WITH PUMP AND PROBE
        ii = 0  # the total number of exposures
        noise = 1000  # initial noise level
        elapsed_time = 0  # sets the duration to 0

        if mode == 'measure':
            self.ncore = 20
        elif mode == 'freerun':
            self.ncore = 2

        # initialize the end condition
        end_condition = False

        # get the starting time
        start_time = time.time()

        while (self.stop_request == False) and (end_condition == False):
            '''
            While not told to stop and while the end condition remains false.
            Causes the measurement cycle to repeat until the time (or noise) limit is met.
            '''
            D = self.pump_probe_core(self.ncore)  # the core data acquisition routine
    
            if ii == 0:  # first interation
                CA = D[0]
                CB = D[1]
            else:  # all following iterations
                CA = (CA * ii + D[0]) / (ii + 1) # running average
                CB = (CB * ii + D[1]) / (ii + 1)

            # finds the average of CA from beginning to inz (what is inz?)
            nza = np.mean(CA[:inz])  # if pumponly, use the old self.nza
            DTA = (CA - nza)/nza

            if self.ps.sel_chon[1]:
                ''' if using both channels, also repeat for channel B '''
                nzb = np.mean(CB[:inz])
                DTB = (CB - nzb)/nzb
                DTD = DTA - DTB
            else:  # zero out all arrays
                nzb = 0.0 * nza
                DTB = 0.0 * DTA
                DTD = 0.0 * DTA

            self.output_box.write_to_output_box(f'Current exposure: {ii * 20 * self.ps.nseg}') # idk what the 20 is for
            
            if self.ps.sel_chon[1]:
                noise = np.std(DTB[:inz])
            else:
                noise = np.std(DTA[:inz])
            self.output_box.write_to_output_box(f'noise of DT/T trace= {noise}')

            self.output_box.write_to_output_box(f'absolute peak value: {np.amax(D[0])-np.mean(D[0][:inz])} V')
            if self.ps.sel_chon[1]:
                finalDT = np.mean(DTB[2*inz:])
                self.output_box.write_to_output_box(f'final DT/T: {finalDT}')

            self.plot_right_graph(CB) # self.tA is defined in this function
            
            num_points_to_smooth = 5
            CB_conv = np.convolve(CB,np.ones(num_points_to_smooth),mode='same')/num_points_to_smooth
            self.right_graph.draw_plot(self.tA, CB_conv, color='tab:orange') # smoothed voltage reading

            self.left_graph.clear_lines()
            self.left_graph.draw_plot(self.tA,D[0],color='black')      # trigger channel
            self.left_graph.draw_plot(self.tA,D[1],color='tab:blue') # most recent voltage reading
            
            if ii==0:
                # adjust the y-axis limits so that the data is always within the visible range
                min = np.min(CB)
                max = np.max(CB)
                self.right_graph.axes.set_ylim(min-0.05*min,max+0.05*max)
                self.right_graph.refresh()

            qtw.QApplication.processEvents()
            if mode == 'freerun':
                CA = CA * 0.0
                if self.ps.sel_chon[1]:
                    CB = CB * 0.0
                ii = 0

            if self.doing_ratio_checking == True:
                elapsed_time = time.time()-start_time
                if elapsed_time > 3:
                    break
        
            end_condition = self.check_end_condition(start_time, noise)

            if self.abort_request == True:
                self.turn_off_devices()
                end_condition = True
            # end of while loop
            ii += 1

        if self.flag_TATPC:
            return DTA, CB, nza, nzb, CB_conv  # if user desires TA on CHA and TPC on CHB
        else:
            return DTA, DTD, nza, nzb, CB_conv

    def pump_probe_core(self, ns: int) -> tuple[np.ndarray, np.ndarray]:
        """
        Calls to the picoscope to begin gathering data. The picoscope runs a block and it is averaged together.
        
        # Inputs:
        - ns: the number of segments to average together
        
        # Returns:
        A tuple containing the trigger channel data and the averaged data.
        """

        if self.ps.sel_chon[1]:
            self.output_box.write_to_output_box('CHB is ON')
        else:
            self.output_box.write_to_output_box('CHB is OFF')

        for kk in range(ns):
            self.ps.flag_newdata = False
            self.ps.run_block()
            while not self.ps.flag_newdata:  # wait until driver supplies data
                qtw.QApplication.processEvents()
                time.sleep(0.01)
            
            if self.ps.sel_rapid_block and self.ps.sel_chon[1]:
                ADCA = np.mean(self.ps.ADCA, 0)
                ADCB = np.mean(self.ps.ADCB, 0)
            elif not self.ps.sel_rapid_block and self.ps.sel_chon[1]:
                ADCA = self.ps.ADCA
                ADCB = self.ps.ADCB
            elif self.ps.sel_rapid_block and not self.ps.sel_chon[1]:
                ADCA = np.mean(self.ps.ADCA, 0)
                ADCB = None
            else: # if neither then
                ADCA = self.ps.ADCA
                ADCB = None
            
            if kk == 0:
                CA = ADCA
                if self.ps.sel_chon[1]:
                    CB = ADCB
                else:
                    CB = 0.0*CA
            else:
                if self.ps.sel_chon[1]:
                    CB = (CB * kk + ADCB) / (kk + 1) # running average
                else:
                    CA = (CA * kk + ADCA) / (kk + 1)
                    CB = 0.0*CA
            if self.abort_request == True:
                return CA, CB
        return CA, CB

    def find_voltage_range(self):
        """ Finds the range of the voltage trace and adjusts the range on the picoscope accordingly. """

        self.ps.sel_v = 7
        self.ps.offset = 0 # resets the offset to zero

        self.ps.enable_channel(self.ps.data_channel, self.ps.sel_v, self.ps.offset, "Data") # temporarily sets the data channel to a range of +/- 2 V with no offset

        def round_down(value: float, possible_values: list[int]) -> tuple[int, int]:
            """ Rounds a value down to the nearest value in a list of possible values. """
            for i in range(len(possible_values)):
                if (value*1.20/2) <= possible_values[i]: # gives 20% wiggle room
                    return possible_values[i], i
                
            # if the value is greater than all of the possible voltage values, then return the largest one.
            return possible_values[-1], len(possible_values)
        
        possible_voltage_ranges    = [1e-2,2e-2,5e-2,0.1,0.2,0.5,1,2,5,10,20] # V
        wide_offset_voltage_ranges = [                       0.5,1,2,5,10,20] # V
        # calls pump_probe_core to gather data
        CA, CB = self.pump_probe_core(20)
        
        self.plot_right_graph(CB)
        
        # finds the range of the voltage trace
        minCB, maxCB = np.min(CB), np.max(CB)
        voltage_range = maxCB - minCB
        
        # round down to the nearest possible voltage range for the picoscope
        chosen_voltage_range, self.ps.sel_v = round_down(voltage_range, possible_voltage_ranges) # in V
        
        # the offset is how far to push the voltage down by in order to stay within the range of the picoscope.
        self.ps.offset = -minCB-chosen_voltage_range/2 # V
        
        if (self.ps.offset <= -0.25) & (self.ps.sel_v <= 4):
            # this is needed because with narrow voltage ranges, the picoscope is restricted to an offset of +/- 250 mV.
            chosen_voltage_range, temp_v = round_down(voltage_range, wide_offset_voltage_ranges) # in V
            self.ps.sel_v = temp_v + 5
            self.ps.offset = -abs(minCB)-chosen_voltage_range/2 # V

        self.output_box.write_to_output_box(f'Observed Range: ({minCB:.3f}, {maxCB:.3f}) V -> {voltage_range:.3f} V')
        self.output_box.write_to_output_box(f'Voltage Range: +/-{chosen_voltage_range} V ({self.ps.sel_v})')
        self.output_box.write_to_output_box(f'Offset: {self.ps.offset:.3f} V')
        self.output_box.write_to_output_box(f'Perceived Picoscope Range: ({-chosen_voltage_range:.3f}, {chosen_voltage_range:.3f}) V')
        self.output_box.write_to_output_box(f'Effective Picoscope Range: ({-chosen_voltage_range -self.ps.offset:.3f}, {chosen_voltage_range -self.ps.offset:.3f}) V')

        # re-open the data channel on the picoscope to have the appropriate voltage range
        self.ps.enable_channel(self.ps.data_channel, self.ps.sel_v, self.ps.offset, "Data") 

    def check_end_condition(self, start_time: float, noise) -> bool:
        """ Checks if the end condition is satisfied. """
        
        elapsed_time: float = time.time() - start_time
        self.output_box.write_to_output_box(f"Elapsed Time:        {elapsed_time:.3f}")
        self.output_box.write_to_output_box(f"Current Noise Level: {noise:.3e}")
        
        end_time : bool = elapsed_time > self.time_limit_value
        end_noise: bool = noise        < self.target_noise
        
        match self.time_noise:
            case 0: return end_time
            case 1: return end_noise
            case 2: return (end_time | end_noise)
            case _: return False

    def total_runtime(self, start_time: float) -> None:
        # end of for loop
        self.output_box.write_to_output_box("All measurements finished.")
        # measures and prints the time elapsed
        finish_time = time.time()
        measurement_time = finish_time - start_time
        self.output_box.write_to_output_box(f"Total Time: {measurement_time:.2f} seconds")

    def save_dialog(self) -> str:
        file_name = qtw.QFileDialog.getSaveFileName(self, 'Save File', 'D:\\', '*.dat')
        file_name = str(file_name)[2:-11]
        self.output_box.write_to_output_box(file_name)
        return file_name

    def save_data(self, voltage, time, file_name, LED_power) -> None:
        """ Saves the data. """
        bias_power_str = f"{LED_power:07.4f}" # total length of up to 7 with 4 decimal places, as a float
        
        full_file_name = f"{file_name[:-4]}_{bias_power_str}.dat"
        np.savetxt(full_file_name, np.stack((time,voltage),axis=-1), fmt='%.7e', delimiter=',')
        #save_matrix_milano_convention(range(raus.shape[1]), t, raus, file_name, suff=bias_power_str)  # added suffix parameter
        self.output_box.write_to_output_box('Measurement with current ' + bias_power_str + ' Amps is finished and saved.')

    def pump_only(self, inz, nza0, nzb0, DTD, mode='measure') -> tuple:
        """
        the core loop
        - mode = measure: saves and runs infinitel
        - mode = freerun: runs looped
        - mode = pumponly: assumes pump is off, and subtracts reference data from very last pump-probe measurement
        """

        ii = 0  # the total number of exposures
        while not self.stop_request:
            D = self.pump_probe_core(20)  # the core data acquisition routine
            if ii == 0:
                CA = D[0]
                CB = D[1]
            else:
                CA = (CA * ii + D[0]) / (ii + 1)
                CB = (CB * ii + D[1]) / (ii + 1)

            nza = np.mean(CA[:inz])  # if pumponly, use the old self.nza
            Apu = (CA - nza)/nza0
            DTAx = DTD - Apu

            if self.ps.sel_chon[1]:
                nzb = np.mean(CB[:inz])
                Bpu = (CB - nzb)/nzb0
                DTDx = DTD - (Apu - Bpu)

            self.output_box.write_to_output_box(nza0)
            self.output_box.write_to_output_box(nzb0)

            self.output_box.write_to_output_box('Current exposure:' + str(ii * 20 * self.ps.nseg))
            ii += 1

            if self.ps.sel_chon[1]:
                noise = np.std(DTDx[:inz])
            else:
                noise = np.std(DTAx[:inz])
            self.output_box.write_to_output_box('noise = ' + str(noise))
            self.output_box.write_to_output_box('absolute peak value:' + str(np.amax(CA)-nza) + 'V')

            self.plot_waveform(Apu, Bpu)
            if self.ps.sel_chon[1]:
                self.plot_pp(DTDx)
            else:
                self.plot_pp(DTAx)
        return DTAx, DTDx, Apu, Bpu

    def plot_right_graph(self, CB: np.ndarray) -> None:
        """ Calculates the time axis and plots the averaged voltage on the right graph. """
        
        self.tA = np.linspace(self.ps.sel_dt, self.total_time, self.ps.sel_nsamp) # generates the time axis
        self.tA = self.tA - self.tA[CB.argmax()]                                  # centers the time axis over the peak of the data channel
        self.right_graph.clear_lines()
        voc = np.mean(CB[:CB.argmax()//2])
        self.right_graph.draw_baseline(self.tA, voc) # baseline (Voc)
        self.right_graph.draw_plot(self.tA, CB, color='tab:blue') # averaged voltage reading

def main():
    app = qtw.QApplication(sys.argv)
    form = muspp_ui()
    form.showMaximized()
    app.exec()

if __name__ == "__main__":
    main()