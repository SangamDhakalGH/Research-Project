import sys
import os
import time
import locale
import scipy.ndimage
import PyQt6.QtCore as qtc
import PyQt6.QtGui as qtg
import PyQt6.QtWidgets as qtw
import matplotlib as mp
mp.use('QtAgg')
#import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from picoscope import picoscope
from MatplotlibCanvas import MatplotlibCanvas
from femtofunc import binnen, save_matrix_milano_convention
from TLDC2200_pyvisa import generate_current_list, LED_open_connection, LED_close_connection
import MotorControl

mp.rcParams['font.size'] = 8
locale.setlocale(locale.LC_NUMERIC, 'C')

class muspp_ui(qtw.QMainWindow):    
    def __init__(self):
        ''' Setting up the interface. '''
        qtw.QMainWindow.__init__(self)

        layout      = qtw.QHBoxLayout()
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
        # Box. I don't know what purpose this has.
        self.textEdit = qtw.QTextEdit(maximumWidth=300, minimumWidth=100)
        self.textEdit.setText("Fit started...")
        # Minimum for LED current
        self.min_LED_current_label = qtw.QLabel(text='Min Bias Light Current (mA)', maximumWidth=200)
        self.min_LED_current       = qtw.QLineEdit(text='10',                       maximumWidth=100)
        # Maximum for LED current
        self.max_LED_current_label = qtw.QLabel(text='Max Bias Light Current (mA)', maximumWidth=200)
        self.max_LED_current       = qtw.QLineEdit('10000',                         maximumWidth=100)
        # Number of points to set LED at
        self.steps_LED_current_label = qtw.QLabel(text='Number of Bias Values', maximumWidth=150)
        self.steps_LED_current       = qtw.QLineEdit(text='15',                 maximumWidth=100)
        # Drop-down Menu for Lin or Log evenly spaced light intensities
        self.label_LED_spacing = qtw.QLabel(text='Bias Intensity Spacing', maximumWidth=150)
        self.combo_LED_spacing = qtw.QComboBox(maximumWidth=100)
        self.combo_LED_spacing.addItem("Linear",      userData=qtc.QVariant(0))
        self.combo_LED_spacing.addItem("Logarithmic", userData=qtc.QVariant(1))
        # Voltage Combo Box
        self.label_V = qtw.QLabel(text='Y Scale', maximumWidth=100)
        self.combo_V = qtw.QComboBox(maximumWidth=100)
        # loops to add the voltage options (changed from explicit declaration by Steven)
        voltage_options = ['10 mV',  '20 mV', '50 mV', '100 mV', '200 mV',
                           '500 mV', '1 V',     '2 V',    '5 V',   '10 V']
        for i in range(0, 10):
            self.combo_V.addItem(voltage_options[i], userData=qtc.QVariant(i))
        self.combo_V.setCurrentIndex(9)  # sets default ot '10 V'
        # Time Interval Input Box
        self.label_dt    = qtw.QLabel(text='Interval (s)',      maximumWidth=100)
        self.line_dt     = qtw.QLineEdit(text='1e-8',           maximumWidth=100)
        # Number of Samples Input Box
        self.label_nsamp = qtw.QLabel(text='Number of Samples', maximumWidth=100)
        self.line_nsamp  = qtw.QLineEdit(text='5000',           maximumWidth=100)
        # Number of Segments Input Box
        self.label_nseg  = qtw.QLabel(text='Number of Segments',maximumWidth=100)
        self.line_nseg   = qtw.QLineEdit(text='300',            maximumWidth=100)
        # Delta(V)/V ratio box
        self.label_desired_ratio = qtw.QLabel(text='Delta (V/V ratio)', maximumWidth=200)
        self.line_desired_ratio  = qtw.QLineEdit(text='0.05',           maximumWidth=100)
        # Noise or Time Combo Box
        self.label_Time_Noise    = qtw.QLabel(text='Ends when?',        maximumWidth=100)
        self.combo_Time_Noise    = qtw.QComboBox(maximumWidth=150)
        self.combo_Time_Noise.addItem("Time",                  userData=qtc.QVariant(0))
        self.combo_Time_Noise.addItem("Noise",                 userData=qtc.QVariant(1))
        self.combo_Time_Noise.addItem("Whichever Comes First", userData=qtc.QVariant(2))
        # Time Input Box
        self.label_time = qtw.QLabel(text='Time per Measurement (s)', maximumWidth=200)
        self.line_time  = qtw.QLineEdit(text='30',                    maximumWidth=100)
        # Target Noise Input Box
        self.label_target_noise    = qtw.QLabel(text='Target Noise',  maximumWidth=100)
        self.line_target_noise     = qtw.QLineEdit(text='2e-3',       maximumWidth=100)
        # Trigger Channel Combo Box
        self.label_trig_channel    = qtw.QLabel(text='Trigger Channel', maximumWidth=100)
        self.combo_trigger_channel = qtw.QComboBox(maximumWidth=100)
        self.combo_trigger_channel.addItem('CH A', userData=qtc.QVariant(0))
        self.combo_trigger_channel.addItem('CH B', userData=qtc.QVariant(1))
        self.combo_trigger_channel.addItem('EXT',  userData=qtc.QVariant(4))
        self.combo_trigger_channel.setCurrentIndex(2)  # sets default to 'EXT'
        # Trigger Direction Combo Box
        self.label_trigger_direction = qtw.QLabel(text='Trigger Direction', maximumWidth=100)
        self.combo_trigger_direction = qtw.QComboBox(maximumWidth=100)
        self.combo_trigger_direction.addItem('rising',  userData=qtc.QVariant(0))
        self.combo_trigger_direction.addItem('falling', userData=qtc.QVariant(1))
        # Trigger Threshold Input Box
        self.label_trigger_threshold = qtw.QLabel(text='Trig Thresh.(mV)',  maximumWidth=100)
        self.line_trigger_threshold  = qtw.QLineEdit(text='3000',           maximumWidth=100)

        ##### Place GUI Elements on the Left Side #####
        # adds the list of widgets utilizing a loop (was originally adding elements explicitely which caused lots of repitition)
        widget_list = [self.textEdit,                self.min_LED_current_label,   self.min_LED_current,    self.max_LED_current_label, self.max_LED_current,
                       self.steps_LED_current_label, self.steps_LED_current,       self.label_LED_spacing,  self.combo_LED_spacing,     self.label_V,
                       self.combo_V,  self.label_dt, self.line_dt,                 self.label_nsamp,        self.line_nsamp,            self.label_nseg,
                       self.line_nseg,               self.label_desired_ratio,     self.line_desired_ratio, self.label_Time_Noise,      self.combo_Time_Noise,
                       self.label_time,              self.line_time,               self.label_target_noise, self.line_target_noise,     self.label_trig_channel,
                       self.combo_trigger_channel,   self.label_trigger_direction, self.combo_trigger_direction, self.label_trigger_threshold, self.line_trigger_threshold]
        for i in widget_list:
            layout_cmb.addWidget(i)

        layout.addLayout(layout_cmb)

        ##### Define and Place Graphical Elements #####
        layout_graph1 = qtw.QVBoxLayout()
        layout_graph2 = qtw.QVBoxLayout()
        
        self.left_graph   = MatplotlibCanvas(self)
        self.left_toolbar = NavigationToolbar(self.left_graph, self)
        layout_graph1.addWidget(self.left_toolbar)
        layout_graph1.addWidget(self.left_graph)
        
        self.right_graph   = MatplotlibCanvas(self)
        self.right_toolbar = NavigationToolbar(self.right_graph, self)
        self.right_graph.axes.set_ylim(-0.2,1.5)
        self.right_graph.axes.set_ylabel("Voltage (V)")
        layout_graph2.addWidget(self.right_toolbar)
        layout_graph2.addWidget(self.right_graph)
    
        layout.addLayout(layout_graph1)
        layout.addLayout(layout_graph2)
        
        ##### Define GUI Elements on the Right Side #####
        self.psh_open_devices  = qtw.QPushButton(text='Open Devices',          maximumWidth=150)
        self.psh_voc_scan      = qtw.QPushButton(text='Scan Voc Range',        maximumWidth=150, toolTip='Check the Voc across the bias light range in order to determine the range that is worth measuring.')
        self.psh_find_trig_lev = qtw.QPushButton(text='Find Trigger Level',    maximumWidth=150)
        self.psh_close_devices = qtw.QPushButton(text='Close Devices',         maximumWidth=150)
        self.psh_pp         = qtw.QPushButton(  text='Pump-Probe 1',           maximumWidth=150, toolTip='TPV/TPC w/o the ability to control the delta(V)/V ratio.')
        self.psh_pp_filter  = qtw.QPushButton(  text='Pump-Probe 1 w/ Filter', maximumWidth=150, toolTip='TPV w/ the ability to control the delta(V)/V ratio.')
        self.psh_pp_filter_TPC_check = qtw.QCheckBox(text='Retain Filter Pos for TPC?',          toolTip='Saves the positions of the filter in millimeters to have identical conditions for TPC measurments.')
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
                        self.psh_pp,           self.psh_pp_filter, self.psh_pp_filter_TPC_check,
                        self.psh_pp2,          self.psh_pp3,       self.check_pumponly,    self.check_TATPC,
                        self.check_chref,      self.psh_pp_fr,     self.psh_stop,          self.psh_abort]
        for i in widget_list2: layout_psh.addWidget(i)

        layout.addLayout(layout_psh)
        self.widget.setLayout(layout)
        self.setCentralWidget(self.widget)
        self.setWindowTitle("Fitting progress")

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
        self.psh_open_devices.clicked.connect(self.open_devices)
        self.psh_voc_scan.clicked.connect(self.Voc_Scan)
        self.psh_close_devices.clicked.connect(self.close_devices)
        self.psh_find_trig_lev.clicked.connect(self.find_trigger_level)
        self.psh_pp.clicked.connect(       lambda mode="measure": self.pump_probe_one(mode))
        self.psh_pp_filter.clicked.connect(lambda mode='measure': self.pump_probe_one_filter(mode))
        self.psh_pp2.clicked.connect(      lambda mode='measure': self.pump_probe_two(mode))
        self.psh_pp3.clicked.connect(      lambda mode='measure': self.pump_probe_three(mode))
        self.psh_pp_fr.clicked.connect(    lambda mode='freerun': self.pump_probe_freerun(mode)) # add lambda fct to indicate that this is free-run
        self.psh_stop.clicked.connect(self.stop)
        self.psh_abort.clicked.connect(self.abort)
        self.line_trigger_threshold.returnPressed.connect(self.refresh_trig)
        self.combo_trigger_direction.currentIndexChanged.connect(self.refresh_trig)

    def find_trigger_level(self):
        ''' stupid by foot routine to do what the operator would do otherwise:
        scan the trig level and find the range in which it fires
        thb set the level to the middle of this range and hope this is highest slope
        (cannot detect slope because it is not a channel)
        '''

        print('setting trigger level...')

        threshs = np.logspace(0, 4, 200)

        self.setup_experiment()

        self.stop_request = False
        self.read_user_settings()
        self.ps.nseg = 2
        self.ps.sel_dt    = 1e-7
        self.ps.sel_nsamp = 1000
        self.ps.sel_autotrig_ms = 100
        self.ps.sel_threshold   = 1
        self.ps.set_channels()
        self.ps.get_timebase()
        self.ps.set_trigger()
        if self.ps.sel_rapid_block:
            self.ps.set_memory_segments()
            self.ps.set_no_of_captures()

        trigmin = 0
        trigmax = 32767
        for n in threshs:
            self.ps.sel_threshold = n
            self.ps.set_trigger()
            tic = time.time()
            self.ps.flag_newdata = False
            print('before run block')
            self.ps.run_block()
            print('after run block')
            qtw.QApplication.processEvents()
            while not self.ps.flag_newdata:  # wait until driver supplies data
                time.sleep(0.01)
            self.ps.flag_newdata = False
            dt = time.time() - tic
            print('have been waiting for trigger for ' + str(dt) + ' s.')
            if (dt < self.ps.sel_autotrig_ms / 1000.) and (trigmin == 0):
                trigmin = n  # lower trig level has been found
            if (dt > self.ps.sel_autotrig_ms / 1000.) and (trigmin > 0):
                trigmax = n
                break
        triglev = 0.5 * (trigmin + trigmax)
        print(f'setting trigger level  to {triglev} which is average between {trigmin} and {trigmax}.')
        #self.ps.sel_threshold = triglev
        self.ps.set_sel_threshold(triglev)
        self.line_trigger_threshold.setText(str(int(np.rint(triglev))))

    def refresh_trig(self):
        print('calling refresh trig')
        self.ps.sel_threshold = int(str(self.line_trigger_threshold.text()))
        self.ps.sel_direction = self.combo_trigger_direction.currentIndex() + 2
        self.ps.set_trigger()

    def read_user_settings(self):
        print('Reading User Settings...')
        # Voltage Range
        self.ps.sel_v  = self.combo_V.currentIndex()
        # Not sure what this is
        self.ps.sel_dt = float(str(self.line_dt.text()))
        # Number of Samples
        self.ps.sel_nsamp = int(str(self.line_nsamp.text()))
        # Number of Segments
        self.ps.nseg = int(str(self.line_nseg.text()))
        # Target Noise
        self.target_noise = float(str(self.line_target_noise.text()))
        # Pump Only? true/false checkbox
        self.measure_pumponly = self.check_pumponly.isChecked()
        # TA/TPC checkbox
        self.flag_TATPC = self.check_TATPC.isChecked()
        # Trigger Channel
        self.trigger_channel = self.combo_trigger_channel.currentIndex()
        if self.trigger_channel != 2:
            self.ps.sel_trig_ch = self.trigger_channel
        else:
            self.ps.sel_trig_ch = 4
        # Trigger Threshold Value
        self.ps.sel_threshold = int(str(self.line_trigger_threshold.text()))
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

        print("User Settings Read Successfully")

    def open_devices(self):
        ''' Function for Open Devices button. '''
        print("Trying to Open Picoscope...")
        try:
            self.ps = picoscope()
            self.ps.open_device()
        except:
            print("Connection to Picoscope was unsuccessful")
        else:
            print("Connection to Picoscope was successful")
        
        print("Trying to open the LED connection...")
        try:
            self.LED_instr, self.rm = LED_open_connection()
        except:
            print("Connection to LED was unsuccessful")
        else:
            print("Connection to LED was successful")
        
        print("Trying to open connection to motor...")
        try:
            self.motor, self.serial_num = MotorControl.load_dll()
            MotorControl.open_motor_connection(self.motor,self.serial_num)
        except:
            print("Connection to motor was unsuccessful")
        else:
            print("Connection to motor was successful. Please allow the motor to home itself.")

    def stop(self):
        ''' Function for Stop button. '''
        self.stop_request = True

    def abort(self):
        ''' Function for Abort button. '''
        self.stop_request = True
        self.abort_request = True

    def close_devices(self):
        ''' Function for Close Devices button '''
        self.ps.close_device()
        LED_close_connection(self.LED_instr, self.rm)
        MotorControl.close_motor_connection(self.motor,self.serial_num)
        print("All device connections closed.")

    def setup_experiment(self):
        self.fig2_yscale = [-self.ps.range_fctr[self.ps.sel_v],
                            self.ps.range_fctr[self.ps.sel_v]]
        self.stop_request = False
        self.abort_request = False
        self.read_user_settings()
        self.prepare_plots()
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

    def Voc_Scan(self):
        """
        This function sweeps the whole range of the bias light in order to gather the Voc (or Jsc).
        The user then determines the appropriate bias light range to use for `pump_probe`
        """

        self.setup_experiment()
        
        # to save the resulting table
        file_name = qtw.QFileDialog.getSaveFileName(self, 'Save File', 'D:\\', '*.csv')
        file_name = str(file_name[0])
        print(file_name)
        
        # Define the values for the LED
        self.LED_current_list = generate_current_list(self.min_current, self.max_current, self.steps_LED, self.LED_spacing)  # in Amps
        
        # list for the Voc values
        Vocs = [np.nan] * len(self.LED_current_list)

        self.right_graph.axes.cla()
        self.right_graph.axes.set_title("$V_{OC}$ Values")
        self.right_graph.axes.set_xscale('log')
        self.right_graph.axes.set_ylabel('$V_{OC}$')
        self.right_graph.axes.set_xlabel('Bias Light Current (A)')
        self.right_graph.axes.set_xlim(self.min_current/2000,self.max_current/500) 

        for id, LED_power in enumerate(np.flip(self.LED_current_list)):
            self.left_graph.axes.set_title('$V_{OC}$ for LED set to ' + f'{round(LED_power,4)}')
            
            # sets the LED current
            self.LED_instr.write("SOURCE1:CCURENT:CURRENT " + str(LED_power)) # type: ignore
            # turn on LED and wait a short time
            self.LED_instr.write("OUTPUT1:STATE ON") # type: ignore
            print("LED on")
            time.sleep(0.5)
            # output applied current
            print("Applied LED current (Amps):" + self.LED_instr.query("SENSE3:CURRENT:DATA?")) # type: ignore
            
            # starting time
            start_time = time.time()
            elapsed_time = 0
            
            # initial reading
            _, CB = self.pump_probe_core(20)
            self.left_graph.axes.cla()
            self.left_graph.axes.set_ylim(np.min(CB)-0.05, np.max(CB)+0.05)
            self.left_graph.axes.set_ylabel('$V_{OC}$ (V)')
            self.left_graph.axes.set_xlabel('Sample')
            left_line = self.left_graph.axes.plot(CB,color='tab:blue')
            
            Voc = self.find_voc_or_jsc(left_line,start_time,elapsed_time)
            
            print(f'The Voc for {round(LED_power,4)} Amps is {round(Voc,4)} Volts.')
            Vocs[-(id+1)] = Voc
            
            # turn off LED and wait a short time
            self.LED_instr.write("OUTPUT1:STATE OFF") # type: ignore
            print("LED off")

            # update the graph to show a scatter of all Voc values on the right
            self.right_graph.axes.set_ylim(1e-2,np.nanmax(np.array(Vocs))+0.1)
            for artist in self.right_graph.axes.lines:
                artist.remove()
            self.right_graph.axes.scatter(self.LED_current_list, Vocs,color="tab:blue")
            self.right_graph.draw()
            qtw.QApplication.processEvents()
            if self.stop_request:
                return None
            
            self.wait_for_cooling(LED_power,30)
            
        print("All Voc checks complete")

        # a dialog window pops up prompting user to place 50 ohm resistor in line if they want to also measure Jsc
        Jsc_dialog_box = qtw.QMessageBox(self)
        Jsc_dialog_box.setWindowTitle("TPV to TPC Switch")
        Jsc_dialog_box.setIcon(qtw.QMessageBox.Icon.Question)
        Jsc_dialog_box.setStandardButtons(qtw.QMessageBox.StandardButton.Yes | qtw.QMessageBox.StandardButton.No)
        Jsc_dialog_box.setText("If you would like to measure Jsc as well, place the 50 ohm resistor inline and press 'Yes'. \n Otherwise, press 'No'.")
        
        button_clicked = Jsc_dialog_box.exec()

        if button_clicked == qtw.QMessageBox.StandardButton.Yes:
            self.setup_experiment()
            
            # list for the Voc (or Jsc) values
            Jscs = [np.nan] * len(self.LED_current_list)
            self.right_graph.axes.cla()
            self.right_graph.axes.set_title("$J_{SC}$ Values")
            self.right_graph.axes.set_xscale('log')
            self.right_graph.axes.set_ylabel('$J_{SC}$')
            self.right_graph.axes.set_xlabel('Bias Light Current (A)')
            self.right_graph.axes.set_xlim(self.min_current/2000,self.max_current/500) 

            for id, LED_power in enumerate(np.flip(self.LED_current_list)):
                self.left_graph.axes.set_title('$V_{OC}$ for LED set to ' + f'{round(LED_power,4)}')
                
                # sets the LED current
                self.LED_instr.write("SOURCE1:CCURENT:CURRENT " + str(LED_power)) # type: ignore
                # turn on LED and wait a short time
                self.LED_instr.write("OUTPUT1:STATE ON") # type: ignore
                print("LED on")
                time.sleep(0.5)
                # output applied current
                print("Applied LED current (Amps):" + self.LED_instr.query("SENSE3:CURRENT:DATA?")) # type: ignore
                
                # starting time
                start_time = time.time()
                elapsed_time = 0
                
                # initial reading
                _, CB = self.pump_probe_core(20)
                self.left_graph.axes.cla()
                self.left_graph.axes.set_ylim(np.min(CB)-0.01, np.max(CB)+0.01)
                self.left_graph.axes.set_ylabel('$J_{SC}$ (V)')
                self.left_graph.axes.set_xlabel('Sample')
                left_line = self.left_graph.axes.plot(CB,color='tab:blue')
                
                Jsc = self.find_voc_or_jsc(left_line,start_time,elapsed_time)
                
                print(f'The Jsc for {round(LED_power,4)} Amps is {round(Jsc,4)} Volts.')
                Jscs[-(id+1)] = Jsc
                
                # turn off LED and wait a short time
                self.LED_instr.write("OUTPUT1:STATE OFF") # type: ignore
                print("LED off")

                # update the graph to show a scatter of all Voc values on the right
                self.right_graph.axes.set_ylim(-0.1,np.nanmax(np.array(Jscs))+0.1)
                for artist in self.right_graph.axes.lines:
                    artist.remove()
                self.right_graph.axes.scatter(self.LED_current_list, Jscs,color="tab:blue")
                self.right_graph.draw()
                qtw.QApplication.processEvents()
                if self.stop_request:
                    return None
                
                self.wait_for_cooling(LED_power,30)
                
            print("All Voc checks complete")
            # save just the Voc data
            data = pd.DataFrame({'Bias (A)': self.LED_current_list, 'Voc (V)': Vocs, 'Jsc (V)': Jscs})
            data.to_csv(file_name,index=False)
        else:
            # save just the Voc data
            data = pd.DataFrame({'Bias': self.LED_current_list, 'Voc': Vocs})
            data.to_csv(file_name,index=False)

    def wait_for_cooling(self, LED_power,wait_time):
        if LED_power != self.LED_current_list[0]:
            # if not the final light intensity, sleep
            print(f"Waiting for the sample to cool for {wait_time} seconds.")
            delay_begin_time = time.time()
            delay_end_time = delay_begin_time + wait_time
            while time.time() < delay_end_time:
                continue
            print("Waiting complete.")

    def find_voc_or_jsc(self,left_line,start_time,elapsed_time) -> float:
        ii=1
        
        self.left_graph.draw()
        qtw.QApplication.processEvents()
        while (elapsed_time<=5):
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
        
    def pump_probe_three(self, mode='measure'):
        '''wrap around pump probe 
        so pump probe can be written such that it can be called with different time windows
        to be stitched together
        '''

        file_name = qtw.QFileDialog.getSaveFileName(self, 'Save File', self.pfad_open, '*.dat')
        file_name = str(file_name)
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

        self.line_dt.setText(str(self.ps.sel_dt * 10.))
        self.setup_experiment()
        inz = int(self.ps.sel_nsamp * self.ps.sel_trig_pos_rel * 0.8)  # up to here, average is calculated for baseline
        D = self.pump_probe(inz, mode='measure')  # call the main loop
        if self.abort_request: return

        DTA2 = D[0]
        DTD2 = D[1]
        nza2 = D[2]
        nzb2 = D[3]
        t2 = self.t * 1.0

        self.line_dt.setText(str(self.ps.sel_dt * 10.))
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

            self.line_dt.setText(str(self.ps.sel_dt * 0.01))
            self.setup_experiment()
            D = self.pump_only(inz, nza1, nzb1, DTD1, mode='measure')  # call the main loop
            if self.abort_request: return

            DTAx1 = D[0]
            DTDx1 = D[1]
            Apu1 = D[2]
            Bpu1 = D[3]

            self.line_dt.setText(str(self.ps.sel_dt * 10.0))
            self.setup_experiment()
            D = self.pump_only(inz, nza2, nzb2, DTD2, mode='measure')  # call the main loop
            if self.abort_request: return

            DTAx2 = D[0]
            DTDx2 = D[1]
            Apu2 = D[2]
            Bpu2 = D[3]

            self.line_dt.setText(str(self.ps.sel_dt * 10.0))
            self.setup_experiment()
            D = self.pump_only(inz, nza3, nzb3, DTD3, mode='measure')  # call the main loop
            if self.abort_request: return

            DTAx3 = D[0]
            DTDx3 = D[1]
            Apu3 = D[2]
            Bpu3 = D[3]

            self.line_dt.setText(str(self.ps.sel_dt * 0.01))
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

        print(clip2)
        print(clip3)

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
        print(raus.shape)
        print(tall.shape)
        save_matrix_milano_convention(range(raus.shape[1]), tall, raus, file_name, suff='')

        print(' measurement finished.')

    def pump_probe_two(self, mode='measure'):
        '''wrap around pump probe 
        so pump probe can be written such that it can be called with different time windows
        to be stitched together
        '''

        file_name = qtw.QFileDialog.getSaveFileName( self, 'Save File', self.pfad_open, '*.dat')
        file_name = str(file_name)
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

        self.line_dt.setText(str(self.ps.sel_dt * 10.))
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

        self.line_dt.setText(str(self.ps.sel_dt * 0.1))
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

            self.line_dt.setText(str(self.ps.sel_dt * 10.0))
            self.line_nsamp.setText(str(self.ps.sel_nsamp * 10))
            self.setup_experiment()
            inz = int(self.ps.sel_nsamp * self.ps.sel_trig_pos_rel * 0.8)  # up to here, average is calculated for baseline
            D = self.pump_only(inz, nza2, nzb2, DTD2, mode='measure')  # call the main loop
            if self.abort_request: return

            DTAx2 = D[0]
            DTDx2 = D[1]
            Apu2 = D[2]
            Bpu2 = D[3]

            self.line_dt.setText(str(self.ps.sel_dt * 0.1))
            self.line_nsamp.setText(str(int(self.ps.sel_nsamp * 0.1)))

        else:
            DTAx1 = DTA1
            DTDx1 = DTD1
            DTAx2 = DTA2
            DTDx2 = DTD2

            Apu1 = Bpu1 = DTAx1 * 0.0
            Apu2 = Bpu2 = DTAx2 * 0.0

        clip2 = (t2 > t1[-1])

        print(clip2)

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
        print(raus.shape)
        print(tall.shape)
        save_matrix_milano_convention(range(raus.shape[1]), tall, raus, file_name, suff='')

        print(' measurement finished.')

    def pump_probe_freerun(self, mode='freerun'):
        self.mode = mode
        self.setup_experiment()
        inz = int(self.ps.sel_nsamp * self.ps.sel_trig_pos_rel * 0.8)  # up to here, average is calculated for baseline
        # run a block to obtain array size
        D = self.pump_probe(inz, mode='freerun')  # call the main loop

    def pump_probe_one(self, mode='measure'):
        ''' wrap around pump probe 
        so pump probe can be written such that it can be called with different time windows
        to be stitched together (see pump_probe_three)
        '''

        # grab current time to see how long the measurements take
        start_time = time.time()

        self.doing_ratio_checking = False

        file_name = qtw.QFileDialog.getSaveFileName(self, 'Save File', self.pfad_open, '*.dat')
        file_name = str(file_name)[2:-11]

        self.mode = mode
        self.setup_experiment()
        # inz = number of samples * relative position of trigger in collected data * 0.8. Not Sure the purpose of this
        inz = int(self.ps.sel_nsamp * self.ps.sel_trig_pos_rel * 0.8)  # up to here, average is calculated for baseline
        # run a block to obtain array size

        # Define the values for the LED
        self.LED_current_list = generate_current_list(self.min_current, self.max_current, self.steps_LED, self.LED_spacing)  # in Amps
        # this for loop was added by Steven
        for LED_power in np.flip(self.LED_current_list):
            # sets the LED current
            self.LED_instr.write("SOURCE1:CCURENT:CURRENT " + str(LED_power)) # type: ignore
            # turn on LED and wait a short time
            self.LED_instr.write("OUTPUT1:STATE ON") # type: ignore
            print("LED on")
            time.sleep(0.5)
            # output applied current
            print("Applied LED current (Amps):" + self.LED_instr.query("SENSE3:CURRENT:DATA?")) # type: ignore

            # calls the main loop to take measurement
            D = self.pump_probe(inz, mode)

            if self.stop_request:
                continue
            if self.abort_request:
                break

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
            print(raus.shape)

            # save the data
            bias_power_str = f"{LED_power:07.4f}"
            save_matrix_milano_convention(range(raus.shape[1]), self.t, raus, file_name, suff=bias_power_str)  # added suffix parameter
            print('Measurement with current ' + bias_power_str + ' Amps is finished and saved.')

            # turn off LED and wait a short time
            self.LED_instr.write("OUTPUT1:STATE OFF") # type: ignore
            print("LED off")
            self.wait_for_cooling(LED_power,30)
            
        # end of for loop
        print("All measurements finished.")
        # measures and prints the time elapsed
        finish_time = time.time()
        measurement_time = finish_time - start_time
        print("Total Time: " + str(measurement_time))

    def pump_probe_one_filter(self, mode='measure') -> None:
        """
        This function does TPV measurements over a series of bias light intensities while maintaining a ratio of delta(V)/V.
        If desired, the function prompts the user to change the setup to do TPC measurements immediately after TPV.
        """

        # grab current time to see how long the measurements take
        start_time = time.time()

        # initialize list for filter positions
        filter_positions = []

        self.early_return = False

        file_name = qtw.QFileDialog.getSaveFileName(self, 'Save File', self.pfad_open, '*.dat')
        file_name = str(file_name)[2:-11]

        self.mode = mode
        self.setup_experiment()
        # inz = number of samples * relative position of trigger in collected data * 0.8. Not Sure the purpose of this
        inz = int(self.ps.sel_nsamp * self.ps.sel_trig_pos_rel * 0.8)  # up to here, average is calculated for baseline
        # run a block to obtain array size

        # Define the values for the LED
        self.LED_current_list = generate_current_list(self.min_current, self.max_current, self.steps_LED, self.LED_spacing)  # in Amps
        # this for loop was added by Steven
        for LED_power in np.flip(self.LED_current_list):
            self.doing_ratio_checking = True
            final_pass = False

            # sets the LED current
            self.LED_instr.write("SOURCE1:CCURENT:CURRENT " + str(LED_power)) # type: ignore
            # turn on LED and wait a short time
            self.LED_instr.write("OUTPUT1:STATE ON") # type: ignore
            print("LED on")
            time.sleep(0.5)
            # output applied current
            print("Applied LED current (Amps):" + self.LED_instr.query("SENSE3:CURRENT:DATA?")) # type: ignore

            # while loop to check delta(V)/V ratio
            ratio_satisfied = False # initialize break condition for when ratio is met
            while True:
                '''
                Repeats the speicific light intensity until the delta(V)/Voc ratio is met.
                '''

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
                print(raus.shape)

                if final_pass:
                    break
                ratio_satisfied = self.filter_laser_power(DTDx) 
                if ratio_satisfied:
                    self.doing_ratio_checking = False
                    final_pass = True
                
            current_motor_position_mm = MotorControl.get_position(self.motor, self.serial_num)
            filter_positions.append(current_motor_position_mm)

            # save the data
            bias_power_str = f"{LED_power:07.4f}"
            save_matrix_milano_convention(range(raus.shape[1]), self.t, raus, file_name, suff=bias_power_str)  # added suffix parameter
            print('Measurement with current ' + bias_power_str + ' Amps is finished and saved.')

            # turn off LED and wait a short time
            self.LED_instr.write("OUTPUT1:STATE OFF") # type: ignore
            print("LED off")

            if self.early_return==False:
                # move motor small amount in preperation for next light intensity
                current_motor_position_mm = MotorControl.get_position(self.motor, self.serial_num)
                MotorControl.increment_to_new_position(self.motor,self.serial_num,current_motor_position_mm,-1.0)
            else: self.early_return = False # re-arms the boolean variable and moves on. will only be reached when early_return is true
            
            self.wait_for_cooling(LED_power,30)

        # end of for loop
        print("All measurements finished.")
        # measures and prints the time elapsed
        finish_time = time.time()
        measurement_time = finish_time - start_time
        print("Total Time: " + str(measurement_time))

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
        index_peak = np.argmax(voltage_values)
        max_voltage = np.amax(voltage_values)
        
        # take half of the data before that peak and take the median (or mean)
        baseline_voltage = np.mean(voltage_values[:int(index_peak/2)])
        # find the ratio between change in voltage from baseline to peak.
        actual_ratio = (max_voltage - baseline_voltage)/baseline_voltage
        # if ratio is within 10% (for example) of desired ratio then move on to save the data
        minimum_acceptable_ratio = self.voltage_ratio_value-(0.20*self.voltage_ratio_value)
        maximum_acceptable_ratio = self.voltage_ratio_value+(0.20*self.voltage_ratio_value)
        
        print('Baseline Voltage: ' + str(baseline_voltage))
        print('Max Voltage: ' + str(max_voltage))
        print('Acceptable Ratio Range: (' + str(minimum_acceptable_ratio) + ', ' + str(maximum_acceptable_ratio) + ')')
        print('Actual Ratio: ' + str(actual_ratio))
        
        if (minimum_acceptable_ratio < actual_ratio) and (actual_ratio < maximum_acceptable_ratio):
            # if voltage ratio is in the "acceptable region" then end the while loop and move on to saving.
            print("Actual Ratio: " + str(round(actual_ratio,5)) + " is within 10% of the desired ratio of " + str(self.voltage_ratio_value))
            return True # ends the while loop
            
        # get the current position of the motor
        current_motor_position_mm = MotorControl.get_position(self.motor, self.serial_num)
        
        if (actual_ratio < self.voltage_ratio_value):
            # if the current ratio is less than the desired ratio
            # move a small amount to block less laser light
            if (current_motor_position_mm > 53):
                print("Desired ratio is unachievable. Motor is at its extreme position.")
                self.early_return = True
                return True
            
            # when the ratio is far away from the desired ratio, move a larger amount
            increment_mm = 0.5
            wait = 0.75
            if (abs(self.voltage_ratio_value-actual_ratio) >= 0.005):
                increment_mm = 1.0
            if (abs(self.voltage_ratio_value-actual_ratio) >= 0.01):
                increment_mm = 2.0 
                wait = 1.0
            if (abs(self.voltage_ratio_value-actual_ratio) >= 0.02):
                increment_mm = 3.0
                wait = 1.25
            if (abs(self.voltage_ratio_value-actual_ratio) >= 0.03):
                increment_mm = 4.0
                wait = 1.50

            MotorControl.increment_to_new_position(self.motor,self.serial_num,current_motor_position_mm,increment_mm)
            time.sleep(wait)
            return False
        
        if (actual_ratio > self.voltage_ratio_value):
            # if the current ratio is larger than the desired ratio
            if (current_motor_position_mm < 0.05):
                # but do not permit the motor to try going negative relative to its home position
                print("Desired ratio is unachievable. Motor is at its extreme position.")
                self.early_return = True
                return True
            
            increment_mm = -0.5
            wait = 0.75
            if (abs(self.voltage_ratio_value-actual_ratio) >= 0.005):
                increment_mm = -1.0
            if (abs(self.voltage_ratio_value-actual_ratio) >= 0.01):
                increment_mm = -2.0 
                wait = 1.0
            if (abs(self.voltage_ratio_value-actual_ratio) >= 0.02):
                increment_mm = -3.0
                wait = 1.25
            if (abs(self.voltage_ratio_value-actual_ratio) >= 0.03):
                increment_mm = -4.0
                wait = 1.50

            # move a small amount to block more laser light
            #increment_mm = -0.5
            MotorControl.increment_to_new_position(self.motor,self.serial_num,current_motor_position_mm,increment_mm)
            time.sleep(1)
            return False
    
    def pump_probe_one_filter_TPC(self, mode="measure") -> None:
        """
        This function is called at the end of `pump_probe_one_filter` when `retain_filter_positions` is checked as `True`.
        """
        # grab current time to see how long the measurements take
        start_time = time.time()

        # import list of filter positions and reverses it
        filter_positions_mm: list[float] = np.loadtxt("filter_positions.csv").tolist()

        file_name = qtw.QFileDialog.getSaveFileName(self, 'Save File', self.pfad_open, '*.dat')
        file_name = str(file_name)[2:-11]

        self.mode = mode
        self.setup_experiment()
        # inz = number of samples * relative position of trigger in collected data * 0.8. Not Sure the purpose of this
        inz = int(self.ps.sel_nsamp * self.ps.sel_trig_pos_rel * 0.8)  # up to here, average is calculated for baseline
        # run a block to obtain array size

        # Define the values for the LED
        self.LED_current_list = generate_current_list(self.min_current, self.max_current, self.steps_LED, self.LED_spacing)  # in Amps

        # this for loop was added by Steven
        for id, LED_power in enumerate(np.flip(self.LED_current_list)):
            
            # move motor to correct position
            current_position = MotorControl.get_position(self.motor,self.serial_num)
            if type(filter_positions_mm)==float:
                desired_position = filter_positions_mm
            else:
                desired_position = filter_positions_mm[id]
            amount_to_move_filter = desired_position - current_position

            MotorControl.increment_to_new_position(self.motor,self.serial_num,current_position,amount_to_move_filter)
            
            wait_time = 30
            print(f"Waiting for the sample to cool for {wait_time} seconds.")
            delay_begin_time = time.time()
            delay_end_time = delay_begin_time + wait_time
            while time.time() < delay_end_time:
                continue
            print("Waiting complete.")

            # sets the LED current
            self.LED_instr.write("SOURCE1:CCURENT:CURRENT " + str(LED_power)) # type: ignore
            # turn on LED and wait a short time
            self.LED_instr.write("OUTPUT1:STATE ON") # type: ignore
            print("LED on")
            time.sleep(0.5)
            # output applied current
            print("Applied LED current (Amps):" + self.LED_instr.query("SENSE3:CURRENT:DATA?")) # type: ignore

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
            print(raus.shape)

            # save the data
            bias_power_str = f"{LED_power:07.4f}"
            save_matrix_milano_convention(range(raus.shape[1]), self.t, raus, file_name, suff=bias_power_str)  # added suffix parameter
            print('Measurement with current ' + bias_power_str + ' Amps is finished and saved.')

            # turn off LED and wait a short time
            self.LED_instr.write("OUTPUT1:STATE OFF") # type: ignore
            print("LED off")
            time.sleep(0.75)

        # end of for loop
        print("All measurements finished.")
        # measures and prints the time elapsed
        finish_time = time.time()
        measurement_time = finish_time - start_time
        print("Total Time: " + str(measurement_time))

    def pump_probe(self, inz: int, mode='measure') -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        '''
        the core loop
        mode = measure saves and runs infinitely
        mode = freerun runs looped
        '''
        # THE SAMPLE LOOP WITH PUMP AND PROBE
        ii = 0  # the total number of exposures
        noise = 1000  # initial noise level
        elapsed_time = 0  # sets the duration to 0
        peak_CB = 0 # sets initial peak to 0

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
                peak_CB = CB.argmax()
            else:  # all following iterations
                # recalculates running average, accounting for current iteration
                # ((running average * iteration number) + current values) / (iteration number + 1->[bc it starts at 0])
                CA = (CA * ii + D[0]) / (ii + 1)
                
                # peak of the current measurement
                instance_peak = D[1].argmax()
                difference_in_peaks = peak_CB - instance_peak
                CB_old = scipy.ndimage.shift(D[1],difference_in_peaks,mode='nearest')
                CB = (CB * ii + CB_old) / (ii + 1)
                #CB = (CB * ii + D[1]) / (ii + 1)

            # finds the average of CA from beginning to inz (what is inz?)
            nza = np.mean(CA[:inz])  # if pumponly, use the old self.nza
            DTA = (CA - nza)/nza

            if self.ps.sel_chon[1]:
                '''if using both channels, also repeat for channel B'''
                nzb = np.mean(CB[:inz])
                DTB = (CB - nzb)/nzb
                DTD = DTA - DTB
            else:  # zero out all arrays
                nzb = 0.0 * nza
                DTB = 0.0 * DTA
                DTD = 0.0 * DTA

            print('Current exposure:' + str(ii * 20 * self.ps.nseg))
            ii += 1

            if self.ps.sel_chon[1]:
                #noiseA = np.std(CA[:inz])
                #noiseB = np.std(CB[:inz])
                noise = np.std(DTD[:inz])
            else:
                #noiseA = np.std(CA[:inz])
                noise = np.std(DTA[:inz])
            print('noise of DT/T trace= ' + str(noise))
            # print 'absolute noise of CHA: ' + str(noiseA) + ' V'
            # print 'absolute noise of CHB: ' + str(noiseB) + ' V'
            print('absolute peak value:' + str(np.amax(D[0])-np.mean(D[0][:inz])) + 'V')
            if self.ps.sel_chon[1]:
                finalDT = np.mean(DTD[2*inz:])
                print('final DT/T:' + str(finalDT))

            # self.plot_waveform(DTA, DTB)
            if self.flag_TATPC:
                self.plot_pp(CB)
            else:
                if self.ps.sel_chon[1]:
                    self.plot_pp(DTD)
                else:
                    self.plot_pp(DTA)

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
            # end of while loop

        if self.flag_TATPC:
            return DTA, CB, nza, nzb  # if user desires TA on CHA and TPC on CHB
        else:
            return DTA, DTD, nza, nzb

    def pump_probe_core(self, ns) -> tuple[np.ndarray, np.ndarray]:
        # the core sampling routine
        # ns number of run_block calls

        if self.ps.sel_chon[1]:
            print('CHB is ON')
        else:
            print('CHB is OFF')

        for kk in range(ns):
            self.ps.run_block()
            qtw.QApplication.processEvents()
            while not self.ps.flag_newdata:  # wait until driver supplies data
                time.sleep(0.01)
            self.ps.flag_newdata = False
            
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
                CA = (CA * kk + ADCA) / (kk + 1)  # the average voltage signal
                if self.ps.sel_chon[1]:
                    CB = (CB * kk + ADCB) / (kk + 1) # the average voltage signal
                else:
                    CB = 0.0*CA

        # print 'sum of CB:' + str(sum(CB))
        return CA, CB
        # the core loop

    def check_end_condition(self, start_time: float, noise) -> bool:
        end_condition = False
        match self.time_noise:
            case 0: end_condition = self.end_time(start_time)
            case 1: end_condition = self.end_noise(noise)
            case 2: end_condition = (self.end_time(start_time) or self.end_noise(noise))
        return end_condition
    
    def end_time(self, start_time: float) -> bool:
        elapsed_time = time.time() - start_time
        print("Elapsed Time:" + str(elapsed_time))
        return (elapsed_time > self.time_limit_value)
        
    def end_noise(self, noise) -> bool:
        return (noise < self.target_noise)
    
    def pump_only(self, inz, nza0, nzb0, DTD, mode='measure'):
        '''the core loop
        mode = measure saves and runs infinitely
        mode = freerun runs looped
        mode = pumponly assumes pump is off, and subtracts reference data from very last pump-probe measurement
        '''

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

            print(nza0)
            print(nzb0)

            print('Current exposure:' + str(ii * 20 * self.ps.nseg))
            ii += 1

            if self.ps.sel_chon[1]:
                noise = np.std(DTDx[:inz])
            else:
                noise = np.std(DTAx[:inz])
            print('noise = ' + str(noise))
            print('absolute peak value:' + str(np.amax(CA)-nza) + 'V')

            self.plot_waveform(Apu, Bpu)
            if self.ps.sel_chon[1]:
                self.plot_pp(DTDx)
            else:
                self.plot_pp(DTAx)
        return DTAx, DTDx, Apu, Bpu

    def prepare_plots(self):
        self.ttot = self.ps.sel_dt*self.ps.sel_nsamp
        self.t = np.linspace(self.ps.sel_dt, self.ttot, self.ps.sel_nsamp)
        self.t = self.t-self.ps.sel_trig_pos_rel*self.ttot

        if self.mode == 'freerun':
            a2 = binnen(self.t.reshape(-1, 1), self.t, [0], bw=10, bt=0)
            self.t2 = a2[0].reshape(-1,)
        else:
            if len(self.t) > 2*self.plot_nbin:
                a2 = binnen(self.t.reshape(-1, 1), self.t, [0], bw=self.plot_nbin, bt=0)
                self.t2 = a2[0].reshape(-1,)
            else:
                self.t2 = self.t

        self.left_graph.axes.cla()
        # a tuple, hence the comma!
        self.line1, = self.left_graph.axes.plot(self.t2, self.t2*0.0, 'r-')
        # a tuple, hence the comma!
        self.line1b, = self.left_graph.axes.plot(self.t2, self.t2*0.0, 'b-')

        self.left_graph.draw()
        for artist in self.right_graph.fig.gca().lines:
            artist.remove()

        #self.right_graph.axes.cla()

        # a tuple, hence the comma!
        self.line2, = self.right_graph.axes.plot(self.t2, self.t2*0.0, 'r-')
        # if self.ps.sel_chon[1]:
        #    self.line3,=self.ax2.plot(self.t2,self.t2*0.0,'g-') # a tuple, hence the comma!
        #    self.line4,=self.ax2.plot(self.t2,self.t2*0.0,'b-') # a tuple, hence the comma!

        self.right_graph.draw()

    def plot_waveform(self, S1, S2):
        if len(S1) > 2*self.plot_nbin:
            y1 = binnen(S1.reshape(-1, 1), self.t, [0], bw=self.plot_nbin, bt=0)
            y2 = binnen(S2.reshape(-1, 1), self.t, [0], bw=self.plot_nbin, bt=0)
            y1 = y1[0].reshape(-1,)
            y2 = y2[0].reshape(-1,)
        else:
            y1 = S1
            y2 = S2
        self.line1.set_ydata(y1)
        self.line1b.set_ydata(y2)
        self.left_graph.draw()
        qtw.QApplication.processEvents()

    def plot_pp(self, S):
        if self.mode == 'freerun':
            if self.ps.sel_chon[1]:
                y2 = binnen(S.reshape(-1, 1), self.t, [0], bw=10, bt=0)
            else:
                y2 = binnen(S.reshape(-1, 1), self.t, [0], bw=10, bt=0)
            y2 = y2[0].reshape(-1,)

        else:
            if len(S) > 2*self.plot_nbin:
                y2 = binnen(S.reshape(-1, 1), self.t, [0], bw=self.plot_nbin, bt=0)
                y2 = y2[0].reshape(-1,)
            else:
                y2 = S

        self.line2.set_ydata(y2)
        self.right_graph.draw()
        qtw.QApplication.processEvents()

def main():
    app = qtw.QApplication(sys.argv)
    form = muspp_ui()
    form.showMaximized()
    app.exec()

if __name__ == "__main__":
    main()
