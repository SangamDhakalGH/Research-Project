import PyQt6.QtWidgets as qtw

class OutputBox(qtw.QWidget):
    """
    This class is used to display the output of the experiment.
    """
    def __init__(self):
        super().__init__()
        self.layout = qtw.QGridLayout()
        self.setLayout(self.layout)
        self.output_label = qtw.QLabel(text='Experiment Log')
        self.font_1 = self.font()
        self.font_1.setPointSize(14)
        self.output_label.setFont(self.font_1)
        self.save_output_button  = qtw.QPushButton(text='Save Log', maximumWidth=150, clicked=self.save_output_to_file)
        self.clear_output_button = qtw.QPushButton(text='Clear Log',maximumWidth=150, clicked=self.clear_output_box)
        
        self.output_box = qtw.QTextEdit()
        self.output_box.setFontFamily("Consolas")
        self.output_box.setReadOnly(True)
        self.layout.addWidget(self.output_label,0,0,1,1)
        self.layout.addWidget(self.save_output_button,0,1,1,1)
        self.layout.addWidget(self.clear_output_button,0,2,1,1)
        self.layout.addWidget(self.output_box,1,0,1,0)
        
    def save_output_to_file(self):
        """ Saves the output box to a file """
        
        # prompt user to select where to save the file
        file_name = qtw.QFileDialog.getSaveFileName(self, 'Save File', 'D:\\', '*.txt')

        output_contents = self.output_box.toPlainText()
        
        # if file_name is not empty, write to file
        if file_name[0] != '':
            with open(file_name[0], 'w') as f:
                f.write(output_contents)
 
    def clear_output_box(self):
        """ Clears the output box """
        
        self.output_box.clear()
    
    def write_to_output_box(self, text: str):
        print(text)
        self.output_box.setReadOnly(False)
        self.output_box.append(str(text))
        self.output_box.setReadOnly(True)
        qtw.QApplication.processEvents()
