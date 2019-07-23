
# coding=utf-8
import sys
from PyQt4 import QtGui, QtCore
import argparse
# import crop
# import peak_cluster



class MyDialog(QtGui.QDialog):
    def __init__(self):
        super(MyDialog, self).__init__()
        self.gridlayout = QtGui.QGridLayout()
        self.label = QtGui.QLabel("start the task!")
        self.label1 = QtGui.QLabel("waiting for the generation...")
        self.label2 = QtGui.QLabel("task finished!")
        self.cancalButton = QtGui.QPushButton("OK")
        self.gridlayout.addWidget(self.label, 0, 0)
        self.gridlayout.addWidget(self.label1, 1, 0)
        self.gridlayout.addWidget(self.label2, 2, 0)
        self.gridlayout.addWidget(self.cancalButton, 3, 0)
        self.connect(self.cancalButton, QtCore.SIGNAL('clicked()'), self.OnCancel)
        self.setLayout(self.gridlayout)

    def OnCancel(self):
        self.done(0)

class Window(QtGui.QWidget):
    def __init__(self):
        super(Window, self).__init__()
        self.setWindowTitle("Template generation...")
        self.dataset = QtGui.QLabel('Input')
        self.datasetEdit = QtGui.QLineEdit()
        self.output = QtGui.QLabel('Output')
        self.outputEdit = QtGui.QLineEdit()
        self.dbButton = QtGui.QPushButton("Browse")
        self.fileType = QtGui.QLabel('Frames used') 
        self.fileTypeEdit = QtGui.QLineEdit()  
        self.createDatabaseButton = QtGui.QPushButton("OK")
        self.createExitButton = QtGui.QPushButton("Exit")
        gridlayout = QtGui.QGridLayout()
        gridlayout.addWidget(self.dataset, 1, 0)
        gridlayout.addWidget(self.datasetEdit, 1, 1)
        gridlayout.addWidget(self.dbButton, 1, 2)
        gridlayout.addWidget(self.fileType, 2, 0)
        gridlayout.addWidget(self.fileTypeEdit, 2, 1)
        gridlayout.addWidget(self.output, 3, 0)
        gridlayout.addWidget(self.outputEdit, 3, 1)
        gridlayout.addWidget(self.createDatabaseButton, 4, 0)
        gridlayout.addWidget(self.createExitButton, 4, 1)
        self.setLayout(gridlayout)
        self.connect(self.fileTypeEdit, QtCore.SIGNAL('activated(QString)'), self.onActivated)
        self.connect(self.dbButton, QtCore.SIGNAL('clicked()'), self.OnButtonOpen)
        self.connect(self.createDatabaseButton, QtCore.SIGNAL('clicked()'), self.OnButton_createDatabase)
        self.connect(self.createExitButton, QtCore.SIGNAL('clicked()'), QtGui.qApp, QtCore.SLOT('quit()'))

    def OnCancel(self):
        self.done(0)

    def onActivated(self, text):
        self.fileTypeEdit.setItemText(1, text)
        if text == 'mol2':
            dbpath = self.datasetEdit.text()
            pathList = str(dbpath).split('/')
            databaseName = pathList.pop().split('.')[0]
            self.outputEdit.setText(databaseName)


    def OnButtonOpen(self):
        fileName = QtGui.QFileDialog.getOpenFileName(self, 'Open')
        if not fileName.isEmpty():
            self.datasetEdit.setText(fileName) 

    def process_args(self):
        parser = argparse.ArgumentParser()
        parser.add_argument("--prefix", help="pdb result name, default pdb_db_name=myTestDB, default mol2_db_name = mol2fileName", default="myTestDB", type=str)
        parser.add_argument("-path", help="db path", type=str)
        parser.add_argument("--np", help="number of point covering [0,1]", default=50, type=int)
        parser.add_argument("--fix_dx", help="Whether keeping default dx=0.7A or not", default=True, type=bool)
        parser.add_argument("--nmax", help="maximum order of zernike expansion", default=20, type=int)
        parser.add_argument("--qmax", help="maximum q value, for which the intensity to be evaluated", default=0.3, type=float)
        parser.add_argument("--uniform", help="Whether the uniform density is used", default=True, type=bool)
        parser.add_argument("--buildmap", help="Whether xplor map will be constructed or not", default=False, type=bool)
        parser.add_argument("--shift", help="Whether the coordates will be shifted or not", default=False, type=bool)
        parser.add_argument("--coef_out", help="Whether dump zernike moments to pickle files", default=False,type=bool)
        return parser.parse_args()

    def OnButton_createDatabase(self):
        fileTypeEditInput = self.fileTypeEdit.text()
        dbpath = self.datasetEdit.text()
        output = self.outputEdit.text()
        args = self.process_args()
        pathList = str(dbpath).split('/')
          
        dialog = MyDialog()       
        pathList.pop()
        args.path = '/'.join(pathList) + '/'
        args.prefix = str(output)
        print(args)
        crop.run(args)
        dialog.exec_()



app = QtGui.QApplication(sys.argv) # QApplication eats argv in constructor
win = Window()
win.show()
app.exec_()




