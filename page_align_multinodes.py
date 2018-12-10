# coding=utf-8
import sys
from PyQt4 import QtGui, QtCore
from PyQt4.QtGui import *
from PyQt4.QtCore import *
import chi_compute
import argparse
import result_table
import os
import multiprocessing

global args,targetfile
targetfile = os.path.join(os.path.split(sys.path[0])[0], "testRes")

class Window(QtGui.QWidget):
    def __init__(self):
        super(Window, self).__init__()
        self.setWindowTitle("molecules alignment")
        self.initUI()

    def initUI(self):        
    #db box
        self.db = db
        self.dbs = QListWidget()
        self.dbs.addItems(db)
        # self.db_path = QtGui.QLabel('Database Path: ')
        #btn
        self.btn_add = QPushButton('Add')
        self.btn_remove = QPushButton('Remove')
        self.btn_sort = QPushButton('Sort')
        #vertical layout of btn
        self.v_box_db = QVBoxLayout()
        self.v_box_db.addWidget(self.btn_add)
        self.v_box_db.addWidget(self.btn_remove)
        self.v_box_db.addWidget(self.btn_sort)
        self.v_box_db.addStretch(1)
        #horizontal layout of db and btn 
        self.h_box_db = QGridLayout()
        # self.h_box_db.addWidget(self.db_path,0,0)
        self.h_box_db.addWidget(self.dbs,0,1)
        self.h_box_db.addLayout(self.v_box_db,0,2)
        #vertical layout of box_db
        self.box_db = QVBoxLayout()
        self.db_title = QLabel('---Input  Database---')
        self.box_db.addWidget(self.db_title)
        self.box_db.addLayout(self.h_box_db)
        #connect
        self.btn_add.clicked.connect(self.Add)
        self.btn_remove.clicked.connect(self.Remove)
        self.btn_sort.clicked.connect(self.Sort)

    #target box
        self.target_title = QtGui.QLabel('---Input  Target---')
        self.target = QtGui.QLabel('Target Path:')
        self.target_edit = QtGui.QLineEdit()
        self.target_type = QtGui.QLabel('Target Type:')
        self.target_type_edit = QtGui.QComboBox()
        self.target_type_edit.addItem("mol2")
        self.target_type_edit.addItem("pdb")
        self.target_type_edit.addItem("nlm")
        self.target_btn = QtGui.QPushButton("Browse")
        #horizontal layout of target inout
        # self.h_box_path = QHBoxLayout()
        # self.h_box_path.addWidget(self.target)
        # self.h_box_path.addWidget(self.target_edit)
        # self.h_box_path.addWidget(self.target_btn)
        # self.h_box_type = QHBoxLayout()
        # self.h_box_type.addWidget(self.target_type)
        # self.h_box_type.addWidget(self.target_type_edit)
        # #vertical layout of box_target
        # self.box_target = QVBoxLayout()
        # self.box_target.addWidget(self.target_title)
        # self.box_target.addLayout(self.h_box_path)
        # self.box_target.addLayout(self.h_box_type)
        #target box
        self.box_target = QtGui.QGridLayout()
        self.box_target.setSpacing(5)
        self.box_target.addWidget(self.target_title, 0, 0)
        self.box_target.addWidget(self.target, 1, 0)
        self.box_target.addWidget(self.target_edit, 1, 1)
        self.box_target.addWidget(self.target_btn, 1, 2)
        self.box_target.addWidget(self.target_type, 2, 0)
        self.box_target.addWidget(self.target_type_edit, 2, 1)
        #connect
        self.connect(self.target_type_edit, QtCore.SIGNAL('activated(QString)'), self.onActivated)
        self.connect(self.target_btn, QtCore.SIGNAL('clicked()'), self.OnButtonOpen)
    
    #rmax_box
        self.rmax_title = QtGui.QLabel('---Screen By Size ---')
        self.rmax = QtGui.QLabel('Target radius:')
        self.rmax_btn_radius = QPushButton('GET')
        self.rmax_edit = QtGui.QLineEdit()
        self.rmax_diff = QtGui.QLabel('Max Difference:')
        self.rmax_diff_edit = QtGui.QLineEdit()
        self.rmax_btn_ok = QPushButton('OK')
        self.rmax_res_edit = QtGui.QTextEdit()
        #horizontal layout
        self.h_box_rmax = QHBoxLayout()
        self.h_box_rmax.addWidget(self.rmax)
        self.h_box_rmax.addWidget(self.rmax_edit)
        self.h_box_rmax.addWidget(self.rmax_btn_radius)
        self.h_box_rmax.addWidget(self.rmax_diff)
        self.h_box_rmax.addWidget(self.rmax_diff_edit)
        self.h_box_rmax.addWidget(self.rmax_btn_ok)
        #vertical layout of box_fnl
        self.box_rmax = QVBoxLayout()
        self.box_rmax.addWidget(self.rmax_title)
        self.box_rmax.addLayout(self.h_box_rmax)
        self.box_rmax.addWidget(self.rmax_res_edit)
        #connect
        self.connect(self.rmax_btn_radius, QtCore.SIGNAL('clicked()'), self.GetRadius)
        # self.connect(self.rmax_btn_ok, QtCore.SIGNAL('clicked()'), self.OnButtonRmax)
    
    #fnl_box
        self.fnl_title = QtGui.QLabel('---Screen By Fnl ---')
        self.node_num = QtGui.QLabel('Node Number:')
        self.node_num_edit = QtGui.QComboBox()
        self.node_num_edit.addItem("1")
        self.node_num_edit.addItem("2")
        self.process_num = QtGui.QLabel('Process Per Node:')
        self.process_num_edit = QtGui.QComboBox()
        self.process_num_edit.DropDownStyle = 'DropDown'
        self.process_num_edit.addItem("1")
        self.process_num_edit.addItem("2")
        self.process_num_edit.addItem("3")
        self.process_num_edit.addItem("...")
        self.fnl_lable_f = QtGui.QLabel('Top')
        self.fnl_num_edit = QtGui.QLineEdit()
        self.fnl_lable_r = QtGui.QLabel('Result')
        self.fnl_btn_ok = QPushButton('OK')
        self.fnl_res_edit = QtGui.QTextEdit()
        #horizontal layout of node_num and process_num
        self.h_box_setting = QHBoxLayout()
        self.h_box_setting.addWidget(self.node_num)
        self.h_box_setting.addWidget(self.node_num_edit)
        self.h_box_setting.addWidget(self.process_num)
        self.h_box_setting.addWidget(self.process_num_edit)
        self.h_box_setting2 = self.h_box_setting
        #horizontal layout of fnl_number , node_num and process_num
        self.h_box_fnl = QHBoxLayout()
        self.h_box_fnl.addLayout(self.h_box_setting)
        self.h_box_fnl.addWidget(self.fnl_lable_f)
        self.h_box_fnl.addWidget(self.fnl_num_edit)
        self.h_box_fnl.addWidget(self.fnl_lable_r)
        self.h_box_fnl.addWidget(self.fnl_btn_ok)
        #vertical layout of box_fnl
        self.box_fnl = QVBoxLayout()
        self.box_fnl.addWidget(self.fnl_title)
        self.box_fnl.addLayout(self.h_box_fnl)
        self.box_fnl.addWidget(self.fnl_res_edit)
        #connect
        # self.connect(self.fnl_btn_ok, QtCore.SIGNAL('clicked()'), self.OnButtonChi)

    #cnlm_box
        self.cnlm_title = QtGui.QLabel('---Screen By Cnlm ---')
        self.node_num = QtGui.QLabel('Node Number:')
        self.node_num_edit = QtGui.QComboBox()
        self.node_num_edit.addItem("1")
        self.node_num_edit.addItem("2")
        self.process_num = QtGui.QLabel('Process per node:')
        self.process_num_edit = QtGui.QComboBox()
        self.process_num_edit.DropDownStyle = 'DropDown'
        self.process_num_edit.addItem("1")
        self.process_num_edit.addItem("2")
        self.process_num_edit.addItem("3")
        self.process_num_edit.addItem("...")
        self.cnlm_num_label_f = QtGui.QLabel('Top')
        self.cnlm_num_edit = QtGui.QLineEdit()
        self.cnlm_num_label_r = QtGui.QLabel('Result')
        self.cnlm_btn_ok = QPushButton('OK')
        self.cnlm_res_edit = QtGui.QTextEdit()
        #horizontal layout of node_num and process_num
        self.h_box_setting2 = QHBoxLayout()
        self.h_box_setting2.addWidget(self.node_num)
        self.h_box_setting2.addWidget(self.node_num_edit)
        self.h_box_setting2.addWidget(self.process_num)
        self.h_box_setting2.addWidget(self.process_num_edit)
        #horizontal layout of cnlm_number , node_num and process_num
        self.h_box_cnlm = QHBoxLayout()
        self.h_box_cnlm.addLayout(self.h_box_setting2)
        self.h_box_cnlm.addWidget(self.cnlm_num_label_f)
        self.h_box_cnlm.addWidget(self.cnlm_num_edit)
        self.h_box_cnlm.addWidget(self.cnlm_num_label_r)
        self.h_box_cnlm.addWidget(self.cnlm_btn_ok)
        #vertical layout of box_fnl
        self.box_cnlm = QVBoxLayout()
        self.box_cnlm.addWidget(self.cnlm_title)
        self.box_cnlm.addLayout(self.h_box_cnlm)
        self.box_cnlm.addWidget(self.cnlm_res_edit)
        #connect
        # self.connect(self.cnlm_btn_ok, QtCore.SIGNAL('clicked()'), self.OnButtonCC)


    #main box
        self.main_box = QVBoxLayout()
        self.main_box.addLayout(self.box_db)
        self.main_box.addLayout(self.box_target)
        self.main_box.addLayout(self.box_rmax)
        self.main_box.addLayout(self.box_fnl)
        self.main_box.addLayout(self.box_cnlm)
        self.exit_btn = QPushButton('EXIT')
        self.main_box.addWidget(self.exit_btn)
        self.setLayout(self.main_box)
        self.resize(350, 300)
        self.connect(self.exit_btn, QtCore.SIGNAL('clicked()'), QtGui.qApp, QtCore.SLOT('quit()'))  
        

    def process_args(self):
        parser = argparse.ArgumentParser()
        parser.add_argument("-output", help="chi result path", type=str)
        parser.add_argument("-dbpath", help="db path", type=str)
        parser.add_argument("-processnum", help="processnum number", type=int)
        parser.add_argument("-cavitypath", help="cavity path", type=str)
        parser.add_argument("-cavitytype", help="cavity type", type=str)
        parser.add_argument("--nmax", help="nmax", default=20, type=int)
        parser.add_argument("--rmax", help="rmax", type=float)
        return parser.parse_args()

    #定义槽
    def Add(self):
        #添加
        absolute_path = QtGui.QFileDialog.getOpenFileName(self, 'Open file', '.', "")
        if absolute_path:
            self.dbs.addItem(absolute_path)      
         
    def Remove(self):
        #移除
        if QMessageBox.warning(self, 'ok', 'remove this item from the list?', QMessageBox.Ok | QMessageBox.Cancel) == QMessageBox.Ok:
            item_deleted = self.dbs.takeItem(self.dbs.currentRow())
            #将读取的值设置为None
            item_deleted = None
        
    def Sort(self):
        #排序
        self.dbs.sortItems(Qt.AscendingOrder)

    def OnButtonOpen(self):
        absolute_path = QtGui.QFileDialog.getOpenFileName(self, 'Open file', '.', "")
        if absolute_path:
            self.target_edit.setText(absolute_path)

    def onActivated(self, text):
        self.target_type_edit.setItemText(1, text)

    def GetRadius(self):
        self.rmax_edit.setText(10)

    def OnButtonOk(self):
        targetInput = self.target_edit.text()
        targetTypeInput = self.target_type_edit.currentText()
        processNumInput = self.processNumEdit.text()
        fnlInput = self.fnlNumEdit.text()
        cnlmInput = self.cnlmNumEdit.text()
        args = self.process_args()    
        args.cavitypath = str(targetInput)
        args.cavitytype = str(targetTypeInput)
        args.processnum = int(processNumInput)
        args.fnl_round = int(fnlInput)
        args.cnlm_round = int(cnlmInput)
        for i in range(self.addCount):
            datasetEdit_n = "datasetEdit" + str(i)
            dbpathInput = self.datasetEdit_n.text()
            args.dbpath = str(dbpathInput).split(".")[0]
            chi_compute.run(args)
        
        print args
        dialog = MyDialog()
        # chi_compute.run(args)
        fileStr = ''
        for line in open(targetfile, 'r').readlines():
            fileStr+=line
        self.result_table.setText(fileStr)   
        dialog.exec_()

class MyDialog(QtGui.QDialog):
    def __init__(self):
        super(MyDialog, self).__init__()
        self.setWindowTitle("Doing Alignment...")
        self.gridlayout = QtGui.QGridLayout()
        self.label = QtGui.QLabel("start the task!")
        self.label1 = QtGui.QLabel("waiting...")
        self.label2 = QtGui.QLabel("task finished!")
        self.cancalButton = QtGui.QPushButton("OK")
        self.gridlayout.addWidget(self.label, 0, 0)
        self.gridlayout.addWidget(self.label1, 1, 0)
        self.gridlayout.addWidget(self.label2, 2, 0)
        self.gridlayout.addWidget(self.cancalButton, 3, 0)     
        self.setLayout(self.gridlayout)
        self.connect(self.cancalButton, QtCore.SIGNAL('clicked()'), self.OnCancel)

    def OnCancel(self):
        self.done(0)

if __name__ == '__main__':
    app = QtGui.QApplication(sys.argv)
    db = []
    win = Window()
    win.show()
    sys.exit(app.exec_())