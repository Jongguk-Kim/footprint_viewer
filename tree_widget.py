#-*- coding:utf-8 -*-
import sys, os
from PyQt5 import QtCore, QtGui, QtWidgets

from PyQt5.QtWidgets import QMainWindow,QApplication, QTreeWidget, QTreeWidgetItem

def makingFullFilePath_linux(cwd=None, file=None): 
    if isinstance(cwd, type(None)): 
        cwd = os.getcwd()
    if cwd[-1] != "/": cwd+='/'

    if "../" in file: 
        names  =file.strip().split("../")
        cnt = 0 
        for name in names: 
            if "" in name: 
                cnt += 1 
        dirs = cwd.split("/")
        path =""
        for i in range(len(dirs)-cnt): 
            if dirs[i]=='':
                path +='/'
            else:
                path += dirs[i]+'/'
        for name in names: 
            if name != "": 
                path += name 
    elif '/home' == file[:5]: 
        path = file 
    elif '' == file: 
        return cwd
    else: 
        if '/' == file[0]: 
            file = file[1:]
        path = cwd + file 

    return path 

fextension=['.inp', '.ptn']

# class TREEWIDGET(QMainWindow): 
class Ui_Dialog(object):
    def retranslateUi(self, Dialog):
        _translate = QtCore.QCoreApplication.translate
        Dialog.setWindowTitle(_translate("Dialog", "Dialog"))

    def setupUi(self, Dialog, extensions=None, sftp=None):
        self.extensions=extensions
        if isinstance(self.extensions, type(None)): 
            self.extensions=fextension
        super().__init__()
        # self.setGeometry(100, 100, 500, 600)
        # self.setWindowTitle('Directory')

        self.searchedDirectories=[]
        self.searchedFiles=[]

        Dialog.setObjectName("Dialog")
        Dialog.resize(642, 849)
        self.gridLayout = QtWidgets.QGridLayout(Dialog)
        self.gridLayout.setObjectName("gridLayout")
        

        self.retranslateUi(Dialog)
        Dialog.finished['int'].connect(Dialog.close)
        QtCore.QMetaObject.connectSlotsByName(Dialog)

        self.savefile = 'searched_directory.sdt'

        self.dialog = Dialog

        head_item = QtWidgets.QTreeWidgetItem()
        

        self.file_return = False
 
        head_item.setText(0, 'File Name')
        head_item.setText(1, 'Full Path')
        

        self.rootTree = QTreeWidget()
        self.rootTree.setHeaderItem(head_item)
        self.rootTree.header().setVisible(True)
        self.rootTree.setAlternatingRowColors(True)
        self.rootTree.itemCollapsed.connect(self.collapsed_event)
        self.rootTree.itemExpanded.connect(self.expanded_event)
        self.rootTree.itemDoubleClicked.connect(self.doubleclicked_event)
        

        self.gridLayout.addWidget(self.rootTree, 0, 0, 1, 1)


        self.rootColumn = 1
        self.nameColumn = 0
        self.rootTree.setColumnWidth(0, 500)


        # self.setCentralWidget(self.rootTree)
        
        try: 
            self.host, self.user, self.pw= sftp 
        except: 
            self.host='10.82.66.65'
            self.pw =  self.user='h20200155'
        
        self.connection = self.connectFTP()

        if self.connection: 
            if self.user == 'fiper': 
                self.home = '/home/users/fiper/ISLM/ISLM_JobFolder'
            else: 
                self.home = '/home/users/'+self.user
                
            loaded = self.loadData(self.home)
            
            if not loaded: 
               
                homedir = self.sftp.listdir(self.home)
                self.searchedDirectories.append(self.home)
                # print (' ADD :', self.searchedDirectories[-1])
                homedir = sorted(homedir)
                self.items =[]
                cnt = 0 

                for name in homedir: 
                    if '.' != name[0] and '_' != name[0]:
                        cnt += 1 
                        results = self.searchDirectory(self.home+"/"+name)
                        if results:
                            self.items.append(QTreeWidgetItem(self.rootTree))
                            self.items[-1].setText(self.nameColumn, name+"/")
                            self.items[-1].setText(self.rootColumn, self.home)
                            dirs, files = results
                            subs=[]
                            if len(dirs): 
                                self.searchedDirectories.append(self.home+"/"+name)
                                # print (' ADD :', self.searchedDirectories[-1])
                            for dname in dirs: 
                                subs.append(QTreeWidgetItem(self.items[-1]))
                                subs[-1].setText(self.nameColumn, dname+"/")
                                directory =  self.home+"/"+name
                                if directory[-1] == '/': directory=directory[:-1]
                                subs[-1].setText(self.rootColumn,directory)
                            for fn in files: 
                                for ex in self.extensions: 
                                    if ex in fn: 
                                        subs.append(QTreeWidgetItem(self.items[-1]))
                                        subs[-1].setText(self.nameColumn, fn)
                                        directory =  self.home+"/"+name
                                        if directory[-1] == '/': directory=directory[:-1]
                                        subs[-1].setText(self.rootColumn,directory)
                                        self.searchedFiles.append(directory+"/"+fn)
                                        break 
                        else: 
                            for ex in self.extensions: 
                                if ex in name:
                                    self.items.append(QTreeWidgetItem(self.rootTree))
                                    self.items[-1].setText(self.nameColumn, name)
                                    self.items[-1].setText(self.rootColumn, self.home)
                                    self.searchedFiles.append(self.home+"/"+fn)
                                    break 
                try: 
                    fp = open("filename_path", 'r')
                    cDir = fp.readline().strip()
                    fp.close()
                    print(cDir)
                    preDirectory = self.checkHistory(self.home, cDir)
                except: 
                    preDirectory = False 
                preDirectory = False
                if preDirectory: 
                    home = self.home.split("/")
                    crnt = cDir.split("/")
                    n = len(home)
                    N = len(crnt)
                    loopdirectory = self.home +"/"+ crnt[n]
                    for m in range(n+1, N): 
                        loopdirectory += "/"+crnt[m]
                        print ("add dir: %s"%(loopdirectory))
                        self.addItemsToTree(loopdirectory)

                        self.rootTree.itemExpanded

    def addItemsToTree(self, root): 
        homedir = self.sftp.listdir(root)
        self.searchedDirectories.append(root)
        # print (' ADD :', self.searchedDirectories[-1])
        homedir = sorted(homedir)
        self.items =[]
        cnt = 0 

        for name in homedir: 
            if '.' != name[0] and '_' != name[0]:
                cnt += 1 
                results = self.searchDirectory(root+"/"+name)
                print ("   %s"%(name))
                if results:
                    self.items.append(QTreeWidgetItem(self.rootTree))
                    self.items[-1].setText(self.nameColumn, name+"/")
                    self.items[-1].setText(self.rootColumn, root)
                    dirs, files = results
                    subs=[]
                    if len(dirs): 
                        self.searchedDirectories.append(root+"/"+name)
                        print (' ADD DR:', self.searchedDirectories[-1])
                    for dname in dirs: 
                        subs.append(QTreeWidgetItem(self.items[-1]))
                        subs[-1].setText(self.nameColumn, dname+"/")
                        
                        directory =  root+"/"+name
                        if directory[-1] == '/': directory=directory[:-1]
                        subs[-1].setText(self.rootColumn, directory)

                    for fn in files: 
                        for ex in self.extensions: 
                            if ex in fn: 
                                subs.append(QTreeWidgetItem(self.items[-1]))
                                print (' ADD FL :', subs[-1])
                                subs[-1].setText(self.nameColumn, fn)
                                directory =  root+"/"+name
                                if directory[-1] == '/': directory=directory[:-1]
                                subs[-1].setText(self.rootColumn,directory)
                                self.searchedFiles.append(directory+"/"+fn)
                                break 


    def checkHistory(self, orghome, directory):  
        home = orghome.split("/")
        current = directory.split("/")

        N = len(home)
        for i in range(N): 
            if home[i] != current[i]: 
                return orghome 
        else: 
            return True 

    def searchDirectory(self, home):
        files=[]
        dirs=[]
        try: 
            current = self.sftp.listdir(home)
        except : 
            # print ("FILE: %s"%(home+"/"+name))
            return False
        current = sorted(current)
        for sb in current : 
            try: 
                _ = self.sftp.listdir(home+"/"+sb)
                dirs.append(sb)
            except: 
                files.append(sb)

        return [dirs, files]


    def makeTree(self, homedir, rootTree): 
        # print (" > Searching : %s"%(homedir))
        
        for sdir in self.searchedDirectories: 
            if sdir == homedir: 
                return rootTree
        else: 
        
            try: 
                homes = self.sftp.listdir(homedir)
                if homedir[-1] == '/': homedir = homedir[:-1]
                self.searchedDirectories.append(homedir)
                # print (' ADD :', self.searchedDirectories[-1])
                
            except Exception as EX: 
                print (EX, homedir)
                return 
            
            if homedir[-1] == '/': homedir = homedir[:-1]
            # print("###############################################")
            # print ("SEARCHING DIR :", homedir)
            
            homes = sorted(homes)
            items =[]
            for name in homes:
                # if '.' != name[0] and '_' != name[0]:
                    results = self.searchDirectory(homedir+"/"+name)
                    
                    if results:
                        items.append(QTreeWidgetItem(rootTree))
                        items[-1].setText(self.nameColumn, name+"/")
                        items[-1].setText(self.rootColumn, homedir)
                        dirs, files = results
                        subs=[]
                        for dname in dirs: 
                            # print ('dname', dname)
                            if sdir in self.searchedDirectories: 
                                if sdir == homedir+"/"+dname: 
                                    # print (homedir+"/"+dname, ' is exists (directory)')
                                    break 
                            else: 
                                subs.append(QTreeWidgetItem(items[-1]))
                                subs[-1].setText(self.nameColumn, dname+"/")
                                directory =  homedir+"/"+name
                                if directory[-1] == '/': directory=directory[:-1]
                                subs[-1].setText(self.rootColumn, directory)
                                # print ("  Add dir item", dname)
                                # subItems = self.makeTree(homedir+"/"+name +"/"+ dir, subs[-1])
                        for fn in files: 
                            for ex in self.extensions: 
                                if ex in fn:
                                    for sf in self.searchedFiles: 
                                        if sf == homedir+"/"+name+"/"+fn: 
                                            # print (homedir+"/"+name+"/"+fn, ' is exists (file)')
                                            break 
                                    else:
                                        subs.append(QTreeWidgetItem(items[-1]))
                                        subs[-1].setText(self.nameColumn, fn)
                                        subs[-1].setText(self.rootColumn, homedir+"/"+name)
                                        # print('Add File', fn)
                                        break 
                    else: 
                        pass 
            return items
    
    def doubleclicked_event(self, item  ): 
            if item.text(self.nameColumn)[-1] !="/" : 
                # print (item.text(self.rootColumn)+"/"+item.text(self.nameColumn))
                fp=open('filename_path', 'w')
                fp.write("%s\n"%(item.text(self.rootColumn)))
                fp.write("%s\n"%(item.text(self.nameColumn)))
                fp.close()

                self.saveData()
            
                self.dialog.close()
            else: 
                dname = item.text(self.rootColumn)+"/"+item.text(self.nameColumn)
                if dname[-1] == "/": dname = dname[:-1]
                for name in self.searchedDirectories:
                    if dname == name: 
                        # print (" Already exist", dname)
                        break 
                else: 
                    self.makeTree(item.text(self.rootColumn)+"/"+item.text(self.nameColumn), item)

    def saveData(self): 
        fp=open(self.savefile, 'w')
        for fn in self.searchedDirectories: 
            fp.write("%s&\n"%fn)
        for fn in self.searchedFiles: 
            fp.write("%s\n"%fn)
        fp.close()

    def loadData(self, homedir): 
        print('start loading')
        if not os.path.isfile(self.savefile): 
            print ("no file")
            return False 
            
        with open(self.savefile) as F: 
            lines = F.readlines()
        try :
            if homedir != lines[0].strip()[:-1]: 
                print("different directory")
                return False 
        except: 
            print ("no lines")
            return False 

        self.items =[]
        return False  


    def connectFTP(self): 
        # host = '10.82.66.65'
        # pw = self.user = 'h20200155'

        import paramiko as FTP 
        self.ftp = FTP.SSHClient()
        self.ftp.set_missing_host_key_policy(FTP.AutoAddPolicy())
        try: 
            port = 22 
            self.ftp.connect(self.host, username=self.user, password=self.pw, port=port, timeout=3)
            self.sftp = self.ftp.open_sftp()
            return True 
        except Exception as EX:
            print (EX)
            return False 

    def collapsed_event(self, item):
        # print('collapsed_event : ', item.text(0))
        pass 
    def expanded_event(self, item):
        # print('expanded_event : ', item.text(0))
        pass 

if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    Dialog = QtWidgets.QDialog()
    ui = Ui_Dialog()
    ui.setupUi(Dialog)
    Dialog.show()
    sys.exit(app.exec_())
