# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'ImageView.ui'
#
# Created by: PyQt5 UI code generator 5.9.2
#
# WARNING! All changes made in this file will be lost!

## make exe file 
## pip install pyinstaller
## pyinstaller -F *.py 
## pyinstaller .exe --onefile --windowwed *.py 
# from matplotlib import use
# use('TkAgg')
import matplotlib.pyplot as plt
# plt.style.use( 'ggplot') 
from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtWidgets import QListWidget, QListWidgetItem

import math, struct
from os import getcwd
from os.path import isfile
# import numpy as np
import numpy as np 


from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.patches import Arc



TireComponents = [
    'BEAD_R', 'BEAD_L', # 'BD1' ,  Bead
    'BEC' , # Belt Edge Cushion
    'BSW' , # Black SideWall Component
    'BTT' , # Belt rubber component (associated with BT components)
    'CCT' , # Carcass rubber component (associated with C components)
    'SCT' , 
    'NCT' ,
    'CHS' ,
    'CTB' , 'CTR', # Tread Componet
    'ET'  , # Edge Tape
    'FIL' , 'BDF', # Bead Filler
    'HUS' , # Hump Strip
    'JBT' , # Rubber associated with JFC, JEC components
    'L11' , 'IL1', # Inner liner component
    'L12' , 'IL2', # #2 Innerliner
    'LBF' , # Lower Bead Filler
    'UBF' , # Upper Bead Filler
    'RIC' , # Rim Cushion
    'SIR' , # Sidewall Insert Rubber
    'SHW' , # Shoulder Wedge
    'SRTT', # Associated with PK1, PK2, RFM and FLI components
    'SUT' , 'UTR', # SubTread
    'TRW' , # Tread Wing
    
    'WSW' , # While Sidewall
    'C01'  , 'CC1', # Carcass Cord 1 
    'C02'  , 'CC2', # Carcass Cord 2
    'C03'  , 'CC3', # Carcass Cord 3 
    'C04'  , 'CC4', # Carcass Cord 4
    'BT1'  , # Belt 1 
    'BT2'  , # Belt 2
    'BT3'  , # Belt 3
    'BT4'  , # Belt 4
    'JFC1' , 'JFC', # Jointless Full Cap 1
    'JFC2' , # Jointless Full Cap 2
    'JEC1' , 'JEC', # Jointless Edge Cap 1
    'OJEC1', # Overlapped Jointless Edge Cap
    'OJFC1', # Overlapped Jointless Full Cap
    'OJEC2', # Overlapped Jointless Edge Cap
    'OJFC2', # Overlapped Jointless Full Cap
    'PK1'  , # Half Bead Packing
    'PK2'  , # Full Bead Packing
    'RFM'  , # Bead S(RFM)
    'FLI'  , # Bead Flipper
    'CH1'  , 'CH1_R', 'CH1_L',  # Steel Chafer 
    'CH2'  , 'CH2_R', 'CH2_L',  # 1st Nylon Chafer
    'CH3'  , 'CH3_R', 'CH3_L',  # 2nd Nylon Chafer
    'BDC'  , # bead cover 
    'SPC'  ,  ## spiral coil
    'SWS'    # temporary component for Endurance simulation 
]


TreadElset = ['CTB', 'SUT', 'CTR', 'UTR', 'TRW']
ChaferName = ['CH1', 'CH2', 'CH3', 'SCF', 'NCF']
lst_colors = ['black', 'red', 'blue', 'green', 'orange', 'gray','pink', 'brown', 'olive', 'cyan', 'purple']

def readSDB(sdb): 

    CAMBER = [0.0]
    NLB   = []      
    NODE  = []  
    ELM2D = []
    ELM3D = []
    ELB   = []
    DISP  = []
    RR    = []
    RIM_ROAD =[]

    file = open(sdb, 'rb')

    file.seek(0,2); fend = file.tell() #find the end position of binary file
    file.seek(0);   fpos = file.tell() #find the start position of binary file
    # # print fpos, fend
    testfile = open("Node.inp","w")

    
    READ_LENGTH = 0
    
    #Read from start to end
    while file.tell() < fend:
        #distinguish each information from BlockID with 4 byte
        #each information includes title, local coordinates, nodes, elements(2D & 3D), rim, and road
        BlockID  = struct.unpack('i', file.read(4))[0]

        if BlockID == 11:
            BlockLength = struct.unpack('i', file.read(4))[0]
            
            BlockTitle = ''
            for i in range(BlockLength):    
                BlockTitle += str(struct.unpack('c', file.read(1))[0])
            # print (BlockLength, "Block Titie=%s"%(BlockTitle)) 
        elif BlockID == 12:
            #we can check the local coordinate of tire and camber angle by BlockID 12
            LocalCoord = []
            i = 0
            while i < 9:
                LocalCoord.append(struct.unpack('d', file.read(8))[0])
                i = i + 1
            CAMBER[0] = math.atan(LocalCoord[5]/LocalCoord[4])*180/math.pi
            # print CAMBER
            del LocalCoord
        elif BlockID == 21:
            #we can check the total node number(=NodesNUM) and each label & X,Y,Z coordinates by BlockID 21
            NodesNUM = struct.unpack('i', file.read(4))[0]
            i = 0
            while i < NodesNUM:
                Label = struct.unpack('i', file.read(4))[0]
                NLB.append(Label) #for making of abaqus odb file
                i = i + 1
            # print "Total Node Number =", len(NLB)
            testfile.writelines("*NODE,\n")
            i = 0
            while i < NodesNUM:
                X = struct.unpack('d', file.read(8))[0]
                Y = struct.unpack('d', file.read(8))[0]
                Z = struct.unpack('d', file.read(8))[0]
                NODE.append([NLB[i], X, Y, Z, 0.0]) #for making of abaqus odb file
                testfile.writelines(" %8d, %5.6e, %5.6e, %5.6e\n" % (NLB[i], X, Y, Z))
                i = i + 1
        elif BlockID == 31:
            #we can check the total 3d solid element number(=E3DNUM) and each label & nodes by BlockID 31
            E3DNUM = struct.unpack('i', file.read(4))[0]
            i = 0
            while i < E3DNUM:
                ID = struct.unpack('i', file.read(4))[0]
                N1 = struct.unpack('i', file.read(4))[0]
                N2 = struct.unpack('i', file.read(4))[0]
                N3 = struct.unpack('i', file.read(4))[0]
                N4 = struct.unpack('i', file.read(4))[0]
                N5 = struct.unpack('i', file.read(4))[0]
                N6 = struct.unpack('i', file.read(4))[0]
                N7 = struct.unpack('i', file.read(4))[0]
                N8 = struct.unpack('i', file.read(4))[0]
                ELM3D.append([ID, N1, N2, N3, N4, N5, N6, N7, N8])
                i  = i + 1
        elif BlockID == 32:
            #we can check the total 2d element number(=E2DNUM) and each label & nodes by BlockID 32
            E2DNUM = struct.unpack('i', file.read(4))[0]
            testfile.writelines("*ELEMENT, TYPE=M3D4R\n")
            # Surf.count = E2DNUM
            i = 0
            while i < E2DNUM:
                ID = struct.unpack('i', file.read(4))[0]
                N1 = struct.unpack('i', file.read(4))[0]
                N2 = struct.unpack('i', file.read(4))[0]
                N3 = struct.unpack('i', file.read(4))[0]
                N4 = struct.unpack('i', file.read(4))[0]
                ELM2D.append([ID, N1, N2, N3, N4])
                testfile.writelines(" %8d, %8d, %8d, %8d, %8d\n" % (ID, N1, N2, N3, N4))
                i  = i + 1
        elif BlockID == 51 or BlockID == 52:
            #we can check the elset name and each elset number(=ENUM) by BlockID 51&52
            BlockLength = struct.unpack('i', file.read(4))[0]
            Eleset = ''
            i = 0
            while i < BlockLength:
                Eleset = Eleset + str(struct.unpack('c', file.read(1))[0])
                i = i + 1
            BlockFlag   = struct.unpack('i', file.read(4))[0]
            ENUM = struct.unpack('i', file.read(4))[0]
            ESET = []
            i = 0
            while i < ENUM:
                ESET.append(struct.unpack('i', file.read(4))[0])
                i = i + 1
            # print "Finish to Read", Eleset
            del ESET
        elif BlockID == 61:
            #we can check the Rim Information by BlockID 61
            # print "Read Rim Information and Save them!"
            ControlNodeID = struct.unpack('i', file.read(4))[0]
            # print "Rim Control Node =", ControlNodeID
            X = struct.unpack('d', file.read(8))[0]
            Y = struct.unpack('d', file.read(8))[0]
            Z = struct.unpack('d', file.read(8))[0]
            R = struct.unpack('d', file.read(8))[0]
            W = struct.unpack('d', file.read(8))[0]
            G = struct.unpack('i', file.read(4))[0]
            RIM_ROAD.append([ControlNodeID, X, Y, Z])
            # # print "R, W, G =", R, W, G
            i = 0
            while i < G:
                BlockFlag   = struct.unpack('i', file.read(4))[0]
                # # print "BlockFlag = ", BlockFlag
                i = i + 1
            i = 0
            while i < G:
                X1 = struct.unpack('d', file.read(8))[0]
                Y1 = struct.unpack('d', file.read(8))[0]
                X2 = struct.unpack('d', file.read(8))[0]
                Y2 = struct.unpack('d', file.read(8))[0]
                # # print "X1, Y1, X2, Y2 =", X1*1000, Y1*1000, X2*1000, Y2*1000
                i = i + 1
        elif BlockID == 62:
            #we can check the Road Information by BlockID 62
            # print "Read Road Information and Save them!"
            Tref = struct.unpack('d', file.read(8))[0]
            ControlNodeID = struct.unpack('i', file.read(4))[0]
            # print "Road Control Node =", ControlNodeID
            X = struct.unpack('d', file.read(8))[0]
            Y = struct.unpack('d', file.read(8))[0]
            Z = struct.unpack('d', file.read(8))[0]
            R = struct.unpack('d', file.read(8))[0]
            W = struct.unpack('d', file.read(8))[0]
            L = struct.unpack('d', file.read(8))[0]
            RIM_ROAD.append([ControlNodeID, X, Y, Z])
        elif BlockID == 999:
            BlockID     = struct.unpack('i', file.read(4))[0]
            ESP         = ''
            i = 0
            while i < 4:
                ESP     = ESP+ str(struct.unpack('c', file.read(1))[0])
                i = i + 1
            RecordHead  = struct.unpack('i', file.read(4))[0]
            BlockID     = struct.unpack('i', file.read(4))[0]
            SolverInfo  = ''
            i = 0
            while i < 42:
                SolverInfo = SolverInfo+ str(struct.unpack('c', file.read(1))[0])
                i = i + 1
            i = 0
            while i < 9:
                BlockID     = struct.unpack('i', file.read(4))[0]
                i = i + 1
            SimulationType = ''
            i = 0
            while i < 26:
                SimulationType = SimulationType+ str(struct.unpack('c', file.read(1))[0])
                i = i + 1
        else:
            break

    # del NLB, XYZ, ELM2D, ELM3D, DISP
    testfile.close()
    
    for e in ELM3D: 
        if e[7] == e[8]: 
            e[4] = e[5]; e[5] = e[6]; e[6] = e[7]; e[8] = 0; e[7] =0 

    # ELM3D =  np.array(ELM3D)
    # ELM2D = np.array(ELM2D)
    # print ("Node", np.array(NODE))
    # print ("RIM ROAD", RIM_ROAD)

    if len(NODE) > 0: #and len(Surf.surflabel) > 0:
        return NODE, ELM2D, ELM3D, RIM_ROAD
    else:
        return -1

def SDBResult_READ(sdbResult, NODE):
    # NODE, ELB, DISP, RR, HeatUniformGenRate, EnergyLossAccum, VisEnergy, StrainLongTerm, Deformed_RIM_ROAD
    resultsfile = open(sdbResult, 'rb')

    HeatUniformGenRate=[]; HeatHourGlassGenRate=[]
    EnergyLossAccum = []
    VisEnergy = []
    StrainLongTerm = []; StrainLongTermHourglass=[]
    Deformed_RIM_ROAD = []
    Temperature = []
    
    resultsfile.seek(0,2); fend = resultsfile.tell() #find the end position of binary file
    resultsfile.seek(0);   fpos = resultsfile.tell() #find the start position of binary file
    # print fpos, fend
    ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RecordHeaderID = 0
    READ_NUM       = 0
    NodeNUM        = 0
    TREAD_ELMENTNUM= 0
    
    tmpELB = [] 
    

    DISP = []
    while resultsfile.tell() < fend:
        RecordHeaderID     = struct.unpack('i', resultsfile.read(4))[0]
        if RecordHeaderID == 13:
            RecordValue    = struct.unpack('i', resultsfile.read(4))[0]
            # print"RecordHeaderID = ", RecordHeaderID, "Record Value = ", RecordValue
            # print ("RecordHeaderID", RecordHeaderID, ", Record Value", RecordValue)
            OutputStepTime = struct.unpack('d', resultsfile.read(8))[0]
            print ("SDB Result Time=%.2f"%OutputStepTime)
            # printOutputStepTime
            OutputStepID       = struct.unpack('i', resultsfile.read(4))[0]
            # print ("Output step ID", OutputStepID)
            # printOutputStepID
            OutputStepName = ''
            i = 0 
            while i < 26:
                OutputStepName     = OutputStepName + str(struct.unpack('c', resultsfile.read(1))[0])
                i = i + 1
            # print"OutputStepName = ", OutputStepName
            # print ("output step name",OutputStepName )
            OutputStepNo      = struct.unpack('i', resultsfile.read(4))[0]
            # printOutputStepNo
            OutputStepNode    = struct.unpack('i', resultsfile.read(4))[0]
            # printOutputStepNode
            NodeNUM           = struct.unpack('i', resultsfile.read(4))[0]
            # print"NodeNUM = ", NodeNUM
            READ_NUM     = NodeNUM
            i = 0
            while i < READ_NUM:
                NodeID = int(struct.unpack('i', resultsfile.read(4))[0])
                DISP.append([0.0,0.0,0.0]) #for making of abaqus odb file
                i = i + 1
            # print"finish to read the Node Lists!"
        elif RecordHeaderID == 14:
            RecordValue    = struct.unpack('i', resultsfile.read(4))[0]
            # print ("RecordHeaderID = ", RecordHeaderID, "Record Value = ", RecordValue)

            IDATA_TYPE     = struct.unpack('i', resultsfile.read(4))[0]
            IOUT_TYPE      = struct.unpack('i', resultsfile.read(4))[0]
            ICOMP1         = struct.unpack('i', resultsfile.read(4))[0]
            ICOMP2         = struct.unpack('i', resultsfile.read(4))[0]
            ICOMP3         = struct.unpack('i', resultsfile.read(4))[0]
            ICOMP4         = struct.unpack('i', resultsfile.read(4))[0]
            ICOMP5         = struct.unpack('i', resultsfile.read(4))[0]
            ICOMP6         = struct.unpack('i', resultsfile.read(4))[0]
            ICOMP7         = struct.unpack('i', resultsfile.read(4))[0]
            ICOMP8         = struct.unpack('i', resultsfile.read(4))[0]
            IDIRECTION     = struct.unpack('i', resultsfile.read(4))[0]
            ICALC_WARN     = struct.unpack('i', resultsfile.read(4))[0]
            IDMINV         = struct.unpack('i', resultsfile.read(4))[0]
            IDMAXV         = struct.unpack('i', resultsfile.read(4))[0]
            VMIN           = struct.unpack('d', resultsfile.read(8))[0]
            VMAX           = struct.unpack('d', resultsfile.read(8))[0]
            DATA_TITLE     = ''
            i = 0 
            while i < 26:
                DATA_TITLE = DATA_TITLE + str(struct.unpack('c', resultsfile.read(1))[0])
                i = i + 1
            # printDATA_TITLE
            i = 0
            if RecordValue == 102:
                while i < READ_NUM:
                    Value = struct.unpack('d', resultsfile.read(8))[0]
                    DISP[i][0] = Value #for making of abaqus odb file
                    NODE[i][1] = NODE[i][1] + Value
                    i = i + 1
            if RecordValue == 103:
                while i < READ_NUM:
                    Value = struct.unpack('d', resultsfile.read(8))[0]
                    DISP[i][1] = Value #for making of abaqus odb file
                    NODE[i][2] = NODE[i][2] + Value
                    i = i + 1
            if RecordValue == 104:
                while i < READ_NUM:
                    Value = struct.unpack('d', resultsfile.read(8))[0]
                    DISP[i][2] = Value #for making of abaqus odb file
                    NODE[i][3] = NODE[i][3] + Value
                    i = i + 1
            if RecordValue == 3101:
                for i in range(READ_NUM):
                    Value = struct.unpack('d', resultsfile.read(8))[0]
                    NODE[i][4] = NODE[i][4] + Value
            elif RecordValue == 900000054:
                while i < READ_NUM:
                    Value = struct.unpack('d', resultsfile.read(8))[0]
                    i = i + 1
            elif RecordValue == 900001001: # Heat Uniform Generate Rate
                while i < READ_NUM:
                    Value = struct.unpack('d', resultsfile.read(8))[0]
                    HeatUniformGenRate.append(Value)
                    i = i + 1
            elif RecordValue == 900001002: # Heat Uniform Generate Rate
                while i < READ_NUM:
                    Value = struct.unpack('d', resultsfile.read(8))[0]
                    HeatHourGlassGenRate.append(Value)
                    i = i + 1
            elif RecordValue == 900001003: # Energy Uniform Loss Accumulation
                print ("Energy Uniform Loss Accumulation")
                while i < READ_NUM:
                    Value = struct.unpack('d', resultsfile.read(8))[0]
                    EnergyLossAccum.append(Value)
                    i = i + 1
            elif RecordValue == 900001005: # Vis Energy In Uniform Storage
                print ("Vis Energy In Uniform Storage")
                while i < READ_NUM:
                    Value = struct.unpack('d', resultsfile.read(8))[0]
                    VisEnergy.append(Value)
                    ELB.append(tmpELB[i])
                    RR.append([tmpELB[i],Value])
                    i = i + 1
                del tmpELB
                tmpELB = []
            elif RecordValue == 900001007: # Strain Uniform Long Term E
                while i < READ_NUM:
                    Value = struct.unpack('d', resultsfile.read(8))[0]
                    StrainLongTerm.append(Value)
                    i = i + 1
            elif RecordValue == 900001008: # Strain Uniform Long Term E Hourglass
                while i < READ_NUM:
                    Value = struct.unpack('d', resultsfile.read(8))[0]
                    StrainLongTermHourglass.append(Value)
                    i = i + 1
            else:
                 
                while i < READ_NUM:
                    Value = struct.unpack('d', resultsfile.read(8))[0]
                    i = i + 1
        elif RecordHeaderID == 2:
            ELMENTNUM = struct.unpack('i', resultsfile.read(4))[0]
            i = 0
            while i < ELMENTNUM:
                TreadID = struct.unpack('i', resultsfile.read(4))[0]
                tmpELB.append(TreadID)
                i = i + 1
            READ_NUM     = ELMENTNUM
        elif RecordHeaderID == 1:
            ELMENT_NUM   = struct.unpack('i', resultsfile.read(4))[0]
            i = 0
            while i < ELMENT_NUM:
                ID = struct.unpack('i', resultsfile.read(4))[0]
                Deformed_RIM_ROAD.append([ID, 0.0, 0.0, 0.0])
                i = i + 1
            READ_NUM     = ELMENT_NUM
            while resultsfile.tell() < fend:
                RecordHeaderID     = struct.unpack('i', resultsfile.read(4))[0]
                RecordValue    = struct.unpack('i', resultsfile.read(4))[0]
                IDATA_TYPE     = struct.unpack('i', resultsfile.read(4))[0]
                IOUT_TYPE      = struct.unpack('i', resultsfile.read(4))[0]
                ICOMP1         = struct.unpack('i', resultsfile.read(4))[0]
                ICOMP2         = struct.unpack('i', resultsfile.read(4))[0]
                ICOMP3         = struct.unpack('i', resultsfile.read(4))[0]
                ICOMP4         = struct.unpack('i', resultsfile.read(4))[0]
                ICOMP5         = struct.unpack('i', resultsfile.read(4))[0]
                ICOMP6         = struct.unpack('i', resultsfile.read(4))[0]
                ICOMP7         = struct.unpack('i', resultsfile.read(4))[0]
                ICOMP8         = struct.unpack('i', resultsfile.read(4))[0]
                IDIRECTION     = struct.unpack('i', resultsfile.read(4))[0]
                ICALC_WARN     = struct.unpack('i', resultsfile.read(4))[0]
                IDMINV         = struct.unpack('i', resultsfile.read(4))[0]
                IDMAXV         = struct.unpack('i', resultsfile.read(4))[0]
                VMIN           = struct.unpack('d', resultsfile.read(8))[0]
                VMAX           = struct.unpack('d', resultsfile.read(8))[0]
                DATA_TITLE     = ''
                i = 0 
                while i < 26:
                    DATA_TITLE = DATA_TITLE + str(struct.unpack('c', resultsfile.read(1))[0])
                    i = i + 1
                i = 0
                if RecordValue == 102:
                    while i < READ_NUM:
                        Value = struct.unpack('d', resultsfile.read(8))[0]
                        Deformed_RIM_ROAD[i][1] = Value
                        i = i + 1
                if RecordValue == 103:
                    while i < READ_NUM:
                        Value = struct.unpack('d', resultsfile.read(8))[0]
                        Deformed_RIM_ROAD[i][2] = Value
                        i = i + 1
                if RecordValue == 104:
                    while i < READ_NUM:
                        Value = struct.unpack('d', resultsfile.read(8))[0]
                        Deformed_RIM_ROAD[i][3] = Value
                        i = i + 1
                if RecordValue == 106:
                    while i < READ_NUM:
                        Value = struct.unpack('d', resultsfile.read(8))[0]
                        i = i + 1
                if RecordValue == 107:
                    while i < READ_NUM:
                        Value = struct.unpack('d', resultsfile.read(8))[0]
                        i = i + 1
                if RecordValue == 108:
                    while i < READ_NUM:
                        Value = struct.unpack('d', resultsfile.read(8))[0]
                        i = i + 1
                else:
                    while i < READ_NUM:
                        Value = struct.unpack('d', resultsfile.read(8))[0]
                        i = i + 1
        else:
            break

    del tmpELB
    resultsfile.close()
    
    HeatGen = [HeatUniformGenRate, HeatHourGlassGenRate]
    iSED = [StrainLongTerm, StrainLongTermHourglass]

    print ("SDB results are loaded")
    
    
    if len(HeatGen) > 0:
        return np.array(NODE), Deformed_RIM_ROAD, HeatGen, iSED
    else:
        return -1

def LayoutMesh_From_axi(axi="", limit=10000, output="" ):
    print (" AXI FILE : %s"%(axi))
    with open(axi) as I: 
        lines = I.readlines()
    fp = open(output, 'w')
    fp.write("**************************************\n")
    fp.write("** TIRE MESH from AXI to debug P3DM\n")
    fp.write("**************************************\n")

    SectorNodes = 0;     AllNodes = 0;     Node_Sectors = 0; 
    SectorEls = 0; AllElements = 0; EL_Sectors = 0 

    nodes=[]
    elements=[]
    

    cmd=None 
    for line in lines: 
        if "**" in line  :
            continue 
        if "*" in line: 

            if "*NODE" in line.upper(): 
                cmd ="ND"
                fp.write("*NODE, SYSTEM=R\n")
            elif "*ELEMENT" in line.upper() and "M3D4" in line.upper(): 
                cmd = 'RB'
                fp.write("*ELEMENT, TYPE=MGAX1\n")
            elif "*ELEMENT" in line.upper() and "C3D6" in line.upper(): 
                cmd = 'C6'
                fp.write("*ELEMENT, TYPE=CGAX3H\n")
            elif "*ELEMENT" in line.upper() and "C3D8" in line.upper(): 
                cmd = 'C8'
                fp.write("*ELEMENT, TYPE=CGAX4H\n")
            elif "*ELSET," in line.upper() and "ELSET=" in line.upper(): 
                cmd = 'ES'
                fp.write(line)
            elif "*SURFACE," in line.upper() and "NAME=PRESS" in line.upper(): 
                cmd = 'PS'
                fp.write(line)
            elif "*SURFACE," in line.upper() and "NAME=RIC_R" in line.upper(): 
                cmd = 'PR'
                fp.write(line)
            elif "*SURFACE," in line.upper() and "NAME=RIC_L" in line.upper(): 
                cmd = 'PL'
                fp.write(line)
            elif "*SURFACE," in line.upper() and "NAME=TIREBODY" in line.upper(): 
                cmd = 'BD'
                fp.write(line)
            elif "*SURFACE," in line.upper() and ("NAME=TIE_M" in line.upper() or ("NAME=" in line.upper() and "_TIE" in line.upper())): 
                cmd = 'TM'
                fp.write(line)
            elif "*SURFACE," in line.upper() and ("NAME=TIE_S" in line.upper() or ("NAME=" in line.upper() and "_TIE" in line.upper())):  
                cmd = 'TS'
                fp.write(line)
            elif "*TIE," in line.upper() : 
                cmd = 'TD'
                fp.write(line)
            elif "NIDOFFSET" in line.upper(): 
                cmd = None 
            else: 
                cmd = None 

            # print ("%s, %s"%(line.strip(), cmd))

        else: 
            # SectorNodes = 0;     AllNodes = 0;     Node_Sectors = 0; 
            # SectorEls = 0; AllElements = 0; EL_Sectors = 0 

            if cmd=="ND" : 
                data = line.split(",")
                if int(data[0]) < limit: 
                    fp.write("%s, %s, %s, %s\n"%(data[0], data[3].strip(), data[2], data[1]))
                    SectorNodes += 1 
                AllNodes += 1 
                nodes.append([int(data[0].strip()), float(data[1].strip()), float(data[2].strip()), float(data[3].strip())])
            elif cmd =="RB": 
                data = line.split(",")
                if int(data[0]) < limit:
                    fp.write("%s, %s, %s\n"%(data[0], data[1], data[2].strip()))
                    SectorEls+= 1
                AllElements += 1  
                elements.append([int(data[0].strip()), int(data[1].strip()), int(data[2].strip()), int(data[3].strip()), int(data[4].strip()), 0, 0, 0, 0])  
            elif cmd =="C6": 
                data = line.split(",")
                if int(data[0]) < limit:
                    fp.write("%s, %s, %s, %s\n"%(data[0], data[4], data[5].strip(), data[6].strip()))
                    SectorEls+= 1
                AllElements += 1  

                elements.append([int(data[0].strip()), int(data[1].strip()), int(data[2].strip()), int(data[3].strip()), int(data[4].strip()), int(data[5].strip()), int(data[6].strip()), 0, 0])
            elif cmd =="C8": 
                data = line.split(",")
                if int(data[0]) < limit:
                    fp.write("%s, %s, %s, %s, %s\n"%(data[0], data[5], data[6].strip(), data[7].strip(), data[8].strip()))
                    SectorEls+= 1
                AllElements += 1  
                elements.append([int(data[0].strip()), int(data[1].strip()), int(data[2].strip()), int(data[3].strip()), int(data[4].strip()), int(data[5].strip()), int(data[6].strip()), int(data[7].strip()), int(data[8].strip())])
            elif cmd =="ES": 
                data = line.split(",")
                cnt = 0 
                for d in data: 
                    d = d.strip()
                    if d != "": 
                        if int(d) < limit: 
                            cnt += 1
                            fp.write("%s,"%(d))
                if cnt > 0: fp.write("\n")

            elif cmd =="PS": 
                data = line.split(",")
                if int(data[0]) < limit: 
                    if "S3" in data[1]: fp.write("%s, S1\n"%(data[0]))
                    if "S4" in data[1]: fp.write("%s, S2\n"%(data[0]))
                    if "S5" in data[1]: fp.write("%s, S3\n"%(data[0]))
                    if "S6" in data[1]: fp.write("%s, S4\n"%(data[0]))
            elif cmd =="PR": 
                data = line.split(",")
                if int(data[0]) < limit: 
                    if "S3" in data[1]: fp.write("%s, S1\n"%(data[0]))
                    if "S4" in data[1]: fp.write("%s, S2\n"%(data[0]))
                    if "S5" in data[1]: fp.write("%s, S3\n"%(data[0]))
                    if "S6" in data[1]: fp.write("%s, S4\n"%(data[0]))
            elif cmd =="PL": 
                data = line.split(",")
                if int(data[0]) < limit: 
                    if "S3" in data[1]: fp.write("%s, S1\n"%(data[0]))
                    if "S4" in data[1]: fp.write("%s, S2\n"%(data[0]))
                    if "S5" in data[1]: fp.write("%s, S3\n"%(data[0]))
                    if "S6" in data[1]: fp.write("%s, S4\n"%(data[0]))
            elif cmd =="BD": 
                data = line.split(",")
                if int(data[0]) < limit: 
                    if "S3" in data[1]: fp.write("%s, S1\n"%(data[0]))
                    if "S4" in data[1]: fp.write("%s, S2\n"%(data[0]))
                    if "S5" in data[1]: fp.write("%s, S3\n"%(data[0]))
                    if "S6" in data[1]: fp.write("%s, S4\n"%(data[0]))
            elif cmd =="TM": 
                data = line.split(",")
                if int(data[0]) < limit: 
                    if "S3" in data[1]: fp.write("%s, S1\n"%(data[0]))
                    if "S4" in data[1]: fp.write("%s, S2\n"%(data[0]))
                    if "S5" in data[1]: fp.write("%s, S3\n"%(data[0]))
                    if "S6" in data[1]: fp.write("%s, S4\n"%(data[0]))
            elif cmd =="TS": 
                data = line.split(",")
                if int(data[0]) < limit: 
                    if "S3" in data[1]: fp.write("%s, S1\n"%(data[0]))
                    if "S4" in data[1]: fp.write("%s, S2\n"%(data[0]))
                    if "S5" in data[1]: fp.write("%s, S3\n"%(data[0]))
                    if "S6" in data[1]: fp.write("%s, S4\n"%(data[0]))
            elif cmd =="TD": 
                fp.write(line)
            else: 
                continue 
            
    fp.close()

    Node_Sectors = AllNodes / SectorNodes 
    EL_Sectors = AllElements / SectorEls 

    print ("\n ** Sector No. check ")
    print (" All Nodes = %d\n Nodes per Sector = %d\n Sectors = %.1f"%(AllNodes, SectorNodes, Node_Sectors))
    print (" All ELs   = %d\n ELs per Sector   = %d\n Sectors = %.1f"%(AllElements, SectorEls, EL_Sectors))

    nodes = np.array(nodes)
    elements = np.array(elements)


    idx = np.where(nodes[:,0]<limit)[0]
    sectornodes = nodes[idx]

    idx = np.where(elements[:,0]<limit)[0]
    s1els = elements[idx]

    sectornodes = nodes 
    return sectornodes 

def PatternMesh_from_trd(trd="", limit=10000, output="", treadNo=10**7 ):

    print (" writing temporary file", output)
    with open(trd) as I: 
        lines = I.readlines()
    fp = open(output, 'w')
    fp.write("**************************************\n")
    fp.write("** TIRE MESH from TRD to debug P3DM\n")
    fp.write("**************************************\n")

    SectorNodes = 0;     AllNodes = 0;     Node_Sectors = 0; 
    SectorEls = 0; AllElements = 0; EL_Sectors = 0 

    nodes=[]
    elements=[]

    nodeIds=[]

    cmd=None 
    for line in lines: 
        if "**" in line  :
            continue 
        if "*" in line: 

            if "*NODE" in line.upper(): 
                cmd ="ND"
                fp.write("*NODE, SYSTEM=R\n")
            elif "*ELEMENT" in line.upper() and "M3D4" in line.upper(): 
                cmd = 'RB'
                fp.write("*ELEMENT, TYPE=MGAX1\n")
            elif "*ELEMENT" in line.upper() and "C3D6" in line.upper(): 
                cmd = 'C6'
                fp.write("*ELEMENT, TYPE=CGAX3H\n")
            elif "*ELEMENT" in line.upper() and "C3D8" in line.upper(): 
                cmd = 'C8'
                fp.write("*ELEMENT, TYPE=CGAX4H\n")
            elif "*ELSET," in line.upper() and "ELSET=" in line.upper(): 
                cmd = 'ES'
                fp.write(line)
            elif "*SURFACE," in line.upper() and "NAME=XTRD" in line.upper(): 
                cmd = 'PS'
                fp.write(line)
            elif "*SURFACE," in line.upper() and "NAME=YTIE" in line.upper(): 
                cmd = 'PR'
                fp.write(line)
            elif "*SURFACE," in line.upper() and "NAME=CONT" in line.upper(): 
                cmd = 'PL'
                fp.write(line)
            elif "*SURFACE," in line.upper() and "NAME=TIREBODY" in line.upper(): 
                cmd = 'BD'
                fp.write(line)
            elif "*SURFACE," in line.upper() and ("NAME=TIE_M" in line.upper() or ("NAME=" in line.upper() and "_TIE" in line.upper())): 
                cmd = 'TM'
                fp.write(line)
            elif "*SURFACE," in line.upper() and ("NAME=TIE_S" in line.upper() or ("NAME=" in line.upper() and "_TIE" in line.upper())):  
                cmd = 'TS'
                fp.write(line)
            elif "*TIE," in line.upper() : 
                cmd = 'TD'
                fp.write(line)
            elif "NIDOFFSET" in line.upper(): 
                cmd = None 
            elif "TREADPTN_NIDSTART_NIDOFFSET_EIDSTART_EIDOFFSET" in line.upper(): 
                cmd = None
                data = line.split("=")[1]
                limit = data.split(",")[1]
                limit = int(limit)
                treadNo = int(data.split(",")[0])
            else: 
                cmd = None 

            # print ("%s, %s"%(line.strip(), cmd))

        else: 

            if cmd=="ND" : 
                data = line.split(",")
                if int(data[0]) < limit*2+treadNo: 
                    fp.write("%s, %s, %s, %s\n"%(data[0], data[3].strip(), data[2], data[1].strip()))
                    SectorNodes += 1 
                AllNodes += 1 
                nodes.append([int(data[0].strip()), float(data[1].strip()), float(data[2].strip()), float(data[3].strip())])
                nodeIds.append(int(data[0].strip()))
            elif cmd =="RB": 
                data = line.split(",")
                if int(data[0]) < limit+treadNo:
                    fp.write("%s, %s, %s\n"%(data[0], data[1], data[2].strip()))
                    SectorEls+= 1
                AllElements += 1  
                elements.append([int(data[0].strip()), int(data[1].strip()), int(data[2].strip()), int(data[3].strip()), int(data[4].strip()), 0, 0, 0, 0])  
            elif cmd =="C6": 
                data = line.split(",")
                if int(data[0]) < limit+treadNo:
                    fp.write("%s, %s, %s, %s\n"%(data[0], data[4], data[5].strip(), data[6].strip()))
                    SectorEls+= 1
                AllElements += 1  

                elements.append([int(data[0].strip()), int(data[1].strip()), int(data[2].strip()), int(data[3].strip()), int(data[4].strip()), int(data[5].strip()), int(data[6].strip()), 0, 0])
            elif cmd =="C8": 
                data = line.split(",")
                if int(data[0]) < limit+treadNo:
                    fp.write("%s, %s, %s, %s, %s\n"%(data[0], data[5], data[6].strip(), data[7].strip(), data[8].strip()))
                    SectorEls+= 1
                AllElements += 1  
                elements.append([int(data[0].strip()), int(data[1].strip()), int(data[2].strip()), int(data[3].strip()), int(data[4].strip()), int(data[5].strip()), int(data[6].strip()), int(data[7].strip()), int(data[8].strip())])
            elif cmd =="ES": 
                data = line.split(",")
                cnt = 0 
                for d in data: 
                    d = d.strip()
                    if d !='': 
                        if int(d) < limit+treadNo: 
                            cnt += 1
                            fp.write("%s,"%(d))
                if cnt > 0: fp.write("\n")

            elif cmd =="PS": 
                data = line.split(",")
                if int(data[0]) < limit+treadNo: 
                    if "S3" in data[1]: fp.write("%s, S1\n"%(data[0]))
                    if "S4" in data[1]: fp.write("%s, S2\n"%(data[0]))
                    if "S5" in data[1]: fp.write("%s, S3\n"%(data[0]))
                    if "S6" in data[1]: fp.write("%s, S4\n"%(data[0]))
            elif cmd =="PR": 
                data = line.split(",")
                if int(data[0]) < limit+treadNo: 
                    if "S3" in data[1]: fp.write("%s, S1\n"%(data[0]))
                    if "S4" in data[1]: fp.write("%s, S2\n"%(data[0]))
                    if "S5" in data[1]: fp.write("%s, S3\n"%(data[0]))
                    if "S6" in data[1]: fp.write("%s, S4\n"%(data[0]))
            elif cmd =="PL": 
                data = line.split(",")
                if int(data[0]) < limit+treadNo: 
                    if "S3" in data[1]: fp.write("%s, S1\n"%(data[0]))
                    if "S4" in data[1]: fp.write("%s, S2\n"%(data[0]))
                    if "S5" in data[1]: fp.write("%s, S3\n"%(data[0]))
                    if "S6" in data[1]: fp.write("%s, S4\n"%(data[0]))
            elif cmd =="BD": 
                data = line.split(",")
                if int(data[0]) < limit+treadNo: 
                    if "S3" in data[1]: fp.write("%s, S1\n"%(data[0]))
                    if "S4" in data[1]: fp.write("%s, S2\n"%(data[0]))
                    if "S5" in data[1]: fp.write("%s, S3\n"%(data[0]))
                    if "S6" in data[1]: fp.write("%s, S4\n"%(data[0]))
            elif cmd =="TM": 
                data = line.split(",")
                if int(data[0]) < limit+treadNo: 
                    if "S3" in data[1]: fp.write("%s, S1\n"%(data[0]))
                    if "S4" in data[1]: fp.write("%s, S2\n"%(data[0]))
                    if "S5" in data[1]: fp.write("%s, S3\n"%(data[0]))
                    if "S6" in data[1]: fp.write("%s, S4\n"%(data[0]))
            elif cmd =="TS": 
                data = line.split(",")
                if int(data[0]) < limit+treadNo: 
                    if "S3" in data[1]: fp.write("%s, S1\n"%(data[0]))
                    if "S4" in data[1]: fp.write("%s, S2\n"%(data[0]))
                    if "S5" in data[1]: fp.write("%s, S3\n"%(data[0]))
                    if "S6" in data[1]: fp.write("%s, S4\n"%(data[0]))
            elif cmd =="TD": 
                fp.write(line)
            
            else: 
                continue 
            
    fp.close()

    trdOffset = findout_offset(nodeIds, step=limit, shift=-treadNo)
    print ("# Tread Offset=", trdOffset)

    

    Node_Sectors = AllNodes / SectorNodes 
    EL_Sectors = AllElements / SectorEls 


    print ("\n ** Sector No. check ")
    print (" All Nodes = %d\n Nodes per Sector = %d\n Sectors = %.1f"%(AllNodes, SectorNodes, Node_Sectors))
    print (" All ELs   = %d\n ELs per Sector   = %d\n Sectors = %.1f"%(AllElements, SectorEls, EL_Sectors))

    if Node_Sectors >= 180: ptn = False 
    else: ptn = True 

    

    nodes = np.array(nodes)
    elements = np.array(elements)

    write_Pattern3d_2dMesh(offset=limit, output=output, npn=nodes, els=elements, surfs=None)

    idx = np.where(nodes[:,0]<limit+treadNo)[0]
    s1nodes = nodes[idx]

    idx = np.where(elements[:,0]<limit+treadNo)[0]
    s1els = elements[idx]

    ix = np.where(nodes[:,0]<limit+treadNo)[0]
    pitchnodes = nodes[ix]

    pitchnodes = nodes 
    return pitchnodes, Node_Sectors 

def meshgrid_in_Quadrilateral(xs, ys, num=10, endpoint=True, vs=None): 
    # xs=[x1, x2, x3, x4]; ys=[y1, y2, y3, y4]
    s = np.linspace(-1, 1, num)
    t = np.linspace(-1, 1, num)
    if not endpoint: 
        s = np.delete(s, 0, axis=0); s=np.delete(s, -1, axis=0)
        t = np.delete(s, 0, axis=0); t=np.delete(s, -1, axis=0)
    if isinstance(vs, type(None)): 
        px=[]; py=[]
        for m in s: 
            tx=[]; ty=[]
            for n in t: 
                tx.append( 0.25*((1-m)*(1-n)*xs[0] + (1+m)*(1-n)*xs[1] + (1+m)*(1+n)*xs[2] + (1-m)*(1+n)*xs[3]) )
                ty.append( 0.25*((1-m)*(1-n)*ys[0] + (1+m)*(1-n)*ys[1] + (1+m)*(1+n)*ys[2] + (1-m)*(1+n)*ys[3]) )
            px.append(tx)
            py.append(ty)

        return np.array(px), np.array(py)
    else: 
        px=[]; py=[]; vy=[]
        for m in s: 
            tx=[]; ty=[]; tv=[]
            for n in t: 
                tx.append( 0.25*((1-m)*(1-n)*xs[0] + (1+m)*(1-n)*xs[1] + (1+m)*(1+n)*xs[2] + (1-m)*(1+n)*xs[3]) )
                ty.append( 0.25*((1-m)*(1-n)*ys[0] + (1+m)*(1-n)*ys[1] + (1+m)*(1+n)*ys[2] + (1-m)*(1+n)*ys[3]) )
                tv.append( 0.25*((1-m)*(1-n)*vs[0] + (1+m)*(1-n)*vs[1] + (1+m)*(1+n)*vs[2] + (1-m)*(1+n)*vs[3]) )
            px.append(tx)
            py.append(ty)
            vy.append(tv)

        return np.array(px), np.array(py), np.array(vy)

def Area_polygon(nodes, xy=23): 
    ii = int(xy/10)
    jj = int(xy%10)

    x=[]; y=[]
    for n in nodes: 
        x.append(n[2])
        y.append(n[3])
    x.append(x[0])
    y.append(y[0])
    
    A=np.zeros(3)
    N = len(nodes)
    for i in range(N): 
        s = x[i] * y[i + 1] - x[i + 1] * y[i]
        A[0] += s
        A[1] += (x[i] + x[i + 1]) * s
        A[2] += (y[i] + y[i + 1]) * s
    
    if round(A, 10) ==0: 
        return A 
    A[0] /= 2.0
    A[1] /=(A[0]*6)
    A[2] /=(A[0]*6)

    return A 

def ElementNodalArea(e, npn, dist=False): 

    ix = np.where(npn[:,0]==e[1])[0][0]; n1=npn[ix]
    ix = np.where(npn[:,0]==e[2])[0][0]; n2=npn[ix]
    ix = np.where(npn[:,0]==e[3])[0][0]; n3=npn[ix]
    n5=[0, 0, (n1[2]+n2[2])/2.0, (n1[3]+n2[3])/2.0]
    n6=[0, 0, (n2[2]+n3[2])/2.0, (n2[3]+n3[3])/2.0]
    if e[4] > 0: 
        ix = np.where(npn[:,0]==e[4])[0][0]; n4=npn[ix]
        n7=[0, 0, (n3[2]+n4[2])/2.0, (n3[3]+n4[3])/2.0]
        n8=[0, 0, (n4[2]+n1[2])/2.0, (n4[3]+n1[3])/2.0]
        n9=[0, 0, (n1[2]+n2[2]+n3[2]+n4[2])/4.0, (n1[3]+n2[3]+n3[3]+n4[3])/4.0]

        a1 = Area_polygon([n1,n5,n9,n8])[0]
        a2 = Area_polygon([n2,n6,n9,n5])[0]
        a3 = Area_polygon([n3,n7,n9,n6])[0]
        a4 = Area_polygon([n4,n8,n9,n7])[0]

        if dist: 
            d1 = math.sqrt( (n1[2]-n9[2])**2 + (n1[3]-n9[3])**2  )
            d2 = math.sqrt( (n2[2]-n9[2])**2 + (n2[3]-n9[3])**2  )
            d3 = math.sqrt( (n3[2]-n9[2])**2 + (n3[3]-n9[3])**2  )
            d4 = math.sqrt( (n4[2]-n9[2])**2 + (n4[3]-n9[3])**2  )
    else: 
        n7=[0, 0, (n3[2]+n1[2])/2.0, (n3[3]+n1[3])/2.0]
        n9=[0, 0, (n1[2]+n2[2]+n3[2])/3.0, (n1[3]+n2[3]+n3[3])/3.0]

        a1 = Area_polygon([n1,n5,n9,n7])[0]
        a2 = Area_polygon([n2,n6,n9,n5])[0]
        a3 = Area_polygon([n3,n7,n9,n6])[0]
        a4 = 0
        if dist: 
            d1 = math.sqrt( (n1[2]-n9[2])**2 + (n1[3]-n9[3])**2  )
            d2 = math.sqrt( (n2[2]-n9[2])**2 + (n2[3]-n9[3])**2  )
            d3 = math.sqrt( (n3[2]-n9[2])**2 + (n3[3]-n9[3])**2  )
            d4 = 0
    if dist: 
        return [a1, a2, a3, a4], [d1, d2, d3, d4]
    return [a1, a2, a3, a4]

def nodes_Master_Slave_Tie(ties, surfs, els, name=None): 
    nmaster=[]; nslaves=[]
    
    for tie in ties: 
        if name: 
            if name != tie[0]: continue 

        for surf in surfs: 
            if surf[0] == tie[1][1]: 
                mst = surf 
                break 
        for surf in surfs: 
            if surf[0] == tie[1][0]: 
                slvs = surf[1:] 
                break

        for e in els: 
            if e[0] == mst[1][0]: 
                if mst[1][1] =='S1': 
                    nmaster.append(e[1]); nmaster.append(e[2])
                elif mst[1][1] =='S2': 
                    nmaster.append(e[2]); nmaster.append(e[3])
                else: 
                    if e[4]: 
                        if mst[1][1] =='S3': 
                            nmaster.append(e[3]); nmaster.append(e[4])
                        else:
                            nmaster.append(e[1]); nmaster.append(e[4])
                    else: 
                        nmaster.append(e[1]); nmaster.append(e[3])

                break 
        for slv in slvs: 
            for e in els: 
                if e[0] == slv[0]: 
                    if slv[1] == 'S1': 
                        nslaves.append(e[1])
                        nslaves.append(e[2])
                    elif slv[1] == 'S2': 
                        nslaves.append(e[2])
                        nslaves.append(e[3])
                    else: 
                        if e[4]: 
                            if slv[1] == 'S3': 
                                nslaves.append(e[3])
                                nslaves.append(e[4])
                            else: 
                                nslaves.append(e[4])
                                nslaves.append(e[1])
                        else: 
                            nslaves.append(e[3])
                            nslaves.append(e[1])

    nmst = np.array(nmaster)
    nslv = np.array(nslaves)
    return np.unique(nslv), np.unique(nmst)
            
def findout_offset(numbers, step=10000, shift=0):
    nn =[]
    for n in numbers: 
        nn.append(n+shift)
    numbers = np.array(nn)
    imax = np.max(numbers) 
    # print (numbers, imax, shift)

    if imax < step: 
        return step
    offset=0
    steps = int(imax / step) 
    nlast = imax - steps*step 

    idxs = np.where(numbers>steps*step)[0]
    counting = len(idxs)

    # print (" Last counting", counting, steps, nlast, steps*step )

    for i in range(1, 100): 
        ix1 = np.where(numbers[:]>step*(i-1))[0]
        ix2 = np.where(numbers[:]<=step*i)[0]
        ix = np.intersect1d(ix1, ix2)
        # if i==1: print ("~%d : counting %d"%(step*i, len(ix)))
        
        if len(ix) == counting: 
            offset = i *step 
            break 
    # print (i)
    return offset 

def write_Pattern3d_2dMesh(offset=10000, output='trd.tmp', npn=None, els=None, surfs=None): 
    treadNo = 10**7 
    limit = offset 

    nodeIds = npn[:,0]
    limit = findout_offset(nodeIds, step=limit, shift=-treadNo)

    print ("# TREAD OFFSET ", limit)

    fp = open(output, 'w')
    fp.write("**************************************\n")
    fp.write("** TIRE MESH from TRD to debug P3DM\n")
    fp.write("**************************************\n")
    fp.write("*NODE, SYSTEM=R\n")
    for n in npn: 
        if n[0] < limit*2+treadNo: 
            r =  math.sqrt( n[3]**2 + n[1]**2) 
            fp.write("%10d, %12.7f, %12.7f, %12.7f\n"%(n[0],r, n[2], 0))

    fp.write("*ELEMENT, TYPE=CGAX4H\n")
    c3 =False 
    cnt = 0 
    for el in els: 
        cnt += 1 
        if el[0] < limit+treadNo and el[7] > 0:
            if cnt %2 : 
                fp.write("%10d, %10d, %10d, %10d, %10d\n"%(el[0], el[1], el[2], el[6], el[5]))
            else: 
                fp.write("%10d, %10d, %10d, %10d, %10d\n"%(el[0]+limit, el[2], el[3], el[7], el[6]))
        else: 
            c3 = True    
    if c3: 
        for el in els: 
            if el[0] < limit+treadNo and el[7] == 0:
                fp.write("%10d, %10d, %10d, %10d, %10d\n"%(el[0], el[1], el[2], el[5], el[4]))
                # fp.write("%10d, %10d, %10d, %10d, %10d\n"%(el[0]+limit, el[2], el[3], el[6], el[5]))
        

# ElementCenterValueToInnerValues
def Distribution_From_CenterValue(npn, elements, elsets, vmin=0.0, vmax=1.0, tie=None, surface=None): 
    # els =[]
    # for el in elements: 


    ########################################################
    ## Check if the node is a node on slave edge or not
    ## then calculate nodal value by distance
    nodes_slave, nodes_master = nodes_Master_Slave_Tie(tie, surface, npel)
    ########################################################
    ## Nodal Area for Calculation for all elements 
    el_area =[]
    for e in npel: 
        areas, dists = ElementNodalArea(e, npn, dist=True)

    for elset in elsets: ## elset =[name, e1, e2,...]
        neg=[]; pos=[]
        memb = False 
        for i in range(1, len(elset)): 
            ix = np.where(npel[:,0]== elset[i])[0][0]; e = npel[ix]
            if e[3] ==0: 
                memb = True
                break 
            ix = np.where(npn[:,0]==e[1])[0][0]; n1=npn[ix]
            if npn[ix][2] >=0: pos.append(e)
            else: neg.append(e)
        if memb : 
            continue 
        
        esets =[pos, neg]
        for eset in esets: 
            px=[]; py=[]; pv=[]
            cnt = 0 
            for i in eset: 
                ix = np.where(npel[:,0]== i)[0][0]; e = npel[ix]

                ndx1 = np.where(nodes_slave==e[1])[0] 
                ndx2 = np.where(nodes_slave==e[2])[0] 
                ndx3 = np.where(nodes_slave==e[3])[0] 
                if e[4]>0: 
                    ndx4 = np.where(nodes_slave==e[4])[0] 
                    if len(ndx1) or len(ndx2) or len(ndx3) or len(ndx4): 
                        slavenode = True 
                    else: 
                        slavenode = False 
                else: 
                    if len(ndx1) or len(ndx2) or len(ndx3) : 
                        slavenode = True 
                    else: 
                        slavenode = False

                ix = np.where(npn[:,0]==e[1])[0][0]; n1=npn[ix]
                ix = np.where(npn[:,0]==e[2])[0][0]; n2=npn[ix]
                ix = np.where(npn[:,0]==e[3])[0][0]; n3=npn[ix]
                if e[4]>0: 
                    ix = np.where(npn[:,0]==e[4])[0][0]; n4=npn[ix]
                    xs =[n1[2], n2[2], n3[2], n4[2]]
                    ys =[n1[3], n2[3], n3[3], n4[3]]
                else: 
                    xs =[n1[2], n2[2], n3[2], n3[2]]
                    ys =[n1[3], n2[3], n3[3], n3[3]]
                    
                if slavenode:
                    nodalArea = nodalAreaByDistance(e, npn, elset) 
                    nodalArea = [0, 0, 0, 0] ## << nodal area by distance 
                else: 
                    nodalArea = ElementNodalArea(e, npn) 
                sumArea = 0 
                for na in nodalArea: 
                    sumArea += na 
                

                
                vs1 = average (sum(nodalArea_at_node * element_center_value ) / sumArea )  ## 식 정교화 필요 

                vs =[vs1 ]  ## how to?? 

                if cnt ==0: 
                    px, py, pv = meshgrid(xs, ys, vs=vs)
                    cnt += 1
                else: 
                    tx, ty, tv = meshgrid(xs, ys, vs=vs)

                    px = np.append(px, tx, axis=0)
                    py = np.append(py, ty, axis=0)
                    pv = np.append(pv, tv, axis=0)

            

def CenterValueToInnerValues(N, E, NR=0.98, G= 0.15E-3, **args):
    for key, value in args.items():
        if key == 'nr':
            NR = value
        if key == 'pointgap' or key == 'Pointgap' or key == 'pointdist' or key == 'pointdistance' or key == 'g':
            G = value
            
            
    Values=[]
    EN = len(E.Element)
    
    MasterEdge, SlaveEdge = E.MasterSlaveEdge(N)
    NodeSlave = ListNodeOnlySlave(MasterEdge, SlaveEdge)
    ## NodeSlave =[[SlaveNodeID, Element Name of the node, Slave EL No(ID), Slave Edge N1, N2, Master Element Name, Master EL No, Master Edge N1, N2], ..]
    for i in range(EN): 
        if E.Element[i][6] == 3 or E.Element[i][6] == 4:
            if len(E.Element[i]) == 11: 
                AlreadyDone = 0
                break
            elif len(E.Element[i]) == 19: 
                AlreadyDone = 1
                break
            else: 
                print ("CHECK ELEMENT Member.. (Element Value to Node Value)")
                return 0
        
    
    if AlreadyDone == 0 :
        for i in range(EN): 
            if E.Element[i][6] == 3 or E.Element[i][6] == 4:
                nN = E.Element[i][6]
                NodalArea = ElementNodalArea(E.Element[i], N)
                for j in range(3):
                    E.Element[i].append(NodalArea[j])
                if nN == 4:
                    E.Element[i].append(NodalArea[3])
                else:
                    E.Element[i].append('')
            else:
                E.Element[i].append('')
                E.Element[i].append('')
                E.Element[i].append('')
                E.Element[i].append('')
    else:
        # print 'Data is Filled already', E.Element[i]
        for i in range(EN): 
            if E.Element[i][6] == 3 or E.Element[i][6] == 4:
                nN = E.Element[i][6]
                NodalArea = ElementNodalArea(E.Element[i], N)
                for j in range(3):
                    E.Element[i][11+j]= NodalArea[j]
                if nN == 4:
                    E.Element[i][14]= NodalArea[3]
                else:
                    E.Element[i][14] = ''
                    
    # print 'Nodal Area, ', E.Element[i]  
    E.Save()
    EL4 = ELEMENT()
    NodesEL4=NODE()
    NoEL4 = 10000
    nc = 0
    for i in range(EN):
        # print 'i', i, E.Element[i]
        # print E.Element[i][6]
        if E.Element[i][6] == 3 or E.Element[i][6] == 4:
            nN = E.Element[i][6]
            
            for j in range(1, nN+1):
                nValue=[]
                SumArea = 0.0
                
                #######################################################
                ## Check if the node is a node on slave edge or not 
                ##  then calculate nodal value with length 
                #######################################################
                isSlave = 0
                sN = len(NodeSlave)
                TieShareEL=[]
                for k in range(sN): 
                    if E.Element[i][j] == NodeSlave[k][0]  and NodeSlave[k][1] == NodeSlave[k][5] and NodeSlave[k][1] == E.Element[i][5]:
                        isSlave = 1
                        MasterEL = NodeSlave[k][6]
                        TieShareEL.append(NodeSlave[k][2])
                        # print 'SlaveNODE', NodeSlave[k]
                #######################################################
                ## if the node is a slave node... 
                if isSlave ==1:
                    TieShareEL.append(MasterEL)
                    # print 'Calculation Slave', TieShareEL, 'at', E.Element[i][j], 'of', E.Element[i][0]
                    counting = 0
                    sN = len(TieShareEL)
                    SumLength = 0.0
                    for k in range(sN):
                        for m in range(EN):
                            if TieShareEL[k] == E.Element[m][0]: 
                                counting += 1
                                N1 = N.NodeByID(E.Element[i][j])
                                dist = math.sqrt((N1[2]-E.Element[m][8])*(N1[2]-E.Element[m][8]) + (N1[3]-E.Element[m][9])*(N1[3]-E.Element[m][9]))
                                nValue.append([dist, E.Element[m][10]])
                                SumLength += dist
                    if len(nValue) == 1:
                        if AlreadyDone == 0:
                            E.Element[i].append(nValue[0][1])
                        else: 
                            E.Element[i][14+j]= nValue[0][1]
                    else:
                        sN = len(nValue)
                        NodalValue = 0.0
                        for k in range(sN):
                            NodalValue += (1 - nValue[k][0] / SumLength) * nValue[k][1]
                        if AlreadyDone == 0:
                            E.Element[i].append(NodalValue)
                        else: 
                            E.Element[i][14+j]= NodalValue
                ## if the node is NOT a slave node... 
                else:
                    nValue.append([E.Element[i][10+j], E.Element[i][10]])
                    SumArea += E.Element[i][10+j]
                    for k in range(EN):
                        if i != k and  (E.Element[k][6] == 3 or E.Element[k][6] == 4) and (E.Element[i][5] == E.Element[k][5] ):
                            nNk = E.Element[k][6]
                            for m in range(1, nNk+1):
                                if E.Element[i][j] == E.Element[k][m]:
                                    nValue.append([E.Element[k][10+m], E.Element[k][10]])
                                    SumArea += E.Element[k][10+m]
                    
                    vN = len(nValue)
                    NodalValue = 0.0
                    for k in range(vN):
                        NodalValue += nValue[k][0] / SumArea * nValue[k][1]
                    if AlreadyDone == 0:
                        E.Element[i].append(NodalValue)
                    else:
                        E.Element[i][14+j]= NodalValue
            if nN == 3:
                if AlreadyDone == 0:
                    E.Element[i].append('')
                else:
                    E.Element[i][18]= ''
                
            ###### DONE : The End of Nodal Value Calculation at all nodes of a Element
             
            if NR == 0.3:
                for m in range(len(V)):
                    if math.isnan(V[m]) == False:
                        Values.append([X[m], Y[m], V[m]])

            else:
                # print '* ', E.Element[i]
                N1 = N.NodeByID(E.Element[i][1])
                N2 = N.NodeByID(E.Element[i][2])
                N3 = N.NodeByID(E.Element[i][3])
                if E.Element[i][6] == 3:
                    # print '*** [3]', E.Element[i]
                    X = [N1[2], N2[2], N3[2], N3[2]]
                    Y = [N1[3], N2[3], N3[3], N3[3]]
                    V = [E.Element[i][15], E.Element[i][16], E.Element[i][17]]
                    # Stx = (X[0], X[1], X[2])
                    # Sty = (Y[0], Y[1], Y[2])
                    # Stv = (V[0], V[1], V[2])
                    # print "6", Stx, Sty, Stv
                    # results = _islm.ValuesOfPointsInElement(Stx, Sty, Stv, NR, G)

                    px, py, pv = meshgrid_in_Quadrilateral(X, Y, vs=V)

                    del(Stx)
                    del(Sty)
                    del(Stv)
                    Nn = len(results)
                    for m in range(int(Nn/3)):
                        if math.isnan(results[m * 3 + 2]) == False:
                            Values.append([results[m*3], results[m*3+1], results[m*3+2]])

                elif E.Element[i][6] == 4:
                    N4 = N.NodeByID(E.Element[i][4]) 
                    
                    NodesEL4.Add( [NoEL4+nc*4 + 1, N1[1], N1[2], N1[3], E.Element[i][15] ])
                    NodesEL4.Add( [NoEL4+nc*4 + 2, N2[1], N2[2], N2[3], E.Element[i][16] ])
                    NodesEL4.Add( [NoEL4+nc*4 + 3, N3[1], N3[2], N3[3], E.Element[i][17] ])
                                                                                                   
                    NodesEL4.Add( [NoEL4+nc*4 + 4, N4[1], N4[2], N4[3], E.Element[i][18] ])
                    
                    EL4.Add( [ NoEL4+nc+1, NoEL4+nc*4 + 1, NoEL4+nc*4 + 2, NoEL4+nc*4 + 3, NoEL4+nc*4 + 4, E.Element[i][5], E.Element[i][6], E.Element[i][7]    ] )
                    
                    nc += 1
                    
        else:
            if AlreadyDone == 0:
                E.Element[i].append('')
                E.Element[i].append('')
                E.Element[i].append('')
                E.Element[i].append('') 
                    
    
    xyv = Innerpoints(G, EL4, NodesEL4, XY=23)
    I = len(xyv.Node)
    for i in range(I):
            Values.append( [ xyv.Node[i][2], xyv.Node[i][3], xyv.Node[i][4] ]) 
    
    return Values


def Printlist(plist, head=5, tail=5, all=0): 
    if not isinstance(plist, list): 
        print ("## the instance is not a list")
        return 
    N = len(plist)
    cnt = 0 
    if all >= 0: 
        for i, pl in enumerate(plist): 
            
            if i < head: print ("%d"%(i), pl)
            elif i > N-tail-1: print ("%d"%(i), pl)
            else: 
                if cnt ==0: print (" ''''' ")
                cnt += 1
    else: 
        for pl in plist: 
            print ("%d"%(i), pl)

def lstArea(NodeList, Node, XY=23, **args):   ## Calculate Area of a polygon
    errorimage = 1
    for key, value in args.items():
        if key == 'xy' :   XY = int(value)
        if key == 'error': errorimage= int(value)
    
    if len(Node) > 0:
        node = np.array(Node) 
    
        ii = int(XY/10)
        jj = int(XY)%10
        
        n = len(NodeList)
        x = [];       y = []
        
        for nd in NodeList: 
            ix = np.where(node[:,0]==nd)[0][0]
            Ni = node[ix]
            x.append(Ni[ii])
            y.append(Ni[jj])
    
        x.append(x[0])
        y.append(y[0])
        
        A = [0.0, 0.0, 0.0]

        for i in range(n):
            s = x[i] * y[i + 1] - x[i + 1] * y[i]
            A[0] += s
            A[1] += (x[i] + x[i + 1]) * s
            A[2] += (y[i] + y[i + 1]) * s

        A[0] = A[0] / 2.0
        try:
            A[1] = A[1] / A[0] / 6
            A[2] = A[2] / A[0] / 6
        except:
            if errorimage > 0: 
                print ("!! Error to calculate Area ", A)
                pNode = NODE()
                for i in range(n):
                    ix = np.where(node[:,0]==nd)[0][0]
                    Ni = node[ix]
                    pNode.Add(Ni)
                pNode.Image("Error_Area_"+str(pNode.Node[0][0])+".png", size=1.5)
            return [0.0, 0.0, 0.0]

        if A[0] < 0:
            A[0] = -A[0]
            # print 'Negative Area Calculation! '
        return A
    else:
        print ("Length of Node is 0")
        return [0.0, 0.0, 0.0]
def FindEdge(elements):
    # elements=[Element No, 1st Node, 2nd Node, 3rd Node, 4th Node, Mat_name, '']
    # edge element [node1, node2, ElsetName, FaceID, elementNo, Tie Definition No]
    i = 0
    edges = []
    while i < len(elements):
        if elements[i][6] == 4:
            edges.append([elements[i][1], elements[i][2], elements[i][5], 'S1', elements[i][0], -1])
            edges.append([elements[i][2], elements[i][3], elements[i][5], 'S2', elements[i][0], -1])
            edges.append([elements[i][3], elements[i][4], elements[i][5], 'S3', elements[i][0], -1])
            edges.append([elements[i][4], elements[i][1], elements[i][5], 'S4', elements[i][0], -1])
        elif elements[i][6] == 3:
            edges.append([elements[i][1], elements[i][2], elements[i][5], 'S1', elements[i][0], -1])
            edges.append([elements[i][2], elements[i][3], elements[i][5], 'S2', elements[i][0], -1])
            edges.append([elements[i][3], elements[i][1], elements[i][5], 'S3', elements[i][0], -1])
        i += 1
    return edges
def FindFreeEdge(edges):
    freeEdge = []
    npedge = []
    for e in edges: 
        npedge.append([e[0], e[1]])
    npedge = np.array(npedge)

    free = []
    cnt = 0 
    for ed in npedge: 
        idx1 = np.where(npedge[:,:]==int(ed[0]))[0]
        idx2 = np.where(npedge[:,:]==int(ed[1]))[0]
        idx = np.intersect1d(idx1, idx2) 
            
        if len(idx) == 1: 
            free.append(edges[cnt])
        cnt += 1 
    # print ("NO of Free Edge=", len(free))
    return free
def FreeEdge(edge):
    FEdge = EDGE()
    
    for i in range(len(edge.Edge)):
        if edge.Edge[i][5] == -1:
            j = i + 1
            count = 0
            while j < len(edge.Edge):
                if edge.Edge[j][5] == -1:
                    if edge.Edge[i][0] == edge.Edge[j][0] and edge.Edge[i][1] == edge.Edge[j][1]:
                        count += 1
                        edge.Edge[i][5] = -2
                        edge.Edge[j][5] = -2
                        break
                    elif edge.Edge[i][0] == edge.Edge[j][1] and edge.Edge[i][1] == edge.Edge[j][0]:
                        count += 1
                        edge.Edge[i][5] = -2
                        edge.Edge[j][5] = -2
                        break
                j += 1
            if count == 0:
                edge.Edge[i][5] = 0
                FEdge.Add(edge.Edge[i])
    return FEdge
def OuterEdge(FreeEdge, Node, Element):

    # FreeEdge.Image(Node, "EDGEIMAGE")
    
    N = len(FreeEdge.Edge)
    
    MinY = 9.9E20
    
    cNodes = [0]
    for i in range(N):
        N1 = Node.NodeByID(FreeEdge.Edge[i][0])
        N2 = Node.NodeByID(FreeEdge.Edge[i][1])
        if N1[3] < MinY:
            MinY = N1[3]
            cNodes[0] = N1[0]
        if N2[3] < MinY:
            MinY = N2[3]
            cNodes[0] = N2[0]
    if cNodes[0] == 0:
        cNodes[0] = Node.NodeIDByCoordinate('z', 0.0, closest=1)

    MAX = 10000   ## max iteration for searching  error
    ShareNodePos = []
    #    connectedEdge = []
    outEdge = EDGE()

    ## Find a 1st surround edge (IL at the center)
    low = 9.9E20
    i = 0
    savei = 0
    while i < len(cNodes):
        j = 0
        while j < len(Node.Node):
            if cNodes[i] == Node.Node[j][0]:
                if Node.Node[j][3] < low:
                    low = Node.Node[j][3]
                    savei = j
            j += 1
        i += 1

    i = 0
    while i < len(FreeEdge.Edge):
        if Node.Node[savei][0] == FreeEdge.Edge[i][0]:
            break
        i += 1

    ## End of 1st Outer Edge finding (IL1)
    # print (i, len(FreeEdge.Edge))
    # print (FreeEdge.Edge[i])
    FreeEdge.Edge[i][5] = 1
    outEdge.Add(FreeEdge.Edge[i])
    iFirstNode = FreeEdge.Edge[i][0]

    count = 0

    #    i=  # i is no matter how big, because i is redefined when next edge is found

    while i < len(FreeEdge.Edge):
        count += 1
        if count > MAX:
            print ('[INPUT] CANNOT FIND OUTER EDGES IN THE MODEL')
            del (outEdge)
            outEdge = EDGE()
            return outEdge
        j = 0
        while j < len(FreeEdge.Edge):
            if i != j:
                if FreeEdge.Edge[i][1] == FreeEdge.Edge[j][0]:
                    # print ('edge[i][1], [j][0] ', FreeEdge.Edge[i], FreeEdge.Edge[j], 'i=', i)
                    ShareNodePos.append(j)
                    # print (ShareNodePos, FreeEdge.Edge[ShareNodePos[0]][0])
            j = j + 1
        if len(ShareNodePos) != 0:
            if FreeEdge.Edge[ShareNodePos[0]][0] == iFirstNode:
                break
        else:
            print ('[INPUT] CANNOT FIND CONNECTED FREE EDGE. CHECK TIE CONDITION')
            del (outEdge)
            outEdge = EDGE()
            return outEdge
        # print ('sharenodePos count = ', len(ShareNodePos))

        if len(ShareNodePos) == 1:
            FreeEdge.Edge[ShareNodePos[0]][5] = 1
            outEdge.Add(FreeEdge.Edge[ShareNodePos[0]])
            i = ShareNodePos[0]
            del ShareNodePos
            ShareNodePos = []
        else:
            if FreeEdge.Edge[i][4] == FreeEdge.Edge[ShareNodePos[0]][4]:
                tmpPos = ShareNodePos[1]
            else:
                SHARE = ShareEdge(FreeEdge.Edge[i][4], FreeEdge.Edge[ShareNodePos[1]][4], Element)
                if SHARE == 1:
                    tmpPos = ShareNodePos[0]
                else:
                    tmpPos = ShareNodePos[1]

                    #######################################################
                    nfe1 = 0; nfe2 = 0
                    for fe in FreeEdge.Edge:
                        if fe[4] == FreeEdge.Edge[tmpPos][4]:
                            # print (fe)
                            nfe1 += 1
                        if fe[4] == FreeEdge.Edge[ShareNodePos[0]][4]:
                            # print (fe)
                            nfe2 += 1
                    # print ("nfe=", nfe, FreeEdge[tmpPos])
                    if nfe1 < nfe2:
                        tmpPos = ShareNodePos[0]
                    elif nfe1 == nfe2:
                        tienode = FreeEdge.Edge[tmpPos][0]
                        nc = 0
                        for fe in FreeEdge.Edge:
                            if fe[4] == FreeEdge.Edge[tmpPos][4] and fe[1] == tienode: 
                                nc += 1
                                break
                        if nc == 0:   tmpPos = ShareNodePos[0]
                    ########################################################

            FreeEdge.Edge[tmpPos][5] = 1
            outEdge.Add(FreeEdge.Edge[tmpPos])
            i = tmpPos
            del ShareNodePos
            ShareNodePos = []
            
    return outEdge
def ShareEdge(m, n, Elements):
    p = ElementShape(m, Elements)
    q = ElementShape(n, Elements)
    lst = []
    if type(lst) != type(Elements): 
        N = len(Elements.Element)
        for i in range(N):
            if m == Elements.Element[i][0]:
                k = i
            if n == Elements.Element[i][0]:
                l = i

        count = 0
        for i in range(1, p+1):
            for j in range(1, q+1):
                if Elements.Element[k][i] == Elements.Element[l][j]:
                    count += 1
    else: 
        for i, el in enumerate(Elements): 
            if m == el[0]: k = i
            if n == el[0]: l = i
        count = 0 
        for i in range(1, p+1):
            for j in range(1, q+1):
                if Elements[k][i] == Elements[l][j]:
                    count += 1

    if count >= 2:
        return 1  # Edge shared
    else:
        return 0
def ElementShape(k, Elements):
    # k = element No.
    lst = []
    if type(lst) != type(Elements): 
        for el in Elements.Element: 
            if k == el[0]: return el[6]
    else: 
        for el in Elements: 
            if k == el[0]: 
                return el[6]

    print (k, 'Element was not found')
    return 0
def NextEdge(edge, alledges): 
    next = -1
    for i, ed in enumerate(alledges): 
        if ed[0] == edge[1]: 
            next = i
            break
    # if next == -1: print (f"no matching to Next edge       {edge}")
    return next 
def PreviousEdge(edge, alledges): 
    next = -1
    for i, ed in enumerate(alledges): 
        if ed[1] == edge[0]: 
            next = i
            break
    # if next == -1: print (f"no matching to Previous edge   {edge}")
    return next 
def Contact_relation_2Elements(e1, e2): 
    m1=[]; m2=[]
    for i in range(1, 5):
        if e1[i] ==0: continue 
        for j in range(1, 5): 
            if e2[j] == 0: continue 
            if e1[i] == e2[j] : 
                if i == 1: m1.append([e1[i], 0])
                else: m1.append([e1[i], i])
                m2.append([e2[j], j])

    if len(m1) == 0: 
        return None, 0 
    elif len(m1) ==1: 
        if m1[0][1] == 0: 
            pos = 1
        else: pos = m1[0][1] 
        return 'Point', pos 
    else:
        if int(e1[4]) != 0: 
            if m1[0][1] + m1[1][1] == 2 : face = 1
            if m1[0][1] + m1[1][1] == 5 : face = 2
            if m1[0][1] + m1[1][1] == 7 : face = 3
            if m1[0][1] + m1[1][1] == 4 : face = 4
        else: 
            if m1[0][1] + m1[1][1] == 2 : face = 1
            if m1[0][1] + m1[1][1] == 5 : face = 2
            if m1[0][1] + m1[1][1] == 3 : face = 3
        return 'Edge', face 
def Element2D_NextElement(el, solids, nodes, start=1, next=2): 
    ## solids = [No, N1, N2, N3, N4, No_of_nodes(3 or 4)]
    x =2; y=3 
    if el[5] == 4: 
        nf = start + next 
        if nf > 4 : nf -= 4  
        nfnode = [el[nf], el[nf+1]]
        if nf + 1 == 5 : nfnode = [ el[nf], el[1] ] 

        ix1 = np.where(solids[:, 1:5] == nfnode[0])[0] 
        ix2 = np.where(solids[:, 1:5] == nfnode[1])[0] 
        ix  = np.intersect1d(ix1, ix2)

        if len(ix) == 1: 
            n_solid = []
        elif len(ix) == 2: 
            if solids[ix[0]][0] == el[0] : n_solid = solids[ix[1]]
            else:                          n_solid = solids[ix[0]]
        else: 
            print ("!!!! ERROR. No FOUND NEXT 2D ELEMENT. ")
            sys.exit()
    elif el[5] ==3: 

        nf = start + next 
        if nf > 3 : nf -= 3  
        nfnode = [el[nf], el[nf+1]]
        if nf + 1 == 4  : nfnode = [ el[nf], el[1] ] 

        ix1 = np.where(solids[:, 1:5] == nfnode[0])[0] 
        ix2 = np.where(solids[:, 1:5] == nfnode[1])[0] 
        ix  = np.intersect1d(ix1, ix2)

        if len(ix) == 1: 
            n_solid = []
        elif len(ix) == 2: 
            if solids[ix[0]][0] == el[0] : n_solid = solids[ix[1]]
            else:                          n_solid = solids[ix[0]]
        else: 
            print ("!!!! ERROR. No FOUND NEXT 2D ELEMENT. ")
            sys.exit()
        
    return n_solid 
def ElementDuplicationCheck(AllElements):

    dup = 0 
    els = []
    for el in AllElements: 
        els.append([el[0], el[1], el[2], el[3], el[4], el[6]]) 
    els = np.array(els)    

    N0 = len(els)
    eln = els[:,0] 
    eln = np.unique(eln)

    if N0 > len(eln): 
        print ("##  %d elements are duplicates"%(N0-len(eln)))
        dup =  1 
        return dup 

    idx = np.where(els[:,5] == 2)[0]
    el2 = els[idx]
    if len(idx) > 0: 
        el2 = els[idx]
        for el in el2: 
            ix1 = np.where(el2[:, 1:3] == el[1])[0]
            ix2 = np.where(el2[:, 1:3] == el[2])[0]
            ix = np.intersect1d(ix1, ix2)
            if len(ix) > 1: 
                print ('## Error! Rebar ' + str(el[0]) + ' is overlapped.')
                dup =  1
    
    idx = np.where(els[:,5] == 3)[0]
    if len(idx) > 0: 
        el3 = els[idx]
        for el in el3: 
            ix1 = np.where(el3[:, 1:4] == el[1])[0]
            ix2 = np.where(el3[:, 1:4] == el[2])[0]
            ix3 = np.where(el3[:, 1:4] == el[3])[0]

            ix = np.intersect1d(ix1, ix2)
            ix = np.intersect1d(ix, ix3)
            if len(ix) > 1: 
                print ('## Error! CGAX3H ' + str(el[0]) + ' is overlapped.')
                dup =  1
    idx = np.where(els[:,5] == 4)[0]
    if len(idx) > 0: 
        el4 = els[idx]
        for el in el4: 
            ix1 = np.where(el4[:, 1:5] == el[1])[0]
            ix2 = np.where(el4[:, 1:5] == el[2])[0]
            ix3 = np.where(el4[:, 1:5] == el[3])[0]
            ix4 = np.where(el4[:, 1:5] == el[4])[0]

            ix = np.intersect1d(ix1, ix2)
            ix = np.intersect1d(ix, ix3)
            ix = np.intersect1d(ix, ix4)
            if len(ix) > 1: 
                print ('## Error! CGAX4H ' + str(el[0]) + ' is overlapped.')
                dup =  1
    return dup 

def ChaferDivide(Elements, ChaferName, Elset, Node):

    asymmetric =0
    I = len(Elements)
    J = len(Node)
    for i in range(I):
        if Elements[i][5]=="BT2":
            for j in range(J):
                if Node[j][0] == Elements[i][1]:
                    if Node[j][2] == 0:
                        relement=Elements[i][0]
                if Node[j][0] == Elements[i][2]:
                    if Node[j][2] == 0:
                        lelement=Elements[i][0]
    
    Offset = abs(relement - lelement )
    if Offset != 2000:
        for j in range(J):
            if Node[j][2] > 0.0 :
                Offset = Node[j][0]
                break
        # print ("Asymmetric Model %d"%(Offset))
        asymmetric =1

                    
    for i in range(len(Elements)):
        if Elements[i][5] == "BEAD_R" and Elements[i][8] < 0.0 :
            Elements[i][5] = 'BEAD_L'

    if asymmetric ==1:
        Offset = 0
    return Elements, Elset, Offset
def FindSolidElementBetweenMembrane(m1, m2, Elements):
    # Data types of m1, m2 are string
    between = []
    Elm1 = []
    Elm2 = []

    for i in range(len(Elements)):
        if Elements[i][5] == m1 or Elements[i][5] == m2:
            for j in range(len(Elements)):
                if j != i and Elements[j][6] == 3:
                    for k in range(1, 3):
                        if k == 1:
                            m = 2
                        else:
                            m = 1
                        for l in range(1, 4):
                            n = l + 1
                            if n > 3:
                                n = l - 3

                            if Elements[i][k] == Elements[j][l] and Elements[i][m] == Elements[j][n]:
                                if Elements[i][5] == m1:
                                    Elm1.append(Elements[j])
                                else:
                                    Elm2.append(Elements[j])
                                break

                elif j != i and Elements[j][6] == 4:
                    for k in range(1, 3):
                        if k == 1:
                            m = 2
                        else:
                            m = 1
                        for l in range(1, 5):
                            n = l + 1
                            if n > 4:
                                n = l - 4

                            if Elements[i][k] == Elements[j][l] and Elements[i][m] == Elements[j][n]:
                                if Elements[i][5] == m1:
                                    Elm1.append(Elements[j])
                                else:
                                    Elm2.append(Elements[j])
                                break

    for i in range(len(Elm1)):
        for j in range(len(Elm2)):
            if Elm1[i][0] == Elm2[j][0]:
                between.append(Elm2[j][0])
                break

    m1f=[]
    m2f=[]
    for i in range(len(Elm1)):
        match=0
        for j in range(len(between)):
            if Elm1[i][0] == between[j]:
                match = 1
                break
        if match == 0:
            m1f.append(Elm1[i])

    for i in range(len(Elm2)):
        match = 0
        for j in range(len(between)):
            if Elm2[i][0] == between[j]:
                match = 1
                break
        if match == 0:
            m2f.append(Elm2[i])
    ########################################################
    for i in range(len(m2f)):
        for j in range(len(m1f)):
            match = 0
            for m in range(1, 5):
                for n in range(1, 5):
                    if m2f[i][m] == m1f[j][n] and m1f[j][n] != '':
                        match += 1
            if match == 2:
                between.append(m2f[i][0])
                between.append(m1f[j][0])
                # print "Between Appended"
                break


    return between
def existEL(SetName, Elsets):
    Exist =0
    for i in range(len(Elsets)):
        if SetName in Elsets[i][0]:
            Exist = 1
            break
    return Exist
def BeadWidth(Element, Node):
    BDCore = []

    for i in range(len(Element)):
        if Element[i][5] == 'BEAD_L' or Element[i][5] == 'BD1':
            BDCore.append(Element[i])
    if len(BDCore) == 0:
        for i in range(len(Element)):
            if Element[i][5] == 'BEAD_R' or Element[i][5] == 'BD1':
                BDCore.append(Element[i])
    BDEdge = FindEdge(BDCore)
    BDFree = FindFreeEdge(BDEdge)

    nodes = []
    for i in range(len(BDFree)):
        if i == 0:
            nodes.append(BDFree[i][0])
            nodes.append(BDFree[i][1])
        else:
            N = len(nodes)
            match1 = 0
            match2 = 0
            for j in range(N):
                if BDFree[i][0] == nodes[j]:
                    match1 = 1
                if BDFree[i][1] == nodes[j]:
                    match2 = 1
            if match1 == 0:
                nodes.append(BDFree[i][0])
            if match2 == 0:
                nodes.append(BDFree[i][1])
    center = lstArea(nodes, Node)
    min = 100000.0
    max = -100000.0
    npn = np.array(Node)
    for nd in nodes: 
        ix = np.where(npn[:,0]==nd)[0][0]
        C = npn[ix]
        if abs(C[2]) > max:
            max = abs(C[2])
        if abs(C[2]) < min:
            min = abs(C[2])

    if abs(center[1]) > abs(max) or abs(center[1]) < abs(min):
        center[1] = (abs(max)+abs(min))/2.0

    return abs(min), abs(max), abs(center[1])
def ElementCheck (Elements, Nodes):
    # rebar Connectivity
    # Element Shape 
    
    tmpEL = []
    for i in range(len(Elements)):
        tmpEL.append(Elements[i][0])
        tmpEL.append(Elements[i][1])
        tmpEL.append(Elements[i][2])
        if Elements[i][3] == '':
            tmpEL.append(0)
        else:
            tmpEL.append(Elements[i][3])
        if Elements[i][4] == '':
            tmpEL.append(0)
        else:
            tmpEL.append(Elements[i][4])
        tmpEL.append(Elements[i][6])

    tmpND = []
    tmpCoord = []
    for i in range(len(Nodes)):
        tmpND.append(Nodes[i][0])
        tmpCoord.append(Nodes[i][1])
        tmpCoord.append(Nodes[i][2])
        tmpCoord.append(Nodes[i][3])
    
    tupleEL =tuple(tmpEL)
    tupleND =tuple(tmpND)
    tupleCo = tuple(tmpCoord)

    try: 
        Results  = _islm.SolidRebarCheck(tupleEL, 6, tupleND, tupleCo, 3, 5)

        Message=[]
        if Results > 100000:
            if Results % 10 == 1:
                Message.append(["ERROR::PRE::[INPUT] More than 1 Sold Element is distorted."])
                Results = 0
            if (Results/10) % 10 == 1:
                Message.append(["ERROR::PRE::[INPUT] More than 1 Rebar is disconnected."])
                Results = 0
            if (Results/100) % 10 == 1:
                Message.append(["ERROR::PRE::[INPUT] Some Rebar Elements are Defined Twice or More."])
                Results = 0
            if (Results/1000) % 10 == 1:
                Message.append(["ERROR::PRE::[INPUT] Some CGAX3H Elements are Defined Twice or More."])
                Results = 0
            if (Results/10000) % 10 == 1:
                Message.append(["ERROR::PRE::[INPUT] Some CGAX4H Elements are Defined Twice or More."])
                Results = 0
        else:
            Results = 1
    except: 
        Message=[]
        Results = 1 

    return Results, Message
def TieSurface(Elements, Nodes):
    AllEdges=FindEdge(Elements)
    FreeEdges=FindFreeEdge(AllEdges)

    Nodes = np.array(Nodes)
    idx = np.where(Nodes[:,3] > 0.005)[0]
    Nodes = Nodes[idx]
    CenterNodes = FindCenterNodes(Nodes)

    OutEdges=FindOutEdges(FreeEdges, CenterNodes, Nodes, Elements)

    InnerFreeEdges = FindInnerFreeEdges(FreeEdges, OutEdges)
    ND =NODE()
    for n in Nodes:
        ND.Add(n)
    MasterEdges, SlaveEdges, Tie_error = DivideInnerFreeEdgesToMasterSlave(InnerFreeEdges, ND)
    return MasterEdges, SlaveEdges, OutEdges,CenterNodes, FreeEdges, AllEdges, Tie_error
def FindCenterNodes(Nodes):
    centers=[]
    for i in range(len(Nodes)):
        if Nodes[i][2] == 0:
            centers.append(Nodes[i][0])
    return centers
def FindOutEdges(FreeEdge, CenterNodes, Nodes, Elements):
    # node of free edge is shared with another tie node
    # ShareNodePos -> 1 or 2
    # for fe in FreeEdge:
        # if fe[2] == "SUT":        print fe
    # for i in range(I):
    MAX = 10000
    ShareNodePos = []
    outEdge = []

    npn = np.array(Nodes)
    idxs = np.where(npn[:,2] > 0)[0]
    posNodes = npn[idxs]
    zMin = np.min(posNodes[:,3])
    idx = np.where(posNodes[:,3]==zMin)[0]
    if len(idx) > 0: 
        rightToe = posNodes[idx[0]]
        if len(idx)>1: 
            for ix in idx: 
                print("### Outer surface right bottom : ", posNodes[ix])
        
    for k, ix in enumerate(FreeEdge):
        if ix[0] == rightToe[0]: 
            i = k 
            break 
    FreeEdge[i][5] = 1
    outEdge.append(FreeEdge[i])
    iFirstNode = FreeEdge[i][0]

    count = 0

    #    i=  # i is no matter how big, because i is redefined when next edge is found

    while i < len(FreeEdge):
        count += 1
        if count > MAX:
            print ('[INPUT] CANNOT FIND OUTER EDGES IN THE MODEL (too much iteration)')
            outEdge = []
            return outEdge
        j = 0
        while j < len(FreeEdge):
            if i != j:
                if FreeEdge[i][1] == FreeEdge[j][0]:
                    # print 'edge[i][1], [j][0] ', FreeEdge.edge[i], FreeEdge.edge[j], 'i=', i
                    ShareNodePos.append(j)
            j = j + 1
        #        print ('**', ShareNodePos)

        #        ShareNodePos=FindNextEdge()
        #        print (ShareNodePos, FreeEdge[ShareNodePos[0]][0])
        if len(ShareNodePos) != 0:
            if FreeEdge[ShareNodePos[0]][0] == iFirstNode:
                break
        else:
            print ('[INPUT] CANNOT FIND CONNECTED FREE EDGE. CHECK TIE CONDITION')
            outEdge=[]
            return outEdge
        # print 'sharenodePos count = ', len(ShareNodePos)
        if len(ShareNodePos) == 1:
            FreeEdge[ShareNodePos[0]][5] = 1
            outEdge.append(FreeEdge[ShareNodePos[0]])
            # print ("1,", FreeEdge[ShareNodePos[0]])
            i = ShareNodePos[0]

            del ShareNodePos
            ShareNodePos = []
        else:
            if FreeEdge[i][4] == FreeEdge[ShareNodePos[0]][4]:
                tmpPos = ShareNodePos[1]
                # print ("passed here")
            else:
                SHARE = ShareEdge(FreeEdge[i][4], FreeEdge[ShareNodePos[1]][4], Elements)
                if SHARE ==1:
                    tmpPos = ShareNodePos[0]
                else:
                    tmpPos = ShareNodePos[1]
                    nfe1 = 0; nfe2 = 0
                    for fe in FreeEdge:
                        if fe[4] == FreeEdge[tmpPos][4]:
                            # print (fe)
                            nfe1 += 1
                        if fe[4] == FreeEdge[ShareNodePos[0]][4]:
                            # print (fe)
                            nfe2 += 1
                    # print ("nfe=", nfe, FreeEdge[tmpPos])
                    if nfe1 < nfe2:
                        tmpPos = ShareNodePos[0]
                    elif nfe1 == nfe2:
                        tienode = FreeEdge[tmpPos][0]
                        nc = 0
                        for fe in FreeEdge:
                            if fe[4] == FreeEdge[tmpPos][4] and fe[1] == tienode: 
                                nc += 1
                                break
                        if nc == 0:   tmpPos = ShareNodePos[0]
                    

            FreeEdge[tmpPos][5] = 1
            outEdge.append(FreeEdge[tmpPos])
            # print ("2, ", FreeEdge[ShareNodePos[0]], FreeEdge[ShareNodePos[1]], "-", FreeEdge[tmpPos])
            i = tmpPos
            del ShareNodePos
            ShareNodePos = []
            
    # print ("#### len", len(outEdge))
    return outEdge
def FindInnerFreeEdges(Oedge, medge): 
    residuals=[]
    for ed1 in Oedge: 
        mch = 0 
        for ed2 in medge: 
            if ed1[0] == ed2[0] and ed1[1] == ed2[1] : 
                mch = 1 
                break 
        if mch == 0: 
            residuals.append(ed1)
    return residuals 
def DivideInnerFreeEdgesToMasterSlave(edges, node_class): 
    nodes = node_class 

    masters=[]
    i = 0 
    while i < len(edges): 
        con1 = 0;         con2 = 0 
        c1e = [];         c2e = []
        
        N01 = nodes.NodeByID(edges[i][0])
        N02 = nodes.NodeByID(edges[i][1])
        ML = NodeDistance(N01, N02) 
        Ly=[N01[2], N02[2]]; Lz=[N01[3], N02[3]]
        MinY = min(Ly); MaxY = max(Ly)
        MinZ = min(Lz); MaxZ = max(Lz)
        tslave = []
        for j in range(len(edges)): 
            
            if i == j : continue 
            
            if edges[i][0] == edges[j][0] or edges[i][0] == edges[j][1] :
                N1 = nodes.NodeByID(edges[j][0])
                N2 = nodes.NodeByID(edges[j][1])
                SL = NodeDistance(N1, N2)  
                if ML < SL: 
                    continue 
                if edges[i][0] == edges[j][0]: 
                    if N2[2] >= MinY and N2[2] <= MaxY and N2[3] >= MinZ and N2[3] <= MaxZ: 
                        dist = DistanceFromLineToNode2D(N2, [N01, N02], onlydist=1)
                        if dist < .10E-03 : 
                            con1 = 1 
                            c1e = edges[j]
                            tslave.append(edges[j])
                else: 
                    if N1[2] >= MinY and N1[2] <= MaxY and N1[3] >= MinZ and N1[3] <= MaxZ: 
                        dist = DistanceFromLineToNode2D(N1, [N01, N02], onlydist=1)
                        if dist < .10E-03 : 
                            con1 = 1 
                            c1e = edges[j]
                            tslave.append(edges[j])

            if edges[i][1] == edges[j][0] or edges[i][1] == edges[j][1] : 
                N1 = nodes.NodeByID(edges[j][0])
                N2 = nodes.NodeByID(edges[j][1])
                SL = NodeDistance(N1, N2)  
                if ML < SL: 
                    continue 

                if edges[i][1] == edges[j][0]: 
                    if N2[2] >= MinY and N2[2] <= MaxY and N2[3] >= MinZ and N2[3] <= MaxZ: 
                        dist = DistanceFromLineToNode2D(N2, [N01, N02], onlydist=1)
                        if dist < .10E-03 : 
                            con2 = 1 
                            c2e = edges[j]
                            tslave.append(edges[j])
                else: 
                    if N1[2] >= MinY and N1[2] <= MaxY and N1[3] >= MinZ and N1[3] <= MaxZ: 
                        dist = DistanceFromLineToNode2D(N1, [N01, N02], onlydist=1)
                        if dist < .10E-03 : 
                            con2 = 1 
                            c2e = edges[j]
                            tslave.append(edges[j])

            if con1 == 1 and con2 ==1: 
                break 
        
        if con1 ==1 or con2 == 1: 
            masters.append([edges[i], tslave])
        i+=1 

    excluding = []
    isError = 0 
    TIE_ERROR = []
    for e in masters: 
        
        excluding.append(e[0])
        excluding.append(e[1][0])
        if len(e[1]) <2: 
            print ("## Error to find Tie Master surface (%d)"%(e[0][4]))
            TIE_ERROR.append(e[0][4])
            isError = 1
            continue 
            
        excluding.append(e[1][1])
    if isError == 1: 
        master_edge=[];  slave_edge=[]
        return master_edge, slave_edge, TIE_ERROR

    print ("* Tie No. in layout mesh =%d"%(len(masters)))
    master_edge = []
    slave_edge = []
    for i, ed in enumerate(masters): 
        master_edge.append(ed[0])
        s_temp=[]
        s_temp.append(ed[1][0])
        
        s_temp.append(ed[1][1])
        if ed[1][0][0] == ed[1][1][1] or ed[1][0][1] == ed[1][1][0] : 
            slave_edge.append(s_temp)
            # print ("* Slave Edges %2d, No=%d"%(i, len(s_temp)))
            continue 

        nexts = ConnectedEdge(ed[1][0], edges, exclude=excluding)
        if len(nexts) ==0: 
            # print (len(master_edge))
            # print ("--")
            # print (len(slave_edge))
            # print ("**")
            # print (TIE_ERROR)
            return master_edge, slave_edge, TIE_ERROR
        
        s_temp.append(nexts[0])
        # print (ed[0], ":", nexts)
        if nexts[0][0] != ed[1][1][1] and nexts[0][1] != ed[1][1][0] : 
            excluding.append(nexts[0]) 
            nexts = ConnectedEdge(nexts[0], edges, exclude=excluding)
            while len(nexts): 
                s_temp.append(nexts[0])
                excluding.append(nexts[0]) 
                nexts = ConnectedEdge(nexts[0], edges, exclude=excluding)
        
            # print (ed[0], ":::", nexts)
        slave_edge.append(s_temp)
        # print ("**Slave Edges %2d, No=%d"%(i, len(s_temp)))

    return master_edge, slave_edge, TIE_ERROR
def NodeDistance(N1, N2, xy=0): 
    if xy > 10: 
        x = int(xy/10); y=int(xy%10)
        return math.sqrt((N2[x] - N1[x])**2 + (N2[y] - N1[y])**2)
    else: 
        return math.sqrt((N2[1]-N1[1])*(N2[1]-N1[1]) + (N2[2]-N1[2])*(N2[2]-N1[2]) + (N2[3]-N1[3])*(N2[3]-N1[3]))
def DistanceFromLineToNode2D(N0, nodes=[], xy=12, onlydist=0):
    x = int(xy/10)
    y = int(xy%10)

    N1=nodes[0]
    N2=nodes[1]
    if len(nodes) ==2: 
        if round(N2[x]-N1[x], 6) !=0: 
            a = (N2[y]-N1[y])/(N2[x]-N1[x])
            A = -a
            C = a * N1[x] - N1[y]

            ## intersection position : N 
            cx = (-a * (-a*N1[x] + N1[y]) +     (N0[x] + a * N0[y]) )/ (1 + a*a)
            cy = (     (-a*N1[x] + N1[y]) + a * (N0[x] + a * N0[y]) )/ (1 + a*a)
            N=[-1, 0, 0, 0]
            N[x] = cx
            N[y] = cy
            distance = abs(A*N0[x]+N0[y]+C) / sqrt(A*A+1)
        else: 
            distance = abs(N0[x] - N1[x])
            N=[-1, 0, 0, 0]
            N[x] = N1[x]
            N[y] = N0[y]
        if onlydist ==1: 
            return distance
        else: 
            return distance, N 
def Surfaces(OutEdges, Node, OffsetLeftRight, TreadElset, AllElements):
    npn = np.array(Node)

    Press=[]; RimContact=[]; TreadToRoad=[]
    EOffset = NOffset = OffsetLeftRight
    Offset=[EOffset, NOffset]

    low = 100000000.0
    startNode = 0
    nextedge = []
    edgeNo = 0
    i = 0
    opposite = 0
    tmpY=0
    while i < len(Node):
        if Node[i][3] < low and Node[i][0] > Offset[0] and Node[i][2] > 0:
            low = Node[i][3]
            startNode = Node[i][0]
            tmpY = Node[i][2]
        i += 1
        if i > 100000:
            print ("[INPUT] Cannot Find the 1st Node for Pressure")
            return Press, RimContact, TreadToRoad

    for i in range(len(Node)):
        if low == Node[i][3] and tmpY == -Node[i][2]:
            opposite = Node[i][0]
            break

    i = 0
    while i < len(OutEdges):
        if OutEdges[i][0] == startNode:
            edgeNo = i
            break
        i += 1
        if i > 100000:
            print ("[INPUT] Cannot Find the 1st Edge for Pressure")
            return Press, RimContact, TreadToRoad

    
    # print ('**', OutEdges[i]); count=0
    method = 1
    if method ==1 : 
        for edge in OutEdges: 
            if "IL" in edge[2] or "L11" in edge[2] or "HUS" in edge[2] or "RIC" in edge[2]:  
                Press.append(edge)

        MAXY = 0
        MINY = 100000000.0
        YS=[]
        for i in range(len(AllElements)):
            if AllElements[i][5] == 'BEAD_R' or AllElements[i][5] == 'BEAD_L' or AllElements[i][5] == 'BD1':
                # print ("BD1 Elements", AllElements[i])
                for j in [1, 2, 3, 4]:
                    if AllElements[i][j] != '' and  AllElements[i][j] != 0:
                        ix = np.where(npn[:,0]==AllElements[i][j])[0][0]
                        YS.append(abs(npn[ix][2]))
        YS = np.array(YS)
        MINY = np.min(YS)

        i = 0
        while i < len(Press): 
            if Press[i][2] == "RIC" or Press[i][2] == "HUS": 
                ix = np.where(npn[:,0]==Press[i][0])[0][0]; n1=npn[ix]
                ix = np.where(npn[:,0]==Press[i][1])[0][0]; n2=npn[ix]

                if (abs(n1[2]) > MINY or abs(n2[2]) > MINY) : 
                    del(Press[i])
                    continue 
            i += 1



    elif method ==0 :  
        i = edgeNo
        Press.append(OutEdges[edgeNo])
        count = 0
        while OutEdges[i][1] != opposite:
            print(OutEdges[i])
            nextedge = FindNextEdge(i, OutEdges)
            i = nextedge[0]
            # print ('**', OutEdges[i])
            Press.append(OutEdges[i])
            if count > 2 and Press[len(Press)-2][4]  == OutEdges[i][4]:
                break
            if count > 1000:  # in case of infinite loop!!
                break
            count+=1
        #    print ('No of press', len(Press))
        # ADD Bead Base Edge as Pressure Surface
        #    print (AllElements[0])

        MAXY = 0
        MINY = 100000000.0
        YS=[]
        for i in range(len(AllElements)):
            if AllElements[i][5] == 'BEAD_R' or AllElements[i][5] == 'BEAD_L' or AllElements[i][5] == 'BD1':
                # print ("BD1 Elements", AllElements[i])
                for j in [1, 2, 3, 4]:
                    if AllElements[i][j] != '' and  AllElements[i][j] != 0:
                        ix = np.where(npn[:,0]==AllElements[i][j])[0][0]
                        YS.append(abs(npn[ix][2]))
        YS = np.array(YS)
        MAXY = np.max(YS)
        MINY = np.min(YS)
        AVGY = (MAXY + MINY) / 2.0
        # print ("find No", AVGY, MAXY, MINY)

        iNext = nextedge[0]
        nextedge = FindNextEdge(iNext, OutEdges)
        # ValueY = find_z(OutEdges[nextedge[0]][0], Node)
        ix = np.where(npn[:,0]==OutEdges[nextedge[0]][0])[0][0]
        ValueY = npn[ix][2]
        c=0
        # while abs(ValueY) < AVGY:
        while abs(ValueY) < MINY:
            # print ('C=', c, OutEdges[nextedge[0]])
            Press.append(OutEdges[nextedge[0]])
            iNext = nextedge[0]
            nextedge = FindNextEdge(iNext, OutEdges)
            ix = np.where(npn[:,0]==OutEdges[nextedge[0]][1])[0][0]
            ValueY = npn[ix][2]
            c += 1
            if c > 100000:
                print ('[INPUT] Cannot Find the next Pressure Edge (Right)')
                return Press, RimContact, TreadToRoad


        previousedge = FindPreviousEdge(edgeNo, OutEdges)
        # ValueY = find_z(OutEdges[previousedge[0]][1], Node)
        ix = np.where(npn[:,0]==OutEdges[previousedge[0]][1])[0][0]
        ValueY = npn[ix][2]
        # print ("AVG=%.2f, MIN=%.2f, Current Y=%.2f"%(AVGY*1000, MINY*1000, ValueY*1000))
        c=0
        # while abs(ValueY) < AVGY:
        while abs(ValueY) < MINY:
            # print (OutEdges[previousedge[0]])
            # print ("AVG=%.2f, MIN=%.2f, Current Y=%.2f"%(AVGY*1000, MINY*1000, ValueY*1000))
            Press.append(OutEdges[previousedge[0]])
            iNext = previousedge[0]
            previousedge = FindPreviousEdge(iNext, OutEdges)
            ix = np.where(npn[:,0]==OutEdges[previousedge[0]][0])[0][0]
            ValueY = npn[ix][2]
            c += 1
            if c > 100000:
                print ('[INPUT] Cannot Find the next Pressure Edge (Left)')
                return Press, RimContact, TreadToRoad

    #    print ('No of press', len(Press))

    if len(Press) < 1:
        print ('[INPUT] No Surface was created for Inner Pressure')
        return Press, RimContact, TreadToRoad
        # logline.append(['ERROR::PRE::[INPUT] No Surface was created for Inner Pressure\n'])
    else:
        print ('* All Edges for Pressure are searched.')
        # logline.append(['* All Surfaces for Pressure are searched. \n'])

    i = 0
    while i < len(OutEdges):
        if OutEdges[i][2] == 'HUS' or OutEdges[i][2] == 'RIC':
            ipress = 0
            if ipress == 0:
                RimContact.append(OutEdges[i])
                # print "Rim Contact", OutEdges[i]
        i += 1
        if i > 100000:
            print ('[INPUT] Cannot Find the Next Outer Edges ')
            return Press, RimContact, TreadToRoad

    #############################################
    ## ADD 5 edges of BSW
    #############################################
    NoOfAddingEdge = 5
    for i in range(len(RimContact)):
        for j in range(len(OutEdges)):
            if RimContact[i] == OutEdges[j]:
                break
        ne = FindNextEdge(j, OutEdges)
        m=ne[0]
        if OutEdges[m][2] == 'BSW':
            RimContact.append(OutEdges[m])
            for j in range(NoOfAddingEdge-1):
                ne= FindNextEdge(m, OutEdges)
                n=ne[0]
                if OutEdges[n][2] == 'BSW':
                    RimContact.append(OutEdges[n])
                m = n
        ne = FindPreviousEdge(j, OutEdges)
        m = ne[0]
        if OutEdges[m][2] == 'BSW':
            RimContact.append(OutEdges[m])
            for j in range(NoOfAddingEdge-1):
                ne= FindPreviousEdge(m, OutEdges)
                n = ne[0]
                if OutEdges[n][2] == 'BSW':
                    RimContact.append(OutEdges[n])
                m = n
    ###############################################

    if len(RimContact) < 1:
        print ('ERROR::PRE::[INPUT] No Surface was created for Rim Contact  ')
        return Press, RimContact, TreadToRoad
        # logline.append(['ERROR::PRE::[INPUT] No Surface was created for Rim Contact\n'])
    else:
        print ('* All Edges for Rim Contact are searched.')
        # logline.append(['* All Surfaces for Rim Contact are searched. \n'])
    
     ## npn .. node number... , AllElements
    ix = np.where(npn[:,2] <0.1E-03)[0]
    ix1 = np.where(npn[:,2] >-0.1E-03)[0]
    ix = np.intersect1d(ix, ix1)
    mincenter = np.min(npn[ix,3]) - 5.0E-03

    for edge in OutEdges: 
        if edge[2] == "BSW" or edge[2] == "CTR" or edge[2] == "CTB" or edge[2] == "SUT" or edge[2] == "UTR" or edge[2] == "TRW" or edge[2] == "BTT" : 
            ix = np.where(npn[:,0] == edge[0])[0][0]; n1 = npn[ix]
            ix = np.where(npn[:,0] == edge[1])[0][0]; n2 = npn[ix]
            if n1[3] >= mincenter or n2[3] >= mincenter: 
                 TreadToRoad.append(edge)

    if len(TreadToRoad) < 1:
        print ('[INPUT] No Surface was created for Road Contact')
        return Press, RimContact, TreadToRoad
        # logline.append(['ERROR::PRE::[INPUT] No Surface was created for Road Contact\n'])
    else:
        print ('* All Edges for Road Contact are searched.')
        # logline.append(['* All Surfaces for Road Contact are searched. \n'])

    return Press, RimContact, TreadToRoad

def FindNextEdge(refEdge, Edges):
    tmpi = refEdge;
    connected = []
    for m in range(len(Edges)):
        if tmpi != m:
            if Edges[tmpi][1] == Edges[m][0]:
                connected.append(m)
    return connected
def FindPreviousEdge(refEdge, Edges):
    tmpi = refEdge;
    connected = []
    for m in range(len(Edges)):
        if tmpi != m:
            if Edges[tmpi][0] == Edges[m][1]:
                connected.append(m)
    return connected
def ConnectedEdge(edge, edges, exclude=[]): 
    con = []
    for e in edges: 
        if e[0] == edge[0] and e[1] == edge[1]: 
            continue 
        if e[1] == edge[0] or e[0] == edge[0]: 
            exc= 0 
            for ex in exclude: 
                if ex[0] == e[0] and ex[1] == e[1] : exc = 1
            if exc == 0:  con.append(e)
        if e[1] == edge[1] or e[0] == edge[1]: 
            exc= 0 
            for ex in exclude: 
                if ex[0] == e[0] and ex[1] == e[1] : exc = 1
            if exc == 0:  con.append(e)

    return con 
class StdoutRedirect(QtCore.QObject):
    printOccur = QtCore.pyqtSignal(str, str, name="print")
 
    def __init__(self, *param):
        QtCore.QObject.__init__(self, None)
        self.daemon = True
        self.sysstdout = sys.stdout.write
        self.sysstderr = sys.stderr.write
 
    def stop(self):
        sys.stdout.write = self.sysstdout
        sys.stderr.write = self.sysstderr
 
    def start(self):
        sys.stdout.write = self.write
        sys.stderr.write = lambda msg : self.write(msg, color="red")
 
    def write(self, s, color="black"):
        sys.stdout.flush()
        self.printOccur.emit(s, color)


def writeworkingdirectory(readfile, dfile='pdir.dir'): 
    cwd=''
    drs = readfile.split("/")
    for i, dr in enumerate(drs): 
        cwd += dr + '/'
        if i == len(drs) -2 : break 
    f= open(dfile, "w")
    f.write(cwd)
    f.close()

    return cwd 

class EDGE:
    def __init__(self):
        self.Edge = []

    def Help(self):
        print ("*********************************************************************************")
        print ("EDGE : Node1, Node2, Elset_Name, FacdID, Element_No, D")
        print (" D : -1= Edge, 0 = Free Edge, -2 = not Free Edge, 1= outer edges, Above 1(2~) = Tie No")
        print ("***************************************************************************")
        print ("** Related Function *******************************************************")
        print ("** - self.Add([]) -> Add a component")
        print ("** - self.Print() -> Print all edges on the screen")
        print ("** - self.Save(file_name) -> save all edges in a file")
        print ("** - edge.Image(Node, Image_File_Name, ImageSize(dpi)) -> save a image")
        print ("***************************************************************************")

    def Add(self, edge):
        self.Edge.append(edge)
    def Sort(self, item=0, reverse=False):
        ## item < 0 sort by connection
        ## item = N sort by Nth number
        if item < 0:
            starts = []
            if reverse == False:
                for i, iedge in enumerate(self.Edge):
                    start = 1
                    for j, jedge in enumerate(self.Edge):
                        if i != j and iedge[0] == jedge[1]:
                            start = 0
                    if start ==1:
                        starts.append(i)

                sortedEdge=EDGE()
                for n in starts:
                    sortedEdge.Add(self.Edge[n])
                    for i, iedge in enumerate(self.Edge):
                        for j, jedge in enumerate(self.Edge):
                            if sortedEdge.Edge[len(sortedEdge.Edge)-1][1] == jedge[0]:
                                sortedEdge.Add(jedge)
                                break
            else:
                for i, iedge in enumerate(self.Edge):
                    start = 1
                    for j, jedge in enumerate(self.Edge):
                        if i != j and iedge[1] == jedge[0]:
                            start = 0
                    if start ==1:
                        starts.append(i)

                sortedEdge=EDGE()
                for n in starts:
                    sortedEdge.Add(self.Edge[n])
                    for i, iedge in enumerate(self.Edge):
                        for j, jedge in enumerate(self.Edge):
                            if sortedEdge.Edge[len(sortedEdge.Edge)-1][0] == jedge [1]:
                                sortedEdge.Add(jedge)
                                break
            for i, edge in enumerate(sortedEdge.Edge):
                self.Edge[i]=edge

        else:
            sortedEdge = EDGE()
            for i, edge in enumerate(self.Edge):
                sortedEdge.Add(edge)
                if i == 0:
                    continue
                else:
                    I = len(sortedEdge.Edge)
                    for j, sedge in enumerate(sortedEdge.Edge):
                        if reverse == True:
                            if sedge[item] < edge[item]:
                                del(sortedEdge.Edge[I-1])
                                sortedEdge.Edge.insert(j, edge)
                                I = j 
                                break
                        else:
                            if sedge[item] > edge[item]:
                                del(sortedEdge.Edge[I-1])
                                sortedEdge.Edge.insert(j, edge)
                                I = j 
                                break


        for i, edge in enumerate(sortedEdge.Edge):
            self.Edge[i] = edge
        del(sortedEdge)    
    
class NODE:
    def __init__(self):
        self.Node = []

    def Add(self, d):
        self.Node.append(d)

    
    def NodeByID(self, n, SORT=0, **args):
        for key, value in args.items():
            if key == 'sort':
                SORT=int(value)
        N = len(self.Node)
        if SORT ==1:
            sorted(self.Node, key=lambda x:x[0])
        if N>100:
            if self.Node[int(N / 4)][0] > n:
                k1 = 0
                k2 = int(N / 4)
            elif self.Node[int(N / 2)][0] > n:
                k1 = int(N / 4)
                k2 = int(N / 2)
            elif self.Node[int(3 * N / 4)][0] > n:
                k1 = int(N / 2)
                k2 = 3 * int(N / 4)
            else:
                k1 = 3 * int(N / 4)
                k2 = N
            for i in range(k1, k2):
                if self.Node[i][0] == n:
                    return self.Node[i]
        for i in range(N):
            if self.Node[i][0] == n:
                return self.Node[i]
        # print ("Cannot Find Node (%d)"%(n))
        NullList = [0, 0.0, 0.0, 0.0]
        return NullList
    def Sort(self, item=0, reverse=False):
        tmpNode = NODE()
        try:
            arr = self.Node[:, item]
        except:
            npNode = np.array(self.Node)
            arr = npNode[:, item]
        if reverse == False: args = np.argsort(arr)
        else:                args = np.argsort(arr)[::-1]
        for nd in self.Node:    tmpNode.Add(nd) 
        sortedlist = []
        for i, arg in enumerate(args):
            self.Node[i] = tmpNode.Node[int(arg)]
            sortedlist.append(tmpNode.Node[int(arg)])
        del(tmpNode)
        sortedlist = np.array(sortedlist)
        return sortedlist

    def DeleteDuplicate(self):
        npary = np.array(self.Node)
        uniques = np.unique(npary)
        N = len(uniques)
        for i, nd in enumerate(self.Node): 
            if i < N :             nd = uniques[i]
        i = N 
        while i < len(self.Node): del(self.Node[i])
    def NodeIDByCoordinate(self, PO, v, closest=0, **args):
    
        N = len(self.Node)
        
        if closest != 0:
            min = 1000000000.0
            
            if PO == 'x' or PO == 'X':
                for i in range(N):
                    if abs(self.Node[i][1]-v) < min:
                        min = self.Node[i][1]
                        ClosestNode = self.Node[i][0]
            elif PO == 'y' or PO == 'Y':
                for i in range(N):
                    if abs(self.Node[i][2]-v) < min:
                        min = self.Node[i][2]
                        ClosestNode = self.Node[i][0]
            elif PO == 'z' or PO == 'Z':
                for i in range(N):
                    if abs(self.Node[i][3]-v) < min:
                        min = self.Node[i][3]
                        ClosestNode = self.Node[i][0]
            else:
                print ("* Check INPUT x/y/z - you input %s"%(PO))
            return ClosestNode
        else: 
            IDs = []
            if PO == 'x' or PO == 'X':
                for i in range(N):
                    if self.Node[i][1] == v:
                        IDs.append(self.Node[i][0])
            elif PO == 'y' or PO == 'Y':
                for i in range(N):
                    if self.Node[i][2] == v:
                        IDs.append(self.Node[i][0])
            elif PO == 'z' or PO == 'Z':
                for i in range(N):
                    if self.Node[i][3] == v:
                        IDs.append(self.Node[i][0])
            else:
                print ("* Check INPUT x/y/z - you input", PO)
            if len(IDs) == 0:
                print ("* Matching Node (", PO, ":", v,")was not found!!")
            return IDs
    
class ELSET:
    def __init__(self):
        self.Elset = []

    def Help(self):
        print ("*************************************************************")
        print ("** [[Elset1, E11, E12, ..], [Elset2, E21, E22, ...], ...]")
        print ("*************************************************************")
        print ("** Related Functions")
        print ("** self.AddName(name) -> Add a new ELSET Name")
        print ("** self.AddNumber(n, name) -> Add a member ")

    def AddName(self, name):
        exist = 0
        for i in range(len(self.Elset)):
            if self.Elset[i][0] == name:
                exist = 1
                break
        if exist == 0:
            self.Elset.append([name])

    def AddNumber(self, n, name):
        for i in range(len(self.Elset)):
            if self.Elset[i][0] == name:
                self.Elset[i].append(n)
    def Add(self, n, name): 
        ## check if the elset name exists or not. 
        fm = -1 
        for i, eset in enumerate(self.Elset): 
            if eset[0].upper() == name.upper(): 
                fm = i 
                break 
        if fm == -1: 
            self.Elset.append([name, n])
        else: 
            f = 0 
            for i, en in enumerate(self.Elset[fm]): 
                if i ==0: continue 
                if n == en: 
                    f = 1 
                    break 
            if f ==0: 
                self.Elset[fm].append(n)

    def Print(self):
        for eset in self.Elset:
            print (eset)

class ELEMENT:
    def __init__(self):
        self.Element = []
        
    def Help(self):
        print ("*************************************************************************************")
        print ("** Element_ID, Node1, Node2, Node3, Node4, Elset_Name, How_Many_Nodes,")
        print ("**             Area/Length, Center_Y(Lateral), Center_Z(Vertical)")
        print ("** Related Functions *****************************************************")
        print ("** - self.Add([]) => Add a ITEM[]")
        print ("** - self.Delete(n) => delete a Element ID 'n'")
        print ("** - self.AddItem(n, v) => Add a item 'v' to Element ID 'n'")
        print ("** - self.DeleteItem(n, j) => delete 'j'th item from Element ID 'n'")
        print ("** - self.Print() => print Element on the screen")
        print ("** - self.Save(File_Name) => Save Element in a file")
        print ("** - self.ElementByID(n) => return Element List of Element ID 'n'")
        print ("** - self.ElementsBySetname('Name') => return a list of element IDs with the same ELSET Name")
        print ("** - self.ElsetNames() => return a list with Elset Names of the class")
        print ("** - self.ChangeSetNameByElset(name1, name2) => Change Elset name from name1 to name2")
        print ("** - self.ChangeSetNameByID(n, Name) => Change Elset name of Element ID 'n' to Name")
        print ("** - self.Image(Node(Class), Image_file_Name, Image_Size)")
        print ("** - self.Elset('Elset_Name') --> create a Element Class of 'Elset_Name'")
        print ("** - self.AllEdge() -> create a EDGE Class from this Element class")
        print ("** - self.FreeEdge() -> return a Edge Class with free edges")
        print ("** - self.OuterEdge(Node) -> Return a Edge class with Outer edges")
        print ("** - self.TieEdge(Node) --> Return a Edge Class with Tie edges")
        print ("*************************************************************************************")

    def Add(self, e):
        self.Element.append(e)
    def Delete(self, n):
        N = len(self.Element)
        for i in range(N):
            if self.Element[i][0] == n:
                del (self.Element[i])
                break
                
    def Sort(self, item=0, reverse=False):
        sortedElement = ELEMENT()
        for i, element in enumerate(self.Element):
            sortedElement.Add(element)
            if i == 0:
                continue
            else:
                I = len(sortedElement.Element)
                for j, selement in enumerate(sortedElement.Element):
                    if reverse == True:
                        if selement[item] < element[item]:
                            del(sortedElement.Element[I-1])
                            sortedElement.Element.insert(j, element)
                            I = j 
                            break
                    else:
                        if selement[item] > element[item]:
                            del(sortedElement.Element[I-1])
                            sortedElement.Element.insert(j, element)
                            I = j 
                            break
        for i, element in enumerate(sortedElement.Element):
            self.Element[i] = element
        del(sortedElement)

    def Tuple(self):
        I = len(self.Element)
        J = len(self.Element[0])
        LST = []
        for i in range(I):
            for j in range(J):
                if j == 5:
                    LST.append(0)
                if self.Element[i][6] == 2 and j == 3:
                    LST.append(0)
                if self.Element[i][6] == 2 and j == 4:
                    LST.append(0)
                if self.Element[i][6] == 3 and j == 4:
                    LST.append(0)
                LST.append(self.Element[i][j])
        
        tuplelst = tuple(LST)
        
        return tuplelst, J
    
    def ElsetToEdge(self, SetName):
        SetElement = self.Elset(SetName)
        
        mEdge=EDGE()
        N = len(SetElement.Element)
        if N == 0: 
            return mEdge
        if SetElement.Element[0][6] == 2:
            for i in range(N): 
                mEdge.Add([SetElement.Element[i][1], SetElement.Element[i][2], SetElement.Element[i][5], 'S0', SetElement.Element[i][0], SetElement.Element[i][7]])
                
        else:
            for i in range(N):
                mEdge.Add([SetElement.Element[i][1], SetElement.Element[i][2], SetElement.Element[i][5], 'S1', SetElement.Element[i][0], 0])
                mEdge.Add([SetElement.Element[i][2], SetElement.Element[i][3], SetElement.Element[i][5], 'S2', SetElement.Element[i][0], 0])
                if SetElement.Element[i][6] == 3: 
                    mEdge.Add([SetElement.Element[i][3], SetElement.Element[i][1], SetElement.Element[i][5], 'S3', SetElement.Element[i][0], 0])
                if SetElement.Element[i][6] == 4: 
                    mEdge.Add([SetElement.Element[i][3], SetElement.Element[i][4], SetElement.Element[i][5], 'S3', SetElement.Element[i][0], 0])
                    mEdge.Add([SetElement.Element[i][4], SetElement.Element[i][1], SetElement.Element[i][5], 'S4', SetElement.Element[i][0], 0])
        return mEdge
        
    def DeleteDuplicate(self):
        i = 0
        while i < len(self.Element):
            j = i + 1
            while j < len(self.Element):
                if self.Element[i][0] == self.Element[j][0]:
                    del (self.Element[j])
                    j += -1
                j += 1
                if j >= len(self.Element):
                    break
            i += 1
            if i >= len(self.Element):
                break
    def Combine(self, element):
        N=len(element.Element)
        for i in range(N): 
            self.Add(element.Element[i])
            
    
    def Print(self):
        N = len(self.Element)
        J = len(self.Element[0])
        for i in range(N):
            line = ''
            for j in range(J):
                if j != J-1:
                    line += str(self.Element[i][j]) + ", "
                else:
                    line += str(self.Element[i][j]) 
            print (line)

    def Save(self, file='ELEMENT.txt'):
        N = len(self.Element)
        f = open(file, "w")
        fline = []
        for i in range(N):
            text = str(self.Element[i]) + '\n'
            fline.append([text])
        f.writelines('%s' % str(item[0]) for item in fline)
        f.close()

    def AddItem(self, n, d):
        N = len(self.Element)
        for i in range(N):
            if self.Element[i][0] == n:
                self.Element[i].append(d)
                break

    def DeleteItem(self, n, j):
        N = len(self.Element)
        for i in range(N):
            if self.Element[i][0] == n:
                del (self.Element[i][j])
                break


    def ElementByID(self, n):
        N = len(self.Element)
        # self.Element.sort()
        if self.Element[int(N / 4)][0] > n:
            k1 = 0
            k2 = int(N / 4)
        elif self.Element[int(N / 2)][0] > n:
            k1 = int(N / 4)
            k2 = int(N / 2)
        elif self.Element[3 * int(N / 4)][0] > n:
            k1 = int(N / 2)
            k2 = 3 * int(N / 4)
        else:
            k1 = 3 * int(N / 4)
            k2 = N
        for i in range(k1, k2):
            if self.Element[i][0] == n:
                return self.Element[i]
        for i in range(N):
            if self.Element[i][0] == n:
                return self.Element[i]
        print ("Cannot Find Element (%d)"%(n))
        NullList = [0, 0, 0, 0, 0, '', 0, 0.0, 0.0, 0.0]
        return NullList

    def ElementBySetname(self, Name):
        SetList = []
        N = len(self.Element)
        for i in range(N):
            if self.Element[i][5] == Name:
                SetList.append(self.Element[i][0])
        return SetList

    def Elset(self, SetName):
        ESet = ELEMENT()
        N = len(self.Element)
        for i in range(N):
            if self.Element[i][5] == SetName:
                ESet.Add(self.Element[i])
        return ESet
    def Nodes(self, **args):

        Node = NODE()
        for key, value in args.items():
            if key == 'Node' or key == 'node':
                Node = value


        I = len(self.Element)
        NL = []
        
        for i in range(I):
            for j in range(1, self.Element[i][6]+1):
                NL.append(self.Element[i][j])
        
        i=0
        while i < len(NL):
            j = i+1
            while j < len(NL):
                if NL[i] == NL[j]:
                    del(NL[j])
                    j -= 1
                j+=1
            i+= 1
        
        if len(Node.Node) == 0:
            # print ("# Node IDs in ELSET : Node_ID=[]")
            return NL

        I = len(NL)
        NC = NODE()
        for i in range(I):
            N = Node.NodeByID(NL[i])
            NC.Add(N)
        # print ("# Nodes in ELSET : Node=NODE()")
        return NC

    def ElsetNames(self):
        Names = []
        N = len(self.Element)

        for i in range(N):
            isSet = 0
            M = len(Names)
            for j in range(M):
                if Names[j] == self.Element[i][5]:
                    isSet = 1
                    break
            if isSet == 0:
                Names.append(self.Element[i][5])
        Names.sort()
        return Names
    def NodesInElement(self, node):
        eNode = NODE()
        N = len(self.Element)
        for i in range(N):
            for j in range(1, self.Element[i][6]+1): 
                eNode.Add(node.NodeByID(self.Element[i][j]))
        eNode.DeleteDuplicate()
        return eNode
                
    def AllEdge(self):
        Name = EDGE()
        N = len(self.Element)
        for i in range(N):
            if self.Element[i][6] == 4:
                Name.Add([self.Element[i][1], self.Element[i][2], self.Element[i][5], 'S1', self.Element[i][0], -1])
                Name.Add([self.Element[i][2], self.Element[i][3], self.Element[i][5], 'S2', self.Element[i][0], -1])
                Name.Add([self.Element[i][3], self.Element[i][4], self.Element[i][5], 'S3', self.Element[i][0], -1])
                Name.Add([self.Element[i][4], self.Element[i][1], self.Element[i][5], 'S4', self.Element[i][0], -1])
            elif self.Element[i][6] == 3:
                Name.Add([self.Element[i][1], self.Element[i][2], self.Element[i][5], 'S1', self.Element[i][0], -1])
                Name.Add([self.Element[i][2], self.Element[i][3], self.Element[i][5], 'S2', self.Element[i][0], -1])
                Name.Add([self.Element[i][3], self.Element[i][1], self.Element[i][5], 'S3', self.Element[i][0], -1])
        return Name

    def FreeEdge(self):
        edges = self.AllEdge()
        freeEdge = FreeEdge(edges)
        return freeEdge

    def OuterEdge(self, Node):
        FEdges = self.FreeEdge()
        # print(len(FEdges.Edge))
        OEdges = OuterEdge(FEdges, Node, self)
        return OEdges

def Mesh2DInformation(InpFileName, comments=True):
    

    with open(InpFileName) as INP:
        lines = INP.readlines()

    Node = NODE()
    Element = ELEMENT()
    Elset = ELSET()
    surface =[]
    isurface =[]

    surface =[]
    tie = []
    itie =[]

    Comments = []
    iComment = 0

    x = 2; y=3
    xy = 23
    spt = ''

    rims = []
    rimtemp = []

    TIRE_SIZE = ''
    TIRE_PATTERN =''
    TIRE_RIM_WIDTH = ''
    TIRE_RIM_DIA = ''
    TIRE_GD=''
    TIRE_TDW=''
    TIRE_OD =''
    TIRE_BSD =""
    TIRE_BT_LIFT=""
    TIRE_ReBELT=''
    TIRE_GROUP=''

    for line in lines:
        if "**" in line:
            if "SIZE" in line.upper() and not "MESH" in line.upper():  TIRE_SIZE= line.strip()
            if "PATTERN" in line.upper(): TIRE_PATTERN= line.strip()
            if "LAYOUT RIM WIDTH" in line.upper() and ":" in line.upper() : TIRE_RIM_WIDTH= line.strip()
            if "RIM DIA" in line.upper() and ":" in line.upper() : TIRE_RIM_DIA= line.strip()
            if "GROOVE DEPTH" in line.upper(): TIRE_GD= line.strip()
            if "TREAD DESIGN WIDTH" in line.upper(): TIRE_TDW = line.strip() 
            if "CAVITY_OD" in line.upper(): TIRE_OD= line.strip()
            if "BEAD SET DISTANCE" in line.upper(): TIRE_BSD= line.strip()
            if "BELT LIFT RATIO" in line.upper(): TIRE_BT_LIFT= line.strip()
            if "REINFORCEMENT BELT" in line.upper(): TIRE_ReBELT= line.strip()
            if "CLASS_CODE" in line.upper(): TIRE_GROUP= line.strip()
            # print(line.strip())
        elif "\n" == line: 
            continue 

        else:
            iComment += 1
            
            if "*" in line: 
                # print(line.strip())

                word = list(line.split(','))

                command = list(word[0].split('*'))
                if "*HEADING" in line.upper(): 
                    spt = 'HD'
                if "*STEP" in line.upper():
                    break  
                elif "*NODE" in line.upper() and not "OUTPUT" in line.upper(): 
                    spt = 'ND'
                elif "*ELEMENT" in line.upper() and "TYPE" in line.upper() : 

                    ename = ''
                    
                    words = line.split(",")
                    for word in words: 
                        if "TYPE" in word.upper(): 
                            EL = word.split("=")[1].strip()
                        if "ELSET" in word.upper(): 
                            ename = word.split("=")[1].strip()
                    
                    if EL == 'MGAX1':
                        spt = 'M1'
                    elif 'CGAX3'  in EL or 'CAX3'  in EL: 
                        spt = 'C3'
                    elif 'CGAX4' in EL or 'CAX4'  in EL: 
                        spt = 'C4'
                    else:
                        spt = 'NN'

                    # print ("$$", line.strip(), word, EL, spt)
                    
                elif "*SURFACE" in line.upper(): 
                    if "INTERACTION" in line.upper() or "BEHAVIOR" in line.upper() or "PROPERTY" in line.upper(): 
                        spt = ''
                    elif "SEGMENT" in line: 
                        spt = 'SF_SEG'
                        if len(rimtemp) > 0: rims.append(rimtemp)
                        rimtemp = []
                    else:
                        try: 
                            spt = 'SF'
                            name = word[2].split('=')[1].strip()
                            if len(isurface) > 0: surface.append(isurface)
                            isurface = [name]
                        except:
                            spt = ''
                            if comments: print ("* Error to reasd surface:", line.strip())
                elif "*TIE" in line.upper(): 
                    try: 
                        spt = 'TI'
                        for wd in word: 
                            if 'NAME' in wd.upper(): 
                                name = wd.split('=')[1].strip()
                                break 
                        if len(itie) > 0: tie.append(itie)
                        itie = [name]
                    except:
                        spt = ''
                        if comments: print ("* Error to reasd TIE:", line.strip())
                elif "*ELSET" in line.upper(): 
                    words = line.split(",")
                    try:
                        for wd in words: 
                            if "ELSET" in wd.upper() and not "*" in wd: 
                                name = wd.split('=')[1].strip()
                                Elset.AddName(name)
                        spt = 'ES'
                    except:
                        spt = ''
                        if comments: print("* Error to read elset:", line.strip())

                elif "*NSET" in line.upper(): 
                    
                    try:
                        spt = 'NS'
                        name = word[1].split('=')[1].strip()
                    except:
                        spt = ''
                        if comments: print("* Error to read Nset", line.strip())
                    

                elif "*INCLUDE" in line.upper(): 
                    spt = 'INC'
                    words = line.split(",")
                    incfile = ''
                    for word in word: 
                        if "INPUT" in word.upper(): 
                            incfile = word.split("=")[1].strip()
                            break 
                    if isfile(incfile): 
                        Node, sel, sset, t_rims = MeshInclude(incfile, Node, Elset)
                    else: 
                        list_file = InpFileName.split("/")
                        fN = len(list_file)
                        workingdirectory=''
                        for dm, dr in enumerate(list_file): 
                            if dm!= fN -1: 
                                workingdirectory+=dr + "/"
                        incfile = workingdirectory + incfile 
                        if isfile(incfile): Node, sel, sset, t_rims = MeshInclude(incfile, Node, Elset)
                    if isfile(incfile): 
                        if comments: print ("Including ", incfile)
                        for se in sel.Element: 
                            Element.Add(se)
                        if len(t_rims) > 0: 
                            for  rm in t_rims: 
                                rims.append(rm)
                else:
                    spt = ''
            else:
                word = list(line.split(','))
                for i, w in enumerate(word): 
                    word[i] =w.strip()

                if spt == 'HD':
                    pass
                if spt == 'ND':
                    try: 
                        if len(word) ==4: 
                            Node.Add([int(word[0]), float(word[3]), float(word[2]), float(word[1])])
                        elif len(word) == 5: ## 5th item : temperature 
                            Node.Add([int(word[0]), float(word[3]), float(word[2]), float(word[1]), float(word[4])])
                    except:
                        if comments: print ("* Error to read node:", line.strip())

                    
                if spt == 'M1':
                    # Element   [EL No,                  N1,          N2,  N3, N4,'elset Name', N,  Area/length, CenterX, CenterY]
                    try: 
                        N1 = Node.NodeByID(int(word[1]))
                        N2 = Node.NodeByID(int(word[2]))
                        C = [[N1[x], N1[y]], [N2[x], N2[y]]]
                        Element.Add([int(word[0]), int(word[1]), int(word[2]), 0, 0, ename, 2, math.sqrt(math.pow(N1[2] - N2[2], 2) + math.pow(N1[3] - N2[3], 2)+ math.pow(N1[1] - N2[1], 2)), (N1[2] + N2[2]) / 2.0, (N1[3] + N2[3]) / 2.0, C, 0])
                    except:
                        # if comments: print ("* Error to read MGAX1:", line.strip())
                        pass 
                if spt == 'C3':
                    try:
                        A, C = Area([int(word[1]), int(word[2]), int(word[3])], Node)
                        Element.Add([int(word[0]), int(word[1]), int(word[2]), int(word[3]), 0, ename, 3, A[0], A[1], A[2], C, 0])
                    except:
                        # if comments: print ("* Error to read CGAX3:", line.strip())
                        pass 
                if spt == 'C4':
                    # print (line.strip())
                    # try:
                        # print (line.strip())
                        A, C = Area([int(word[1]), int(word[2]), int(word[3]), int(word[4])], Node)
                        # print(A, C)
                        Element.Add([int(word[0]), int(word[1]), int(word[2]), int(word[3]), int(word[4]), ename, 4, A[0], A[1], A[2], C, 0])
                        # print (Element.Element[-1])
                    # except:
                    #     # if comments: print ("* Error to read CGAX4:", line.strip())
                    #     pass 

                if spt == 'NS':
                    pass
                if spt == 'ES':
                    try: 
                        for wd in word: 
                            if wd != "":  Elset.AddNumber(int(wd), name)
                    except: 
                        if comments: 
                            print ("* Discarding the Elset %s"%(name))
                            print (" %s"%(line.strip()))


                if spt == 'SF':
                    try: 
                        isurface.append([int(word[0]), word[1].strip()])
                    except: 
                        if comments: 
                            print ("* Discarding the Elset %s"%(name))
                            print (" %s"%(line.strip()))
                if spt == 'SF_SEG': 
                    rimtemp.append(word)

                if spt == 'TI':
                    try:
                        itie.append([word[0].strip(), word[1].strip()])
                    except:
                        if comments: print ("* Error to read Surface:", line.strip())
                else:
                    pass
    # print ("STEP 1")
    
    
    # if len(Element.Element) == 0: 
    #     if comments: print ("# NO DATA in FILE!")
    #     return Node, Element, Elset, surface, tie, xy, rims
    if len(rimtemp) > 0: rims.append(rimtemp)
    rimtemp = [] 
    if len(isurface) > 0: surface.append(isurface)
    if len(itie) > 0: tie.append(itie)

    
    # print ("STEP 2")
    i = 0
    while i < len(surface):
        for j in range(i+1, len(surface)):
            if j < len(surface) : 
                if surface[i][0] == surface[j][0]: 
                    for k in range(1, len(surface[j])):
                        surface[i].append(surface[j][k])
                    del(surface[j])
        i += 1
        if i > 10000: break 
    # print ("STEP 3")
    ename = ''
    for eset in Elset.Elset: 
        f = 0 
        for cmpn in TireComponents: 
            if eset[0].upper() == cmpn: 
                f = 1
                break 
        if f ==1: 
            for i, es in enumerate(eset): 
                if i == 0: 
                    ename = es
                for j, el in enumerate(Element.Element): 
                    if el[0] == es: 
                        Element.Element[j][5] = ename
                        Element.Element[j][11] = 1
    # print ("STEP *")
    for el in Element.Element: 
        if el[11] ==0: 
            Elset.Add(el[0], el[5])
        el = [el[0], el[1], el[2], el[3], el[4], el[5],el[6],el[7],el[8], el[9], el[10]]
    # print ("STEP **")
    x = 0; y=0
    x0=0; y0=0; z0=0
    xp=0; yp=0; zp=0
    xn=0; yn=0; zn=0

    iz = 3; ix=1

    # 
    # ixs = np.where(npn[:,3]<=0.001)[0]
    # N1 = len(ixs)
    # ixs = np.where(npn[:,1]<=0.001)[0]
    # N2 = len(ixs)

    # if N1 < N2 : 
    #     ix = 3; iz = 1 
    # else: 
    #     ix = 1; iz = 3

    for i, nd in enumerate(Node.Node):
        if nd[1] != 0: x0 = 1
        if nd[2] != 0: y0 = 1
        if nd[3] != 0: z0 = 1
        if nd[1] > 0: xp = 1
        if nd[2] > 0: yp = 1
        if nd[3] > 0: zp = 1
        if nd[1] < 0: xn = -1
        if nd[2] < 0: yn = -1
        if nd[3] < 0: zn = -1
        Node.Node[i][iz] = math.sqrt(nd[1]**2 + nd[3]**2)
        Node.Node[i][ix] = 0.0
    
    for i, el in enumerate(Element.Element):
        if el[6] == 3: 
            A, C = Area([el[1], el[2], el[3]], Node, XY=xy)
            # Element.Add([el[0], el[1], el[2], el[3], '', el[5], 3, A[0], A[1], A[2], C])
            el = [el[0], el[1], el[2], el[3], '', el[5], 4, A[0], A[1], A[2], C]
        elif el[6] ==4: 
            A, C = Area([el[1], el[2], el[3], el[4]], Node, XY=xy)
            el = [el[0], el[1], el[2], el[3], el[4], el[5], 4, A[0], A[1], A[2], C]
        else:
            N1 = Node.NodeByID(el[1])
            N2 = Node.NodeByID(el[2])
            C = [[N1[2], N1[3]], [N2[2], N2[3]]]
            # Element.Add([el[0], el[1], el[2], '', '', el[5], 2, el[7], (N1[2] + N2[2]) / 2.0, (N1[3] + N2[3]) / 2.0, C])
            el = [el[0], el[1], el[2], '', '', el[5], 2, el[7], (N1[2] + N2[2]) / 2.0, (N1[3] + N2[3]) / 2.0, C]
    
    # print ("STEP *****")
    
    if xy != 23: 
        for i, el in enumerate(Element.Element):
            if el[6] == 3: 
                A, C = Area([el[1], el[2], el[3]], Node, XY=xy)
                Element.Add([el[0], el[1], el[2], el[3], '', el[5], 3, A[0], A[1], A[2], C])
            elif el[6] ==4: 
                A, C = Area([el[1], el[2], el[3], el[4]], Node, XY=xy)
                Element.Element[i] = [el[0], el[1], el[2], el[3], el[4], el[5], 4, A[0], A[1], A[2], C]
            else:
                N1 = Node.NodeByID(el[1])
                N2 = Node.NodeByID(el[2])
                C = [[N1[int(x/10)], N1[y]], [N2[int(x/10)], N2[y]]]
                Element.Add([el[0], el[1], el[2], '', '', el[5], 2, el[7], (N1[int(x/10)] + N2[int(x/10)]) / 2.0, (N1[y] + N2[y]) / 2.0, C])
    # print ("STEP *****")
    if comments: 
        if TIRE_GROUP != '': print (TIRE_GROUP)
        if TIRE_SIZE !="": print (TIRE_SIZE)
        if TIRE_PATTERN !="": print (TIRE_PATTERN)
        if TIRE_OD !="": print (TIRE_OD) 
        if TIRE_TDW !="": print (TIRE_TDW)
        if TIRE_GD !="": print (TIRE_GD)

        if TIRE_RIM_WIDTH !="": print (TIRE_RIM_WIDTH)
        if TIRE_RIM_DIA !="": print (TIRE_RIM_DIA)

        if TIRE_BSD !="": print (TIRE_BSD)
        if TIRE_BT_LIFT !="": print (TIRE_BT_LIFT)
        if TIRE_ReBELT != '': print(TIRE_ReBELT)


    td = []
    elements = []
    for el in Element.Element: 
        elements.append([el[0], el[1], el[2], el[3], el[4]])
        if el[5] == "CTB" or el[5] == "SUT" : 
            td.append(el)
    try: 
        elements = np.array(elements)
        elno = elements[:,0]
        elno=np.unique(elno)
        nds = elements[:,1:5]
        nds = nds.flatten()
        nds = np.unique(nds)
        nds = np.delete(nds, 0)
        el2 = np.where(elements[:,3] ==0)[0]
        el4 = np.where(elements[:,4] > 0)[0]
    

        if comments: 
            print ("\n**    NO. of EL=%d, Node=%d"%(len(elements), len(nds))) 
            print ("**    EL 2N=%d, 3N=%d, 4N=%d "%(len(el2), len(elements) - len(el2) - len(el4), len(el4)))
            print ("**    ELEMENTs CAP/SUB=%d (else=%d)"%(len(td), len(elements)-len(td)))
    except: 
        pass 

    npn = np.array(Node.Node)
    i = 0 
    while i < len(Element.Element): 
        ix = np.where(npn[:,0] == Element.Element[i][1])[0]
        if not len(ix): 
            del(Element.Element[i])
            continue 
        ix = np.where(npn[:,0] == Element.Element[i][2])[0]
        if not len(ix): 
            del(Element.Element[i])
            continue 
        if Element.Element[i][3] : 
            ix = np.where(npn[:,0] == Element.Element[i][3])[0]
            if not len(ix): 
                del(Element.Element[i])
                continue 
        if Element.Element[i][4] : 
            ix = np.where(npn[:,0] == Element.Element[i][4])[0]
            if not len(ix): 
                del(Element.Element[i])
                continue

        i += 1


    return Node, Element, Elset, surface, tie, xy, rims
def MeshInclude(incfile, Node, Elset, xy=23): 

    with open(incfile) as INP:
        lines = INP.readlines()
    x = int(xy/10)
    y = int(xy%10)

    cmd = ""
    ename = ''
    esetname = ''
    
    Element = ELEMENT() 
    rims=[]
    rimtemp = []
    for line in lines: 
        if "**" in line: 
            continue 
        if "*" in line: 
            if "*NODE" in line.upper(): 
                cmd = "nd"
            if "*ELEMENT" in line.upper(): 
                if "MGAX1" in line.upper():  cmd = 'el2'
                if "CGAX3" in line.upper():  cmd = 'el3'
                if "CGAX4" in line.upper():  cmd = 'el4'
                ename = ''
                if "ELSET" in line.upper(): 
                    list_texts = line.split(",")
                    for txt in list_texts: 
                        if "ELSET" in txt.upper() :
                            ename = txt.split("=")[1].strip() 
            if "*ELSET" in line.upper(): 
                cmd="es"
                word = line.split(",")
                esetname = word[1].split('=')[1].strip()
                # if esetname != "BETWEEN_BELTS" and esetname != "BD1" and esetname != "BetweenBelts":
                Elset.AddName(esetname)
            if "*SURFACE" in line.upper(): 
                if "INTERACTION" in line.upper() or "BEHAVIOR" in line.upper() or "PROPERTY" in line.upper(): 
                    cmd = ''
                elif "SEGMENT" in line: 
                    cmd = 'SF_SEG'
                    if len(rimtemp) > 0: rims.append(rimtemp)
                    rimtemp = []
                else:
                    cmd = ''
        else: 
            word = list(line.split(','))
            for i, w in enumerate(word): 
                word[i] =w.strip()

            if cmd =="nd": 
                try: tx = float(word[3])
                except: tx = 0.0
                try: ty = float(word[2])
                except: ty = 0.0
                try: tz = float(word[1])
                except: tz = 0.0
                Node.add([int(word[0]), tx, ty, tz])

            if cmd == 'el2':
                # Element   [EL No,                  N1,          N2,  N3, N4,'elset Name', N,  Area/length, CenterX, CenterY]
                N1 = Node.NodeByID(int(word[1]))
                N2 = Node.NodeByID(int(word[2]))
                C = [[N1[x], N1[y]], [N2[x], N2[y]]]
                Element.Add([int(word[0]), int(word[1]), int(word[2]), '', '', ename, 2, math.sqrt(math.pow(N1[2] - N2[2], 2) + math.pow(N1[3] - N2[3], 2)+ math.pow(N1[1] - N2[1], 2)), (N1[2] + N2[2]) / 2.0, (N1[3] + N2[3]) / 2.0, C, 0])
            if cmd == 'el3':
                A, C = Area([int(word[1]), int(word[2]), int(word[3])], Node)
                Element.Add([int(word[0]), int(word[1]), int(word[2]), int(word[3]), '', ename, 3, A[0], A[1], A[2], C, 0])
            if cmd == 'el4':
                A, C = Area([int(word[1]), int(word[2]), int(word[3]), int(word[4])], Node)
                Element.Add([int(word[0]), int(word[1]), int(word[2]), int(word[3]), int(word[4]), ename, 4, A[0], A[1], A[2], C, 0])
            if cmd == 'es': 
                try: 
                    for i in range(len(word)):
                        Elset.AddNumber(int(word[i]), esetname)
                        # print (word[i], "is added to", esetname)
                except: 
                    pass
            if cmd == 'SF_SEG': 
                rimtemp.append(word)

    if len(rimtemp) > 0: rims.append(rimtemp)
    
    return Node, Element, Elset, rims 

def RegenenerateMesh(savefile, meshfile, allnodes, Element, Elset, sdb=False): 

    Element, Elset, LROffset = ChaferDivide(Element.Element, ChaferName, Elset.Elset, allnodes.Node)
    OffsetLeftRight= LROffset

    exist3bt = existEL('BT3', Elset)
    exist4bt = existEL('BT4', Elset)

    BetweenBelts = ['BETWEEN_BELTS']
    Between = []
    if exist3bt == 1:
        Between = FindSolidElementBetweenMembrane('BT1', 'BT2', Element)
        BetweenBelts = BetweenBelts + Between
        print ('* Elset "BETWEEN_BELTS" added (%d).' % (len(Between)))
        Between = FindSolidElementBetweenMembrane('BT2', 'BT3', Element)
        BetweenBelts = BetweenBelts + Between
        print ('* Elset "BETWEEN_BELTS" added (%d).' %(len(Between)))
    else:
        Between = FindSolidElementBetweenMembrane('BT1', 'BT2', Element)
        BetweenBelts = BetweenBelts + Between
        print ('* Elset "BETWEEN_BELTS" added (%d).' %(len(Between)))
    if exist4bt == 1:
        Between = FindSolidElementBetweenMembrane('BT3', 'BT4', Element)
        BetweenBelts = BetweenBelts + Between
        print ('* Elset "BETWEEN_BELTS" added (%d).' %(len(Between)))

    if len(BetweenBelts)>0:
        for i, eset in enumerate(Elset): 
            if eset[0] == 'BETWEEN_BELTS': 
                del(Elset[i])
                break 
        Elset.append(BetweenBelts)
    
    del (BetweenBelts)
    del (Between)

    GeneralElement = ELEMENT()
    sws =  ELEMENT()
    for e in Element:
        if e[5] != "SWS":        GeneralElement.Add(e)
        else:                    sws.Add(e)
    e, message = ElementCheck(GeneralElement.Element, allnodes.Node)

    if len(sws.Element)>0:
        i=0
        while i < len(Element):
            if Element[i][5] == 'SWS': 
                del(Element[i])
                i -=1
            i+= 1

        for i, es in enumerate(Elset):
            if es[0]  == 'SWS': 
                ElsetSWS = es 
                del(Elset[i])

    if e == 0:
        for i in range(len(message)):
            print (message[i])
        return "TIE ERROR"
    
    BDmin, BDMax, Center = BeadWidth(Element, allnodes.Node)
    beadcorefile = "bead.tmp"
    bd=open(beadcorefile, "w")
    bd.writelines("%6.3f" %((BDMax-BDmin)*1000))
    print ("* Bead Core width=%.2f"%((BDMax-BDmin)*1000))
    bd.close()

    MasterEdge, SlaveEdge, OutEdges, CenterNodes, FreeEdges, AllEdges, TieError = TieSurface(Element, allnodes.Node)
    if len(TieError) > 0: return TieError

    if len(MasterEdge) != len(SlaveEdge):
        print ("## Errot to fine Tie : Master(%dEA)/Slave(%dEA)"%(len(MasterEdge), len(SlaveEdge)))
        return "TIE Master/Slave Number is not matched"

    Nslave = []
    for i, me in enumerate(MasterEdge): 
        me[5] = i+1 
        # print (me)
        temp = []
        for se in SlaveEdge[i]: 
            se[5] = i+1 
            temp.append(se)
            # print (" >", se)
        Nslave.append(temp)

    PressureSurface, RimContactSurface, TreadToRoadSurface = Surfaces(OutEdges, allnodes.Node, OffsetLeftRight, TreadElset, Element)

    duplel = ElementDuplicationCheck(Element)

    Write2DFile(savefile, allnodes.Node, Element, Elset, TreadToRoadSurface, PressureSurface, RimContactSurface, MasterEdge, Nslave, OffsetLeftRight, CenterNodes)
    with open(meshfile) as mesh: 
        lines = mesh.readlines()

    for line in lines: 
        if "*INCLUDE" in line: 
            incfile = line.split(",")[1]
            incfile = incfile.split("=")[1].strip()
            if isfile(incfile): 
                with open(incfile) as mesh: 
                    lines = mesh.readlines()
            else: 
                list_file = meshfile.split("/")
                fN = len(list_file)
                workingdirectory=''
                for dm, dr in enumerate(list_file): 
                    if dm!= fN -1: 
                        workingdirectory+=dr + "/"
                incfile = workingdirectory + incfile 

                if isfile(incfile): 
                    with open(incfile) as mesh: 
                        lines = mesh.readlines()
            if isfile(incfile): 
                f = open(savefile, 'a')
                for line in lines: 
                    f.write(line)
                f.close()

    return ""
                    
# RegenenerateMesh(savefile, self.Node.Node, Element, Elset, TreadToRoadSurface, PressureSurface, RimContactSurface, MasterEdge, Nslave, OffsetLeftRight, CenterNodes, Comments=comments)
def Write2DFile(FileName, Node, AllElements, Elset, TreadToRoad, Press, RimContact, MasterEdges, SlaveEdges, Offset, CenterNodes, Comments="", addelset=[]):

    f = open(FileName, 'w')

    fline = []
    for i in range(len(Comments)):
        fline.append([Comments[i]])
    fline.append(['*NODE, SYSTEM=R, NSET=ALLNODES\n'])

    i = 0
    while i < len(Node):
        fline.append(['%10d, %15.6E, %15.6E, %15.6E\n' % (Node[i][0], Node[i][3], Node[i][2], Node[i][1])])
        i += 1
    i = 0
    fline.append(['*ELEMENT, TYPE=MGAX1, ELSET=ALLELSET\n'])
    while i < len(AllElements):
        if AllElements[i][6] == 2:
            fline.append(['%10d, %10d, %10d\n' % (AllElements[i][0], AllElements[i][1], AllElements[i][2])])
        i += 1
    i = 0

    fline.append(['*ELEMENT, TYPE=CGAX3H, ELSET=ALLELSET\n'])
    while i < len(AllElements):
        if AllElements[i][6] == 3:
            fline.append(['%10d, %10d, %10d, %10d\n' % (AllElements[i][0], AllElements[i][1], AllElements[i][2], AllElements[i][3])])
        i += 1
    i = 0
    fline.append(['*ELEMENT, TYPE=CGAX4H, ELSET=ALLELSET\n'])
    while i < len(AllElements):
        if AllElements[i][6] == 4:
            fline.append(['%10d, %10d, %10d, %10d, %10d\n' % (AllElements[i][0], AllElements[i][1], AllElements[i][2], AllElements[i][3], AllElements[i][4])])
        i += 1
    isCH1=0
    isCH2=0
    isBDr=0
    isBDl=0
    for i in range(len(Elset)):
        fline.append(["*ELSET, ELSET=%s\n" %(Elset[i][0])])
        if 'CH1' in Elset[i][0]:
            isCH1=1
        if 'CH2' in Elset[i][0]:
            isCH2=1
        if 'BEAD_R' in Elset[i][0]:
            isBDr =1
        if 'BEAD_L' in Elset[i][0]:
            isBDl =1


        k = 0
        for j in range(1, len(Elset[i])):
            if ((k + 1) % 10 != 0):
                if (k +2) == len(Elset[i]):
                    fline.append(['%8d\n' % (Elset[i][j])])
                else:
                    fline.append(['%8d,' % (Elset[i][j])])
            else:
                fline.append(['%8d\n' % (Elset[i][j])])
            k += 1
    if isCH1 == 1:
        fline.append(['*ELSET,  ELSET=CH1\n'])
        fline.append([' CH1_R, CH1_L\n'])
    if isCH2 == 1:
        fline.append(['*ELSET,  ELSET=CH2\n'])
        fline.append([' CH2_R, CH2_L\n'])

    for i in range(len(addelset)):
        fline.append(["*ELSET, ELSET=%s\n" %(addelset[i][0])])
        k = 0
        for j in range(1, len(addelset[i])):
            if ((k + 1) % 10 != 0):
                if (k +2) == len(addelset[i]):
                    fline.append(['%8d\n' % (addelset[i][j])])
                else:
                    fline.append(['%8d,' % (addelset[i][j])])
            else:
                fline.append(['%8d\n' % (addelset[i][j])])
            k += 1


    i = 0
    fline.append(['*SURFACE, TYPE=ELEMENT, NAME=CONT\n'])
    while i < len(TreadToRoad):
        fline.append(['%6d, %s\n' % (TreadToRoad[i][4], TreadToRoad[i][3])])
        i += 1

    i = 0
    fline.append(['*SURFACE, TYPE=ELEMENT, NAME=PRESS\n'])
    while i < len(Press):
        fline.append(['%6d, %s\n' % (Press[i][4], Press[i][3])])
        i += 1

    i = 0
    fline.append(['*SURFACE, TYPE=ELEMENT, NAME=RIC_R\n'])
    i=0
    while i < len(RimContact):
        for nd in Node:
            if RimContact[i][0] == nd[0]: 
                if nd[2]<0: 
                    fline.append(['%6d, %s\n' % (RimContact[i][4], RimContact[i][3])])
                    break
        i += 1

    i = 0
    fline.append(['*SURFACE, TYPE=ELEMENT, NAME=RIC_L\n'])
    while i < len(RimContact):
        for nd in Node:
            if RimContact[i][0] == nd[0]: 
                if nd[2]>0: 
                    fline.append(['%6d, %s\n' % (RimContact[i][4], RimContact[i][3])])
                    break
        i += 1


    cnt = 0 
    for mst in MasterEdges: 
        cnt += 1
        fline.append(['*SURFACE, TYPE=ELEMENT, NAME=Tie_m' + str(cnt) + '\n'])
        fline.append(['%6d, %s\n' % (mst[4], mst[3])])
        fline.append(['*SURFACE, TYPE=ELEMENT, NAME=Tie_s' + str(cnt) + '\n'])

        for slt in SlaveEdges[cnt-1]: 
            fline.append(['%6d, %s\n' % (slt[4], slt[3])])
        
        fline.append(['*TIE, NAME=TIE_' + str(cnt) + '\n'])
        fline.append(['Tie_s' + str(cnt) + ', ' + 'Tie_m' + str(cnt) + '\n'])
        
    ###########################################################

    if isBDr ==1 and isBDl ==1:
        fline.append(['*ELSET, ELSET=BD1\n BEAD_R, BEAD_L\n'])
    elif isBDr == 1 :
        fline.append(['*ELSET, ELSET=BD1\n BEAD_R\n'])
    elif isBDl ==1:
        fline.append(['*ELSET, ELSET=BD1\n BEAD_L\n'])

    Bdr = [0, 0]
    for i in range(len(AllElements)):
        if AllElements[i][5] == 'BEAD_L':
            Bdr[0] = AllElements[i][1]
        if AllElements[i][5] == 'BEAD_R':
            Bdr[1] = AllElements[i][1]

    if Bdr[0] != 0:
        fline.append(['\n*NSET, NSET=BD_L\n'])
        fline.append([str(Bdr[0])])
    if Bdr[1] != 0:
        fline.append(['\n*NSET, NSET=BD_R\n'])
        fline.append([str(Bdr[1])])


    if Bdr[0] != 0 and Bdr[1] != 0:
        fline.append(['\n*NSET, NSET=BDR\n BD_R, BD_L\n'])
    elif Bdr[1] != 0:
        fline.append(['\n*NSET, NSET=BDR\n BD_R\n'])
    elif Bdr[0] != 0 :
        fline.append(['\n*NSET, NSET=BDR\n BD_L\n'])
    else:
        pass

    fline.append(['*NSET, NSET=CENTER\n'])

    i = 0
    while i < len(CenterNodes):
        if ((i + 1) % 10 == 0):
            fline.append(['%8d\n' % (CenterNodes[i])])
        else:
            fline.append(['%8d,' % (CenterNodes[i])])
        i += 1

    f.writelines('%s' % str(item[0]) for item in fline)

    f.close()

def Edge_ElementBoundary(list_element): 
    ## the group of 3 or 4 node elements 

    edges=EDGE()
    for el in list_element: 
        # edges.append([el[1], el[2], 1, el[0], 0])
        if el[6] == 3: 
            edges.Add([el[1], el[2], 1, el[0], 0])
            edges.Add([el[2], el[3], 2, el[0], 0])
            edges.Add([el[3], el[1], 3, el[0], 0])
        elif el[6] == 4: 
            edges.Add([el[1], el[2], 1, el[0], 0])
            edges.Add([el[2], el[3], 2, el[0], 0])
            edges.Add([el[3], el[4], 3, el[0], 0])
            edges.Add([el[4], el[1], 4, el[0], 0])

    if len(edges.Edge) ==0: 
        # print ("No edges to return")
        return edges
    # print ("No of all edges =%d"%(len(edges.Edge)))

    nedges = np.array(edges.Edge)

    bnd = EDGE()
    for i, ed in enumerate(nedges): 
        if ed[4] == 0: 
            ix1 = np.where(nedges[:,1] == ed[0])[0]
            ix2 = np.where(nedges[:,0] == ed[1])[0]
            ix = np.intersect1d(ix1, ix2)
            if len(ix) == 0: 
                bnd.Add(ed) 
                nedges[i][4] = 1 
            elif len(ix) == 1: 
                nedges[i][4] = -1 
                nedges[ix[0]][4] = -1 
    return bnd 

def Area(NodeList, Node, XY=23, **args):   ## Calculate Area of a polygon
    errorimage = 1
    for key, value in args.items():
        if key == 'xy' :   XY = int(value)
        if key == 'error': errorimage= int(value)
    
    if len(Node.Node) > 0:
    
        ii = int(XY/10)
        jj = int(XY)%10
        
        n = len(NodeList)
        x = []
        y = []
        coordinates = []
        for i in range(n):
            Ni = Node.NodeByID(NodeList[i])
            x.append(Ni[ii])
            y.append(Ni[jj])
            coordinates.append([Ni[ii], Ni[jj]])

        x.append(x[0])
        y.append(y[0])
        A = [0.0, 0.0, 0.0]

        for i in range(n):
            s = x[i] * y[i + 1] - x[i + 1] * y[i]
            A[0] += s
            A[1] += (x[i] + x[i + 1]) * s
            A[2] += (y[i] + y[i + 1]) * s

        A[0] = A[0] / 2.0
        try:
            A[1] = A[1] / A[0] / 6
            A[2] = A[2] / A[0] / 6
        except:
            if errorimage > 0: 
                pass 

            return [0.0, 0.0, 0.0], None 

        if A[0] < 0:
            A[0] = -A[0]
            # print 'Negative Area Calculation! '
        return A, coordinates
    else:
        print ("Length of Node is 0")
        print (NodeList)
        print (Node.Node[0])
        return [0.0, 0.0, 0.0], None 

def Allelsets(filename): 
    with open (filename) as  INP: 
        lines = INP.readlines()

    cmd = "p"
    elsets=[]
    elset =[]
    for line in lines: 
        if "**" in line: 
            continue 
        if "*" in line: 
            if "*ELSET" in line: 
                if len(elset) > 0: 
                    elsets.append(elset)
                    # print (elset)
                    elset = []

                cmd = 'set'
                word = line.split(",")[1]
                name = word.split("=")[1]
                name = name.strip() 
                elset=[name]
            else: 
                if len(elset) > 0: 
                    elsets.append(elset)
                    # print (elset)
                    elset = []
                cmd = 'p'
        else: 
            if cmd == 'set':
                data = line.split(",")
                for d in data: 
                    try: 
                        n = int(d.strip())
                    except:
                        n = 0 
                    if n > 0: 
                        elset.append(n)

    return elsets 

def Generate_nodes_for_OE(element, node, fname=None, outeredge=None):

    npn = np.array(node.Node)

    tp = open(fname, 'w') 

    onode=[]
    for ed in outeredge.Edge: 
        onode.append(ed[0])
        onode.append(ed[1])
    onode = np.array(onode)
    onode = np.unique(onode)

    tp.write("*NODE, NSET=OUTER\n")
    cx = 1; cz=3   ## 
    for on in onode: 
        ix =np.where(npn[:,0]==on)[0][0]
        n = npn[ix]
        tp.write("%10d, %12.6f, %12.6f, %12.6f\n"%(n[0], n[cx], n[2], n[cz]))

    memb=['C01'  , 'CC1', # Carcass Cord 1 
    'C02'  , 'CC2', # Carcass Cord 2
    'C03'  , 'CC3', # Carcass Cord 3 
    'C04'  , 'CC4', # Carcass Cord 4
    'BT1'  , # Belt 1 
    'BT2'  , # Belt 2
    'BT3'  , # Belt 3
    'BT4'  , # Belt 4
    'JFC1' , 'JFC', # Jointless Full Cap 1
    'JFC2' , # Jointless Full Cap 2
    'JEC1' , 'JEC', # Jointless Edge Cap 1
    'OJEC1', # Overlapped Jointless Edge Cap
    'OJFC1', # Overlapped Jointless Full Cap
    'OJEC2', # Overlapped Jointless Edge Cap
    'OJFC2', # Overlapped Jointless Full Cap
    'PK1'  , # Half Bead Packing
    'PK2'  , # Full Bead Packing
    'RFM'  , # Bead S(RFM)
    'FLI'  , # Bead Flipper
    'CH1'  , 'CH1_R', 'CH1_L',  # Steel Chafer 
    'CH2'  , 'CH2_R', 'CH2_L',  # 1st Nylon Chafer
    'CH3'  , 'CH3_R', 'CH3_L',  # 2nd Nylon Chafer
    # 'BDC'  , # bead cover 
    'SPC'  ,  ## spiral coil
    ]

    for mem in memb: 
        el = element.Elset(mem)
        if len(el.Element): 
            tp.write("*NODE, NSET=%s\n"%(mem))
            el_nodes = el.Nodes(node=node)
            for n in el_nodes.Node: 
                tp.write("%10d, %12.6f, %12.6f, %12.6f\n"%(n[0], n[cx], n[2], n[cz]))

    tp.close()


def LayoutToProfile(element, node, output='edge', color='darkgray', lw=0.5, counting=1, meshfile=None, nodeout=False): 
    if isinstance(node, list): npn = np.array(node)
    else:                      npn = np.array(node.Node)

    

    outer0 = element.OuterEdge(node)

    if meshfile and nodeout: 
        Generate_nodes_for_OE(element, node, outeredge=outer0, fname=meshfile[:-4]+"-node_outer.inp")

        
    for el in element.Element: 
        if el[6] ==2: 
            outer0.Add([el[1], el[2], el[5], 0, el[0], 0])
    BDR = element.Elset('BEAD_R')
    BDL = element.Elset('BEAD_L')


    r=[]
    mr = 0 
    if len(BDR.Element) > 0:  
        r=[]
        eBDR = BDR.FreeEdge()
        for ed in eBDR.Edge: 
            outer0.Add(ed)
            ix = np.where(npn[:,0]==ed[0])[0][0]; r.append(npn[ix][3])
            ix = np.where(npn[:,0]==ed[1])[0][0]; r.append(npn[ix][3])
    if len(BDL.Element) > 0:  
        eBDL = BDL.FreeEdge()
        for ed in eBDL.Edge: 
            outer0.Add(ed)
            ix = np.where(npn[:,0]==ed[0])[0][0]; r.append(npn[ix][3])
            ix = np.where(npn[:,0]==ed[1])[0][0]; r.append(npn[ix][3])

    if len(r) > 0: 
        r = np.array(r)

        mr = np.min(r) *1000
        
        if   mr < 170.0:      inch = 13
        elif mr < 185.0:    inch = 14
        elif mr < 198.0:    inch = 15
        elif mr < 210.0:    inch = 16
        elif mr < 224.0:    inch = 17
        elif mr < 230.0:    inch = 17.5
        elif mr < 238.0:    inch = 18
        elif mr < 247.5:    inch = 19
        elif mr < 255.0:    inch = 19.5
        elif mr < 265.0:    inch = 20
        elif mr < 275.0:    inch = 21
        elif mr < 286.0:    inch = 22
        elif mr < 293.0:    inch = 22.5
        elif mr < 305.0:    inch = 23
        elif mr < 315.0:    inch = 24
        else:               inch = round(mr/25.4*2, 1)

        if counting > 0: print ("* Estimated Rim Dia.=", inch, "inch", "\n  Min. BDC R=%.2fmm(D=%.1finch)"%(mr, mr*2/25.4))


    if output =='edge':    
        return outer0, mr  


    elif output == 'line': 
        # polygon = plt.Polygon([[n1[x], n1[y]], [n2[x], n2[y]]], color='black', lw=0.5)
        
        x = 2; y=3 
        polygons = []
        for ed in outer0.Edge: 
            ix = np.where(npn[:,0]==ed[0])[0][0]; n1=npn[ix]
            ix = np.where(npn[:,0]==ed[1])[0][0]; n2=npn[ix]
            # polygon = plt.Polygon([[n1[x], n1[y]], [n2[x], n2[y]]], color=color, lw=lw)
            polygons.append([[[n1[x], n1[y]], [n2[x], n2[y]]], color, lw])

        return polygons, mr 


def Color(elsetname):                ## Define Color set 
    c = ''
    if elsetname == 'CTR' or elsetname == 'CTB':
        c = 'darkgray'
    elif elsetname == 'UTR' or elsetname == 'SUT':
        c = 'lightpink'
    elif elsetname == 'CC1' or elsetname == 'C01' or elsetname == 'C02':
        c = 'lightsalmon'
    elif elsetname == 'CCT':
        c = 'purple'
    elif elsetname == 'BTT':
        c = 'steelblue'
    elif elsetname == 'FIL' or elsetname == 'LBF':
        c = 'green'
    elif elsetname == 'UBF':
        c = 'lightpink'
    elif elsetname == 'IL1' or elsetname == 'L11':
        c = 'y'
    elif elsetname == 'BSW':
        c = 'yellowgreen'
    elif elsetname == 'HUS':
        c = 'steelblue'
    elif elsetname == 'RIC':
        c = 'darkgray'
    elif elsetname == 'SHW':
        c = 'darkcyan'
    elif elsetname == 'BD1':
        c = 'dimgray'
    elif elsetname == 'BDC':
        c = 'black'
    elif elsetname == 'MEMB':
        c = 'black'
    elif elsetname == 'DOT':
        c = 'red'
    elif elsetname == 'PRESS':
        c = 'blue'
    elif elsetname == 'RIM':
        c = 'red'
    elif elsetname == 'TDBASE':
        c = 'aqua'
    elif elsetname == 'TDROAD':
        c = 'coral'
    elif elsetname == 'BDTOP':
        c = 'gray'
    else:
        c = 'silver'
    return c

class Ui_MainWindow(object):

    def __init__(self):
        super().__init__()
        self.setupUi(MainWindow)

    def setupUi(self, MainWindow):
        MainWindow.setObjectName("MainWindow")
        MainWindow.resize(1164, 803)

        icon = QtGui.QIcon()
        icon.addPixmap(QtGui.QPixmap("mesh.ico"), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        MainWindow.setWindowIcon(icon)
        
        self.centralwidget = QtWidgets.QWidget(MainWindow)
        self.centralwidget.setObjectName("centralwidget")
        self.gridLayout = QtWidgets.QGridLayout(self.centralwidget)
        self.gridLayout.setObjectName("gridLayout")
        self.ImageType = QtWidgets.QLabel(self.centralwidget)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.ImageType.sizePolicy().hasHeightForWidth())
        self.ImageType.setSizePolicy(sizePolicy)
        self.ImageType.setMaximumSize(QtCore.QSize(150, 30))
        self.ImageType.setFrameShape(QtWidgets.QFrame.Box)
        self.ImageType.setAlignment(QtCore.Qt.AlignLeading|QtCore.Qt.AlignLeft|QtCore.Qt.AlignVCenter)
        self.ImageType.setObjectName("ImageType")
        self.gridLayout.addWidget(self.ImageType, 2, 1, 1, 1)
        
        self.ImageSize = QtWidgets.QLabel(self.centralwidget)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.ImageSize.sizePolicy().hasHeightForWidth())
        self.ImageSize.setSizePolicy(sizePolicy)
        self.ImageSize.setMaximumSize(QtCore.QSize(300, 30))
        self.ImageSize.setFrameShape(QtWidgets.QFrame.Box)
        self.ImageSize.setAlignment(QtCore.Qt.AlignLeading|QtCore.Qt.AlignLeft|QtCore.Qt.AlignVCenter)
        self.ImageSize.setObjectName("ImageSize")
        self.gridLayout.addWidget(self.ImageSize, 0, 1, 1, 1)
        self.horizontalLayout_4 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_4.setObjectName("horizontalLayout_4")

        self.pushButton_elsetSeq = QtWidgets.QPushButton(self.centralwidget)
        self.pushButton_elsetSeq.setObjectName("pushButton_elsetSeq")
        self.horizontalLayout_4.addWidget(self.pushButton_elsetSeq)
        
        
        self.pushButton_SurfaceSeq = QtWidgets.QPushButton(self.centralwidget)
        self.pushButton_SurfaceSeq.setObjectName("pushButton_SurfaceSeq")
        self.horizontalLayout_4.addWidget(self.pushButton_SurfaceSeq)
        self.gridLayout.addLayout(self.horizontalLayout_4, 7, 1, 1, 2)
        self.verticalLayout = QtWidgets.QVBoxLayout()
        self.verticalLayout.setObjectName("verticalLayout")
        self.Label_fName = QtWidgets.QLabel(self.centralwidget)
        self.Label_fName.setMinimumSize(QtCore.QSize(0, 25))
        self.Label_fName.setMaximumSize(QtCore.QSize(16777215, 25))
        self.Label_fName.setFrameShape(QtWidgets.QFrame.Box)
        self.Label_fName.setText("")
        self.Label_fName.setObjectName("Label_fName")
        self.verticalLayout.addWidget(self.Label_fName)
        self.gridLayout.addLayout(self.verticalLayout, 0, 0, 1, 1)

        self.verticalLayout_5 = QtWidgets.QVBoxLayout()
        self.verticalLayout_5.setObjectName("verticalLayout_5")

        self.pushButton_redraw = QtWidgets.QPushButton(self.centralwidget)
        self.pushButton_redraw.setObjectName("pushButton_redraw")
        # self.gridLayout.addWidget(self.pushButton_redraw, 9, 1, 1, 2)
        self.verticalLayout_5.addWidget(self.pushButton_redraw)

        self.pushButton_search = QtWidgets.QPushButton(self.centralwidget)
        self.pushButton_search.setObjectName("pushButton_search")
        self.verticalLayout_5.addWidget(self.pushButton_search)
        
        
        self.gridLayout.addLayout(self.verticalLayout_5, 5, 1, 1, 2)
        self.verticalLayout_11 = QtWidgets.QVBoxLayout()
        self.verticalLayout_11.setObjectName("verticalLayout_11")
        self.checkBox_Spress = QtWidgets.QCheckBox(self.centralwidget)
        self.checkBox_Spress.setMaximumSize(QtCore.QSize(150, 20))
        self.checkBox_Spress.setObjectName("checkBox_Spress")
        self.verticalLayout_11.addWidget(self.checkBox_Spress)
        self.checkBox_sRicR = QtWidgets.QCheckBox(self.centralwidget)
        self.checkBox_sRicR.setMaximumSize(QtCore.QSize(150, 20))
        self.checkBox_sRicR.setObjectName("checkBox_sRicR")
        self.verticalLayout_11.addWidget(self.checkBox_sRicR)
        self.checkBox_sRicL = QtWidgets.QCheckBox(self.centralwidget)
        self.checkBox_sRicL.setMaximumSize(QtCore.QSize(150, 20))
        self.checkBox_sRicL.setObjectName("checkBox_sRicL")
        self.verticalLayout_11.addWidget(self.checkBox_sRicL)
        self.checkBox_Road = QtWidgets.QCheckBox(self.centralwidget)
        self.checkBox_Road.setMaximumSize(QtCore.QSize(150, 20))
        self.checkBox_Road.setObjectName("checkBox_Road")
        self.verticalLayout_11.addWidget(self.checkBox_Road)
        self.gridLayout.addLayout(self.verticalLayout_11, 4, 2, 1, 1)
        self.verticalLayout_8 = QtWidgets.QVBoxLayout()
        self.verticalLayout_8.setObjectName("verticalLayout_8")

        self.list_widget = QListWidget(self.centralwidget)
        self.list_widget.setGeometry(0,0, 300, 200)
        self.verticalLayout_8.addWidget(self.list_widget)


        self.textBrowser = QtWidgets.QTextBrowser(self.centralwidget)
        self.textBrowser.setObjectName("textBrowser")
        self.verticalLayout_8.addWidget(self.textBrowser)

        self.push_addingPoint = QtWidgets.QPushButton(self.centralwidget)
        self.push_addingPoint.setObjectName("push_addingPoint")
        self.verticalLayout_8.addWidget(self.push_addingPoint)
        self.txt_coordinate_point = QtWidgets.QLineEdit(self.centralwidget)
        self.txt_coordinate_point.setAlignment(QtCore.Qt.AlignCenter)
        self.txt_coordinate_point.setObjectName("txt_coordinate_point")
        self.verticalLayout_8.addWidget(self.txt_coordinate_point)

        self.pushButton_colordepth = QtWidgets.QPushButton(self.centralwidget)
        self.pushButton_colordepth.setObjectName("pushButton_colordepth")
        self.verticalLayout_8.addWidget(self.pushButton_colordepth)
        self.gridLayout.addLayout(self.verticalLayout_8, 10, 1, 1, 2)
        self.horizontalLayout_3 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_3.setObjectName("horizontalLayout_3")
        self.verticalLayout_9 = QtWidgets.QVBoxLayout()
        self.verticalLayout_9.setObjectName("verticalLayout_9")
        self.cdepth_solid = QtWidgets.QLabel(self.centralwidget)
        self.cdepth_solid.setAlignment(QtCore.Qt.AlignCenter)
        self.cdepth_solid.setObjectName("cdepth_solid")
        self.verticalLayout_9.addWidget(self.cdepth_solid)
        self.cdepth_memb = QtWidgets.QLabel(self.centralwidget)
        self.cdepth_memb.setAlignment(QtCore.Qt.AlignCenter)
        self.cdepth_memb.setObjectName("cdepth_memb")
        self.verticalLayout_9.addWidget(self.cdepth_memb)
        self.horizontalLayout_3.addLayout(self.verticalLayout_9)
        self.verticalLayout_10 = QtWidgets.QVBoxLayout()
        self.verticalLayout_10.setObjectName("verticalLayout_10")
        self.Cdepth_solid_val = QtWidgets.QLineEdit(self.centralwidget)
        self.Cdepth_solid_val.setAlignment(QtCore.Qt.AlignCenter)
        self.Cdepth_solid_val.setObjectName("Cdepth_solid_val")
        self.verticalLayout_10.addWidget(self.Cdepth_solid_val)
        self.Cdepth_memb_val = QtWidgets.QLineEdit(self.centralwidget)
        self.Cdepth_memb_val.setAlignment(QtCore.Qt.AlignCenter)
        self.Cdepth_memb_val.setObjectName("Cdepth_memb_val")
        self.verticalLayout_10.addWidget(self.Cdepth_memb_val)
        self.horizontalLayout_3.addLayout(self.verticalLayout_10)


        self.checkBox_EnergyLoss = QtWidgets.QCheckBox(self.centralwidget)
        self.checkBox_EnergyLoss.setMaximumSize(QtCore.QSize(150, 20))
        self.checkBox_EnergyLoss.setObjectName("checkBox_EnergyLoss")
        self.verticalLayout_11.addWidget(self.checkBox_EnergyLoss)

        self.checkBox_SED = QtWidgets.QCheckBox(self.centralwidget)
        self.checkBox_SED.setMaximumSize(QtCore.QSize(150, 20))
        self.checkBox_SED.setObjectName("checkBox_SED")
        self.verticalLayout_11.addWidget(self.checkBox_SED)


        self.gridLayout.addLayout(self.horizontalLayout_3, 11, 1, 1, 2)

        self.horizontalLayout_2 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_2.setObjectName("horizontalLayout_2")

        self.verticalLayout_6 = QtWidgets.QVBoxLayout()
        self.verticalLayout_6.setObjectName("verticalLayout_6")
        self.search_node = QtWidgets.QLabel(self.centralwidget)
        self.search_node.setAlignment(QtCore.Qt.AlignCenter)
        self.search_node.setObjectName("search_node")
        self.verticalLayout_6.addWidget(self.search_node)
        self.search_element = QtWidgets.QLabel(self.centralwidget)
        self.search_element.setAlignment(QtCore.Qt.AlignCenter)
        self.search_element.setObjectName("search_element")
        self.verticalLayout_6.addWidget(self.search_element)
        self.horizontalLayout_2.addLayout(self.verticalLayout_6)
        
        self.verticalLayout_7 = QtWidgets.QVBoxLayout()
        self.verticalLayout_7.setObjectName("verticalLayout_7")
        self.serach_node_id = QtWidgets.QLineEdit(self.centralwidget)
        self.serach_node_id.setAlignment(QtCore.Qt.AlignCenter)
        self.serach_node_id.setObjectName("serach_node_id")
        self.verticalLayout_7.addWidget(self.serach_node_id)
        self.search_el_id = QtWidgets.QLineEdit(self.centralwidget)
        self.search_el_id.setAlignment(QtCore.Qt.AlignCenter)
        self.search_el_id.setObjectName("search_el_id")
        self.verticalLayout_7.addWidget(self.search_el_id)
        self.horizontalLayout_2.addLayout(self.verticalLayout_7)

        self.gridLayout.addLayout(self.horizontalLayout_2, 6, 1, 1, 2)
        
        
        self.verticalLayout_4 = QtWidgets.QVBoxLayout()
        

        self.checkBox_NdNo = QtWidgets.QCheckBox(self.centralwidget)
        self.checkBox_NdNo.setMaximumSize(QtCore.QSize(130, 20))
        self.checkBox_NdNo.setChecked(False)
        self.checkBox_NdNo.setObjectName("checkBox_NdNo")
        self.verticalLayout_4.addWidget(self.checkBox_NdNo)
        self.checkBox_ElNo = QtWidgets.QCheckBox(self.centralwidget)
        self.checkBox_ElNo.setMaximumSize(QtCore.QSize(130, 20))
        self.checkBox_ElNo.setObjectName("checkBox_ElNo")
        self.verticalLayout_4.addWidget(self.checkBox_ElNo)
        self.checkBox_Surface = QtWidgets.QCheckBox(self.centralwidget)
        self.checkBox_Surface.setMaximumSize(QtCore.QSize(130, 20))
        self.checkBox_Surface.setObjectName("checkBox_Surface")
        self.verticalLayout_4.addWidget(self.checkBox_Surface)
        self.checkBox_Tie = QtWidgets.QCheckBox(self.centralwidget)
        self.checkBox_Tie.setMaximumSize(QtCore.QSize(130, 20))
        self.checkBox_Tie.setObjectName("checkBox_Tie")
        self.verticalLayout_4.addWidget(self.checkBox_Tie)

        self.checkBox_Temperaure = QtWidgets.QCheckBox(self.centralwidget)
        self.checkBox_Temperaure.setMaximumSize(QtCore.QSize(150, 20))
        self.checkBox_Temperaure.setObjectName("checkBox_Temperaure")
        self.verticalLayout_4.addWidget(self.checkBox_Temperaure)


        self.gridLayout.addLayout(self.verticalLayout_4, 4, 1, 1, 1)
        self.verticalLayout_2 = QtWidgets.QVBoxLayout()
        self.verticalLayout_2.setObjectName("verticalLayout_2")

        self.gridLayout.addLayout(self.verticalLayout_2, 1, 0, 11, 1)


        self.horizontalLayout_5 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_5.setObjectName("horizontalLayout_5")

        self.checkBox_NodeOut = QtWidgets.QCheckBox(self.centralwidget)
        self.checkBox_NodeOut.setMaximumSize(QtCore.QSize(130, 20))
        self.checkBox_NodeOut.setChecked(False)
        self.checkBox_NodeOut.setObjectName("checkBox_node_out")
        self.horizontalLayout_5.addWidget(self.checkBox_NodeOut)

        self.checkBox_ElsetNode = QtWidgets.QCheckBox(self.centralwidget)
        self.checkBox_ElsetNode.setMaximumSize(QtCore.QSize(130, 20))
        self.checkBox_ElsetNode.setChecked(False)
        # self.checkBox_ElsetNode.setObjectName("checkBox_ElsetNode")
        self.horizontalLayout_5.addWidget(self.checkBox_ElsetNode)

        self.checkBox_ElsetElement = QtWidgets.QCheckBox(self.centralwidget)
        self.checkBox_ElsetElement.setMaximumSize(QtCore.QSize(130, 20))
        self.checkBox_ElsetElement.setChecked(False)
        # self.checkBox_ElsetElement.setObjectName("checkBox_ElsetElement")
        self.horizontalLayout_5.addWidget(self.checkBox_ElsetElement)

        self.gridLayout.addLayout(self.horizontalLayout_5, 8, 1, 1, 2)


        self.horizontalLayout = QtWidgets.QHBoxLayout()
        self.horizontalLayout.setObjectName("horizontalLayout")
        self.pushButton_dotsize = QtWidgets.QPushButton(self.centralwidget)
        self.pushButton_dotsize.setMaximumSize(QtCore.QSize(150, 25))
        self.pushButton_dotsize.setObjectName("pushButton_dotsize")
        self.horizontalLayout.addWidget(self.pushButton_dotsize)
        self.lineEdit = QtWidgets.QLineEdit(self.centralwidget)
        self.lineEdit.setMaximumSize(QtCore.QSize(150, 25))
        self.lineEdit.setAlignment(QtCore.Qt.AlignCenter)
        self.lineEdit.setObjectName("lineEdit")
        self.horizontalLayout.addWidget(self.lineEdit)
        self.gridLayout.addLayout(self.horizontalLayout, 3, 1, 1, 2)
        self.InfoSize = QtWidgets.QLabel(self.centralwidget)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.InfoSize.sizePolicy().hasHeightForWidth())
        self.InfoSize.setSizePolicy(sizePolicy)
        self.InfoSize.setMinimumSize(QtCore.QSize(150, 0))
        self.InfoSize.setMaximumSize(QtCore.QSize(150, 30))
        self.InfoSize.setFrameShape(QtWidgets.QFrame.Box)
        self.InfoSize.setAlignment(QtCore.Qt.AlignLeading|QtCore.Qt.AlignLeft|QtCore.Qt.AlignVCenter)
        self.InfoSize.setObjectName("InfoSize")
        self.gridLayout.addWidget(self.InfoSize, 1, 1, 1, 1)
        self.push_regen_mesh = QtWidgets.QPushButton(self.centralwidget)
        self.push_regen_mesh.setMaximumSize(QtCore.QSize(150, 30))
        self.push_regen_mesh.setObjectName("push_regen_mesh")
        self.gridLayout.addWidget(self.push_regen_mesh, 1, 2, 1, 1)
        self.push_origin = QtWidgets.QPushButton(self.centralwidget)
        self.push_origin.setMaximumSize(QtCore.QSize(150, 30))
        self.push_origin.setObjectName("push_origin")
        self.gridLayout.addWidget(self.push_origin, 2, 2, 1, 1)
        self.push_comparing = QtWidgets.QPushButton(self.centralwidget)
        self.push_comparing.setMaximumSize(QtCore.QSize(150, 30))
        self.push_comparing.setObjectName("push_comparing")
        self.gridLayout.addWidget(self.push_comparing, 0, 2, 1, 1)
        MainWindow.setCentralWidget(self.centralwidget)
        self.menubar = QtWidgets.QMenuBar(MainWindow)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 1164, 21))
        self.menubar.setObjectName("menubar")
        self.menuFILE = QtWidgets.QMenu(self.menubar)
        self.menuFILE.setObjectName("menuFILE")
        MainWindow.setMenuBar(self.menubar)
        self.statusbar = QtWidgets.QStatusBar(MainWindow)
        self.statusbar.setObjectName("statusbar")
        MainWindow.setStatusBar(self.statusbar)
        self.actionOpen = QtWidgets.QAction(MainWindow)
        self.actionOpen.setObjectName("actionOpen")
        self.actionClose = QtWidgets.QAction(MainWindow)
        self.actionClose.setObjectName("actionClose")

        self.actionAxiOpen = QtWidgets.QAction(MainWindow)
        self.actionAxiOpen.setObjectName("actionAxiOpen")

        self.actionTrdOpen = QtWidgets.QAction(MainWindow)
        self.actionTrdOpen.setObjectName("actionTrdOpen")

        self.actionSDBOpen = QtWidgets.QAction(MainWindow)
        self.actionSDBOpen.setObjectName("actionSDBOpen")

        self.menuFILE.addAction(self.actionOpen)
        self.menuFILE.addAction(self.actionClose)
        self.menuFILE.addAction(self.actionAxiOpen)
        self.menuFILE.addAction(self.actionTrdOpen)
        self.menuFILE.addAction(self.actionSDBOpen)
        self.menubar.addAction(self.menuFILE.menuAction())


        #################################################################################
        # self.checkBox_NdNo.stateChanged.connect(self.checkBoxState)
        self.checkBox_NdNo.stateChanged.connect(self.Add_Rem_NodeID)
        self.checkBox_ElNo.stateChanged.connect(self.Add_Rem_ELID)
        self.checkBox_Surface.stateChanged.connect(self.checkBoxState)
        self.checkBox_Tie.stateChanged.connect(self.checkBoxState)
        self.checkBox_sRicR.stateChanged.connect(self.checkBoxState)
        self.checkBox_sRicL.stateChanged.connect(self.checkBoxState)
        self.checkBox_Spress.stateChanged.connect(self.checkBoxState)
        self.checkBox_Road.stateChanged.connect(self.checkBoxState)
        self.checkBox_Temperaure.stateChanged.connect(self.add_results)
        self.checkBox_EnergyLoss.stateChanged.connect(self.add_results)
        self.checkBox_SED.stateChanged.connect(self.add_results)
        self.checkBox_Temperaure.setDisabled(True)
        self.checkBox_EnergyLoss.setDisabled(True)
        self.checkBox_SED.setDisabled(True)

        # self.checkBox_ElsetNode.stateChanged.connect(self.checkElsetNode)
        # self.checkBox_ElsetElement.stateChanged.connect(self.checkElsetElement)

        self.retranslateUi(MainWindow)
        self.actionClose.triggered.connect(MainWindow.close)
        self.pushButton_dotsize.clicked.connect(self.drawDots)
        self.lineEdit.returnPressed.connect(self.drawDots)
        self.pushButton_redraw.clicked.connect(self.Redraw)
        self.push_comparing.clicked.connect(self.AddComparingLayouts)
        # self.pushButton_redraw.clicked.connect(self.DisplayElset.clear)
        # self.pushButton_redraw.clicked.connect(self.DisplaySurface.clear)
        # self.pushButton_redraw.clicked.connect(self.draw2Dmesh)

        self.push_regen_mesh.clicked.connect(self.RegenMesh)
        self.push_origin.clicked.connect(self.ReloadMesh)

        self.list_widget.clearSelection()

        self.list_widget.setSelectionMode(
            QtWidgets.QAbstractItemView.ExtendedSelection
        )
        self.list_widget.itemClicked.connect(self.selectionQt_list)
        # self.list_widget.currentRowChanged.connect(self.addElsetBoundary)

        
        
        self.pushButton_elsetSeq.clicked.connect(self.ConvertToProfile)
        self.pushButton_SurfaceSeq.clicked.connect(self.clearElset)
        

        self.pushButton_search.clicked.connect(self.addSearchObjects)
        self.serach_node_id.returnPressed.connect(self.addSearchObjects)
        self.search_el_id.returnPressed.connect(self.addSearchObjects)
        self.Cdepth_solid_val.returnPressed.connect(self.setColordepth)
        self.Cdepth_memb_val.returnPressed.connect(self.setColordepth)
        self.pushButton_colordepth.clicked.connect(self.setColordepth)

        self.txt_coordinate_point.returnPressed.connect(self.addPoint)
        self.push_addingPoint.clicked.connect(self.addPoint)

        ############################################################# 
        self.figure = myCanvas()
        self.canvas = self.figure.canvas
        self.toolbar = self.figure.toolbar
        self.actionOpen.triggered.connect(self.read2DMeshFile)
        self.actionAxiOpen.triggered.connect(self.callAxi)
        self.actionTrdOpen.triggered.connect(self.callTrd)
        self.actionSDBOpen.triggered.connect(self.callSDB)
        self.verticalLayout_2.addWidget(self.toolbar)
        self.verticalLayout_2.addWidget(self.canvas)
        MainWindow.setCentralWidget(self.centralwidget)  ## QDialog : self.setLayout(self.verticalLayout_2)
        ###### #######################################################

        self.meshfile = ''
        self.meshread = 0
        self.node = NODE()
        self.npn =[]
        self.element = ELEMENT()
        self.elset = ELSET()
        self.surface = []
        self.tie =[]
        self.TieError = []

        self.Init_ND = NODE()
        self.Init_EL =ELEMENT()
        self.Init_elset =[]
        self.Init_surface = []
        self.Init_tie = [] 
        self.CurrentElset=''
        self.ElsetNames=[]
        self.allelsets=[]
        self.npel=[]
        self.rims =[]
        

        self.solidepth = 0.5
        self.membdepth = 0.8

        self.searchnode = []
        self.searchelement = []

        self.xy = 23

        self.countelset=-1; self.elsetshow=0
        self.countsurface=-1; self.surfaceshow=0

        self.layoutcomparing = -1
        self.comparingmesh = ''
        self.layoutcounting = 0  
        self.Profiles=[]
        self.Bead_Min_R = []

        self.axifile =''; self.trdfile=''
        self.tmpaxi='axi.tmp'; self.tmptrd='trd.tmp'
        self.limit_offset = 10000

        self.readSDB=0
        self.Eloss =[]
        self.SED =[]
        self.n3d = []
        self.e3dM=[]
        self.e3dS=[]
        self.rim_center =[]
        self.road =[]
        self.E_node=NODE()
        self.S_node=NODE()

        self.temprature_contour=[]
        self.Eloss_contour=[]
        self.SED_contour=[]
        self.contourDotNum=10


        self.dfile = "pdir.dir"
        cwd =getcwd()
        if isfile(self.dfile) == False: 
            df = cwd + '/' +self.dfile
            self.cwd=writeworkingdirectory(df, dfile=self.dfile)
        else: 
            f = open(self.dfile, 'r')
            line =f.readlines()
            f.close()
            self.cwd =line[0]


        QtCore.QMetaObject.connectSlotsByName(MainWindow)

        self._stdout = StdoutRedirect()
        self._stdout.start()
        self._stdout.printOccur.connect(lambda x : self._append_text(x))


    def selectionQt_list(self): 
        self.figure.Clear_ElsetBoundary()

        items = self.list_widget.selectedItems()

        if self.checkBox_ElsetNode.isChecked() == True: 
            nodeids = 1 
        else: 
            nodeids = 0 
        if self.checkBox_ElsetElement.isChecked() == True: 
            elementids = 1 
        else: 
            elementids = 0 
        
        selected = []
        surf = []
        for name in items:
            if "Surface__" in name.text(): 
                surf.append(name.text()[9:]) 
            for eset in self.allelsets: 
                if eset[0]==name.text(): 
                    selected.append(eset[1:])
                    break 
        selected_memb=[]
        selected_solid=[]
        for elist in selected: 
            for n in elist: 
                for el in self.element.Element: 
                    if n == el[0]: 
                        if el[3] > 0: 
                            selected_solid.append(el)
                        else: 
                            selected_memb.append(el)
                        break 

        selectedsf=[]
        for name in surf: 
            for sf in self.surface: 
                if sf[0] == name: 
                    N = len(sf)
                    for i in range(1, N): 
                        selectedsf.append(sf[i])
        



        nodes = np.array(self.node.Node)
        self.figure.Add_ElementsBoundary(selected_solid, selected_memb, nodes, self.npel, nid=nodeids, eid=elementids, surf=selectedsf)
        

    def addElsetBoundary(self):
        try: 
            row = self.list_widget.currentRow()
            print (" Selected Elset : %s"%(self.ElsetNames[row]), end=" ")
            if self.checkBox_ElsetNode.isChecked() == True: 
                nodeids = 1 
            else: 
                nodeids = 0 
            if self.checkBox_ElsetElement.isChecked() == True: 
                elementids = 1 
            else: 
                elementids = 0 
            
            self.figure.Add_ElsetBoundary(self.ElsetNames[row], self.element, self.node, nodeids, elementids, self.allelsets, self.npel)
        except: 
            print (">> Nothing to show")

    def _append_text(self, msg):
        self.textBrowser.moveCursor(QtGui.QTextCursor.End)
        self.textBrowser.insertPlainText(msg)
        QtWidgets.QApplication.processEvents(QtCore.QEventLoop.ExcludeUserInputEvents)
    def clearElset(self): 
        self.list_widget.clearSelection()
        self.figure.Clear_ElsetBoundary()
    
    def read2DMeshFile(self): 
        self.meshfile, _ = QtWidgets.QFileDialog.getOpenFileName(None, "Select File", self.cwd, "File Open(*.inp *.msh)")
        self.layoutcomparing = -1
        self.textBrowser.clear()
        if self.meshfile and isfile(self.meshfile):
            self.Cdepth_memb_val.setText('0.8')
            self.Cdepth_solid_val.setText("0.5")
            self.call2Dmesh(manualfile=self.meshfile)

    def call2Dmesh(self, manualfile='', pattern=False): 
        self.meshfile = manualfile 
        if self.meshfile and isfile(self.meshfile):
            self.layoutcounting = 1
            
            # print ("READ MESH %s"%(self.meshfile) )
            self.node, self.element, self.elset, self.surface, self.tie, self.xy, self.rims = Mesh2DInformation(self.meshfile)
            self.npn = np.array(self.node.Node)
            
            if len(self.npn[0]) > 4 : self.checkBox_Temperaure.setEnabled(True)
            else: self.checkBox_Temperaure.setEnabled(False)
            self.checkBox_EnergyLoss.setEnabled(False)
            self.checkBox_SED.setEnabled(False)
            
            print ("**** Done reading mesh file")
            # sys.exit()
            self.list_widget.clear()


            if len(self.element.Element) == 0 and len(self.node.Node): 
                self.figure.gplot(self.node, self.element)
                return 

            if len(self.element.Element) ==0: 
                with open(self.meshfile) as MF: 
                    lines = MF.readlines()

                print ("********************")
                print (" Command line in file")
                print ("********************")
                for line in lines: 
                    if "*" in line and not "**" in line: 
                        print (line, end="")
                        
                print ("\n# Error! to read mesh file.")
                return

            # print ("DRAW NORMAL")
        
            self.ElsetNames=[]
            # self.allelsets = Allelsets(self.meshfile)

            self.allelsets  = self.elset.Elset
            npel = []
            for el in self.element.Element:
                el[0] = int(el[0])
                el[1] = int(el[1])
                el[2] = int(el[2])
                if el[3] =="": el[3] = 0 
                if el[4] =="": el[4] = 0 
                el[3] = int(el[3])
                el[4] = int(el[4])
                npel.append([el[0], el[1], el[2], el[3], el[4]])
            self.npel = np.array(npel)

            self.allelsets = sorted(self.allelsets, key=lambda name:name[0])
            for en in self.allelsets:
                if len(en) > 1: 
                    self.ElsetNames.append(en[0])
                    item = QListWidgetItem(en[0])
                    self.list_widget.addItem(item)

            self.surface = sorted(self.surface, key=lambda name:name[0] )
            for sf in self.surface: 
                self.ElsetNames.append("Surface__"+sf[0])
                item = QListWidgetItem("Surface__"+sf[0])
                self.list_widget.addItem(item)
            
            # self.draw2Dmesh(self.meshfile)
            self.Label_fName.setText(self.meshfile)

            self.checkBox_NdNo.setChecked(False)
            self.checkBox_ElNo.setChecked(False)


            self.cwd=writeworkingdirectory(self.meshfile, dfile=self.dfile)

            ## reordering ... 
            temp = []
            for i, sf in enumerate(self.surface):
                temp.append(sf)
                if i == 0:
                    continue
                else:
                    I = len(temp)
                    for j, tmp in enumerate(temp):
                        if len(tmp) < len(sf):
                            del(temp[I-1])
                            temp.insert(j, sf)
                            I = j 
                            break
            for i, sf in enumerate(temp):
                self.surface[i] = sf
            del(temp)

            temp = []
            for i, sf in enumerate(self.elset.Elset):
                temp.append(sf)
                if i == 0:
                    continue
                else:
                    I = len(temp)
                    for j, tmp in enumerate(temp):
                        if len(tmp) < len(sf):
                            del(temp[I-1])
                            temp.insert(j, sf)
                            I = j 
                            break
            for i, sf in enumerate(temp):
                self.elset.Elset[i] = sf
            del(temp)
            


            # text = "## No of Nodes = "+ str(len(self.node.Node)) + ", Elements =" +  str(len(self.element.Element)) + ", Elsets = "+ str(len(self.elset.Elset))
            # print (text)

            self.ImageSize.setText("Nodes = "+ str(len(self.node.Node)))
            self.InfoSize.setText("Elements = "+ str(len(self.element.Element)))
            self.ImageType.setText("Elsets = "+  str(len(self.elset.Elset)))
            # if len(self.node.Node) > 0 and len(self.element.Element) > 0: 
            self.figure.getplotinformation(self.node, self.element, self.elset, self.surface, self.tie,  xy=self.xy)
            self.draw2Dmesh(self.meshfile)

    def callAxi(self): 
        self.axifile, _ = QtWidgets.QFileDialog.getOpenFileName(None, "Select File", self.cwd, "File Open(*.axi)")
        if self.axifile: 
            nodes = LayoutMesh_From_axi(self.axifile, output=self.tmpaxi, limit=self.limit_offset)
            self.call2Dmesh(manualfile = self.tmpaxi)
    def callTrd(self): 
        self.trdfile, _ = QtWidgets.QFileDialog.getOpenFileName(None, "Select File", self.cwd, "File Open(*.trd)")
        if self.trdfile: 
            nodes, sectors = PatternMesh_from_trd(self.trdfile, output=self.tmptrd, limit=self.limit_offset)
            if sectors < 180: 
                ptn = True 
            else: 
                ptn = False 

            self.call2Dmesh(manualfile = self.tmptrd, pattern=ptn)
    def callSDB(self): 
        self.sdbresult, _ = QtWidgets.QFileDialog.getOpenFileName(None, "Select File", self.cwd, "File Open(*.sdb0*)")
        if self.sdbresult: 
            self.temprature_contour=[]
            self.Eloss_contour=[]
            self.SED_contour=[]
            self.sdbmodel=self.sdbresult[:-3] 
            self.cwd=writeworkingdirectory(self.sdbmodel, dfile=self.dfile)
            print ('SDB: %s'%(self.sdbmodel.split("/")[-1]))

            npn, np2d, np3d, rim_road = readSDB(self.sdbmodel) 
            d_npn, deformed_Rim_road, iELD, iSED =SDBResult_READ(self.sdbresult, npn)

            EnergyLoss =[]
            for e, e1, e2 in zip(np3d, iELD[0], iELD[1]): 
                if e1 + e2 < 0: 
                    EnergyLoss.append([e[0], 0.0])
                else: 
                    EnergyLoss.append([e[0], e1 + e2])
            StrainEnergyDensity =[]
            for e, e1, e2 in zip(np3d, iSED[0], iSED[1]): 
                if e1 + e2 < 0: 
                    StrainEnergyDensity.append([e[0], 0.0])
                else: 
                    StrainEnergyDensity.append([e[0], e1 + e2])

            self.rim_center = deformed_Rim_road[0]
            self.road = deformed_Rim_road[1]

            self.Eloss = np.array(EnergyLoss)
            self.SED = np.array(StrainEnergyDensity)
            self.n3d = np.array(d_npn)
            self.e3dM = np.array(np2d)
            self.e3dS = np.array(np3d)

            print (" No. 3D El=%d, Eloss=%d, SED=%d"%(len(self.e3dS), len(self.Eloss ), len(self.SED)))
            self.readSDB =  1 
            self.E_node=NODE()
            self.S_node=NODE()
            self.checkBox_Temperaure.setEnabled(True)
            self.checkBox_EnergyLoss.setEnabled(True)
            self.checkBox_SED.setEnabled(True)


    def RegenMesh(self): 
        self.layoutcomparing = -1
        if self.meshfile =="": 
            print ("\n* Need to load a mesh\n") 
            return 
        if len(self.element.Element) ==0: 
            print ("\n* No mesh information to regenerate\n")
            return 

        filename = self.meshfile.split(".")[0]+"-2d"
        self.savefile, _= QtWidgets.QFileDialog.getSaveFileName(None, "Save files as", filename, "New Mesh(*.inp)")

        if self.savefile:
            self.textBrowser.clear()
            self.Cdepth_memb_val.setText('0.8')
            self.Cdepth_solid_val.setText("0.5")
            self.layoutcounting = 1

             # def WriteRegeneratedMesh(FileName, Node, AllElements, Elset, TreadToRoad, Press, RimContact, MasterEdges, SlaveEdges, Offset, CenterNodes, Comments="", addelset=[]):
            self.TieError  = RegenenerateMesh(self.savefile, self.meshfile, self.node, self.element, self.elset)

            self.node = NODE()
            self.element = ELEMENT()
            self.elset = ELSET()
            self.surface = []
            self.tie =[]

            # Regen = PTN.MESH2D(self.meshfile, self.savefile)

            if isfile(self.savefile) :
                self.node, self.element, self.elset, self.surface, self.tie, self.xy, self.rims = Mesh2DInformation(self.savefile)
                print ("")
            else: 
                self.node, self.element, self.elset, self.surface, self.tie, self.xy, self.rims = Mesh2DInformation(self.meshfile)
                
                print ("\n## Original Mesh reloaded.")
            if len(self.node.Node[0]) > 4 or len(self.temprature_contour): self.checkBox_Temperaure.setEnabled(True)
            else: self.checkBox_Temperaure.setEnabled(False)

            self.list_widget.clear()
            self.ElsetNames=[]
            # self.allelsets = Allelsets(self.savefile)
            self.allelsets  = self.elset.Elset
            
            npel = []
            for el in self.element.Element:
                el[0] = int(el[0])
                el[1] = int(el[1])
                el[2] = int(el[2])
                if el[3] =="": el[3] = 0 
                if el[4] =="": el[4] = 0 
                el[3] = int(el[3])
                el[4] = int(el[4])
                npel.append([el[0], el[1], el[2], el[3], el[4]])
            self.npel = np.array(npel)
            self.allelsets = sorted(self.allelsets, key=lambda name:name[0])
            for en in self.allelsets: 
                self.ElsetNames.append(en[0])
                item = QListWidgetItem(en[0])
                self.list_widget.addItem(item)
            self.surface = sorted(self.surface, key=lambda name:name[0] )
            for sf in self.surface: 
                self.ElsetNames.append("Surface__"+sf[0])
                item = QListWidgetItem("Surface__"+sf[0])
                self.list_widget.addItem(item)
            
            self.checkBox_NdNo.setChecked(False)
            self.checkBox_ElNo.setChecked(False)
            temp = []
            for i, sf in enumerate(self.surface):
                temp.append(sf)
                if i == 0:
                    continue
                else:
                    I = len(temp)
                    for j, tmp in enumerate(temp):
                        if len(tmp) < len(sf):
                            del(temp[I-1])
                            temp.insert(j, sf)
                            I = j 
                            break
            for i, sf in enumerate(temp):
                self.surface[i] = sf
            del(temp)

            temp = []
            for i, sf in enumerate(self.elset.Elset):
                temp.append(sf)
                if i == 0:
                    continue
                else:
                    I = len(temp)
                    for j, tmp in enumerate(temp):
                        if len(tmp) < len(sf):
                            del(temp[I-1])
                            temp.insert(j, sf)
                            I = j 
                            break
            for i, sf in enumerate(temp):
                self.elset.Elset[i] = sf
            del(temp)



            self.figure.getplotinformation(self.node, self.element, self.elset, self.surface, self.tie,  xy=self.xy, add2d=self.TieError)
            self.draw2Dmesh(self.savefile)
            self.Label_fName.setText(self.savefile)

    def ReloadMesh(self): 
        # self.figure = myCanvas()
        # self.canvas = self.figure.canvas
        # self.toolbar = self.figure.toolbar
        self.layoutcomparing = -1
        self.layoutcounting = 1
        self.Cdepth_memb_val.setText('0.8')
        self.Cdepth_solid_val.setText("0.5")
        if self.meshfile =="" or not isfile(self.meshfile) : 
            return 
        self.textBrowser.clear()
        self.node = NODE()
        self.element = ELEMENT()
        self.elset = ELSET()
        self.surface = []
        self.tie =[]
        # self.figure.nodechars=[]
        # self.figure.elchars=[]

        self.node, self.element, self.elset, self.surface, self.tie, self.xy, self.rims = Mesh2DInformation(self.meshfile)
        if len(self.node.Node[0]) > 4 or len(self.temprature_contour): 
            self.checkBox_Temperaure.setEnabled(True)
        else: 
            self.checkBox_Temperaure.setEnabled(False)
        print("")
        if len(self.element.Element) ==0: 
            with open(self.meshfile) as MF: 
                lines = MF.readlines()
            i = 0 
            for line in lines: 
                i += 1 
                if i <30: 
                    print (line)

            print ("# Error to read mesh file.")
            return 
        
        self.list_widget.clear()
        self.ElsetNames=[]
        # self.allelsets = Allelsets(self.meshfile)
        self.allelsets  = self.elset.Elset
        npel = []
        for el in self.element.Element:
            el[0] = int(el[0])
            el[1] = int(el[1])
            el[2] = int(el[2])
            if el[3] =="": el[3] = 0 
            if el[4] =="": el[4] = 0 
            el[3] = int(el[3])
            el[4] = int(el[4])
            npel.append([el[0], el[1], el[2], el[3], el[4]])
        self.npel = np.array(npel)
        self.allelsets = sorted(self.allelsets, key=lambda name:name[0])
        for en in self.allelsets: 
            self.ElsetNames.append(en[0])
            item = QListWidgetItem(en[0])
            self.list_widget.addItem(item)
        self.surface = sorted(self.surface, key=lambda name:name[0] )
        for sf in self.surface: 
            self.ElsetNames.append("Surface__"+sf[0])
            item = QListWidgetItem("Surface__"+sf[0])
            self.list_widget.addItem(item)
        self.textBrowser.clear()
        self.checkBox_NdNo.setChecked(False)
        self.checkBox_ElNo.setChecked(False)
        temp = []
        for i, sf in enumerate(self.surface):
            temp.append(sf)
            if i == 0:
                continue
            else:
                I = len(temp)
                for j, tmp in enumerate(temp):
                    if len(tmp) < len(sf):
                        del(temp[I-1])
                        temp.insert(j, sf)
                        I = j 
                        break
        for i, sf in enumerate(temp):
            self.surface[i] = sf
        del(temp)

        temp = []
        for i, sf in enumerate(self.elset.Elset):
            temp.append(sf)
            if i == 0:
                continue
            else:
                I = len(temp)
                for j, tmp in enumerate(temp):
                    if len(tmp) < len(sf):
                        del(temp[I-1])
                        temp.insert(j, sf)
                        I = j 
                        break
        for i, sf in enumerate(temp):
            self.elset.Elset[i] = sf
        del(temp)

        self.figure.getplotinformation(self.node, self.element, self.elset, self.surface, self.tie,  xy=self.xy)
        self.draw2Dmesh(self.meshfile)
        self.Label_fName.setText(self.meshfile)

    def Redraw(self):
        self.layoutcomparing = -1
        self.layoutcounting = 1
        self.figure.clearplot()
        self.checkBox_NdNo.setChecked(False)
        self.checkBox_ElNo.setChecked(False)

        self.figure.getplotinformation(self.node, self.element, self.elset, self.surface, self.tie,  xy=self.xy)
        self.draw2Dmesh(self.meshfile)
        
    def drawDots(self):
        try: 
            size = float(self.lineEdit.text())
            self.figure.AddNodes(self.node, size)
        except:
            pass 

    def addPoint(self): 
        # try: 
            text = self.txt_coordinate_point.text()
            PT = NODE()
            if ']' in text: 
                pts = text.split("]")
                npts = []
                for pt in pts: 
                    pt = pt.split("[")
                    if len(pt)>1: 
                        pt = pt[1]
                        npts.append(pt)
                
                for pt in npts: 
                    pt = pt.split(",")
                    try: 
                        x = float(pt[0]); y=float(pt[1])
                        if x>100.0 or y>100.0: 
                            x /= 1000; y /=1000
                        PT.Add([0, 0, x, y])
                    except:
                        continue 
            else: 
                pt = text.split(",")
                try: 
                    x = float(pt[0]); y=float(pt[1])
                    if x>100.0 or y>100.0: 
                        x /= 1000; y /=1000
                    PT.Add([0, 0, x, y])
                except:
                    pass 

            size = 5.0
            if len(PT.Node): 
                self.figure.OnlyAddingNode(PT, size)
        # except:
        #     pass 

    def checkBoxState(self):
        self.draw2Dmesh(self.meshfile)

    def sdb_elementValue(self, e_value, vtype='sum', log=False): 
        value_element=[]
        for el in self.element.Element:
            if el[3] ==0: continue  
            # if el[5] == 'CTB' or el[5] == 'CTR' or el[5] == 'SUT' or el[5] == 'UTR': 
            #     td =10**7
            #     # continue 
            # else: td = 0 
            
            sum = 0 
            values=[]
            ix = np.where(e_value[:,0] % self.limit_offset == el[0])[0]
            if log: 
                print ("* Log value ")
                for x in ix: 
                    if e_value[x][1] >0: 
                        sum+=math.log(e_value[x][1], 10)
                        values.append(math.log(e_value[x][1], 10))
            else: 
                for x in ix: 
                    sum+=e_value[x][1]
                    values.append(e_value[x][1])
            if vtype=='sum': 
                eld = sum #* 2*math.pi * el[7] * el[9]
            if vtype == 'max': 
                eld = max(values)

            value_element.append([el[0], el[1], el[2], el[3], el[4], eld]  )
            
        return np.array(value_element)

    def add_results(self): 
        # self.limit_offset = 10000
        # # self.node, self.element, self.elset, self.surface, self.tie, self.xy, self.rims 
        # self.Eloss =[];  self.SED =[];   self.n3d = [];    self.e3dM=[];   self.e3dS=[]
        # self.element.Element = [E_id, N1, N2, N3, N4, 'Name', node_counting, Center_2, Center_3, Area, [[n1x, n1y], ...[n4x, n4y] ], 1]
        # nodes_slave, nodes_master = nodes_Master_Slave_Tie(self.tie, self.surface, self.element.Element)
        size = float(self.lineEdit.text())
        self.np_nodes =  np.array(self.node.Node)
        num = int(float(self.Cdepth_memb_val.text()))
        if num < 10: num  =10 
        if self.readSDB and self.meshfile: 
            if self.checkBox_Temperaure.isChecked():
                
                if not len(self.temprature_contour) or num != self.contourDotNum : 
                    self.contourDotNum = num
                    
                    els=[]
                    for el in self.element.Element: 
                        if el[3] !=0 or el[3] != '': 
                            if el[4] =='': el[4]=0 
                            els.append([el[0], el[1], el[2], el[3], el[4]])
                    nodes=[]
                    tmps=[]
                    for nd in self.node.Node: 
                        ix = np.where(self.n3d[:, 0]%10000==nd[0])[0]
                        if len(ix): 
                            nodes.append([nd[0], nd[1], nd[2], nd[3], self.n3d[ix[0]][4]])
                            tmps.append(self.n3d[ix[0]][4])
                    vmin = min(tmps); vmax=max(tmps)
                    self.np_nodes=np.array(nodes)
                    self.np_els =  np.array(els)
                    print ("** Inner Max Temperature =%.2f"%(vmax))
                    self.temprature_contour = self.figure.plot_temperature(self.np_nodes,  self.np_els, \
                        vmin=vmin*1.05, vmax=vmax*0.95, levels=50, size=size, num=num)
                    
                else: 
                    
                    tns = self.temprature_contour[2].reshape(-1)
                    vmin = np.min(tns)*1.05
                    vmax = np.max(tns)*0.95
                    _, _, _ = self.figure.plot_temperature(self.np_nodes, self.np_els, vmin=vmin*1.05, vmax=vmax*0.95,\
                         levels=50, contour=self.temprature_contour, size=size, num=num)
            
            elif self.checkBox_EnergyLoss.isChecked():
                if len(self.Eloss_contour) == 0 or num != self.contourDotNum:
                    self.Cdepth_solid_val.setText("1.75")
                    vrange = 1.75
                    self.contourDotNum = num 
                    self.ELoss2D =  self.sdb_elementValue(self.Eloss, log=False, vtype='sum')
                    # npEN = np.array(self.E_node.Node, dtype=np.float64)
                    idxs = np.where(self.ELoss2D[:,5]>=0)[0]
                    values = self.ELoss2D[idxs]

                    values = self.ELoss2D 

                    std = np.std(values[:,5])
                    avg =  np.average(values[:,5])
                    print (" E-Loss Standard Deviation=%.3E, AVG=%.3E"%(std, avg))
                    vmin=avg-abs(std*vrange); vmax=avg+abs(std*vrange)
                    if std > avg:
                        std = avg * 0.5  
                    if vmin < 0: vmin = 0 
                    if vmax <= vmin: vmax = vmin + abs(avg*0.01)
                    vmin =0
                    self.Eloss_contour = self.figure.plot_ElementValue(self.np_nodes,  self.ELoss2D, vmin=vmin*1.2, vmax=vmax*0.8, \
                        num=num, levels=100, size=size, col=5, contour=None, cmap='rainbow')
                    self.Cdepth_solid_val.setText(str(vrange))
                    

                else: 
                    vrange = float(self.Cdepth_solid_val.text())
                    idxs = np.where(self.ELoss2D[:,5]>=0)[0]
                    values = self.ELoss2D[idxs]

                    std = np.std(values[:,5])
                    avg =  np.average(values[:,5])
                    # print ("Standard Deviation E-Loss=%.3E, avg=%.3E"%(std, avg))
                    vmin=avg-abs(std*vrange); vmax=avg+abs(std*vrange)
                    if vmin < 0: vmin = 0 
                    if vmax <= vmin: vmax = vmin + abs(avg*0.01)
                    vmin =0
                    try: 
                        _ = self.figure.plot_ElementValue(self.np_nodes,  self.ELoss2D, vmin=vmin*1.2, vmax=vmax*0.8,\
                             num=num, levels=100, size=size, col=5, contour=self.Eloss_contour, cmap='rainbow')
                    except: 
                        print ( " Check scale value ")
                        return 

            elif self.checkBox_SED.isChecked():
                if len(self.SED_contour) == 0 or num != self.contourDotNum:
                    if len(self.Eloss_contour) == 0:
                        self.Cdepth_memb_val.setText('10')
                        num = 10
                    self.contourDotNum = num 
                    self.Cdepth_solid_val.setText("1.75")
                    vrange = 1.75
                    self.SED_2D =  self.sdb_elementValue(self.SED, vtype='max')
                    # npEN = np.array(self.E_node.Node, dtype=np.float64)
                    idxs = np.where(self.SED_2D[:,5]>=0)[0]
                    values = self.SED_2D[idxs]

                    std = np.std(values[:,5])
                    avg =  np.average(values[:,5])
                    print (" SED Standard Deviation =%.3E, AVG=%.3E"%(std, avg))
                    if std > avg:
                        std = avg * 0.5  
                    vmin=avg-abs(std*vrange); vmax=avg+abs(std*vrange)
                    if vmin < 0: vmin = 0 
                    if vmax <= vmin: vmax = avg + abs(avg*0.01)
                    vmin =0
                    self.SED_contour = self.figure.plot_ElementValue(self.np_nodes,  self.SED_2D, vmin=vmin*1.2, vmax=vmax*0.8,\
                         num=num, levels=100, size=size, col=5, contour=None)
                    self.Cdepth_solid_val.setText(str(vrange))
                else: 
                    vrange = float(self.Cdepth_solid_val.text())
                    idxs = np.where(self.SED_2D[:,5]>=0)[0]
                    values = self.SED_2D[idxs]

                    std = np.std(values[:,5])
                    avg =  np.average(values[:,5])
                    # print ("Standard Deviation E-Loss=%.3E, avg=%.3E"%(std, avg))
                    vmin=avg-abs(std*vrange); vmax=avg+abs(std*vrange)
                    if vmin < 0: vmin = 0 
                    if vmax <= vmin: vmax = avg + abs(avg*0.01)
                    vmin =0
                    try: 
                        _ = self.figure.plot_ElementValue(self.np_nodes,  self.SED_2D, vmin=vmin*1.2, vmax=vmax*0.8, \
                            num=num, levels=100, size=size, col=5, contour=self.SED_contour)
                    except: 
                        print ( " Check scale value ")
                        return 

                return 
                
            else: 
                self.Redraw()
                

        elif self.meshfile: 
            if self.checkBox_Temperaure.isChecked():
                size = float(self.lineEdit.text())
                if not len(self.temprature_contour) : 
                    self.np_nodes = np.array(self.node.Node) 
                    
                    els=[]
                    for el in self.element.Element: 
                        if el[3] !=0 or el[3] != '': 
                            if el[4] =='': el[4]=0 
                            els.append([el[0], el[1], el[2], el[3], el[4]])
                    self.np_els =  np.array(els)
                    vmin = min(self.np_nodes[:,4]); vmax=max(self.np_nodes[:,4])
                    
                    print ("** Inner Max Temperature =%.2f"%(vmax))
                    self.temprature_contour = self.figure.plot_temperature(self.np_nodes,  self.np_els, vmin=vmin*1.05, vmax=vmax*0.95, levels=50, size=size)
                    
                else: 
                    tns = self.temprature_contour[2].reshape(-1)
                    vmin = np.min(tns)*1.05
                    vmax = np.max(tns)*0.95
                    _, _, _ = self.figure.plot_temperature(self.np_nodes, self.np_els, vmin=vmin*1.05, vmax=vmax*0.95, levels=50, contour=self.temprature_contour, size=size)
            

    def Add_Rem_NodeID(self): 
        if self.checkBox_NdNo.isChecked(): 
            self.figure.NIDShow(show=1)  
        else: 
            self.figure.NIDShow(show=0)  

    def Add_Rem_ELID(self): 
        if self.checkBox_ElNo.isChecked(): 
            self.figure.EIDShow(show=1) 
        else: 
            self.figure.EIDShow(show=0) 

    def setColordepth(self): 
        try: 
            self.solidepth = float(self.Cdepth_solid_val.text())
            self.membdepth = float(self.Cdepth_memb_val.text())
        
            if self.solidepth>1.0: 
                self.solidepth = 1.0 
                self.Cdepth_solid_val.setText("1.0")
            if self.solidepth<0.00: 
                self.solidepth = 0.05 
                self.Cdepth_solid_val.setText("0.05")
            
            if self.membdepth>1.0: 
                self.membdepth = 1.0 
                self.Cdepth_memb_val.setText("1.0")
            if self.membdepth<0.0:  
                self.membdepth = 0.05
                self.Cdepth_memb_val.setText("0.05")
        
        
            self.draw2Dmesh(self.meshfile)
        except:
            print ("## Inappropriate values are input.")
        
    def AddComparingLayouts(self): 
        ## def LayoutToProfile(element, node, output='edge', color='darkgray', lw=0.5): 
        basic_line_width = 0.5
        if self.layoutcounting ==0: 
            print ("\n* Need to load a mesh first.\n")
            return
        elif self.layoutcounting == 1: 
            self.Bead_Min_R = []
            self.Profiles = []
            profile0, mr0= LayoutToProfile(self.element, self.node, output='line', color=lst_colors[0], lw=basic_line_width, counting=0)
            if mr0 : 
                self.Bead_Min_R.append(mr0)
            self.Profiles.append(profile0)
        CN = len(lst_colors)
        color = lst_colors[self.layoutcounting % CN]
        self.comparingmesh, _ = QtWidgets.QFileDialog.getOpenFileName(None, "Select File", self.cwd, "File Open(*.inp *.msh)")
        if self.comparingmesh: 
            nameadded = self.comparingmesh.split("/")[-1]
            print ("\n## Layout to compare is added ")
            print ("   '%s'\n"%(nameadded))
            self.com_node, self.com_element, self.com_elset, self.com_surface, self.com_tie, self.com_xy, self.rims = Mesh2DInformation(self.comparingmesh)
            print ("")
            if len(self.com_element.Element) ==0 and len(self.com_node.Node):
                self.figure.AddNodes(self.com_node, size=2)
                return 
            elif len(self.com_element.Element) ==0: 
                with open(self.comparingmesh) as CM: 
                    lines = CM.readlines()
                print ("******************************************")
                for line in lines: 
                    if "*" in line and not "**" in line: 
                        print (line, end="")
                
                print ("\n ## Error to load the layout to compare")
                self.figure.Draw_profiles(self.Profiles, self.Bead_Min_R)
                return 
            self.layoutcounting += 1
            profile, mr= LayoutToProfile(self.com_element, self.com_node, output='line', color=color, lw=basic_line_width, counting=self.layoutcounting)
            if mr > 0: 
                self.Bead_Min_R.append(mr)
            self.Profiles.append(profile)
            self.figure.Draw_profiles(self.Profiles, self.Bead_Min_R)
        else: 
            
            self.figure.Draw_profiles(self.Profiles, self.Bead_Min_R)
            ## plotting image 
        
    def Comparing(self):

        self.layoutcomparing = self.layoutcomparing * -1
        if self.layoutcomparing == 1: 
            if len(self.element.Element) ==0: 
                print ("\n* Need to load a mesh first.\n")
                return 

            self.comparingmesh, _ = QtWidgets.QFileDialog.getOpenFileName(None, "Select File", self.cwd, "File Open(*.inp *.msh)")

            if self.comparingmesh and isfile( self.comparingmesh): 
                nameadded = self.comparingmesh.split("/")[-1]
                print ("\n## Layout to compare is added ")
                print ("   '%s'\n"%(nameadded))
                self.com_node, self.com_element, self.com_elset, self.com_surface, self.com_tie, self.com_xy, self.rims = Mesh2DInformation(self.comparingmesh)
                print ("")
                if len(self.com_element.Element) ==0: 
                    with open(self.comparingmesh) as CM: 
                        lines = CM.readlines()
                    print ("******************************************")
                    for line in lines: 
                        if "*" in line and not "**" in line: 
                            print (line, end="")
                    
                    print ("\n ## Error to load the layout to compare")
                    return 

                allouter = self.com_element.OuterEdge(self.com_node)
                comparing_edgeset=[allouter]
                memb = EDGE()
                for el in self.com_element.Element: 
                    if el[6] ==2: 
                        memb.Add([el[1], el[2], el[5], 0, el[0], 0])
                comparing_edgeset.append(memb)

                BD = self.com_element.ElsetToEdge("BD1")
                comparing_edgeset.append(BD)



                allouter = self.element.OuterEdge(self.node)
                initial_edgeset=[allouter]
                memb = EDGE()
                for el in self.element.Element: 
                    if el[6] ==2: 
                        memb.Add([el[1], el[2], el[5], 0, el[0], 0])
                initial_edgeset.append(memb)
                BD = self.element.ElsetToEdge("BD1")
                initial_edgeset.append(BD)

                self.figure.ComparingMode(initial_edgeset, self.node, self.xy, comparing_edgeset, self.com_node, self.com_xy)
            else: 
                try: 
                    
                    del(self.com_node)
                    del(self.com_element)
                    del(self.com_elset)
                    del(self.com_surface)
                    del(self.com_tie)
                    del(self.com_xy)
                    self.draw2Dmesh(self.meshfile) 
                    print ("  Layout to compare is removed")
                except:
                    pass 

    def ConvertToProfile(self): 
        basic_line_width = 0.5
        if self.layoutcounting ==0: 
            print ("\n* Need to load a mesh first.\n")
            return
        self.layoutcomparing = self.layoutcomparing * -1
        if self.layoutcomparing > 0: 
            self.Bead_Min_R = []
            self.Profiles = []
            profile0, mr0= LayoutToProfile(self.element, self.node, output='line', color=lst_colors[0], lw=basic_line_width, counting=0,\
                 meshfile=self.meshfile, nodeout=self.checkBox_NodeOut.isChecked())
            self.Bead_Min_R.append(mr0)
            self.Profiles.append(profile0)
            self.figure.Draw_profiles(self.Profiles, self.Bead_Min_R)
        else: 
            self.Redraw()


    def addSearchObjects(self):
        self.searchnode = []
        self.searchelement = []
        maxnid = 0
        maxeid = 0
        for n in self.node.Node:
            if maxnid < n[0]: 
                maxnid = n[0]
        for n in self.element.Element:
            if maxeid < n[0]: 
                maxeid = n[0]
        text = self.serach_node_id.text()
        if "~" in text or "," in text:
            numbers = text.split(",")
            for nm in numbers:
                nn = nm.split("~")
                if len(nn) == 2:
                    nstart = 0; nend = 0
                    try:               nstart = int(nn[0].strip())
                    except:            nstart = 0
                    try:               nend = int(nn[1].strip())
                    except:            nend = 0
                    
                    if nstart > 0 and nend > 0 : 
                        for k in range(nstart, nend+1) : 
                            self.searchnode.append(k)
                    elif nstart > 0: 
                        for k in range(nstart, maxnid+1) : 
                            self.searchnode.append(k)
                    elif nend > 0: 
                        for k in range(1, nend+1) : 
                            self.searchnode.append(k)
                elif len(nn) == 1:
                    nstart = 0
                    try:         nstart = int(nn[0].strip())
                    except:      nstart = 0
                    if nstart > 0 : self.searchnode.append(nstart)
                elif len(nn) == 1:
                    nstart = 0
                    try:        nstart = int(nn[0].strip())
                    except:     nstart = 0
                    if nstart > 0 : self.searchnode.append(nstart)
        else:
            numbers = text.split(" ")
            if len(numbers): 
                for num in numbers: 
                    try: 
                        self.searchnode.append(int(num))
                    except: 
                        continue 

        text = self.search_el_id.text()
        if "~" in text or "," in text:
            numbers = text.split(",")
            for nm in numbers:
                nn = nm.split("~")
                if len(nn) == 2:
                    nstart = 0; nend = 0
                    try:        nstart = int(nn[0].strip())
                    except:     nstart = 0
                    try:        nend = int(nn[1].strip())
                    except:     nend = 0
                    
                    if nstart > 0 and nend > 0 : 
                        for k in range(nstart, nend+1) : 
                            self.searchelement.append(k)
                    elif nstart > 0: 
                        for k in range(nstart, maxeid+1) : 
                            self.searchelement.append(k)
                    elif nend > 0: 
                        for k in range(1, nend+1) : 
                            self.searchelement.append(k)
                elif len(nn) == 1:
                    nstart = 0
                    try:       nstart = int(nn[0].strip())
                    except:    nstart = 0
                    if nstart > 0 : self.searchelement.append(nstart)
        else:
            numbers = text.split(" ")
            if len(numbers): 
                for num in numbers: 
                    try: 
                        self.searchelement.append(int(num))
                    except: 
                        continue 

        # self.searchnode = int(self.serach_node_id.text())
        # self.searchelement = int(self.search_el_id.text())
        # print (self.searchnode)
        # print ("#####################")
        # print (self.searchelement)
        # self.draw2Dmesh(self.meshfile)
        self.figure.SearchShow(sel=self.searchelement, snode=self.searchnode, node=self.node, element=self.element )

    
    def drawSurface(self): 
        try: 
            self.surfaceshow=1
            self.countsurface +=1
            if self.countsurface >= len(self.surface): self.countsurface -= len(self.surface)
            dsurface = self.surface[self.countsurface]
            delset = self.elset.Elset[self.countelset]
            self.DisplaySurface.setText(dsurface[0])
            self.figure.Displaysets(dsurface, delset, self.element, cdepth=self.solidepth, mdepth=self.membdepth,  xy=self.xy, sshow=self.surfaceshow, eshow=self.elsetshow)
        except:
            pass
    
            
    def drawElset(self): 
        try: 
            self.elsetshow = 1
            self.countelset += 1
            if self.countelset >= len(self.elset.Elset): self.countelset -= len(self.elset.Elset)
            dsurface = self.surface[self.countsurface]
            delset = self.elset.Elset[self.countelset]
            self.DisplayElset.setText(delset[0])
            self.figure.Displaysets(dsurface, delset, self.element, cdepth=self.solidepth, mdepth=self.membdepth,  xy=self.xy, sshow=self.surfaceshow, eshow=self.elsetshow)
        except: 
            pass

    def retranslateUi(self, MainWindow):
        _translate = QtCore.QCoreApplication.translate
        MainWindow.setWindowTitle(_translate("MainWindow", "2D Mesh Viewer"))
        self.ImageSize.setText(_translate("MainWindow", "Nodes"))
        self.InfoSize.setText(_translate("MainWindow", "Elements"))
        self.ImageType.setText(_translate("MainWindow", "Elsets"))
        self.checkBox_NdNo.setText(_translate("MainWindow", "Node No"))
        self.checkBox_ElNo.setText(_translate("MainWindow", "Element No"))
        self.checkBox_NodeOut.setText(_translate("MainWindow", "Save Outer Nodes"))

        self.checkBox_ElsetNode.setText(_translate("MainWindow", "Node No"))
        self.checkBox_ElsetElement.setText(_translate("MainWindow", "Element No"))


        self.checkBox_Surface.setText(_translate("MainWindow", "Surface"))
        self.checkBox_Tie.setText(_translate("MainWindow", "Tie"))
        self.pushButton_dotsize.setText(_translate("MainWindow", "Node size"))
        self.pushButton_dotsize.setShortcut(_translate("MainWindow", "Ctrl+N"))
        self.pushButton_dotsize.setToolTip("Ctrl+N")
        self.lineEdit.setText(_translate("MainWindow", "0"))
        
        
        self.search_node.setText(_translate("MainWindow", "NODE"))
        self.search_element.setText(_translate("MainWindow", "        ELEMENT        "))
        self.serach_node_id.setText(_translate("MainWindow", "0"))
        self.search_el_id.setText(_translate("MainWindow", "0"))
        
        self.pushButton_colordepth.setText(_translate("MainWindow", "set Transparency (0~1)"))
        self.pushButton_colordepth.setShortcut(_translate("MainWindow", "Ctrl+T"))            ## transparency setting 
        self.pushButton_colordepth.setToolTip("Ctrl+T")

        self.push_addingPoint.setText(_translate("MainWindow", "Add Point"))
        self.txt_coordinate_point.setText(_translate("MainWindow", "[x, y],"))
        
        self.cdepth_solid.setText(_translate("MainWindow", "Solid"))
        self.cdepth_memb.setText(_translate("MainWindow", "            Membrane            "))
        self.Cdepth_solid_val.setText(_translate("MainWindow", "0.5"))
        self.Cdepth_memb_val.setText(_translate("MainWindow", "0.8"))
        
        self.pushButton_redraw.setText(_translate("MainWindow", "Redraw"))
        self.pushButton_redraw.setShortcut(_translate("MainWindow", "Ctrl+R"))                ## re-draw original image 
        self.pushButton_redraw.setToolTip("Ctrl+R")

        self.checkBox_Spress.setText(_translate("MainWindow", "Surface Pressure"))
        self.checkBox_sRicR.setText(_translate("MainWindow", "Surface RIC_R"))
        self.checkBox_sRicL.setText(_translate("MainWindow", "Surface RIC_L"))
        self.checkBox_Road.setText(_translate("MainWindow", "Surface Road Contact"))
        self.pushButton_elsetSeq.setText(_translate("MainWindow", "view Mode"))
        self.pushButton_elsetSeq.setShortcut(_translate("MainWindow", "Ctrl+M"))          ## Write mesh 
        self.pushButton_elsetSeq.setToolTip("Ctrl+M")
        
        self.pushButton_SurfaceSeq.setText(_translate("MainWindow", "Delete sets"))
        self.pushButton_SurfaceSeq.setShortcut(_translate("MainWindow", "Ctrl+D"))    ## delete sets from image 
        self.pushButton_SurfaceSeq.setToolTip("Ctrl+D")

        self.push_regen_mesh.setText(_translate("MainWindow", "Save as "))
        self.push_regen_mesh.setShortcut(_translate("MainWindow", "Ctrl+S"))          ## Write mesh 
        self.push_regen_mesh.setToolTip("Ctrl+S")

        self.pushButton_search.setText(_translate("MainWindow", "Find nodes/elements"))
        self.pushButton_search.setShortcut(_translate("MainWindow", "Ctrl+F"))  
        self.pushButton_search.setToolTip("Ctrl+F, sample : 12, 1~20")                 ## Searching 


        self.push_origin.setText(_translate("MainWindow", "reLoad"))
        self.push_origin.setShortcut(_translate("MainWindow", "Ctrl+L"))               ## again read original 
        self.push_origin.setToolTip("Ctrl+L")
        self.push_origin.setDisabled(True)

        self.push_comparing.setText(_translate("MainWindow", "add Profile to compare"))
        self.push_comparing.setShortcut(_translate("MainWindow", "Ctrl+P"))             ## profile add 
        self.push_comparing.setToolTip("Ctrl+P")

        self.menuFILE.setTitle(_translate("MainWindow", "FILE"))
        self.actionOpen.setText(_translate("MainWindow", "Open"))
        self.actionOpen.setShortcut(_translate("MainWindow", "Ctrl+O"))               ## open file 
        self.actionAxiOpen.setText(_translate("MainWindow", "Axi 3D Open"))
        self.actionAxiOpen.setShortcut(_translate("MainWindow", "Shift+a"))
        self.actionTrdOpen.setText(_translate("MainWindow", "Trd 3D Open"))
        self.actionTrdOpen.setShortcut(_translate("MainWindow", "Shift+t"))
        self.actionSDBOpen.setText(_translate("MainWindow", "SDB Result Open"))
        self.actionSDBOpen.setShortcut(_translate("MainWindow", "Shift+O"))   
        self.actionClose.setText(_translate("MainWindow", "Close"))
        self.actionClose.setShortcut(_translate("MainWindow", "Ctrl+Q"))              ## Quit program 


        #### Adding SDB results 
        self.checkBox_Temperaure.setText(_translate("MainWindow", "Temperature"))
        self.checkBox_EnergyLoss.setText(_translate("MainWindow", "E-Loss(Sum)"))
        self.checkBox_SED.setText(_translate("MainWindow", "SED(Max.)"))

        
    def draw2Dmesh(self, fileName, add2d=[] ):
        # print (self.meshfile)
        self.layoutcomparing = -1
        if self.meshfile !="":
            self.surfaceshow = 0 
            self.elsetshow = 0
            self.countelset=-1
            self.countsurface=-1
            
            size = float(self.lineEdit.text())
            if self.checkBox_NdNo.isChecked() == True:           nd = 1
            else:                                                  nd = 0
            if self.checkBox_ElNo.isChecked() == True:           el = 1
            else:                                                  el = 0
            if self.checkBox_Surface.isChecked() == True:        sf = 1
            else:                                                  sf = 0
            if self.checkBox_Tie.isChecked() == True:            ti = 1
            else:                                                  ti = 0 
            if self.checkBox_sRicR.isChecked() == True:          srr = 1
            else:                                                  srr = 0 
            if self.checkBox_sRicL.isChecked() == True:          srl = 1
            else:                                                  srl= 0 
            if self.checkBox_Spress.isChecked() == True:         spr = 1
            else:                                                  spr = 0 
            if self.checkBox_Road.isChecked() == True:         srd = 1
            else:                                               srd = 0

            self.figure.gplot( self.node, self.element, self.elset, self.surface, self.tie,  xy=self.xy, \
                size=size, ni=nd, ei=el, si=sf, ti=ti, srr=srr, srl=srl, spr=spr, srd=srd, cdepth=self.solidepth, \
                mdepth=self.membdepth, snode=self.searchnode, sel=self.searchelement, rims=self.rims)
    
class myCanvas(FigureCanvas):
    def __init__(self, parent=None, *args, **kwargs):
        self.figure = plt.figure()
        FigureCanvas.__init__(self, self.figure)
        self.setParent(parent)
        self.canvas = FigureCanvas(self.figure)
        self.toolbar = NavigationToolbar(self.canvas,self)

        self.nid=[];         self.x=[];        self.y=[];        self.polygons=[]
        self.lines=[]
        self.eid=[];        self.ey=[];         self.ex=[];         self.ecolor=[]
        self.ric =[]; self.press=[]; self.otheredge=[]
        self.mtie =[]; self.stie =[]
        self.ric_r=[]; self.ric_l=[]; self.cont=[]

        self.ax =[]
        self.pei =0; self.pni =0


        self.p1=[]
        self.p2=[]
        self.distance = 0.0 
        self.clicked = 0 
        self.sumLength = 0.0
        
        self.dots=[]
        self.lines=[]
        self.chars=[]
        self.achars=[]
        self.cline=[]
        self.llen=[]
        self.fontsize =  8

        self.mclick=0
        self.circle=[]

        self.points = []
        self.cxs=[]
        self.cys=[]
        self.lxs=[]
        self.lys=[]

        self.snap_mode = 1 

        self.nodechars=[]
        self.elchars = []
        self.searchchars=[]
        self.searchpolygons=[]
        self.id_poly =[]
        self.searchdots =[]
        self.settexts=[]

        self.LineEL = []
        
        self.dots=[]

        self.temperature_mode = False 
        self.temperature_node =[]

    def Area(self, ix=[], iy=[]): 
        x =[]; y=[]
        for px, py in zip(ix, iy):
            x.append(px)
            y.append(py)
        x.append(ix[0]); y.append(iy[0])

        A = [0.0, 0.0, 0.0]

        n = len(x)-1

        for i in range(n):
            s = x[i] * y[i + 1] - x[i + 1] * y[i]
            A[0] += s
            A[1] += (x[i] + x[i + 1]) * s
            A[2] += (y[i] + y[i + 1]) * s

        A[0] = A[0] / 2.0


        return A[0]
    def zoom(self, event, base_scale=1.2): 
        
        cur_xlim = self.ax.get_xlim()
        cur_ylim = self.ax.get_ylim()

        xdata = event.xdata # get event x location
        ydata = event.ydata # get event y location

        if event.button == 'down':
            # deal with zoom in
            scale_factor = 1 / base_scale
        elif event.button == 'up':
            # deal with zoom out
            scale_factor = base_scale
        else:
            # deal with something that should never happen
            scale_factor = 1
            print (event.button)

        new_width = (cur_xlim[1] - cur_xlim[0]) * scale_factor
        new_height = (cur_ylim[1] - cur_ylim[0]) * scale_factor
        try: 
            relx = (cur_xlim[1] - xdata)/(cur_xlim[1] - cur_xlim[0])
            rely = (cur_ylim[1] - ydata)/(cur_ylim[1] - cur_ylim[0])

            #and set limits
            plt.xlim([xdata - new_width * (1-relx), xdata + new_width * (relx)])
            plt.ylim([ydata - new_height * (1-rely), ydata + new_height * (rely)])
        except:
            pass 
        self.figure.canvas.draw_idle()
    def onclick(self, event): 
        if event.dblclick: 
            if event.button == 1: 
                if self.temperature_mode: 
                    self.snap_mode = 1
                else: 
                    self.snap_mode *= -1 
                    if self.snap_mode == 1: print (" >> Snap Mode ON")
                    else: print (" >> Snap Mode OFF")
                

    def onReleased(self, event): 
        
        dotsize = 2 
        if event.button == 2: 
            self.mclick += 1

            if self.snap_mode == 1 : #and self.mclick < 4: 
                ind = []
                indx1 = np.where(self.points[:,2]>=event.xdata-0.01)[0]
                indx2 = np.where(self.points[:,2]<=event.xdata+0.01)[0]
                indx = np.intersect1d(indx1, indx2)
                if len(indx) > 0: 
                    indy1 = np.where(self.points[:,3]>=event.ydata-0.01)[0]
                    indy2 = np.where(self.points[:,3]<=event.ydata+0.01)[0]
                    indy = np.intersect1d(indy1, indy2)
                    if len(indy) > 0: 
                        ind = np.intersect1d(indx, indy) 
                if len(ind) > 0 : 
                    mn = []
                    for ix in ind: 
                        l = math.sqrt( (event.xdata - self.points[ix][2])**2 + (event.ydata - self.points[ix][3])**2 )
                        mn.append([ix, l])
                    mn = np.array(mn)
                    lmin = np.min(mn[:,1]) 
                    lx = -1
                    for l in mn: 
                        if l[1] == lmin: 
                            lx = int(l[0])
                            break 
                    if lx >=0: 
                        tx = self.points[lx][2]
                        ty = self.points[lx][3]
                        # print ("B2", self.points[lx])
                else: 
                    self.mclick -= 1
                    return 
            else:
                tx = event.xdata
                ty = event.ydata 


            if self.mclick == 4: 
                self.cxs=[]
                self.cys=[]
                self.mclick=0
                for dot in self.dots:
                    dot.remove()
                self.dots=[]

                for char in self.chars: 
                    char.set_visible(False)

                for cl in self.circle:
                    cl.remove()

                self.circle=[]
            
            elif self.mclick ==3:
                
                sameposition = 0 
                for p1, p2 in zip(self.cxs, self.cys): 
                    if p1 == tx and p2 == ty: 
                        sameposition = 1 
                        break 
                if sameposition == 1: 
                    self.mclick -= 1 
                    return 
                else:  
                    self.cxs.append(tx)
                    self.cys.append(ty)

                    d, = plt.plot(tx, ty, 'o', color='gray', markersize=dotsize)
                    self.dots.append(d)
                    

                    x1 = self.cxs[0]; x2=self.cxs[1]; x3=self.cxs[2]
                    y1 = self.cys[0]; y2=self.cys[1]; y3=self.cys[2]
                    A = x1*(y2-y3) - y1 *(x2-x3) + x2*y3 - x3*y2
                    B = (x1*x1 + y1*y1)*(y3-y2) +(x2**2 + y2**2)*(y1-y3) + (x3**2+y3**2)*(y2-y1)
                    C = (x1**2 + y1**2)*(x2-x3)+(x2**2+y2**2)*(x3-x1) + (x3*x3 + y3*y3)*(x1-x2)
                    D = (x1*x1 + y1*y1)*(x3*y2-x2*y3)+(x2*x2+y2*y2)*(x1*y3-x3*y1)+(x3*x3+y3*y3)*(x2*y1-x1*y2)
                    SQRT = B*B + C*C - 4*A*D  

                    # print ("A", A)
                    # print ("B", B)
                    # print ("C", C)
                    # print ("D", D)
                    # print ("S", SQRT)

                    if A == 0 or SQRT < 0.0: 
                        print (" The 3 nodes cannot make a circle.")
                        print (" (%10.3E, %10.3E), (%10.3E, %10.3E)\n (%10.3E, %10.3E)"%(self.cxs[0]*1000,self.cys[0]*1000, self.cxs[1]*1000,self.cys[1]*1000, self.cxs[2]*1000,self.cys[2]*1000))

                        for dot in self.dots:
                            dot.remove()
                        self.dots=[]
                        self.cxs=[]
                        self.cys=[]
                        self.mclick=0

                    else: 
                        cx = -B/A/2.0
                        cy = -C/A/2.0
                        
                        R = math.sqrt(SQRT) / 2/abs(A)

                        self.cxs.append(tx)
                        self.cys.append(ty)
                        # d, = plt.plot(cx, cy, 'o', color='red')
                        # self.dots.append(d)

                         
                        ch = plt.text((self.cxs[0]+self.cxs[1])/2.0, (self.cys[0]+self.cys[1])/2.0, "R="+str(round(R*1000, 2)), size=self.fontsize, color='black')
                        self.chars.append(ch)

                        if R > 1E10: 
                            print (" The circle cannot be drawn.")
                            print (" Center = %.3E, %.3E"%(-B/A*500, -C/A*500))
                            print (" Radius = %.3E"%(math.sqrt(SQRT) /abs(A)*500))
                        else: 
                            crcl = plt.Circle((cx, cy), R, color='orange', fill=False, linewidth=dotsize/3)
                            self.ax.add_artist(crcl)
                            self.circle.append(crcl)
                    self.mclick = 0 
                    self.cxs=[]
                    self.cys=[]
                        

            else:
                sameposition = 0 
                for p1, p2 in zip(self.cxs, self.cys): 
                    if p1 == tx and p2 == ty: 
                        sameposition = 1 
                        break 
                if sameposition == 1: 
                    self.mclick -= 1 
                else: 
                    self.cxs.append(tx)
                    self.cys.append(ty)

                    d, = plt.plot(tx, ty, 'o', color='gray', markersize=dotsize)
                    self.dots.append(d)

            self.plot_current_display_range()
            # current_xlim=self.ax.get_xlim()
            # current_ylim=self.ax.get_ylim()
            # plt.xlim(current_xlim[0], current_xlim[1])
            # plt.ylim(current_ylim[0], current_ylim[1])
            # self.figure.canvas.draw_idle()

        elif event.button == 1: 
            self.clicked =0  
            self.mclick = 0
            self.sumLength = 0.0
            
            self.cxs=[]
            self.cys=[]
            self.lxs=[]
            self.lys=[]

            for dot in self.dots: 
                dot.remove()
            for line in self.lines: 
                line.remove()
            for char in self.chars: 
                char.set_visible(False)
            for char in self.achars: 
                char.set_visible(False)
            for cl in self.cline: 
                cl.remove()
            for ch in self.llen:
                ch.set_visible(False)
            for cl in self.circle:
                cl.remove()

            self.circle=[]

            self.dots=[]
            self.lines=[]
            self.chars=[]
            self.achars=[]

            self.cline=[]
            self.llen=[]

            self.plot_current_display_range()

        elif event.button ==3:
            self.clicked += 1

            if self.snap_mode == 1: 
                if self.temperature_mode: 
                    ind = []
                    indx1 = np.where(self.temperature_node[:,2]>=event.xdata-0.01)[0]
                    indx2 = np.where(self.temperature_node[:,2]<=event.xdata+0.01)[0]
                    indx = np.intersect1d(indx1, indx2)
                    if len(indx) > 0: 
                        indy1 = np.where(self.temperature_node[:,3]>=event.ydata-0.01)[0]
                        indy2 = np.where(self.temperature_node[:,3]<=event.ydata+0.01)[0]
                        indy = np.intersect1d(indy1, indy2)
                        if len(indy) > 0: 
                            ind = np.intersect1d(indx, indy) 
                    if len(ind) > 0 : 
                        mn = []
                        for ix in ind: 
                            l = math.sqrt( (event.xdata - self.temperature_node[ix][2])**2 + (event.ydata - self.temperature_node[ix][3])**2 )
                            mn.append([ix, l])
                        mn = np.array(mn)
                        lmin = np.min(mn[:,1]) 
                        lx = -1
                        for l in mn: 
                            if l[1] == lmin: 
                                lx = int(l[0])
                                break 
                        if lx >=0: 
                            tx = self.temperature_node[lx][2]
                            ty = self.temperature_node[lx][3]

                            for dot in self.dots: 
                                dot.remove()
                            for char in self.chars: 
                                char.set_visible(False)
                            self.chars=[]
                            self.dots=[]
                            
                            d, = plt.plot(tx, ty, 'o', color='black')
                            self.dots.append(d)
                            shift = 0.002
                            ch = plt.text(tx+shift, ty+shift, 'Temperature=%.2f'%(self.temperature_node[lx][4]), size=self.fontsize)
                            self.chars.append(ch)

                            self.plot_current_display_range()

                            # current_xlim=self.ax.get_xlim()
                            # current_ylim=self.ax.get_ylim()
                            # plt.xlim(current_xlim[0], current_xlim[1])
                            # plt.ylim(current_ylim[0], current_ylim[1])
                            
                            # self.figure.canvas.draw_idle()
                
                    return 
                else: 
                    ind = []
                    indx1 = np.where(self.points[:,2]>=event.xdata-0.01)[0]
                    indx2 = np.where(self.points[:,2]<=event.xdata+0.01)[0]
                    indx = np.intersect1d(indx1, indx2)
                    if len(indx) > 0: 
                        indy1 = np.where(self.points[:,3]>=event.ydata-0.01)[0]
                        indy2 = np.where(self.points[:,3]<=event.ydata+0.01)[0]
                        indy = np.intersect1d(indy1, indy2)
                        if len(indy) > 0: 
                            ind = np.intersect1d(indx, indy) 
                    if len(ind) > 0 : 
                        mn = []
                        for ix in ind: 
                            l = math.sqrt( (event.xdata - self.points[ix][2])**2 + (event.ydata - self.points[ix][3])**2 )
                            mn.append([ix, l])
                        mn = np.array(mn)
                        lmin = np.min(mn[:,1]) 
                        lx = -1
                        for l in mn: 
                            if l[1] == lmin: 
                                lx = int(l[0])
                                break 
                        if lx >=0: 
                            tx = self.points[lx][2]
                            ty = self.points[lx][3]
                            # print ("B3", self.points[lx])
                    else: 
                        self.clicked -= 1
                        return 

                    
            else: 
                tx = event.xdata; ty = event.ydata

            prev = 0 
            for xs,ys in zip(self.lxs, self.lys): 
                if tx == xs and ty == ys: 
                    prev = 1 
                    break 
            if prev == 1: 
                self.clicked -= 1
                return 

            self.lxs.append(tx)
            self.lys.append(ty)
            
            d, = plt.plot(tx, ty, 'o', color='red', markersize=dotsize)
            self.dots.append(d)
            # print ("%.6f, %.6f"%(tx, ty))
            N = len(self.lxs)-1
            if N> 0: 
                self.distance = round( math.sqrt((self.lxs[N]-self.lxs[N-1])**2 + (self.lys[N]-self.lys[N-1])**2 ) *1000, 2)
                self.sumLength += self.distance
                print (" Total Length =%.3fmm"%(self.sumLength))
                ch = plt.text((self.lxs[N]+self.lxs[N-1])/2.0, (self.lys[N]+self.lys[N-1])/2.0, str(self.distance), size=self.fontsize)
                self.chars.append(ch)

                ln, = plt.plot([self.lxs[N-1], self.lxs[N]],[self.lys[N-1], self.lys[N]], color='red', linewidth=dotsize/3)
                self.lines.append(ln)

                if self.clicked > 2: 

                    sx = 0; sy=0
                    for x, y in zip(self.lxs, self.lys): 
                        sx += x
                        sy += y
                    cx = sx/(float(N)+1)
                    cy = sy/(float(N)+1)

                    area = self.Area(ix=self.lxs, iy=self.lys)
                    for achar in self.achars: 
                        achar.set_visible(False)

                    ach= plt.text(cx, cy, "A="+ str(round(area*1_000_000, 1)), color='gray', size=self.fontsize)
                    self.achars.append(ach)

                    for char in self.llen:
                        char.set_visible(False)

                    for line in self.cline: 
                        line.remove()
                    self.cline=[]
                    
                    self.distance = round( math.sqrt((self.lxs[N]-self.lxs[0])**2 + (self.lys[N]-self.lys[0])**2 ) *1000, 2)
                    # ch = plt.text((self.lxs[N]+self.lxs[0])/2.0, (self.lys[N]+self.lys[0])/2.0, str(self.distance), color='gray', size=self.fontsize)
                    # self.llen.append(ch)
                    ln, = plt.plot([self.lxs[0], self.lxs[N]],[self.lys[0], self.lys[N]], color='gray', linestyle="--", linewidth=dotsize/3 )
                    self.cline.append(ln)
            # current_xlim=self.ax.get_xlim()
            # current_ylim=self.ax.get_ylim()
            # plt.xlim(current_xlim[0], current_xlim[1])
            # plt.ylim(current_ylim[0], current_ylim[1])
            # self.figure.canvas.draw_idle()
            self.plot_current_display_range()


    def ComparingMode(self, edge0, node0, xy0, edge1, node1, xy1): 
        self.figure.clear()
        self.ax = self.figure.add_subplot(111)
        self.ax.axis("equal")
        self.figure.canvas.mpl_connect('button_release_event', self.onReleased)
        self.figure.canvas.mpl_connect('scroll_event',self.zoom)

        try: 
            for nch in self.nodechars: 
                nch.set_visible(False)
            for nch in self.elchars: 
                nch.set_visible(False)
        except:
            pass 

        allx = []; ally=[]
        nd0 = np.array(node0.Node) 
        x = int(xy0/10)
        y = int(xy0%10)
        for eds in edge0: 
            for ed in eds.Edge: 
                ix = np.where(nd0[:,0] == ed[0])[0][0]; n1 = nd0[ix]
                ix = np.where(nd0[:,0] == ed[1])[0][0]; n2 = nd0[ix]
                polygon = plt.Polygon([[n1[x], n1[y]], [n2[x], n2[y]]], color='black', lw=0.5)
                self.ax.add_patch(polygon)
        
        nd1 = np.array(node1.Node) 
        x = int(xy1/10)
        y = int(xy1%10)
        for eds in edge1: 
            for ed in eds.Edge: 
                ix = np.where(nd1[:,0] == ed[0])[0][0]; n1 = nd1[ix]
                ix = np.where(nd1[:,0] == ed[1])[0][0]; n2 = nd1[ix]
                polygon = plt.Polygon([[n1[x], n1[y]], [n2[x], n2[y]]], color='red', lw=0.5)
                self.ax.add_patch(polygon)
                
        n0x = nd0[:,x]
        n0y = nd0[:,y]
        
        n1x = nd1[:,x]
        n1y = nd1[:,y]
        
        nx = np.concatenate((n0x, n1x))
        ny = np.concatenate((n0y, n1y))
        plt.xlim(np.min(nx)-0.01, np.max(nx)+0.01)
        plt.ylim(np.min(ny)-0.01, np.max(ny)+0.01)

        self.figure.tight_layout()
        self.figure.canvas.draw()
        
    def Draw_profiles(self, profiles, bdr): 
        
        for nch in self.nodechars: 
            nch.set_visible(False)
        for nch in self.elchars: 
            nch.set_visible(False)

        self.elchars=[]
        self.nodechars=[]

        self.figure.clear()
        self.ax = self.figure.add_subplot(111)
        self.ax.axis("equal")
        xs = []; ys=[]
        for profile in profiles: 
            if len(profile) : 
                for pg in profile: 
                    polygon = plt.Polygon(pg[0], color=pg[1], lw=pg[2])
                    xs.append(pg[0][0][0]); xs.append(pg[0][1][0])
                    ys.append(pg[0][0][1]); ys.append(pg[0][1][1])
                    self.ax.add_patch(polygon)
        
        xs = np.array(xs); ys=np.array(ys)
        xmn = np.min(xs); xmx = np.max(xs)
        ymn = np.min(ys); ymx = np.max(ys)
        plt.xlim(xmn-0.01, xmx+0.01)
        plt.ylim(ymn-0.01, ymx+0.01)
        self.figure.tight_layout()

        self.figure.canvas.mpl_connect('button_release_event', self.onReleased)
        self.figure.canvas.mpl_connect('scroll_event',self.zoom)
        self.figure.canvas.draw()

    def AddNodes(self, nodes, size, color=0, vmin=0, vmax=0 ): 
        for dot in self.dots: 
            dot.remove()
        self.dots = []
        if size > 0.0: 
            Xs=[]; Ys=[]
            npn = np.array(nodes.Node)
            Xs = npn[:,2]; Ys = npn[:,3]
            if color == 0: 
                nd = self.ax.scatter(Xs, Ys, s=size, c='gray' )    
            else: 
                try: 
                    Values = npn[:,4]
                    nd = self.ax.scatter(Xs, Ys, s=size, c=Values, cmap='jet', vmin=vmin, vmax=vmax )    
                except: 
                    nd = self.ax.scatter(Xs, Ys, s=size, c='gray' )    

            self.dots.append(nd)

        self.plot_current_display_range()
        # current_xlim=self.ax.get_xlim()
        # current_ylim=self.ax.get_ylim()
        # plt.xlim(current_xlim[0], current_xlim[1])
        # plt.ylim(current_ylim[0], current_ylim[1])
        # self.figure.canvas.draw_idle()

    def OnlyAddingNode(self, nodes, size, color=0, vmin=0, vmax=0 ): 
        # self.dots = []
        if size > 0.0: 
            Xs=[]; Ys=[]
            npn = np.array(nodes.Node)
            Xs = npn[:,2]; Ys = npn[:,3]
            if color == 0: 
                nd = self.ax.scatter(Xs, Ys, s=size, c='gray' )    
            else: 
                try: 
                    Values = npn[:,4]
                    nd = self.ax.scatter(Xs, Ys, s=size, c=Values, cmap='jet', vmin=vmin, vmax=vmax )    
                except: 
                    nd = self.ax.scatter(Xs, Ys, s=size, c='gray' )    

            self.dots.append(nd)

        self.plot_current_display_range()
        # current_xlim=self.ax.get_xlim()
        # current_ylim=self.ax.get_ylim()
        # plt.xlim(current_xlim[0], current_xlim[1])
        # plt.ylim(current_ylim[0], current_ylim[1])
        # self.figure.canvas.draw_idle()
            
    def Clear_ElsetBoundary(self): 
        try: 
            for line in self.LineEL: 
                line.remove()

            for txt in self.settexts: 
                txt.set_visible(False)
            self.plot_current_display_range()
            # current_xlim=self.ax.get_xlim()
            # current_ylim=self.ax.get_ylim()
            # plt.xlim(current_xlim[0], current_xlim[1])
            # plt.ylim(current_ylim[0], current_ylim[1])
            # self.figure.canvas.draw_idle()
        except: 
            pass 
    
    def Add_ElementsBoundary(self, solid, membrane, nodes, npel, nid=0, eid=0, surf=[]): 
        ## npel.append([el[0], el[1], el[2], el[3], el[4]])
        try: 
            for line in self.LineEL: 
                line.remove()
            for txt in self.settexts: 
                txt.set_visible(False)
        except: 
            pass
        self.LineEL=[]

        if nid ==1 or eid == 1: ids = 1 
        else: ids = 0 
        Edges = EDGE()
        elementcolor='blue'
        nodecolor = 'darkorange'
        if len(membrane) > 0: 
            for num in membrane: 
                ix = np.where(npel[:,0]== num[0])[0][0]
                el = npel[ix]
                Edges.Add([el[1], el[2], 'membrane', 0, el[0], -1])

                if ids == 1: 
                    ix = np.where(nodes[:,0] == el[1])[0][0]; n1 = nodes[ix]
                    ix = np.where(nodes[:,0] == el[2])[0][0]; n2 = nodes[ix]
                    if eid == 1: 
                        ch = self.ax.text((n1[2]+n2[2])/2, (n1[3]+n2[3])/2, str(int(el[0])), size=8, color=elementcolor )
                        self.settexts.append(ch)

                    if nid == 1: 
                        ch = self.ax.text(n1[2], n1[3], str(int(n1[0])), size=8,color=nodecolor )
                        self.settexts.append(ch)
                        ch = self.ax.text(n2[2], n2[3], str(int(n2[0])), size=8,color=nodecolor )
                        self.settexts.append(ch)


        if len(solid) > 0: 
            eset = ELEMENT()
            for num in solid: 
                ix = np.where(npel[:,0] == num[0])[0][0]
                if npel[ix][4] > 0:   el = [npel[ix][0], npel[ix][1], npel[ix][2], npel[ix][3], npel[ix][4], 'solid', 4]
                else: el = [npel[ix][0], npel[ix][1], npel[ix][2], npel[ix][3], npel[ix][4], 'solid', 3]
                eset.Add(el)

                if ids == 1: 
                    ix = np.where(nodes[:,0] == el[1])[0][0]; n1 = nodes[ix]
                    ix = np.where(nodes[:,0] == el[2])[0][0]; n2 = nodes[ix]
                    ix = np.where(nodes[:,0] == el[3])[0][0]; n3 = nodes[ix]
                    if el[6] == 4: ix = np.where(nodes[:,0] == el[4])[0][0]; n4 = nodes[ix]
                    if eid ==1: 
                        ch = self.ax.text((n1[2]+n2[2]+n3[2])/3, (n1[3]+n2[3]+n3[3])/3, str(int(el[0])), size=8, color=elementcolor )
                        self.settexts.append(ch)
                    if nid == 1: 
                        ch = self.ax.text(n1[2], n1[3], str(int(n1[0])), size=8,color=nodecolor )
                        self.settexts.append(ch)
                        ch = self.ax.text(n2[2], n2[3], str(int(n2[0])), size=8,color=nodecolor )
                        self.settexts.append(ch)
                        ch = self.ax.text(n3[2], n3[3], str(int(n3[0])), size=8,color=nodecolor )
                        self.settexts.append(ch)
                        if el[6] == 4: 
                            ch = self.ax.text(n4[2], n4[3], str(int(n4[0])), size=8,color=nodecolor )
                            self.settexts.append(ch)


            sEdges = Edge_ElementBoundary(eset.Element)


            for ed in sEdges.Edge: 
                Edges.Add(ed)

        if len(surf)>0: 
            for sf in surf: 
                try:
                    idx = np.where(npel[:,0] == sf[0])[0][0]
                except:
                    print ("No elements: %d"%(sf[0]))
                    continue 
                el = npel[idx]
                if el[4] ==0: 
                    if sf[1] == 'S1': 
                        n1d = el[1]; n2d = el[2]
                    if sf[1] == 'S2': 
                        n1d = el[2]; n2d = el[3]
                    if sf[1] == 'S3': 
                        n1d = el[3]; n2d = el[1]
                else: 
                    if sf[1] == 'S1': 
                        n1d = el[1]; n2d = el[2]
                    if sf[1] == 'S2': 
                        n1d = el[2]; n2d = el[3]
                    if sf[1] == 'S3': 
                        n1d = el[3]; n2d = el[4]
                    if sf[1] == 'S4': 
                        n1d = el[4]; n2d = el[1]

                Edges.Add([n1d, n2d, 'surface', 0, el[0], -1])

                if ids == 1: 
                    ix = np.where(nodes[:,0] == n1d)[0][0]; n1 = nodes[ix]
                    ix = np.where(nodes[:,0] == n2d)[0][0]; n2 = nodes[ix]
                    if eid == 1: 
                        ch = self.ax.text((n1[2]+n2[2])/2, (n1[3]+n2[3])/2, str(int(el[0])), size=8, color=elementcolor )
                        self.settexts.append(ch)

                    if nid == 1: 
                        ch = self.ax.text(n1[2], n1[3], str(int(n1[0])), size=8,color=nodecolor )
                        self.settexts.append(ch)
                        ch = self.ax.text(n2[2], n2[3], str(int(n2[0])), size=8,color=nodecolor )
                        self.settexts.append(ch)


        try: 
            for edge in Edges.Edge: 
                ix = np.where(nodes[:,0] == edge[0])[0][0]; n1 = nodes[ix]
                ix = np.where(nodes[:,0] == edge[1])[0][0]; n2 = nodes[ix]
                X=[n1[2], n2[2]]
                Y=[n1[3], n2[3]]
                if edge[2] == 'membrane':    line, = self.ax.plot(X, Y, color="black", linewidth=2.0)
                elif edge[2] == 'surface':   line, = self.ax.plot(X, Y, color="blue", linewidth=2.0)
                else:                        line, = self.ax.plot(X, Y, color="green", linewidth=2.0)
                self.LineEL.append(line) 

        except Exception : pass 

        self.plot_current_display_range()
        # current_xlim=self.ax.get_xlim()
        # current_ylim=self.ax.get_ylim()
        # plt.xlim(current_xlim[0], current_xlim[1])
        # plt.ylim(current_ylim[0], current_ylim[1])
        # self.figure.canvas.draw_idle()


    def Add_ElsetBoundary(self, name, elements, nodes, nid, eid,  allelsets, npel): 
        
        try: 
            for line in self.LineEL: 
                line.remove()
            for txt in self.settexts: 
                txt.set_visible(False)
        except: 
            pass 

        searching =[]
        for elset in allelsets: 
            if name == elset[0]: 
                searching = elset[1:]
                
                break 
        
        self.LineEL=[]
        solid=-1
        for el in elements.Element: 
            if el[0]==searching[0]: 
                if el[6] == 2: solid =0 
                else: solid = 1 
                break 

        nodes = np.array(nodes.Node)

        if solid == 1: 
            elset =ELEMENT()
            for num in searching: 
                ix = np.where(npel[:,0]== num)[0][0]
                el = elements.Element[ix]
                elset.Add(el)

                ix = np.where(nodes[:,0] == el[1])[0][0]; n1 = nodes[ix]
                ix = np.where(nodes[:,0] == el[2])[0][0]; n2 = nodes[ix]
                ix = np.where(nodes[:,0] == el[3])[0][0]; n3 = nodes[ix]
                if el[6] == 4: 
                    ix = np.where(nodes[:,0] == el[4])[0][0]; n4 = nodes[ix]

                if eid ==1: 
                    ch = self.ax.text((n1[2]+n2[2]+n3[2])/3, (n1[3]+n2[3]+n3[3])/3, str(int(el[0])), size=8, color='orange' )
                    self.settexts.append(ch)
                if nid == 1: 
                    ch = self.ax.text(n1[2], n1[3], str(int(n1[0])), size=8,color='gray' )
                    self.settexts.append(ch)
                    ch = self.ax.text(n2[2], n2[3], str(int(n2[0])), size=8,color='gray' )
                    self.settexts.append(ch)
                    ch = self.ax.text(n3[2], n3[3], str(int(n3[0])), size=8,color='gray' )
                    self.settexts.append(ch)
                    if el[6] == 4: 
                        ch = self.ax.text(n4[2], n4[3], str(int(n4[0])), size=8,color='gray' )
                        self.settexts.append(ch)


            print ("(%dEA)"%(len(elset.Element)))
            Edges = Edge_ElementBoundary(elset.Element)
        else: 
            Edges = EDGE()
            for num in searching: 
                ix = np.where(npel[:,0]== num)[0][0]
                el = elements.Element[ix]
                Edges.Add([el[1], el[2], el[5], 0, el[0], -1])
            
            print ("(%dEA)"%(len(Edges.Edge)))

        for edge in Edges.Edge: 
            ix = np.where(nodes[:,0] == edge[0])[0][0]; n1 = nodes[ix]
            ix = np.where(nodes[:,0] == edge[1])[0][0]; n2 = nodes[ix]
            X=[n1[2], n2[2]]
            Y=[n1[3], n2[3]]
            line, = self.ax.plot(X, Y, color="green", linewidth=2.0)
            if solid==0 and eid ==1: 
                ch = self.ax.text((n1[2]+n2[2])/2, (n1[3]+n2[3])/2, str(int(edge[4])), size=8, color='orange' )
                self.settexts.append(ch)
            if solid ==0 and nid == 1: 
                ch = self.ax.text(n1[2], n1[3], str(int(n1[0])), size=8,color='gray' )
                self.settexts.append(ch)
                ch = self.ax.text(n2[2], n2[3], str(int(n2[0])), size=8,color='gray' )
                self.settexts.append(ch)


            self.LineEL.append(line) 

        self.plot_current_display_range()
        # current_xlim=self.ax.get_xlim()
        # current_ylim=self.ax.get_ylim()
        # plt.xlim(current_xlim[0], current_xlim[1])
        # plt.ylim(current_ylim[0], current_ylim[1])
        # self.figure.canvas.draw_idle()


    def getplotinformation(self, node, element, elset, surface, tie, xy=23, add2d=[]):
        X = int(xy/10)
        Y = int(xy%10)

        self.nid=[];         self.x=[];         self.y=[]
        self.polygons=[];    self.eid=[];       self.ex = [];      self.ey=[];         self.ecolor=[]
        self.ric_r=[];       self.ric_l=[]; self.cont=[]

        Mcolor = 'red'

             
        for nd in node.Node:
            self.nid.append(nd[0])
            self.x.append(nd[X])
            self.y.append(nd[Y])

        for el in element.Element:
            try: 
                if el[6] == 3: 
                    x1 = el[10][0][0]; y1 = el[10][0][1]
                    x2 = el[10][1][0]; y2 = el[10][1][1]
                    x3 = el[10][2][0]; y3 = el[10][2][1]
                    icolor = Color(el[5])
                    # polygon = plt.Polygon([[x1, y1], [x2, y2], [x3, y3]], color=icolor, alpha=cdepth, lw=MeshLineWidth)
                    polygon = [[[x1, y1], [x2, y2], [x3, y3]], icolor]
                    self.polygons.append(polygon)
                    self.eid.append(el[0])
                    self.ex.append((x1 + x2 + x3 )/3)
                    self.ey.append((y1 + y2 + y3 )/3)
                    
                elif el[6] == 4: 
                    x1 = el[10][0][0]; y1 = el[10][0][1]
                    x2 = el[10][1][0]; y2 = el[10][1][1]
                    x3 = el[10][2][0]; y3 = el[10][2][1]
                    x4 = el[10][3][0]; y4 = el[10][3][1]
                    icolor = Color(el[5])
                    # polygon = plt.Polygon([[x1, y1], [x2, y2], [x3, y3], [x4, y4]], color=icolor, alpha=cdepth, lw=MeshLineWidth)
                    polygon = [[[x1, y1], [x2, y2], [x3, y3], [x4, y4]], icolor]
                    self.polygons.append(polygon)
                    self.eid.append(el[0])
                    self.ex.append((x1 + x2 + x3 + x4)/4)
                    self.ey.append((y1 + y2 + y3 + y4)/4)
                else:
                    x1 = el[10][0][0]; y1 = el[10][0][1]
                    x2 = el[10][1][0]; y2 = el[10][1][1]
                    polygon = [[[x1, y1], [x2, y2]], Mcolor]
                    self.polygons.append(polygon)
                    self.eid.append(el[0])
                    self.ex.append((x1 + x2 )/2)
                    self.ey.append((y1 + y2 )/2)
            except Exception as EX: 
                print (EX)
                continue 
                

        if len(add2d) > 0: 
            icolor = 'black'
            for en in add2d: 
                for el in element.Element: 
                    if en == el[0]: 
                        if el[6] == 3: 
                            x1 = el[10][0][0]; y1 = el[10][0][1]
                            x2 = el[10][1][0]; y2 = el[10][1][1]
                            x3 = el[10][2][0]; y3 = el[10][2][1]
                            polygon = [[[x1, y1], [x2, y2], [x3, y3]], icolor]
                            self.polygons.append(polygon)
                        elif el[6] == 4: 
                            x1 = el[10][0][0]; y1 = el[10][0][1]
                            x2 = el[10][1][0]; y2 = el[10][1][1]
                            x3 = el[10][2][0]; y3 = el[10][2][1]
                            x4 = el[10][3][0]; y4 = el[10][3][1]
                            polygon = [[[x1, y1], [x2, y2], [x3, y3], [x4, y4]], icolor]
                            self.polygons.append(polygon)
                        else: 
                            x1 = el[10][0][0]; y1 = el[10][0][1]
                            x2 = el[10][1][0]; y2 = el[10][1][1]
                            polygon = [[[x1, y1], [x2, y2]], icolor]
                            self.polygons.append(polygon)
                        break 

        ## surface =[[surfac name, [elname, face]], [], ... ]
        self.ric = []
        self.press =[]
        self.otheredge =[]
        N1 = 0; N2=0
        for sf in surface: 
            for i, en in enumerate(sf): 
                if i ==0: continue
                for el in element.Element:
                    if el[0] == en[0]:
                        if el[6] == 3: 
                            if en[1] == "S1": 
                                N1 = el[1]; N2=el[2]
                            if en[1] == 'S2': 
                                N1 = el[2]; N2=el[3]
                            if en[1] == 'S3': 
                                N1 = el[3]; N2=el[1]
                        elif el[6] == 4: 
                            if en[1] == "S1": 
                                N1 = el[1]; N2=el[2]
                            if en[1] == 'S2': 
                                N1 = el[2]; N2=el[3]
                            if en[1] == 'S3': 
                                N1 = el[3]; N2=el[4]
                            if en[1] == 'S4': 
                                N1 = el[4]; N2=el[1]
                        else:
                            pass
                        break
                # print (en, N1, N2)
                s1 = node.NodeByID(N1)
                s2 = node.NodeByID(N2)
                if s1[3] ==0 or s2 ==0: 
                    continue 
                x1 = s1[X];                y1 = s1[Y]
                x2 = s2[X];                y2 = s2[Y]

                if 'RIC' in sf[0] : self.ric.append([[x1, y1], [x2, y2]])
                elif 'PRESS' in sf[0] : self.press.append([[x1, y1], [x2, y2]])
                else:  self.otheredge.append([[x1, y1], [x2, y2]])

                if 'RIC_R' in sf[0]: self.ric_r.append([[x1, y1], [x2, y2]])
                if 'RIC_L' in sf[0]: self.ric_l.append([[x1, y1], [x2, y2]])
                if 'CONT' in sf[0]: self.cont.append([[x1, y1], [x2, y2]])
        self.stie = []
        self.mtie = []
        for t in tie:
            for sf in surface:
                if sf[0] == t[1][0]: ## slave tie
                    for i, en in enumerate(sf): 
                        if i ==0: continue
                        for el in element.Element:
                            if el[0] == en[0]:
                                if el[6] == 3: 
                                    if en[1] == "S1": 
                                        N1 = el[1]; N2=el[2]
                                    if en[1] == 'S2': 
                                        N1 = el[2]; N2=el[3]
                                    if en[1] == 'S3': 
                                        N1 = el[3]; N2=el[1]
                                elif el[6] == 4: 
                                    if en[1] == "S1": 
                                        N1 = el[1]; N2=el[2]
                                    if en[1] == 'S2': 
                                        N1 = el[2]; N2=el[3]
                                    if en[1] == 'S3': 
                                        N1 = el[3]; N2=el[4]
                                    if en[1] == 'S4': 
                                        N1 = el[4]; N2=el[1]
                                else:
                                    pass
                                break
                        s1 = node.NodeByID(N1)
                        s2 = node.NodeByID(N2)
                        if s1[3] ==0 or s2 ==0:  continue 

                        x1 = s1[X];                y1 = s1[Y]
                        x2 = s2[X];                y2 = s2[Y]
                        self.stie.append([[x1, y1], [x2, y2]])


                if sf[0] == t[1][1]:  ## master tie
                    for i, en in enumerate(sf): 
                        if i ==0: continue
                        for el in element.Element:
                            if el[0] == en[0]:
                                if el[6] == 3: 
                                    if en[1] == "S1": 
                                        N1 = el[1]; N2=el[2]
                                    if en[1] == 'S2': 
                                        N1 = el[2]; N2=el[3]
                                    if en[1] == 'S3': 
                                        N1 = el[3]; N2=el[1]
                                elif el[6] == 4: 
                                    if en[1] == "S1": 
                                        N1 = el[1]; N2=el[2]
                                    if en[1] == 'S2': 
                                        N1 = el[2]; N2=el[3]
                                    if en[1] == 'S3': 
                                        N1 = el[3]; N2=el[4]
                                    if en[1] == 'S4': 
                                        N1 = el[4]; N2=el[1]
                                else:
                                    pass
                                break
                        s1 = node.NodeByID(N1)
                        s2 = node.NodeByID(N2)
                        if s1[3] ==0 or s2 ==0:  continue 
                        x1 = s1[X];                y1 = s1[Y]
                        x2 = s2[X];                y2 = s2[Y]
                        self.mtie.append([[x1, y1], [x2, y2]])
    def SearchShow(self, sel=[], snode=[], node=[], element=[], size=0.1):
        X = 2; Y=3 
        nodecolor = 'gray'
        dotcolor = 'black'
        MeshLineWidth = 0.3
        MembWidth = 0.5
        Mcolor = 'red'
        # memtextcolor = 'purple'
        soltextcolor = 'steelblue'
        textsize = 8
        tcolors = ['orange', 'red', 'blue', 'green', 'aqua'] # Slave Tie colors 


        for ch in self.searchchars: 
            ch.set_visible(False)
            # ch.remove()
        
        for py in self.id_poly: 
            py.remove()
            
        try:
            self.searchdots.remove()
        except:
            pass

        self.plot_current_display_range() 
        # self.figure.canvas.draw_idle()

        self.id_poly = []
        self.searchdots= []
        self.searchchars=[]

        

        if len(sel) > 0: 
            spolygon = []
            for el in element.Element:
                for se in sel: 
                    if se == el[0]: 
                        if el[6] == 3: 
                            x1 = el[10][0][0]; y1 = el[10][0][1]
                            x2 = el[10][1][0]; y2 = el[10][1][1]
                            x3 = el[10][2][0]; y3 = el[10][2][1]
                            # icolor = Color(el[5])
                            spolygon.append([[[x1, y1], [x2, y2], [x3, y3]], el[0], el[8], el[9]])
                            
                        elif el[6] == 4: 
                            x1 = el[10][0][0]; y1 = el[10][0][1]
                            x2 = el[10][1][0]; y2 = el[10][1][1]
                            x3 = el[10][2][0]; y3 = el[10][2][1]
                            x4 = el[10][3][0]; y4 = el[10][3][1]
                            # icolor = Color(el[5])
                            spolygon.append([[[x1, y1], [x2, y2], [x3, y3], [x4, y4]], el[0], el[8], el[9]])
                        else:
                            x1 = el[10][0][0]; y1 = el[10][0][1]
                            x2 = el[10][1][0]; y2 = el[10][1][1]
                            spolygon.append([[[x1, y1], [x2, y2]], el[0], el[8], el[9]])

                        break 

        if len(snode)>0: 
            searchednode = []
            searchx =[]; searchy=[]
            for nd in node.Node:
                for sn in snode: 
                    if sn == nd[0]: 
                        searchednode.append(nd[0])
                        searchx.append(nd[X])
                        searchy.append(nd[Y])
                        break

        if len(sel) > 0 and len(spolygon)>0: 
            for i, splg in enumerate(spolygon): 
                if len(splg[0]) == 2: 
                    polygon = plt.Polygon(splg[0], color='black', alpha=1.0, lw=MembWidth)

                else:
                    polygon = plt.Polygon(splg[0], color='black', alpha=1.0, lw=MeshLineWidth)
                self.id_poly.append(polygon)
                self.ax.add_patch(polygon)
                ch = self.ax.text(splg[2], splg[3], splg[1], color='red', size=textsize)
                self.searchchars.append(ch)

        if len(snode) > 0 and len(searchednode) > 0: 
            if size == 0:         isize = 20.0
            else:                 isize = size * 10
            self.searchdots = self.ax.scatter(searchx, searchy, s=isize, c='red')
            for i, d in enumerate(searchednode): 
                ch=self.ax.text(searchx[i], searchy[i], d, color='green', size=textsize)
                self.searchchars.append(ch)

        self.plot_current_display_range()
        # self.figure.canvas.draw_idle()
    def NIDShow(self, show=0): 
        if show ==1: 
            nodecolor = 'gray'
            textsize = 8
            for i, d in enumerate(self.nid): 
                nch = self.ax.text(self.x[i], self.y[i], d, color=nodecolor, size=textsize)
                self.nodechars.append(nch)
        else: 
            for nch in self.nodechars: 
                nch.set_visible(False)
        # current_xlim=self.ax.get_xlim()
        # current_ylim=self.ax.get_ylim()
        # plt.xlim(current_xlim[0], current_xlim[1])
        # plt.ylim(current_ylim[0], current_ylim[1])
        # self.figure.canvas.draw_idle()
        self.plot_current_display_range()
    def EIDShow(self, show=0): 
        if show ==1: 
            soltextcolor = 'steelblue'
            textsize = 8
            for i, d in enumerate(self.eid): 
                nch = self.ax.text(self.ex[i], self.ey[i], d, color=soltextcolor, size=textsize )
                self.elchars.append(nch)
        else: 
            for nch in self.elchars: 
                nch.set_visible(False)
        # current_xlim=self.ax.get_xlim()
        # current_ylim=self.ax.get_ylim()
        # plt.xlim(current_xlim[0], current_xlim[1])
        # plt.ylim(current_ylim[0], current_ylim[1])
        # self.figure.canvas.draw_idle()
        self.plot_current_display_range()
    def Displaysets(self, dsurface, delset, element, cdepth=0.5, mdepth=0.8,  xy=23, sshow=1, eshow=1):
        X = int(xy/10)
        Y = int(xy%10) 
        MeshLineWidth = 0.3
        MembWidth = 0.5
        Mcolor = 'red'
        soltextcolor = 'steelblue'

        for nch in self.nodechars: 
            nch.set_visible(False)
        for nch in self.elchars: 
            nch.set_visible(False)
        self.elchars=[]
        self.nodechars=[]

 

        xmin = 1000.0; xmax = -1000.0
        ymin = 1000.0; ymax = -1000.0
        self.figure.clear()
        axe = self.figure.add_subplot(111)
        axe.axis("equal")
        if sshow ==1: 
            for i, sf in enumerate(dsurface): 
                if i == 0: continue
                for el in element.Element:
                    if sf[0] == el[0]: 
                        if el[6] ==3: 
                            if sf[1] =="S1": shade = [[el[10][0][0], el[10][0][1]], [el[10][1][0], el[10][1][1]]]
                            if sf[1] =="S2": shade = [[el[10][1][0], el[10][1][1]], [el[10][2][0], el[10][2][1]]]
                            if sf[1] =="S3": shade = [[el[10][2][0], el[10][2][1]], [el[10][0][0], el[10][0][1]]]
                        elif el[6] ==4: 
                            if sf[1] =="S1": shade = [[el[10][0][0], el[10][0][1]], [el[10][1][0], el[10][1][1]]]
                            if sf[1] =="S2": shade = [[el[10][1][0], el[10][1][1]], [el[10][2][0], el[10][2][1]]]
                            if sf[1] =="S3": shade = [[el[10][2][0], el[10][2][1]], [el[10][3][0], el[10][3][1]]]
                            if sf[1] =="S4": shade = [[el[10][3][0], el[10][3][1]], [el[10][0][0], el[10][0][1]]]

                        polygon = plt.Polygon(shade, color="black", alpha=mdepth, lw=MembWidth*3)
                        axe.add_patch(polygon)

                        for k in range(el[6]): 
                            if xmax < el[10][k][0]: xmax = el[10][k][0]
                            if xmin > el[10][k][0]: xmin = el[10][k][0]
                            if ymax < el[10][k][1]: ymax = el[10][k][1]
                            if ymin > el[10][k][1]: ymin = el[10][k][1]

                        break

        if eshow ==1: 
            for i, eset in enumerate(delset): 
                if i == 0 : continue 
                for el in element.Element:
                    if eset == el[0]: 
                        if el[6] == 2: 
                            shade = [[el[10][0][0], el[10][0][1]], [el[10][1][0], el[10][1][1]]]
                            icolor = Color(el[5])
                            polygon = plt.Polygon(shade, color=icolor, alpha=mdepth, lw=MembWidth*3 )
                            axe.add_patch(polygon)
                        elif el[6] == 3: 
                            shade = [[el[10][0][0], el[10][0][1]], [el[10][1][0], el[10][1][1]], [el[10][2][0], el[10][2][1]]]
                            icolor = Color(el[5])
                            polygon = plt.Polygon(shade, color=icolor, alpha=cdepth, lw=MeshLineWidth )
                            axe.add_patch(polygon)
                        elif el[6] == 4: 
                            shade = [[el[10][0][0], el[10][0][1]], [el[10][1][0], el[10][1][1]], [el[10][2][0], el[10][2][1]], [el[10][3][0], el[10][3][1]]]
                            icolor = Color(el[5])
                            polygon = plt.Polygon(shade, color=icolor, alpha=cdepth, lw=MeshLineWidth )
                            axe.add_patch(polygon)
                        for k in range(el[6]): 
                            if xmax < el[10][k][0]: xmax = el[10][k][0]
                            if xmin > el[10][k][0]: xmin = el[10][k][0]
                            if ymax < el[10][k][1]: ymax = el[10][k][1]
                            if ymin > el[10][k][1]: ymin = el[10][k][1]
                        break
        if eshow ==0 and sshow ==1: 
            for poly in self.polygons:
                if len(poly[0]) > 2:
                    polygon = plt.Polygon(poly[0], color=poly[1], alpha=cdepth, lw=MeshLineWidth)
                    axe.add_patch(polygon)
                    for k in range(len(poly[0])): 
                        if xmax < poly[0][k][0]: xmax = poly[0][k][0]
                        if xmin > poly[0][k][0]: xmin = poly[0][k][0]
                        if ymax < poly[0][k][1]: ymax = poly[0][k][1]
                        if ymin > poly[0][k][1]: ymin = poly[0][k][1]
            for poly in self.polygons:
                if len(poly[0]) == 2:
                    polygon = plt.Polygon(poly[0], color=poly[1], alpha=mdepth, lw=MembWidth )
                    axe.add_patch(polygon)
                    for k in range(len(poly[0])): 
                        if xmax < poly[0][k][0]: xmax = poly[0][k][0]
                        if xmin > poly[0][k][0]: xmin = poly[0][k][0]
                        if ymax < poly[0][k][1]: ymax = poly[0][k][1]
                        if ymin > poly[0][k][1]: ymin = poly[0][k][1]

        plt.xlim(xmin-0.01, xmax+0.01)
        plt.ylim(ymin-0.01, ymax+0.01)
        self.figure.canvas.draw()
    def gplot(self,  node, element, elset=None, surface=None, tie=None, xy=23, size=0.1, ni=0, ei=0, si=0, ti=0, srr=0, srl=0, spr=0, srd=0, cdepth=0.5,  mdepth=0.8, snode=[], sel=[], rims=[]):
        X = int(xy/10)
        Y = int(xy%10)

        nodecolor = 'gray'
        dotcolor = 'black'
        MeshLineWidth = 0.3
        MembWidth = 0.5
        Mcolor = 'red'
        # memtextcolor = 'purple'
        soltextcolor = 'steelblue'
        textsize = 8
        tcolors = ['orange', 'red', 'blue', 'green', 'aqua'] # Slave Tie colors 

        for nch in self.nodechars: 
            nch.set_visible(False)
        for nch in self.elchars: 
            nch.set_visible(False)

        self.elchars=[]
        self.nodechars=[]

        self.temperature_mode = False 

        # print (" ** ", len(node.Node),  len(element.Element))
        if len(node.Node) == 0 or len(element.Element) == 0: 
            self.figure.clear()
            self.ax = self.figure.add_subplot(111)
            self.ax.axis("equal")
            self.figure.canvas.mpl_connect('button_release_event', self.onReleased)
            self.figure.canvas.draw()
        elif len(node.Node) and len(element.Element) == 0: 
            self.figure.clear()
            self.ax = self.figure.add_subplot(111)
            self.ax.axis("equal")
            npn = np.array(node.Node)

            self.ax.scatter(npn[:,2], npn[:,3]) 

            self.points = np.array(node.Node)

            self.figure.canvas.mpl_connect('button_release_event', self.onReleased)
            self.figure.canvas.mpl_connect('scroll_event',self.zoom)
            self.figure.canvas.mpl_connect('button_press_event', self.onclick)
            self.figure.canvas.draw()


        else:
            if self.pei ==0 and ei ==1: 
                for i, d in enumerate(self.eid):
                    self.ax.text(self.ex[i], self.ey[i], d, color=soltextcolor, size=textsize )

                self.pei = ei
                if ni ==1: 
                    for i, d in enumerate(self.nid): 
                        self.ax.text(self.x[i], self.y[i], d, color=nodecolor, size=textsize)
                self.figure.canvas.mpl_connect('button_release_event', self.onReleased)
                # self.figure.canvas.draw()
                # self.figure.canvas.draw_idle()
                self.plot_current_display_range()
            elif self.pni ==0 and ni ==1: 
                for i, d in enumerate(self.nid): 
                    self.ax.text(self.x[i], self.y[i], d, color=nodecolor, size=textsize)
                self.pni = ni 
                if ei ==1: 
                    for i, d in enumerate(self.eid):
                        self.ax.text(self.ex[i], self.ey[i], d, color=soltextcolor, size=textsize )
                self.figure.canvas.mpl_connect('button_release_event', self.onReleased)
                # self.figure.canvas.draw_idle()
                self.plot_current_display_range()
            else:
                self.figure.clear()
                self.ax = self.figure.add_subplot(111)
                self.ax.axis("equal")

                for poly in self.polygons:
                    if len(poly[0]) > 2:
                        polygon = plt.Polygon(poly[0], color=poly[1], alpha=cdepth, lw=MeshLineWidth)
                        self.ax.add_patch(polygon)
                for poly in self.polygons:
                    if len(poly[0]) == 2:
                        polygon = plt.Polygon(poly[0], color=poly[1], alpha=mdepth, lw=MembWidth )
                        self.ax.add_patch(polygon)
                self.ax.scatter(self.x, self.y, s=size, c=dotcolor)

                if si ==1:
                    sdepth = mdepth
                    if sdepth == 0: sdepth = 0.1
                    swidth = MembWidth
                    if swidth ==0.0 : swidth = 0.1
                    presscolor = 'blue'

                    for psf in self.press:
                        polygon = plt.Polygon(psf, color=presscolor, alpha=sdepth, lw=swidth*3.0)
                        self.ax.add_patch(polygon)
                    riccolor = 'green'
                    for psf in self.ric:
                        polygon = plt.Polygon(psf, color=riccolor, alpha=sdepth, lw=swidth*1.5)
                        self.ax.add_patch(polygon)
                    othercolor = 'black'
                    for psf in self.otheredge:
                        polygon = plt.Polygon(psf, color=othercolor, alpha=sdepth, lw=swidth)
                        self.ax.add_patch(polygon)

                if ti ==1: 
                    
                    tdepth = 0.8
                    tdepth = mdepth
                    if tdepth == 0: tdepth = 0.1
                    swidth = MembWidth
                    if swidth ==0.0 : swidth = 0.1
                    mtiecolor = 'white'
                    for psf in self.mtie:
                        polygon = plt.Polygon(psf, color=mtiecolor, alpha=tdepth, lw=swidth*10.0)
                        self.ax.add_patch(polygon)
                        
                    
                    for i, psf in enumerate(self.stie):
                        stiecolor = tcolors[i%4]
                        polygon = plt.Polygon(psf, color=stiecolor, alpha=tdepth, lw=swidth*2)
                        self.ax.add_patch(polygon)

                if srr ==1:
                    sdepth = mdepth
                    if sdepth == 0: sdepth = 0.1
                    swidth = MembWidth
                    if swidth ==0.0 : swidth = 0.1
                    presscolor = 'green'

                    for psf in self.ric_r:
                        polygon = plt.Polygon(psf, color=presscolor, alpha=sdepth, lw=swidth*3.0)
                        self.ax.add_patch(polygon)

                if srl ==1:
                    sdepth = mdepth
                    if sdepth == 0: sdepth = 0.1
                    swidth = MembWidth
                    if swidth ==0.0 : swidth = 0.1
                    presscolor = 'green'

                    for psf in self.ric_l:
                        polygon = plt.Polygon(psf, color=presscolor, alpha=sdepth, lw=swidth*3.0)
                        self.ax.add_patch(polygon)

                if spr ==1:
                    sdepth = mdepth
                    if sdepth == 0: sdepth = 0.1
                    swidth = MembWidth
                    if swidth ==0.0 : swidth = 0.1
                    presscolor = 'blue'

                    for psf in self.press:
                        polygon = plt.Polygon(psf, color=presscolor, alpha=sdepth, lw=swidth*3.0)
                        self.ax.add_patch(polygon)
                if srd ==1:
                    sdepth = mdepth
                    if sdepth == 0: sdepth = 0.1
                    swidth = MembWidth
                    if swidth ==0.0 : swidth = 0.1
                    presscolor = 'black'

                    for psf in self.cont:
                        polygon = plt.Polygon(psf, color=presscolor, alpha=sdepth, lw=swidth*3.0)
                        self.ax.add_patch(polygon)

                if self.pei ==1 and ei ==1: 
                    for i, d in enumerate(self.eid):
                        self.ax.text(self.ex[i], self.ey[i], d, color=soltextcolor, size=textsize )
                if self.pni ==1 and ni ==1: 
                    for i, d in enumerate(self.nid): 
                        self.ax.text(self.x[i], self.y[i], d, color=nodecolor, size=textsize)

                if len(rims) > 0: 
                    for rim in rims: 
                        for word in rim: 
                            if 'START' in word[0].upper(): 
                                sx = float(word[2])
                                sy = float(word[1])

                            if 'CIRCL' in word[0].upper(): 
                                ex = float(word[2])
                                ey = float(word[1])
                                cx = float(word[4])
                                cy = float(word[3])
                                D = math.sqrt((cx-ex)**2 + (cy-ey)**2)*2

                                n1 = [0, 0, sx, sy]
                                n2 = [0, 0, cx, cy]
                                n3 = [0, 0, ex, ey]
                                angle = Angle_3nodes(n1, n2, n3, xy=23)
                                
                                cn = [0, 0, cx, cy+ 0.1]
                                if n1[2] > n3[2] : 
                                    startangle = Angle_3nodes(n1, n2, cn, xy=23)  * 180.0 / 3.141592  
                                    sp = sx 
                                else: 
                                    startangle = Angle_3nodes(n3, n2, cn, xy=23)  * 180.0 / 3.141592  
                                    sp = ex 
                                if cx < sp: 
                                    if startangle >=90: 
                                        startangle = 450.0 - startangle 
                                    else: 
                                        startangle = 90.0 - startangle 
                                else: 
                                    startangle += 90.0

                                if sp > 0: 
                                    if ey < cy: 
                                        startangle -= angle*180.0/3.141592
                                if sp < 0: 
                                    if sy < cy: 
                                        startangle -= angle*180.0/3.141592
                                
                                
                                arc = Arc((cx, cy), D, D, theta1=0, theta2=angle*180.0/3.141592, angle=startangle, color='black',  linewidth=0.2)
                                self.ax.add_patch(arc)
                                sx = ex
                                sy = ey 

                                # print ("CIRCL", cx,cy, ex,ey)

                            if 'LINE' in word[0].upper(): 
                                ex = float(word[2])
                                ey = float(word[1])

                                self.ax.plot([sx,ex], [sy, ey], color='black', linewidth=0.2 ) 
                                # print ("Line", cx,cy, ex,ey)

                                sx = ex
                                sy = ey 


                self.pei = ei
                self.pni = ni
                # allaxes = self.figure.axes
                # print(allaxes)
                self.figure.tight_layout()

                self.points = np.array(node.Node)

                self.figure.canvas.mpl_connect('button_release_event', self.onReleased)
                self.figure.canvas.mpl_connect('scroll_event',self.zoom)
                self.figure.canvas.mpl_connect('button_press_event', self.onclick)
                self.figure.canvas.draw()

    def plot_temperature(self, npn, els, vmin=25, vmax=100, levels=10, contour=None, size=0, num=10):
        self.clearplot()

        self.ax = self.figure.add_subplot(111)
        self.ax.axis("equal")
        
        if isinstance(contour, type(None)):
            self.temperature_node = npn  
            cnt = 0 
            for el in els: 
                if el[4]>0: 
                    ix = np.where(npn[:,0]==el[1])[0][0]; n1=npn[ix]
                    ix = np.where(npn[:,0]==el[2])[0][0]; n2=npn[ix]
                    ix = np.where(npn[:,0]==el[3])[0][0]; n3=npn[ix]
                    ix = np.where(npn[:,0]==el[4])[0][0]; n4=npn[ix]
                elif el[3] > 0: 
                    ix = np.where(npn[:,0]==el[1])[0][0]; n1=npn[ix]
                    ix = np.where(npn[:,0]==el[2])[0][0]; n2=npn[ix]
                    ix = np.where(npn[:,0]==el[3])[0][0]; n3=npn[ix]
                    n4 = n3 
                if el[3] != 0: 
                    len1=math.sqrt((n1[2]-n2[2])**2 + (n1[3]-n2[3])**2)
                    len2=math.sqrt((n2[2]-n3[2])**2 + (n2[3]-n3[3])**2)
                    len3=math.sqrt((n3[2]-n4[2])**2 + (n3[3]-n4[3])**2)
                    len4=math.sqrt((n4[2]-n1[2])**2 + (n4[3]-n1[3])**2)
                    if len1 > 0.1 or len2 > 0.1 or len3 > 0.1 or len4 > 0.1:   
                        print (" %d, %d, %d, %d"%(n1[0], n2[0], n3[0], n4[0]))
                        print (" %.2f, %.2f, %.2f, %.2f"%(len1*1000, len2*1000, len3*1000, len4*1000))

                    xs =[n1[2], n2[2], n3[2], n4[2]]
                    ys =[n1[3], n2[3], n3[3], n4[3]]
                    vs = [n1[4], n2[4], n3[4], n4[4]]
                    if cnt ==0: 
                        tpx, tpy, tpv =  meshgrid_in_Quadrilateral (xs, ys, vs=vs, num=num)
                        cnt = 1
                    else: 
                        px, py, pv =  meshgrid_in_Quadrilateral (xs, ys, vs=vs, num=num)
                        xrange = [min(xs), max(xs)]
                        yrange = [min(ys), max(ys)]
                        xmin=np.min(px); xmax=np.max(px)
                        ymin=np.min(py); ymax=np.max(py)
                        draw = True  
                        if xrange[0] < xmin or xrange[1] > xmax: 
                            print (" x check %d"%(el[0]))
                            draw = False
                        if yrange[0] < ymin or yrange[1] > ymax: 
                            print (" y check %d"%(el[0]))
                            draw = False
                        if draw: 
                            tpx = np.append(tpx, px,axis=0)
                            tpy = np.append(tpy, py,axis=0)
                            tpv = np.append(tpv, pv,axis=0)

        else: 
            tpx = contour[0]
            tpy = contour[1]
            tpv = contour[2]

        for el in els:    
            if el[3] ==0: 
                ix = np.where(npn[:,0]==el[1])[0][0]; n1=npn[ix]
                ix = np.where(npn[:,0]==el[2])[0][0]; n2=npn[ix]
                X=[n1[2], n2[2]]
                Y=[n1[3], n2[3]]
                line, = self.ax.plot(X, Y, color="lightgray", linewidth=0.5)
        
        # cp = self.ax.contourf(tpx, tpy, tpv, cmap='rainbow', levels=levels, vmin=vmin, vmax=vmax)
        if size==0: size=1.0
        self.ax.scatter(tpx, tpy, c=tpv, cmap='rainbow', vmin=vmin, vmax=vmax, s=size)
        # print ("* Plotting temeprature Contour")
        self.figure.tight_layout()
        self.figure.canvas.draw()
        self.temperature_mode = True 
        

        return [tpx, tpy, tpv]

    def plot_ElementValue(self, npn, els, vmin=25, vmax=100, num=5, levels=10, contour=None, size=0, col=5, type='flat', cmap='rainbow'):
        self.clearplot()

        self.ax = self.figure.add_subplot(111)
        self.ax.axis("equal")
        # if isinstance(contour, type(None)):
        if contour == None : 
             
            cnt = 0 
            for el in els: 
                if el[4]>0: 
                    ix = np.where(npn[:,0]==el[1])[0][0]; n1=npn[ix]
                    ix = np.where(npn[:,0]==el[2])[0][0]; n2=npn[ix]
                    ix = np.where(npn[:,0]==el[3])[0][0]; n3=npn[ix]
                    ix = np.where(npn[:,0]==el[4])[0][0]; n4=npn[ix]
                elif el[3] > 0: 
                    ix = np.where(npn[:,0]==el[1])[0][0]; n1=npn[ix]
                    ix = np.where(npn[:,0]==el[2])[0][0]; n2=npn[ix]
                    ix = np.where(npn[:,0]==el[3])[0][0]; n3=npn[ix]
                    n4 = n3 
                if el[3] != 0: 
                    len1=math.sqrt((n1[2]-n2[2])**2 + (n1[3]-n2[3])**2)
                    len2=math.sqrt((n2[2]-n3[2])**2 + (n2[3]-n3[3])**2)
                    len3=math.sqrt((n3[2]-n4[2])**2 + (n3[3]-n4[3])**2)
                    len4=math.sqrt((n4[2]-n1[2])**2 + (n4[3]-n1[3])**2)
                    if len1 > 0.1 or len2 > 0.1 or len3 > 0.1 or len4 > 0.1:   
                        print (" %d, %d, %d, %d"%(n1[0], n2[0], n3[0], n4[0]))
                        print (" %.2f, %.2f, %.2f, %.2f"%(len1*1000, len2*1000, len3*1000, len4*1000))

                    xs =[n1[2], n2[2], n3[2], n4[2]]
                    ys =[n1[3], n2[3], n3[3], n4[3]]

                    if type=='flat': 
                        vs = [el[col], el[col], el[col], el[col]]
                    if cnt ==0: 
                        tpx, tpy, tpv =  meshgrid_in_Quadrilateral (xs, ys, vs=vs, num=num)
                        cnt = 1
                    else: 
                        px, py, pv =  meshgrid_in_Quadrilateral (xs, ys, vs=vs, num=num)
                        tpx = np.append(tpx, px,axis=0)
                        tpy = np.append(tpy, py,axis=0)
                        tpv = np.append(tpv, pv,axis=0)

        else: 
            tpx = contour[0]
            tpy = contour[1]
            tpv = contour[2]

        for el in els:    
            if el[3] ==0: 
                ix = np.where(npn[:,0]==el[1])[0][0]; n1=npn[ix]
                ix = np.where(npn[:,0]==el[2])[0][0]; n2=npn[ix]
                X=[n1[2], n2[2]]
                Y=[n1[3], n2[3]]
                line, = self.ax.plot(X, Y, color="black", linewidth=1.)
        
        # cp = self.ax.contourf(tpx, tpy, tpv, cmap='rainbow', levels=levels, vmin=vmin, vmax=vmax)
        if size==0: size=1.0
        self.ax.scatter(tpx, tpy, c=tpv, cmap=cmap, vmin=vmin, vmax=vmax, s=size)
        self.figure.tight_layout()
        self.figure.canvas.draw()

        return [tpx, tpy, tpv]

    def clearplot(self):
        self.figure.clear()
        self.temperature_mode = False 
        self.nid=[]
        self.x=[]
        self.y=[]
        self.polygons=[]
        self.lines=[]

        self.eid=[]
        self.ey=[]
        self.ex=[]
        self.ecolor=[]

        for nch in self.nodechars: 
            nch.set_visible(False)
        for nch in self.elchars: 
            nch.set_visible(False)

        self.elchars =[]
        self.nodechars = []

        for ch in self.searchchars: 
            ch.set_visible(False)
            # ch.remove()
        
        for py in self.id_poly: 
            py.remove()
            
        try:
            self.searchdots.remove()
        except:
            pass 
        # self.figure.canvas.draw_idle()
        self.plot_current_display_range()

        self.id_poly = []
        self.searchdots= []
        self.searchchars=[]

    def plot_current_display_range(self): 
        current_xlim=self.ax.get_xlim()
        current_ylim=self.ax.get_ylim()
        plt.xlim(current_xlim[0], current_xlim[1])
        plt.ylim(current_ylim[0], current_ylim[1])
        self.figure.canvas.draw_idle()
         
def Angle_3nodes(n1=[], n2=[], n3=[], xy=0): ## n2 : mid node 
    v1 = [n1[1]-n2[1], n1[2]-n2[2], n1[3]-n2[3] ]
    v2 = [n3[1]-n2[1], n3[2]-n2[2], n3[3]-n2[3] ]
    
    if xy ==0: 
        cos = round((v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2]) /  math.sqrt(v1[0]**2 + v1[1]**2 + v1[2]**2) / math.sqrt(v2[0]**2 + v2[1]**2 + v2[2]**2), 9)
        angle = math.acos(cos)
    else: 
        x = int(xy/10)-1;     y = int(xy%10)-1
        cos = round((v1[x]*v2[x] + v1[y]*v2[y] ) /  math.sqrt(v1[x]**2 + v1[y]**2 ) / math.sqrt(v2[x]**2 + v2[y]**2), 9)
        angle = math.acos(cos)

    return angle
    


if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    MainWindow = QtWidgets.QMainWindow()
    ui = Ui_MainWindow()
    ui.setupUi(MainWindow)
    MainWindow.show()
    sys.exit(app.exec_())



 ########################################
 ## Basic Functions 
 ########################################


