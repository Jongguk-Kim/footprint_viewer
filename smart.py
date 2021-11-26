import struct, math 

import numpy as np 
from os.path import isfile 
from layout import NODE, SURFACE 
from basicfunctions import Angle_Between_Vectors

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

def readSDBResults(sdbmodel, sdbresult, btm=True, sectors=240, offset=10**4, rotating=0): 
    npn, np2d, np3d, rim_road = readSDB(sdbmodel) 
    d_npn, deformed_Rim_road, iELD, iSED =SDBResult_READ(sdbresult, npn)
    
    idx = np.where(d_npn[:,0]<10**7)[0]
    bdn = d_npn[idx]

    if btm: 
        zmax = np.max(bdn[:,3])
        ix = np.where(bdn[:,3]==zmax)[0][0]
        tn = bdn[ix]
        tsector = int(tn[0]/offset)
        bsector = tsector + int(sectors/2)
        if bsector >= sectors: bsector -= sectors 
        ix1 = np.where(bdn[:,0]>offset*bsector)[0]
        ix2 = np.where(bdn[:,0]<offset*bsector+offset)[0]
        idx = np.intersect1d(ix1, ix2)
        btm_sector_node = bdn[idx]
        btm = []
        if not rotating: 
            for bm in btm_sector_node: 
                nid = int(bm[0])%offset 
                btm.append([nid, bm[1], bm[2], bm[3]])
        else: 
            angle = math.radians(rotating)
            for bm in btm_sector_node:
                nid = int(bm[0])%offset 

                r = math.sqrt(bm[1]**2 + bm[3]**2)

                x = math.cos(angle)*bm[2] - math.sin(angle)*r
                y = math.sin(angle)*bm[2] + math.cos(angle)*r

                btm.append([nid, 0.0, x, y])

        fp=open('deformed.tmp', 'w')
        for n in btm: 
            fp.write("%d, %.7f, %7f, %7f\n"%(n[0], n[1], n[2], n[3]))
        fp.close()

        return np.array(btm)

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

def ResultSfric(model="", result="", sfric='', deformed=1, ht=8.0E-03, sfht=1.0E-03):

    """
    Extracting the results from HK-SMART sfric file 
    sfric should be defined as a class of SFRIC -> sfric = SFRIC()
    sfric's member : self.Node, self.Surface, self.Rim, self.Road
    if deformed ==0: this reads only model file. 
    Example :
            sfric = SFRIC()
            ResultSfric(model=osfric, result=rsfric, sfric=sfric, deformed=1)
    """
    sfricfile = model 
    sfricresultfile = result 
    if sfricfile !="" or sfricresultfile !="":
        if  sfricfile == "" : sfricfile = sfricresultfile[:-3]
        ## definition variables 
        EN_SDB_TITLE=11         
        EN_SDB_LOCAL_COORD = 12        
        EN_SDB_NODE = 21
        EN_SDB_ELEM_3D = 31      
        EN_SDB_ELEM_2D = 32
        EN_SDB_NSET    = 41
        EN_SDB_ELSET3D = 51       
        EN_SDB_ELSET2D = 52
        EN_SDB_RIM     = 61    
        EN_SDB_ROAD = 62
        EN_SDB_IMPACTOR = 63
        EN_SDB_END     = 999

        EN_RESULT_TITLE = 13
        EN_RESULT=14
        EN_RESULT_ERR=1101
        EN_RESULT_DIS=101
        EN_RESULT_DIS_X=102 
        EN_RESULT_DIS_Y=103 
        EN_RESULT_DIS_Z=104
        EN_RESULT_ELEMENT = 2   
        EN_RESULT_RIGID=1
        RESULT_Temperature=3101

        RESULT_NodalArea = 900000001
        RESULT_NodalAreaPJTRatio=900000002
        RESULT_FricEnergyRateSum = 900000011
        RESULT_FricEnergyRateX=900000012
        RESULT_FricEnergyRateY=900000013
        RESULT_FricEnergyRateZ=900000014
        RESULT_SlipVelocityMag=900000021
        RESULT_SlipVelocityX=900000022
        RESULT_SlipVelocityY=900000023
        RESULT_SlipVelocityZ=900000024
        RESULT_SlipFlag=900000030
        RESULT_FricCoefficient = 900000031
        RESULT_ContactShearForceMag = 900000041
        RESULT_ContactShearForceX = 900000042
        RESULT_ContactShearForceY=900000043
        RESULT_ContactForceNormal=900000044
        RESULT_ContactShearStressMag = 900000051
        RESULT_ContactShearStressX = 900000052
        RESULT_ContactShearStressY = 900000053
        RESULT_ContactStressNormal = 900000054
        RESULT_FricEnergyAccumSum=900000061
        RESULT_FricEnergyAccumX = 900000062
        RESULT_FricEnergyAccumY = 900000063
        RESULT_FricEnergyAccumZ = 900000064
        RESULT_WEAR = 900328631
    else: 
        return 0

    if isfile(sfricfile) == False: 
        print ("No sfric model file")
        return 0
    file = open(sfricfile, 'rb')
    file.seek(0,2); fend = file.tell() #find the end position of binary file
    file.seek(0);   fpos = file.tell() #find the start position of binary file

    ### SFRIC FILE
    READ_LENGTH = 0
    # CAMBER=[0.0]
    NLB = []
    

    while file.tell() < fend:
        BlockID     = struct.unpack('i', file.read(4))[0]

        if BlockID == EN_SDB_TITLE:
            BlockLength = struct.unpack('i', file.read(4))[0]
            BlockTitle = ''
            for i in range(BlockLength): 
                BlockTitle += str(struct.unpack('c', file.read(1))[0])

        elif BlockID == EN_SDB_LOCAL_COORD:
            LocalCoord = []
            for i in range(9):
                LocalCoord.append(struct.unpack('d', file.read(8))[0])
                # CAMBER[0] = math.atan(LocalCoord[5]/LocalCoord[4])*180/math.pi

        elif BlockID == EN_SDB_NODE:
            NodesNUM = struct.unpack('i', file.read(4))[0]
            for i in range(NodesNUM):
                Label = struct.unpack('i', file.read(4))[0]
                NLB.append(Label)
            for i in range(NodesNUM):
                X = struct.unpack('d', file.read(8))[0]
                Y = struct.unpack('d', file.read(8))[0]
                Z = struct.unpack('d', file.read(8))[0]
                sfric.Node.Add([NLB[i], X, Y, Z, 0.0, 0.0])
            # print ("NODE NUM", NodesNUM)

        elif BlockID == EN_SDB_ELEM_3D:
            E3DNUM = struct.unpack('i', file.read(4))[0]
            E3D = []
            for i in range(E3DNUM): 
                ID = struct.unpack('i', file.read(4))[0]
                N1 = struct.unpack('i', file.read(4))[0]
                N2 = struct.unpack('i', file.read(4))[0]
                N3 = struct.unpack('i', file.read(4))[0]
                N4 = struct.unpack('i', file.read(4))[0]
                N5 = struct.unpack('i', file.read(4))[0]
                N6 = struct.unpack('i', file.read(4))[0]
                N7 = struct.unpack('i', file.read(4))[0]
                N8 = struct.unpack('i', file.read(4))[0]
                E3D.append([ID, N1, N2, N3, N4, N5, N6, N7, N8])

        elif BlockID == EN_SDB_ELEM_2D:
            E2DNUM = struct.unpack('i', file.read(4))[0]
            # Surf.count = E2DNUM
            tsurface = []
            for i in range(E2DNUM):
                ID = struct.unpack('i', file.read(4))[0]
                N1 = struct.unpack('i', file.read(4))[0]
                N2 = struct.unpack('i', file.read(4))[0]
                N3 = struct.unpack('i', file.read(4))[0]
                N4 = struct.unpack('i', file.read(4))[0]
                tsurface.append([ID, N1, N2, N3, N4])

        elif BlockID == EN_SDB_ELSET3D or BlockID == EN_SDB_ELSET2D:
            BlockLength = struct.unpack('i', file.read(4))[0]
            name = ''
            for i in range(BlockLength): 
                name += str(struct.unpack('c', file.read(1))[0])
            
            BlockFlag   = struct.unpack('i', file.read(4))[0]
            ENUM = struct.unpack('i', file.read(4))[0]
            ESET = []
            for i in range(ENUM):
                ESET.append(struct.unpack('i', file.read(4))[0])
            del(ESET)

        elif BlockID == EN_SDB_RIM:
            ControlNodeID = struct.unpack('i', file.read(4))[0]
            X = struct.unpack('d', file.read(8))[0]
            Y = struct.unpack('d', file.read(8))[0]
            Z = struct.unpack('d', file.read(8))[0]
            R = struct.unpack('d', file.read(8))[0]
            W = struct.unpack('d', file.read(8))[0]
            G = struct.unpack('i', file.read(4))[0]
            sfric.Rim.Add([ControlNodeID, X, Y, Z])

            for i in range(G):
                BlockFlag   = struct.unpack('i', file.read(4))[0]

            for i in range(G):
                X1 = struct.unpack('d', file.read(8))[0]
                Y1 = struct.unpack('d', file.read(8))[0]
                X2 = struct.unpack('d', file.read(8))[0]
                Y2 = struct.unpack('d', file.read(8))[0]
        
        elif BlockID == EN_SDB_ROAD:
            Tref = struct.unpack('d', file.read(8))[0]
            ControlNodeID = struct.unpack('i', file.read(4))[0]
            X = struct.unpack('d', file.read(8))[0]
            Y = struct.unpack('d', file.read(8))[0]
            Z = struct.unpack('d', file.read(8))[0]
            R = struct.unpack('d', file.read(8))[0]
            W = struct.unpack('d', file.read(8))[0]
            L = struct.unpack('d', file.read(8))[0]
            sfric.Road.Add([ControlNodeID, X, Y, Z])

        elif BlockID == EN_SDB_END:
            BlockID     = struct.unpack('i', file.read(4))[0]
            ESP         = ''

            for i in range(4):
                ESP += str(struct.unpack('c', file.read(1))[0])

            RecordHead  = struct.unpack('i', file.read(4))[0]
            BlockID     = struct.unpack('i', file.read(4))[0]
            SolverInfo  = ''

            for i in range(42): 
                SolverInfo += str(struct.unpack('c', file.read(1))[0])
            for i in range(9):
                BlockID     = struct.unpack('i', file.read(4))[0] 
            
            SimulationType = ''
            for i in range(26):
                SimulationType += str(struct.unpack('c', file.read(1))[0])

        else:
            break

    file.close()

    if len(sfric.Node.Node) == 0 : 
        print ("There is no information in 'sfric'")
        return 0
    if not deformed : return 0
    if isfile(sfricresultfile) == False: 
        print ("No sfric result file")
        return 0

    initNodes = np.array(sfric.Node.Node)

    if sfricresultfile != "":
        resultsfile = open(sfricresultfile, 'rb')
        resultsfile.seek(0,2); fend = resultsfile.tell() #find the end position of binary file
        resultsfile.seek(0);   fpos = resultsfile.tell() #find the start position of binary file
    else:
        return 0

    DISP=[]
    DeformedRigid=NODE()
    while resultsfile.tell() < fend:
        RecordHeaderID     = struct.unpack('i', resultsfile.read(4))[0]
        # print ("**", RecordHeaderID)
        if RecordHeaderID == EN_RESULT_TITLE:
            RecordValue    = struct.unpack('i', resultsfile.read(4))[0]
            OutputStepTime = struct.unpack('d', resultsfile.read(8))[0]
            OutputStepID       = struct.unpack('i', resultsfile.read(4))[0]
            OutputStepName = ''
            for i in range(26): 
                OutputStepName += str(struct.unpack('c', resultsfile.read(1))[0])
            OutputStepNo      = struct.unpack('i', resultsfile.read(4))[0]
            OutputStepNode    = struct.unpack('i', resultsfile.read(4))[0]
            NodeNUM           = struct.unpack('i', resultsfile.read(4))[0]

            for i in range(NodeNUM): 
                NodeID = int(struct.unpack('i', resultsfile.read(4))[0])
                DISP.append([NodeID, 0.0,0.0,0.0])

        elif RecordHeaderID == EN_RESULT:
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
            for i in range(26):
                DATA_TITLE += str(struct.unpack('c', resultsfile.read(1))[0])
            if RecordValue == EN_RESULT_DIS_X:
                for i in range(NodeNUM):
                    Value = struct.unpack('d', resultsfile.read(8))[0]
                    if math.isnan(Value) == False: sfric.Node.Node[i][1] += Value
                    else:
                        print ("Error in sfric result file (Del X = not a number).")
                        return 0
            elif RecordValue == EN_RESULT_DIS_Y:
                for i in range(NodeNUM):
                    Value = struct.unpack('d', resultsfile.read(8))[0]
                    sfric.Node.Node[i][2] += Value

            elif RecordValue == EN_RESULT_DIS_Z:
                for i in range(NodeNUM):
                    Value = struct.unpack('d', resultsfile.read(8))[0]
                    sfric.Node.Node[i][3] += Value

            elif RecordValue == RESULT_ContactStressNormal:
                for i in range(NodeNUM):
                    Value = struct.unpack('d', resultsfile.read(8))[0]
                    sfric.Node.Node[i][4] += Value
            elif RecordValue == RESULT_ContactForceNormal:
                for i in range(NodeNUM):
                    Value = struct.unpack('d', resultsfile.read(8))[0]
                    sfric.Node.Node[i][5] += Value
            else:
                for i in range(NodeNUM):
                    Value = struct.unpack('d', resultsfile.read(8))[0]

        elif RecordHeaderID == EN_RESULT_ELEMENT:
            TREAD_ELMENTNUM = struct.unpack('i', resultsfile.read(4))[0]
            for i in range(TREAD_ELMENTNUM): 
                TreadID = struct.unpack('i', resultsfile.read(4))[0]
        elif RecordHeaderID == EN_RESULT_RIGID:
            ELMENT_NUM   = struct.unpack('i', resultsfile.read(4))[0]
            rd = "" 

            for i in range(ELMENT_NUM): 
                ID = struct.unpack('i', resultsfile.read(4))[0]
                DeformedRigid.Add([ID, 0.0, 0.0, 0.0])

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
                for i in range(26): 
                    DATA_TITLE += str(struct.unpack('c', resultsfile.read(1))[0])

                if RecordValue == EN_RESULT_DIS_X:
                    for i in range(ELMENT_NUM): 
                        Value = struct.unpack('d', resultsfile.read(8))[0]
                        DeformedRigid.Node[i][1] = Value
                elif RecordValue == EN_RESULT_DIS_Y:
                    for i in range(ELMENT_NUM): 
                        Value = struct.unpack('d', resultsfile.read(8))[0]
                        DeformedRigid.Node[i][2] = Value
                elif RecordValue == EN_RESULT_DIS_Z:
                    for i in range(ELMENT_NUM): 
                        Value = struct.unpack('d', resultsfile.read(8))[0]
                        DeformedRigid.Node[i][3] = Value
                else:
                    for i in range(ELMENT_NUM): 
                        Value = struct.unpack('d', resultsfile.read(8))[0]
        else: 
            break 

    resultsfile.close()

    for rdg in DeformedRigid.Node:
        if rdg[0] == sfric.Rim.Node[0][0]:  
            sfric.Rim.Node[0][1] += rdg[1]
            sfric.Rim.Node[0][2] += rdg[2]
            sfric.Rim.Node[0][3] += rdg[3]
        if rdg[0] == sfric.Road.Node[0][0]:  
            sfric.Road.Node[0][1] += rdg[1]
            sfric.Road.Node[0][2] += rdg[2]
            sfric.Road.Node[0][3] += rdg[3]

    minz  = sfric.Node.Node[0][3]
    for nd in sfric.Node.Node:
        if minz > nd[3]: 
            minz = nd[3] 
    for nd in sfric.Node.Node:
        if nd[3] < minz + ht:
            sfric.pNode.Add(nd)   ## bottom nodes... 

    sfric.Node.Node = np.array(sfric.Node.Node)
    minz = np.min(sfric.Node.Node[:,3])
    ix = np.where (sfric.Node.Node[:,3] < minz + ht)[0]
    sfric.pNode.Node = sfric.Node.Node[ix]
    ndx = sfric.pNode.Node[:,0]

    srf = np.array(tsurface)

    bodyoffset = -1 
    sf1000 = srf[:1000]
    for s in sf1000: 
        OF1 = s[1] - s[4]; OF2 = s[2] - s[3] 
        if OF1 == OF2 and OF1%1000  == 0 and OF1 > 2000 and OF1 <10001:
            bodyoffset = abs(OF1) 

    nsurf = np.max (srf[:,1])
    if nsurf > 10**7 and  bodyoffset == 10000: 
        ix = np.where(srf[:,1]>=10**7)[0]
        srf = srf[ix]
    elif nsurf > 10**6 and bodyoffset < 10000 and bodyoffset > 2000: 
        ix = np.where(srf[:,1]>=10**6)[0]
        srf = srf[ix]

    print ("BODY OFFSET %d"%(bodyoffset))
    for sf in srf:
        coord = []
        ex =0
        for i in range(1, 5):
            re = np.where(ndx == sf[i])
            if len(re[0]) >0 and sfric.pNode.Node[re[0][0]][3] < minz+sfht:
                coord.append(sfric.pNode.Node[re[0][0]])
            else:
                ex =1
                break
        if ex ==0:
            # sfric.Surface.Add([sf[0], sf[1], sf[2], sf[3], sf[4], coord])
            sfric.Surface.Surface.append([sf[0], sf[1], sf[2], sf[3], sf[4], coord])
    print ("No. of contact surface %d"%(len(sfric.Surface.Surface)))
    
    return initNodes

def lateralShift_fromSFRIC(sfric=None, initial=None, offset=10000, trd=10**7): 
    # center =  int(centers[0] - trd )% offset

    npn = np.array(sfric.Node.Node)
    ixs = np.where(npn[:,0]>trd)[0]
    tpn = npn[ixs]
    minHt = np.min(tpn[:,3])
    ixs = np.where(tpn[:,3]<minHt+0.03)[0]

    btms = tpn[ixs]
    ix1 = np.where(btms[:,1]>-0.01)[0]
    ix2 = np.where(btms[:,1]< 0.01)[0]
    ix = np.intersect1d(ix1, ix2)

    ix1 = np.where(btms[:,2]>-0.01)[0]
    ix2 = np.where(btms[:,2]< 0.01)[0]
    iy = np.intersect1d(ix1, ix2)

    ix = np.intersect1d(ix, iy)

    mids = btms[ix]

    up = np.max(mids[:,1])
    dn = np.min(mids[:,1])

    ix = np.where(mids[:,1]==up)[0][0]; up = mids[ix]
    ix = np.where(mids[:,1]==dn)[0][0]; dn = mids[ix]

    v1=[0, 0, 1, 0]
    v2 = [0, up[1]-dn[1], up[2]-dn[2], 0]
    print (v2)

    angle_1 = Angle_Between_Vectors(v1, v2)

    if up[2] > dn[2] : 
        angle_1 = -math.degrees(angle_1)
    else: 
        angle_1 = math.degrees(angle_1)

    
    ix = np.where(initial[:,0] == up[0])[0][0];up1=initial[ix]
    ix = np.where(initial[:,0] == dn[0])[0][0];dn1=initial[ix]

    v2 = [0, up1[1]-dn1[1], up1[2]-dn1[2], 0]
    print (v2)
    angle_2 = Angle_Between_Vectors(v1, v2)
    if up[2] > dn[2] : 
        angle_2 = -math.degrees(angle_2)
    else: 
        angle_2 = math.degrees(angle_2)
    if angle_2 > 90: 
        angle_2 = 180 - angle_2 

    angle = angle_1 - angle_2 

    print (" Rotation angle =%.2f (%.2f, %.2f)"%(angle, angle_1, angle_2))

    print ("%.3f, %.3f, , %.3f, %3f"%(up[2]*1000, up[1]*1000, up1[2]*1000, up1[1]*1000))
    print ("%.3f, %.3f, , %.3f, %3f"%(dn[2]*1000, dn[1]*1000, dn1[2]*1000, dn1[1]*1000))

    shift =  (up1[2]+dn1[2])/2 - (up[2]+dn[2])/2


    return  shift, angle  

def readSMART_Inp(filename): 
    with open(filename) as F: 
        lines = F.readlines()
    
    for line in lines: 
        if "*CONDITION_LOAD" in line: 
            data = line.split("=")[1].strip()
            dt = data.split(",")
            press = float(dt[0].strip())
            load =  float(dt[2].strip())
            speed = float(dt[3].strip())
        if "*CAMBER_ANGLE" in line: 
            data = line.split("=")[1].strip()
            camber = float(data.strip())
        if "*LATERAL_CONTROL" in line: 
            data = line.split("=")[1].strip()
            dt = data.split(",")
            if '1' in dt[0]: 
                slipangle = float(dt[1].strip())
            else: 
                slipangle = 0 
        if "*ROTATION_CONTROL" in line: 
            data = line.split("=")[1].strip()
            dt = data.split(",")
            if '1' in dt[0]: 
                slipratio = float(dt[1].strip())
            else: 
                slipratio = 0  
    return [press, load, speed, camber, slipangle, slipratio]   


class SFRIC: 
    def __init__(self):
        self.Node = NODE()
        self.Surface = SURFACE()
        self.Rim = NODE()
        self.Road = NODE()
        self.pNode = NODE()
        self.cNode = NODE()
