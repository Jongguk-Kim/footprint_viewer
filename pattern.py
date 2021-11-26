import numpy as np 
import math 

from basicfunctions import Angle_Between_Vectors, NormalVector_plane
class PATTERN: 
    def __init__(self, filename,  test=0, start_number=10000000):
        self.GlobalXY=21
        self.Node=[]   ## float : Id, x, y, z
        self.Solid=[]  ## int : id, n1, ..., n8, 6 or 8, center (*10^6 -> make it into integer)
        self.Beam=[];         self.UpFront=[];         self.UpBack=[];         self.LowBack=[];         self.Center=[]

        self.GrooveCenter=[];        self.TranslatedNode=[]
        
        self.Surface=[]
        self.leftprofile =[];        self.rightprofile = []
        self.pitch=[];        self.pitchsequence =[] 
        
        self.NoPitch = 0;        self.ModelGD = 0.0
        self.PatternWidth =0.0;        self.TreadDesignWidth = 0.0
        self.HalfDia = 0.0;        self.diameter = 0.0 
        self.upFrontMaxY = 0.0

        self.Nstart = 10**7;        self.Estart = 10**7
        self.profilescaling = 0.001;        self.pitchscaling = 0.001

        self.surf_pattern_neg_side=[];        self.surf_pattern_pos_side=[];        self.SF_fulldepthgrooveside=[];        self.Edge_bottomSurface = []

        self.IsError = 0;        self.errorcode = 0 
        self.KerfsideSurface = [];        self.edge_top_surface = []

        self.shoulderType = 'R'

        self.GrvUpNode=[] ## nodes on up pitch surface among main groove bottom  
        self.MainGrvEdgeGroup =[]

        ####################################################################################################
        print ("############################################")
        print ("## Reading Pattern Mesh file (*.ptn)")
        print ("############################################")
        print (filename)

        result_reading = self.ReadPtn(filename)
        if result_reading >= 100 or self.TreadDesignWidth==0: 
            self.IsError = 1
            return 
        
        if result_reading ==0 : 
            print ("## Error to read pattern mesh file")
            self.IsError =  1 
            return 
        self.npn = np.array(self.Node)
        self.nps = np.array(self.Solid)

        self.Node_Origin = np.array(self.npn) 
        
        self.pitchlength = self.Pitchlength()
        self.errsolid = [] 
        self.HalfDia = round(self.HalfDia*self.profilescaling, 5)

        if len(self.UpFront) > 0: 
            ix = np.where(self.npn[:,0]== self.UpFront[0][1])[0][0]
            if self.npn[ix][3] != self.HalfDia: 
                print ("# Half diameter=%.2fmm, Guide Ht=%.2fmm"%(self.HalfDia*1000, self.npn[ix][3]*1000))
                print ("  Nodes are shifted (%.4fmm)"%((self.HalfDia - self.npn[ix][3])*1000))
                self.npn[:,3] += self.HalfDia - self.npn[ix][3]

                self.diameter = round(self.npn[ix][3] *2.0, 5)
                
            nds = []
            for nd in self.UpFront:
                nds.append(nd[1]); nds.append(nd[2])
            nds = np.array(nds)
            nds = np.unique(nds)
            upnodes=[]
            for nd in nds:
                ix = np.where(self.npn[:,0]==nd)[0][0]
                upnodes.append(self.npn[ix])
            upnodes = np.array(upnodes)
            self.upFrontMaxY = np.max(upnodes[:, 2])
        else:
            self.upFrontMaxY = self.npn[:,2]

        self.nps, self.Surface, printout, self.errsolid =  Generate_all_surfaces_on_solid(self.npn, self.nps, diameter=self.diameter)
        self.freetop, self.freebottom, self.uncheckedfree, self.Surface = self.Top_Bottom_FreeSurfacesFromAllSurfaces_01(self.Surface, \
            self.npn, radius=self.diameter/2.0, margin=1.0E-03)
        tf =  open('topsurf.tmp', 'w')
        for sf in self.freetop: 
            if sf[10] ==0: 
                tf.write("%d, %d, %d, %d, %d\n"%(sf[0], sf[7], sf[8], sf[9], sf[9]))
            else: 
                tf.write("%d, %d, %d, %d, %d\n"%(sf[0], sf[7], sf[8], sf[9], sf[10]))
        tf.close()
    def Pitchlength(self): 
        pmin = 10000.0
        pmax = -10000.0
        x = int(self.GlobalXY/10)
        y = int(self.GlobalXY%10)
        for bm in self.Center:
            if pmin > bm[3][y-1]: pmin = bm[3][y-1]
            if pmin > bm[4][y-1]: pmin = bm[4][y-1]
            if pmax < bm[3][y-1]: pmax = bm[3][y-1]
            if pmax < bm[4][y-1]: pmax = bm[4][y-1]
        return (pmax-pmin)
    def SurfaceBoundary(self, surface):
        ## surface = [El_id, Face_Id(1~6), type(3 or 4), layer, center X, y, z, n1, n2, n3, n4]
        bndedge=[]
        alledge =[]
        for sf in surface:
            alledge.append([int(sf[7]), int(sf[8]), 0, sf[0]])
            alledge.append([int(sf[8]), int(sf[9]), 0, sf[0]])
            if sf[2] == 3: alledge.append([int(sf[9]), int(sf[7]), 0, sf[0]])
            else:
                alledge.append([int(sf[9]), int(sf[10]), 0, sf[0]])
                alledge.append([int(sf[10]), int(sf[7]), 0, sf[0]])

        npedge = np.array(alledge, dtype=np.int32)
        N = len(npedge)
        for i, eg in enumerate(npedge):
            if eg[2] == -1: continue
            bnd = 1

            ind1 = np.where(npedge[:, 1] == eg[0])
            if len(ind1[0]) > 0 : 
                N = len(ind1[0])
                for j in range(N):
                    if npedge[ind1[0][j]][0] == eg[1]: 
                        npedge[i][2] = -1
                        bnd = 0 
                        break 
            if bnd ==1:
                npedge[i][2] =1
                bndedge.append(npedge[i])
        
        return np.array(bndedge)
    def Top_Bottom_FreeSurfacesFromAllSurfaces_01(self, allSurface, npn, radius=0.0, margin=1.0E-03): 

        # t0 = time.time()

        nodes =[]
        for i, sf in enumerate(allSurface):
            tnode = sf[7:]
            nodes.append(np.sort(tnode))
        nodes = np.array(nodes)
        filter_heightmargin =radius - margin ## 1mm from ht. 
        filter_heightmargin_top =radius - margin/2.0
        
        # print ("ht margin for top=%.3f"%(filter_heightmargin_top*1000))
        # print ("ht margin for btm=%.3f"%(filter_heightmargin*1000))
        # print ("Radius =%.3f"%(radius*1000))

        free = []
        bottom = []
        topfree = []
        for i, sf in enumerate(allSurface): 
            # ifree =1
            # print (" %d, %d, %d"%(sf[0]-10**7, sf[1], sf[3]))
            if sf[3] ==99: continue 
            if sf[2] == 3:
                ind1 = np.where(nodes[:, 3] == nodes[i][3])[0]  ## because nodes[i][0] == 0 
                ind2 = np.where(nodes[:, 1] == nodes[i][1])[0]
                ind3 = np.where(nodes[:, 2] == nodes[i][2])[0]
                ind = np.intersect1d(ind1, ind2, assume_unique=True)
                ind = np.intersect1d(ind,  ind3, assume_unique=True) 
                m = 3                
            else: 
                ind1 = np.where(nodes[:, 0] == nodes[i][0])[0]
                ind2 = np.where(nodes[:, 1] == nodes[i][1])[0]
                ind3 = np.where(nodes[:, 2] == nodes[i][2])[0]
                ind4 = np.where(nodes[:, 3] == nodes[i][3])[0]
                ind = np.intersect1d(ind1, ind2, assume_unique=True)
                ind = np.intersect1d(ind,  ind3, assume_unique=True) 
                ind = np.intersect1d(ind,  ind4, assume_unique=True) 
                m = 0 

            if len(ind) ==2: 
                # allSurface[i][3] = 99  ## 99: not free surface 
                allSurface[ind[0]][3] = 99
                allSurface[ind[1]][3] = 99
            #     ifree = 0
            # if ifree ==1:  ## among free surfaces 
            elif len(ind) ==1: 
                ind= ind[0]
                idx1 = np.where(npn[:,0]==nodes[ind][m])[0][0]
                idx2 = np.where(npn[:,0]==nodes[ind][1])[0][0]
                idx3 = np.where(npn[:,0]==nodes[ind][2])[0][0]
                n1 = npn[idx1]; n2=npn[idx2]; n3 = npn[idx3]

                if allSurface[i][1] == 1 and n1[3] < filter_heightmargin and n2[3] < filter_heightmargin and n3[3] < filter_heightmargin:  ## bottom surface : face = 1
                    allSurface[i][3] = 199  ## bottom surface 
                    bottom.append(sf)
                elif allSurface[i][1] == 2 and n1[3] > filter_heightmargin_top and n2[3] > filter_heightmargin_top and n3[3] > filter_heightmargin_top:  ## top surface : face = 2
                    allSurface[i][3] = 101  ## top surface 
                    topfree.append(sf)
                else: #if ht < filter_heightmargin :
                    allSurface[i][3] = 100  ## free surface 
                    free.append(sf)
        
        # print (" btm surf %d, top=%d, free=%d"%(len(bottom), len(topfree), len(free)))
        # t1 = time.time(); print ("** Top/BTM %.3f "%(t1-t0)); t0 = time.time()
        ## verifying the searching bottom surface #############

        bnd_btm = self.SurfaceBoundary(bottom)
        bnd = np.array(bnd_btm)
        self.Edge_bottomSurface = bnd_btm 

        # print (" the No. of bottom edge =%d"%(len(bnd_btm)))

        bd1 = bnd[0]
        bnd = np.delete(bnd, 0, axis=0)
        # print ("*DEL ", bd1)
        i = 0 
        groups=[]
        group=[bd1]
        added = 0 
        while i < len(bnd): 
            ix = np.where(bnd[:,0]==bd1[1])[0]
            if len(ix) > 0: 
                group.append(bnd[ix[0]])
                # print ("DEL ", bnd[ix[0]])
                bd1 = bnd[ix[0]]

                bnd = np.delete(bnd, ix[0], axis=0)
            else: 
                groups.append(group)
                group=[]
                if len(bnd)> 0: 
                    group=[bnd[0]]
                    bd1 = bnd[0]
                    # print ("*DEL ", bd1)
                    bnd = np.delete(bnd, 0, axis=0)
                    
        if len(group) > 0: 
            groups.append(group)
            group = []

        # t1 = time.time(); print ("** BTM %.3f"%(t1-t0)); t0 = time.time()
        if len(groups) > 1: 
            mg = len(groups[0])
            ibt=0 
            for i, gr in enumerate(groups): 
                if len(gr) > mg: 
                    mg = len(gr)
                    ibt = i 
            
            btm = groups[ibt] 
            btmg = []
            for bt in btm: 
                ix = np.where(npn[:,0] == bt[0])[0][0]; n1 = npn[ix]
                ix = np.where(npn[:,0] == bt[1])[0][0]; n2 = npn[ix]
                if n1[2] != n2[2]: 
                    if n1[2] > n2[2]: btmg.append([bt[0], bt[1], bt[2], bt[3], n2, n1])
                    else:             btmg.append([bt[0], bt[1], bt[2], bt[3], n1, n2]) 
            
            i = 0 
            while i < len(bottom): 
                cn = [0, bottom[i][4], bottom[i][5], bottom[i][6]]
                ix = np.where(npn[:,0] == bottom[i][7])[0][0]; n1 = npn[ix]
                ix = np.where(npn[:,0] == bottom[i][9])[0][0]; n3 = npn[ix]
                f = 0 
                for bt in btmg: 
                    d, P = DistanceFromLineToNode2D(cn, [bt[4], bt[5]], xy=23)
                    if d <= 0.2E-03 and bt[4][2] <= P[2] and P[2] <= bt[5][2] : 
                        f = 1
                        break 
                    if bt[4][2] <= n1[2] and n1[2] <= bt[5][2] and abs(bt[5][3]-bt[4][3]) > abs(bt[5][3]-n1[3]): 
                        f = 1
                        break 
                    if bt[4][2] <= n3[2] and n3[2] <= bt[5][2] and abs(bt[5][3]-bt[4][3]) > abs(bt[5][3]-n3[3]): 
                        f = 1
                        break 

                if f == 0: 
                    free.append(bottom[i])
                    # print ("del", bottom[i][0]-10**7, "f=", bottom[i][1], "dist=%.3f"%(d*1000))
                    del(bottom[i])
                    
                    i -= 1
                i += 1 

        # t1 = time.time(); print ("** FN %.3f"%(t1-t0)); t0 = time.time()
        # t2 = time.time()
        # print ("#####################################")
        # print (" TIME TO SEARCH FREE SURFACE =%.2f"%(t2-t1))
        # print ("#####################################")
        
        return np.array(topfree), np.array(bottom), np.array(free), allSurface 
    def ReadPtn(self, filename, valuereturn=0):

        with open(filename) as PTN: 
            lines = PTN.readlines()
        cmd = ""
        depths=[]
        solidname = ''
        centerbeams = []
        upAft =[]; lwAft = []; UpFwd =[]
        for line in lines:
            if "Regenerated Pattern mesh from P3DM" in line: 
                print ("* This mesh was generated by P3DM.")
                print ("  it is already bended to fit a layout \n")
                return 100

            if "**" in line: 
                continue
            elif "*" in line:
                if "*ELEMENT" in line.upper() and "CGAX4" in line.upper() : 
                    print ("* This mesh may be not a pattern mesh.")
                    print ("  This mesh can not be expanded.")
                    return 101 
                if "PROFILE_SCALING" in line.upper():
                    data = line.split(":")
                    self.profilescaling = float(data[1].strip())
                elif "GROOVE_DEPTH" in line.upper() or ("GROOVE" in line.upper() and 'DEPTH' in line.upper()):
                    data = line.split(":")
                    self.ModelGD = float(data[1].strip()) 
                elif "HALF_DIAMETER" in line.upper():
                    data = line.split(":")
                    self.diameter = float(data[1].strip()) * 2.0 
                    self.HalfDia =  float(data[1].strip())
                elif "CENTER_ANGLE" in line:
                    data = line.split(":")
                    self.centerangle = float(data[1].strip()) 
                elif "PROFILE_LHS" in line.upper():
                    cmd = 'LHS'
                elif "PROFILE_RHS" in line.upper():
                    cmd = 'RHS'
                elif "PITCH_SCALING" in line.upper():
                    data = line.split(":")
                    self.pitchscaling = float(data[1].strip()) 
                elif "TREAD_DESIGN_WIDTH" in line.upper() or ("TREAD" in line.upper() and "DESIGN" in line.upper() and "WIDTH" in line.upper()):
                    data = line.split(":")
                    self.TreadDesignWidth = float(data[1].strip())
                elif "GUIDELINE_TOLERANCE" in line.upper():
                    data = line.split(":")
                    self.guidelinetolerance = float(data[1].strip()) 
                elif "PITCH_DEFINITION_FIRST" in line:
                    cmd = 'P1'
                elif "PITCH_ARRAY_FIRST" in line.upper():
                    data = line.split(",")
                    for dt in data: 
                        if 'EIDSTART' in dt: 
                            self.Estart = int(dt.split("=")[1].strip())
                            self.Estart = 10**7
                        if 'NIDSTART' in dt: 
                            self.Nstart = int(dt.split("=")[1].strip())
                            self.Nstart = 10**7
                        if 'EIDOFFSET' in dt: 
                            self.Eoffset = int(dt.split("=")[1].strip())
                            self.Eoffset = 10000
                        if 'NIDOFFSET' in dt: 
                            self.Noffset = int(dt.split("=")[1].strip())
                            self.Noffset =10000
                        if 'ANGLE' in dt: self.Pangle = float(dt.split("=")[1].strip())
                        if 'DIRECTION' in dt: self.Direction = dt.split("=")[1].strip()
                    cmd = 'PS'
                elif "HEADING" in line.upper():
                    cmd = 'Heading'
                elif "*NODE" in line.upper() and not "FILE" in line.upper() and not 'OUTPUT' in line.upper() :
                    cmd = "ND"
                elif "ELEMENT" in line.upper() and "B31" in line.upper():
                    cmd = "BM"
                elif "ELEMENT" in line.upper() and "C3D8" in line.upper():
                    cmd = "SD8"
                elif "ELEMENT" in line.upper() and "C3D6" in line.upper():
                    cmd = "SD6"
                elif "*ELSET" in line.upper() and "SOLID" in line.upper():
                    if 'generate' in line.lower(): 
                        cmd = 'SOLG'
                        dline = line.split(",")
                        for d in dline: 
                            if "ELSET" in d.upper() and "=" in d: 
                                solidname = d.split("=")[1].strip()
                                break 
                    else:
                        cmd ="SOL"
                        dline = line.split(",")
                        for d in dline: 
                            if "ELSET" in d.upper() and "=" in d: 
                                solidname = d.split("=")[1].strip()
                                break 
                elif "*ELSET" in line.upper() and "CENTER" in line.upper():
                    if 'generate' in line.lower(): 
                        cmd = 'BCENG'
                        dline = line.split(",")
                        for d in dline: 
                            if "ELSET" in d.upper() and "=" in d: 
                                beamname = d.split("=")[1].strip()
                                break
                    else:
                        cmd ="BCEN"
                        dline = line.split(",")
                        for d in dline: 
                            if "ELSET" in d.upper() and "=" in d: 
                                beamname = d.split("=")[1].strip()
                                break 
                elif "ELSET" in line.upper() and "UPFWD" in line.upper():
                    if 'generate' in line.lower(): 
                        cmd = 'BUFG'
                        dline = line.split(",")
                        for d in dline: 
                            if "ELSET" in d.upper() and "=" in d: 
                                beamname = d.split("=")[1].strip()
                                break
                    else:   
                        cmd = "BUF"
                        dline = line.split(",")
                        for d in dline: 
                            if "ELSET" in d.upper() and "=" in d: 
                                beamname = d.split("=")[1].strip()
                                break
                elif "ELSET" in line.upper() and "UPAFT" in line.upper():
                    if 'generate' in line.lower(): 
                        cmd = 'BUAG'
                        dline = line.split(",")
                        for d in dline: 
                            if "ELSET" in d.upper() and "=" in d: 
                                beamname = d.split("=")[1].strip()
                                break
                    else:   
                        cmd = "BUA"
                        dline = line.split(",")
                        for d in dline: 
                            if "ELSET" in d.upper() and "=" in d: 
                                beamname = d.split("=")[1].strip()
                                break
                elif "ELSET" in line.upper() and "LWAFT" in line.upper():
                    if 'generate' in line.lower(): 
                        cmd = 'BLAG'
                        dline = line.split(",")
                        for d in dline: 
                            if "ELSET" in d.upper() and "=" in d: 
                                beamname = d.split("=")[1].strip()
                                break
                    else:
                        cmd = "BLA"
                        dline = line.split(",")
                        for d in dline: 
                            if "ELSET" in d.upper() and "=" in d: 
                                beamname = d.split("=")[1].strip()
                                break
                elif "TREADPTN_NIDSTART_NIDOFFSET_EIDSTART_EIDOFFSET" in line.upper(): 
                    return 0 
                else:
                    cmd =""
            else:
                if cmd =="LHS": 
                    data = line.split(",")
                    PR = round(float(data[0].strip())*self.profilescaling, 6)
                    if PR == 0.0 : PR = 10.0
                    self.leftprofile.append([PR, float(data[1].strip())*self.profilescaling])
                if cmd =="RHS": 
                    data = line.split(",")
                    PR = round(float(data[0].strip())*self.profilescaling, 6)
                    if PR == 0.0 : PR = 10.0
                    self.rightprofile.append([PR, float(data[1].strip())*self.profilescaling])
                if cmd =="P1": 
                    data = line.split(",")
                    self.pitch.append([data[0].strip(), data[4].strip(), data[5].strip(), data[6].strip(), data[7].strip(), data[8].strip()])
                if cmd =="PS":
                    data = line.split(",")
                    if len(data) == 3: 
                        self.pitchsequence.append([int(data[0].strip()), data[1].strip(), int(data[2].strip()) ])
                if cmd =="ND":
                    skip = 0 
                    data = line.split(",")
                    if len(data)>=4: 
                        if float(data[3].strip()) > 0: 
                            self.Node.append([float(int(data[0].strip()) + self.Nstart),  round(float(data[1].strip()) * self.pitchscaling, 7), \
                                        round(float(data[2].strip()) * self.pitchscaling, 7),  round(float(data[3].strip()) * self.pitchscaling, 7)])
                        else: 
                            print ("* All the Z values of the nodes should be positive")
                            print ("  This mesh can not be expanded.")
                            return 101 

                    # except:
                    #     print ("ERROR!! ", line, data)
                    #     return 101
                        

                    # if self.MaxY < float(data[2].strip()) * self.pitchscaling: self.MaxY = float(data[2].strip()) * self.pitchscaling
                    # if self.MinY > float(data[2].strip()) * self.pitchscaling: self.MinY = float(data[2].strip()) * self.pitchscaling
                if cmd == "BM": 
                    data = line.split(",")
                    if len(data) > 1: 
                        self.Beam.append([int(data[0].strip())+self.Estart, int(data[1].strip())+self.Nstart, int(data[2].strip())+self.Nstart])
                if cmd == "SD6": 
                    data = line.split(",")
                    # int(data[0].strip()) + self.Estart
                    self.Solid.append([int(data[0].strip()) + self.Estart, \
                        int(data[1].strip()) + self.Nstart, int(data[3].strip()) + self.Nstart, \
                        int(data[2].strip()) + self.Nstart, int(data[4].strip()) + self.Nstart, \
                        int(data[6].strip()) + self.Nstart, int(data[5].strip()) + self.Nstart, \
                        0, 0, 6])
                if cmd == "SD8": 
                    data = line.split(",")
                    self.Solid.append([int(data[0].strip()) + self.Estart, \
                        int(data[1].strip()) + self.Nstart, int(data[4].strip()) + self.Nstart, \
                        int(data[3].strip()) + self.Nstart, int(data[2].strip()) + self.Nstart, \
                        int(data[5].strip()) + self.Nstart, int(data[8].strip()) + self.Nstart, \
                        int(data[7].strip()) + self.Nstart, int(data[6].strip()) + self.Nstart, 8])                    
                if cmd == "BCEN":
                    data = line.split(",")
                    temp = [beamname]
                    for dt in data:
                        if dt.strip() =="": continue 
                        if dt.strip() !="": bid =  int(dt.strip()) + self.Estart
                        for bm in self.Beam: 
                            if bm[0] == bid: 
                                tmp = [bm[0], bm[1], bm[2]]
                                for nd in self.Node:
                                    if bm[1] == nd[0]:
                                        tmp.append([nd[1], nd[2], nd[3]])
                                        break
                                for nd in self.Node:
                                    if bm[2] == nd[0]:
                                        tmp.append([nd[1], nd[2], nd[3]])
                                        break
                                
                                self.Center.append(tmp)
                                temp.append(tmp)
                    centerbeams.append(temp)

                                
                if cmd == 'BCENG': 
                    data = line.split(",")
                    if len(data) < 3 : continue
                    data[0] = int(data[0].strip())
                    data[1] = int(data[1].strip())
                    data[2] = int(data[2].strip())
                    temp = [beamname]
                    for dn in range(data[0], data[1]+1, data[2]): 
                        bid = dn + self.Estart 
                        for bm in self.Beam: 
                            if bm[0] == bid: 
                                tmp = [bm[0], bm[1], bm[2]]
                                for nd in self.Node:
                                    if bm[1] == nd[0]:
                                        tmp.append([nd[1], nd[2], nd[3]])
                                        break
                                for nd in self.Node:
                                    if bm[2] == nd[0]:
                                        tmp.append([nd[1], nd[2], nd[3]])
                                        break
                                
                                self.Center.append(tmp)
                                temp.append(tmp)
                    centerbeams.append(temp)
                
                if cmd == "BUF":
                    data = line.split(",")
                    temp = [beamname]
                    for dt in data:
                        if dt.strip() =="": continue 
                        bid =  int(dt.strip()) + self.Estart
                        for bm in self.Beam: 
                            if bm[0] == bid: 
                                tmp = [bm[0], bm[1], bm[2]]
                                for nd in self.Node:
                                    if bm[1] == nd[0]:
                                        tmp.append([nd[1], nd[2], nd[3]])
                                        break
                                for nd in self.Node:
                                    if bm[2] == nd[0]:
                                        tmp.append([nd[1], nd[2], nd[3]])
                                        break
                                
                                self.UpFront.append(tmp)
                                temp.append(tmp)
                    UpFwd.append(temp)
                if cmd == 'BUFG': 
                    data = line.split(",")
                    if len(data) < 3 : continue
                    data[0] = int(data[0].strip())
                    data[1] = int(data[1].strip())
                    data[2] = int(data[2].strip())
                    temp = [beamname]
                    for dn in range(data[0], data[1]+1, data[2]): 
                        bid = dn + self.Estart 
                        for bm in self.Beam: 
                            if bm[0] == bid: 
                                tmp = [bm[0], bm[1], bm[2]]
                                for nd in self.Node:
                                    if bm[1] == nd[0]:
                                        tmp.append([nd[1], nd[2], nd[3]])
                                        break
                                for nd in self.Node:
                                    if bm[2] == nd[0]:
                                        tmp.append([nd[1], nd[2], nd[3]])
                                        break
                                
                                self.UpFront.append(tmp)  # # upAft =[]; lwAft = []; UpFwd =[]
                                temp.append(tmp)
                    UpFwd.append(temp)
                if cmd == "BUA":
                    data = line.split(",")
                    temp = [beamname]
                    for dt in data:
                        if dt.strip() =="": continue 
                        bid =  int(dt.strip()) + self.Estart
                        for bm in self.Beam: 
                            if bm[0] == bid: 
                                tmp = [bm[0], bm[1], bm[2]]
                                for nd in self.Node:
                                    if bm[1] == nd[0]:
                                        tmp.append([nd[1], nd[2], nd[3]])
                                        break
                                for nd in self.Node:
                                    if bm[2] == nd[0]:
                                        tmp.append([nd[1], nd[2], nd[3]])
                                        break
                                
                                self.UpBack.append(tmp) # # upAft =[]; lwAft = []; UpFwd =[]
                                temp.append(tmp)
                    upAft.append(temp)
                if cmd == 'BUAG': 
                    data = line.split(",")
                    if len(data) < 3 : continue
                    data[0] = int(data[0].strip())
                    data[1] = int(data[1].strip())
                    data[2] = int(data[2].strip())
                    temp = [beamname]
                    for dn in range(data[0], data[1]+1, data[2]): 
                        bid = dn + self.Estart 
                        for bm in self.Beam: 
                            if bm[0] == bid: 
                                tmp = [bm[0], bm[1], bm[2]]
                                for nd in self.Node:
                                    if bm[1] == nd[0]:
                                        tmp.append([nd[1], nd[2], nd[3]])
                                        break
                                for nd in self.Node:
                                    if bm[2] == nd[0]:
                                        tmp.append([nd[1], nd[2], nd[3]])
                                        break
                                
                                self.UpBack.append(tmp)# # upAft =[]; lwAft = []; UpFwd =[]
                                temp.append(tmp)
                    upAft.append(temp)
                if cmd == "BLA":
                    data = line.split(",")
                    temp = [beamname]
                    for dt in data:
                        if dt.strip() =="": continue 
                        bid =  int(dt.strip()) + self.Estart
                        for bm in self.Beam: 
                            if bm[0] == bid: 
                                tmp = [bm[0], bm[1], bm[2]]
                                for nd in self.Node:
                                    if bm[1] == nd[0]:
                                        tmp.append([nd[1], nd[2], nd[3]])
                                        break
                                for nd in self.Node:
                                    if bm[2] == nd[0]:
                                        tmp.append([nd[1], nd[2], nd[3]])
                                        break
                                
                                self.LowBack.append(tmp)# # upAft =[]; lwAft = []; UpFwd =[]
                                temp.append(tmp)
                    lwAft.append(temp)
                if cmd == 'BLAG': 
                    data = line.split(",")
                    if len(data) < 3 : continue
                    data[0] = int(data[0].strip())
                    data[1] = int(data[1].strip())
                    data[2] = int(data[2].strip())
                    temp = [beamname]
                    for dn in range(data[0], data[1]+1, data[2]): 
                        bid = dn + self.Estart 
                        for bm in self.Beam: 
                            if bm[0] == bid: 
                                tmp = [bm[0], bm[1], bm[2]]
                                for nd in self.Node:
                                    if bm[1] == nd[0]:
                                        tmp.append([nd[1], nd[2], nd[3]])
                                        break
                                for nd in self.Node:
                                    if bm[2] == nd[0]:
                                        tmp.append([nd[1], nd[2], nd[3]])
                                        break
                                
                                self.LowBack.append(tmp)# # upAft =[]; lwAft = []; UpFwd =[]
                                temp.append(tmp)
                    lwAft.append(temp)
                # self.pitch=[P1, ELSET_NAME]
                if cmd == 'SOLG': 
                    temp = []
                    data = line.split(",")
                    if len(data) < 3 : continue
                    data[0] = int(data[0].strip())
                    data[1] = int(data[1].strip())
                    data[2] = int(data[2].strip())
                    for dn in range(data[0], data[1]+1, data[2]): 
                        temp.append(dn + self.Estart)

                    for p in self.pitch: 
                        if p[1] == solidname: 
                            if len(p) == 6: 
                                p.append(temp)
                                # print (" solidname : ", solidname, temp)
                            else: 
                                for t in temp: 
                                    p[-1].append(t)
                            # print (p)

                if cmd == 'SOL':
                    temp = [] 
                    data = line.split(",")
                    for d in data: 
                        d = d.strip()
                        if d!="": 
                            temp.append(int(d)+self.Estart)
                    
                    for p in self.pitch: 
                        
                        if p[1] == solidname:
                            if len(p) == 6: 
                                p.append(temp)
                            #     # p.append(data)
                            #     # print (" solidname", solidname, data)
                            else: 
                                for t in temp: 
                                    p[-1].append(t)

                            # print (p)


        if len(self.Node) > 100000: 
            print ("\n Too many nodes in the mesh (=%d >100,000)"%(len(self.Node)))
            return 0
        if len(self.Solid) == 0: 
            print ("\n No information of pattern mesh\n")
            return 0

        if len(self.pitch) > 1:
            p0=self.pitch[0][0] 
            s0 =self.pitch[0][1] 
            i = 1
            while i < len(self.pitch): 
                if self.pitch[i][0] != p0 and self.pitch[i][1] != s0: 
                    self.CombinePitches(self.pitch, self.pitchsequence, centerbeams, upAft, lwAft, UpFwd)
                    break 
                i += 1 
        
        self.TreadDesignWidth = round(self.TreadDesignWidth/1000, 9)

        tnode = np.array(self.Node)
        zs = tnode[:,3]
        self.diameter = np.max(zs) * 2.0

        # self.diameter = round(self.diameter * self.pitchscaling, 6)
        # self.guidelinetolerance = round(self.guidelinetolerance * self.pitchscaling, 6)
        ix = np.where(tnode[:,3] > self.diameter /2.0 - 0.0005)[0]
        wn = tnode[ix]
        ws = wn[:,2]
        wmin = np.min(ws); wmax=np.max(ws)
        self.PatternWidth = wmax - wmin 
        if self.ModelGD ==0: self.ModelGD = 1.0E-03 
        else: self.ModelGD = self.ModelGD * 0.001

        # if self.TreadDesignWidth == 0: 
        #     self.TreadDesignWidth = self.PatternWidth  - 10.0E-03 
        #     print ("*Design Width was set to 'Total width -10mm'")

        return 1 


def Generate_all_surfaces_on_solid(npn, nps, diameter=0, text=""):  ##  --> GenerateAllSurfaces(self) 
    ## before call this function, 'makenumpyarray()' should be called.

    Solid = []
    Surface=[]
    if diameter == 0: 
        zs = npn[:,3]
        diameter = np.max(zs) * 2 
    R = round(diameter / 2.0, 7)
    R1000 = R*1000
    changed = 0
    img = 0 
    topmargin = R - 0.5E-03
    topmax = 0
    btmmax = 0 
    truncation = 5
    layer = -1

    checkEL = 5841; checkEL1 = 5862

    vertLine = [0, 0, 0, 1]

    for solid in nps:

        
        if solid[9] == 6: N=7     ## solid[9] == 6
        else: N=9                 ## solid[9] == 8

        nodes_coord = []
        elementNodes = []
        angles=[]

        sx = 0.0;   sy = 0.0;      sz = 0.0
        topnode = []
        for i in range(1, N):
            index = np.where(npn[:, 0] == solid[i])
            tx = npn[index[0][0]][1]; ty=npn[index[0][0]][2]; tz =npn[index[0][0]][3]
            sx += tx;                  sy += ty;                 sz += tz 
            if abs(round(tz -R, truncation)) < 0.5E-3 : topnode.append(npn[index[0][0]][0])
            nodes_coord.append([tx,ty,tz])
            idx = index[0][0]
            elementNodes.append(npn[index[0][0]])

        SolidCenter = [round(sx/solid[9], truncation), round(sy/solid[9], truncation), round(sz/solid[9], truncation)] 

        centers = []

        top = 0 
        gap = R*2 
        
        if N == 9: 
            i = 0; j=1; m=2; n=3
            ssx = round((nodes_coord[i][0] + nodes_coord[j][0] + nodes_coord[m][0] + nodes_coord[n][0])  / 4.0, truncation)
            ssy = round((nodes_coord[i][1] + nodes_coord[j][1] + nodes_coord[m][1] + nodes_coord[n][1])  / 4.0, truncation)
            ssz = round((nodes_coord[i][2] + nodes_coord[j][2] + nodes_coord[m][2] + nodes_coord[n][2])  / 4.0, truncation)
            centers.append([ssx, ssy, ssz, 1])
            Normal = NormalVector_plane(elementNodes[i], elementNodes[j], elementNodes[m])
            angles.append(Angle_Between_Vectors(Normal, vertLine)-1.570796)
            if ssz >= topmargin and R-ssz < gap : 
                top =1 
                gap = R-ssz 

            i = 4; j=5; m=6; n=7
            ssx = round((nodes_coord[i][0] + nodes_coord[j][0] + nodes_coord[m][0] + nodes_coord[n][0])  / 4.0, truncation)
            ssy = round((nodes_coord[i][1] + nodes_coord[j][1] + nodes_coord[m][1] + nodes_coord[n][1])  / 4.0, truncation)
            ssz = round((nodes_coord[i][2] + nodes_coord[j][2] + nodes_coord[m][2] + nodes_coord[n][2])  / 4.0, truncation) 
            centers.append([ssx, ssy, ssz, 2])
            Normal = NormalVector_plane(elementNodes[i], elementNodes[j], elementNodes[m])
            angles.append(Angle_Between_Vectors(Normal, vertLine)-1.570796)

            if ssz >= topmargin and R-ssz < gap : 
                top =2 
                gap = R-ssz 

            i = 0; j=1; m=5; n=4
            ssx = round((nodes_coord[i][0] + nodes_coord[j][0] + nodes_coord[m][0] + nodes_coord[n][0])  / 4.0, truncation)
            ssy = round((nodes_coord[i][1] + nodes_coord[j][1] + nodes_coord[m][1] + nodes_coord[n][1])  / 4.0, truncation)
            ssz = round((nodes_coord[i][2] + nodes_coord[j][2] + nodes_coord[m][2] + nodes_coord[n][2])  / 4.0, truncation) 
            centers.append([ssx, ssy, ssz, 3])
            Normal = NormalVector_plane(elementNodes[i], elementNodes[j], elementNodes[m])
            angles.append(Angle_Between_Vectors(Normal, vertLine)-1.570796)
            if ssz >= topmargin and R-ssz < gap : 
                top =3 
                gap = R-ssz 

            i = 1; j=2; m=6; n=5
            ssx = round((nodes_coord[i][0] + nodes_coord[j][0] + nodes_coord[m][0] + nodes_coord[n][0])  / 4.0, truncation)
            ssy = round((nodes_coord[i][1] + nodes_coord[j][1] + nodes_coord[m][1] + nodes_coord[n][1])  / 4.0, truncation)
            ssz = round((nodes_coord[i][2] + nodes_coord[j][2] + nodes_coord[m][2] + nodes_coord[n][2])  / 4.0, truncation) 
            centers.append([ssx, ssy, ssz, 4])
            Normal = NormalVector_plane(elementNodes[i], elementNodes[j], elementNodes[m])
            angles.append(Angle_Between_Vectors(Normal, vertLine)-1.570796)
            if ssz >= topmargin and R-ssz < gap : 
                top =4 
                gap = R-ssz 

            i = 2; j=3; m=7; n=6
            ssx = round((nodes_coord[i][0] + nodes_coord[j][0] + nodes_coord[m][0] + nodes_coord[n][0])  / 4.0, truncation)
            ssy = round((nodes_coord[i][1] + nodes_coord[j][1] + nodes_coord[m][1] + nodes_coord[n][1])  / 4.0, truncation)
            ssz = round((nodes_coord[i][2] + nodes_coord[j][2] + nodes_coord[m][2] + nodes_coord[n][2])  / 4.0, truncation) 
            centers.append([ssx, ssy, ssz, 5])
            Normal = NormalVector_plane(elementNodes[i], elementNodes[j], elementNodes[m])
            angles.append(Angle_Between_Vectors(Normal, vertLine)-1.570796)
            if ssz >=topmargin and R-ssz < gap : 
                top =5 
                gap = R-ssz 

            i = 3; j=0; m=4; n=7
            ssx = round((nodes_coord[i][0] + nodes_coord[j][0] + nodes_coord[m][0] + nodes_coord[n][0])  / 4.0, truncation)
            ssy = round((nodes_coord[i][1] + nodes_coord[j][1] + nodes_coord[m][1] + nodes_coord[n][1])  / 4.0, truncation)
            ssz = round((nodes_coord[i][2] + nodes_coord[j][2] + nodes_coord[m][2] + nodes_coord[n][2])  / 4.0, truncation) 
            centers.append([ssx, ssy, ssz, 6])
            Normal = NormalVector_plane(elementNodes[i], elementNodes[j], elementNodes[m])
            angles.append(Angle_Between_Vectors(Normal, vertLine)-1.570796)
                
            if ssz >=topmargin and R-ssz < gap : 
                top =6 
                gap = R-ssz 
        else:
            i = 0; j=1; m=2
            ssx = round((nodes_coord[i][0] + nodes_coord[j][0] + nodes_coord[m][0] )  / 3.0, truncation)
            ssy = round((nodes_coord[i][1] + nodes_coord[j][1] + nodes_coord[m][1] )  / 3.0, truncation)
            ssz = round((nodes_coord[i][2] + nodes_coord[j][2] + nodes_coord[m][2] )  / 3.0, truncation) 
            centers.append([ssx, ssy, ssz, 1])
            if ssz >= topmargin and R-ssz < gap : 
                top =1 
                gap = R-ssz 

            i = 3; j=4; m=5
            ssx = round((nodes_coord[i][0] + nodes_coord[j][0] + nodes_coord[m][0] )  / 3.0, truncation)
            ssy = round((nodes_coord[i][1] + nodes_coord[j][1] + nodes_coord[m][1] )  / 3.0, truncation)
            ssz = round((nodes_coord[i][2] + nodes_coord[j][2] + nodes_coord[m][2] )  / 3.0, truncation) 
            centers.append([ssx, ssy, ssz, 2])
            if ssz >= topmargin and R-ssz < gap : 
                top =2
                gap = R-ssz 

            i = 0; j=1; m=4; n=3
            ssx = round((nodes_coord[i][0] + nodes_coord[j][0] + nodes_coord[m][0] + nodes_coord[n][0])  / 4.0, truncation)
            ssy = round((nodes_coord[i][1] + nodes_coord[j][1] + nodes_coord[m][1] + nodes_coord[n][1])  / 4.0, truncation)
            ssz = round((nodes_coord[i][2] + nodes_coord[j][2] + nodes_coord[m][2] + nodes_coord[n][2])  / 4.0, truncation) 
            centers.append([ssx, ssy, ssz, 3])
            if ssz >=topmargin and R-ssz < gap : 
                top =3 
                gap = R-ssz 

            i = 1; j=2; m=5; n=4
            ssx = round((nodes_coord[i][0] + nodes_coord[j][0] + nodes_coord[m][0] + nodes_coord[n][0])  / 4.0, truncation)
            ssy = round((nodes_coord[i][1] + nodes_coord[j][1] + nodes_coord[m][1] + nodes_coord[n][1])  / 4.0, truncation)
            ssz = round((nodes_coord[i][2] + nodes_coord[j][2] + nodes_coord[m][2] + nodes_coord[n][2])  / 4.0, truncation) 
            centers.append([ssx, ssy, ssz, 4])
            if ssz >=topmargin and R-ssz < gap : 
                top =4 
                gap = R-ssz  

            i = 2; j=0; m=3; n=5
            ssx = round((nodes_coord[i][0] + nodes_coord[j][0] + nodes_coord[m][0] + nodes_coord[n][0])  / 4.0, truncation)
            ssy = round((nodes_coord[i][1] + nodes_coord[j][1] + nodes_coord[m][1] + nodes_coord[n][1])  / 4.0, truncation)
            ssz = round((nodes_coord[i][2] + nodes_coord[j][2] + nodes_coord[m][2] + nodes_coord[n][2])  / 4.0, truncation) 
            centers.append([ssx, ssy, ssz, 5])
            if ssz >= topmargin and R-ssz < gap : 
                top =5 
                gap = R-ssz 
        
        if top == 0 or N != 9: 
            cnt = 0
            mz = 100000.0
            my1  =centers[0][1]
            my2 = centers[0][1]
            for cz in centers: 
                if cz[1] > my2: my2 = cz[1]
                if cz[1] < my1: my1 = cz[1]
            for i, cz in enumerate(centers): 
                if N == 9: 
                    if cz[2] < mz and abs(angles[i]) > 0.785398 : 
                        mz = cz[2]
                        if my1 < cz[1] and cz[1] < my2: 
                            cnt = i
                else: 
                    if cz[2] < mz  : 
                        mz = cz[2]
                        if my1 < cz[1] and cz[1] < my2: 
                            cnt = i
        elif N==9 : 
            if top == 2  :  cnt = 0 
            elif top == 1 : cnt = 1
            elif  top == 3 :cnt = 4
            elif top == 4 : cnt = 5
            elif top == 5 : cnt = 2
            else :          cnt = 3

        if N ==9: 
            if cnt == 0:   orders = [[0, 1, 4, 3, 2, 5, 8, 7, 6, 9] , [0, 1, 5, 4, 3, 2] ]  # [0, 1, 2, 3, 4, 5]  ## orders = [[el[0], N1[i]],   centers of faces ()]
            elif cnt == 1: orders = [[0, 6, 7, 8, 5, 2, 3, 4, 1, 9] , [1, 0, 3, 4, 5, 2] ]  # [[0,   6, 5, 8, 7,   2, 1, 4, 3, 9] , [1, 0,   2, 5, 4, 3] ]

            elif cnt == 2: orders = [[0, 1, 2, 6, 5, 4, 3, 7, 8, 9] , [2, 4, 0, 3, 1, 5] ]  # [[0,   1, 5, 6, 2,   4, 8, 7, 3, 9] , [2, 4,   5, 1, 3, 0] ]
            elif cnt == 3: orders = [[0, 2, 3, 7, 6, 1, 4, 8, 5, 9] , [3, 5, 0, 4, 1, 2] ]  # [[0,   2, 6, 7, 3,   1, 5, 8, 4, 9] , [3, 5,   2, 1, 4, 0] ]
            elif cnt == 4: orders = [[0, 4, 8, 7, 3, 1, 5, 6, 2, 9]  , [4, 2, 5, 1, 3, 0] ] # [[0,   4, 3, 7, 8,   1, 2, 6, 5, 9]  , [4, 2,   0, 3, 1, 5] ]
            elif cnt == 5: orders = [[0, 1, 5, 8, 4, 2, 6, 7, 3, 9] , [5, 3, 2, 1, 4, 0] ]  # [[0,   1, 4, 8, 5,   2, 3, 7, 6, 9] , [5, 3,   0, 4, 1, 2] ]
            else:  
                orders = [[0, 1, 4, 3, 2, 5, 8, 7, 6, 9] , [0, 1, 5, 4, 3, 2] ]
        else: 
            if cnt == 0: orders =   [[0, 1, 3, 2, 4, 6, 5, 7, 8, 9] , [0, 1, 4, 3, 2] ]    # [[0, 1, 2, 3, 4, 5, 6, 7, 8, 9] , [0, 1, 2, 3, 4] ]
            elif cnt == 1: orders = [[0, 4, 5, 6, 1, 2, 3, 7, 8, 9] , [1, 0, 2, 3, 4] ]    # [[0, 4, 6, 5, 1, 3, 2, 7, 8, 9] , [1, 0, 4, 3, 2] ]

            else: 
                orders =   [[0, 1, 3, 2, 4, 6, 5, 7, 8, 9] , [0, 1, 4, 3, 2] ] 
                errsolid =[solid]
        tmp = []
        for i in orders[0]: 
            tmp.append(solid[i])

        tmp.append(SolidCenter[0])
        tmp.append(SolidCenter[1])
        tmp.append(SolidCenter[2])

        for i in orders[0]: 
            if i ==9: break 
            else: 
                if N == 7 and i > 6: 
                        tmp.append(0.0);                tmp.append(0.0);                tmp.append(0.0)
                else: 
                    tmp.append(nodes_coord[i-1][0])
                    tmp.append(nodes_coord[i-1][1])
                    tmp.append(nodes_coord[i-1][2])

        margin = 0.15E-03
        pl =0 
        
        if N == 9: 
            # layer = -1
            Surface.append([tmp[0], 1, 4, layer, centers[orders[1][0]][0], centers[orders[1][0]][1], centers[orders[1][0]][2], tmp[1],tmp[2],tmp[3],tmp[4]])
            Surface.append([tmp[0],2, 4, layer, centers[orders[1][1]][0], centers[orders[1][1]][1], centers[orders[1][1]][2], tmp[5],tmp[6],tmp[7],tmp[8]])
            
            # layer = -1
            Surface.append([tmp[0],3, 4, layer, centers[orders[1][2]][0], centers[orders[1][2]][1], centers[orders[1][2]][2], tmp[1],tmp[2],tmp[6],tmp[5]])
            Surface.append([tmp[0],4, 4, layer, centers[orders[1][3]][0], centers[orders[1][3]][1], centers[orders[1][3]][2], tmp[2],tmp[3],tmp[7],tmp[6]])
            Surface.append([tmp[0],5, 4, layer, centers[orders[1][4]][0], centers[orders[1][4]][1], centers[orders[1][4]][2], tmp[3],tmp[4],tmp[8],tmp[7]])
            Surface.append([tmp[0],6, 4, layer, centers[orders[1][5]][0], centers[orders[1][5]][1], centers[orders[1][5]][2], tmp[4],tmp[1],tmp[5],tmp[8]])
        else: 
            # layer = -1
            Surface.append([tmp[0],1, 3, layer, centers[orders[1][0]][0], centers[orders[1][0]][1], centers[orders[1][0]][2], tmp[1],tmp[2],tmp[3], 0])
            Surface.append([tmp[0],2, 3, layer, centers[orders[1][1]][0], centers[orders[1][1]][1], centers[orders[1][1]][2], tmp[4],tmp[5],tmp[6],0 ])
            
            # layer = -1
            Surface.append([tmp[0],3, 4, layer, centers[orders[1][2]][0], centers[orders[1][2]][1], centers[orders[1][2]][2], tmp[1],tmp[2],tmp[5],tmp[4]])
            Surface.append([tmp[0],4, 4, layer, centers[orders[1][3]][0], centers[orders[1][3]][1], centers[orders[1][3]][2], tmp[2],tmp[3],tmp[6],tmp[5]])
            Surface.append([tmp[0],5, 4, layer, centers[orders[1][4]][0], centers[orders[1][4]][1], centers[orders[1][4]][2], tmp[3],tmp[1],tmp[4],tmp[6]])

        Solid.append(tmp)
        if cnt >0: 
            changed+=1

    del(centers)
    del(orders)
    nps = np.array(Solid, dtype=np.float64)
    Surface = np.array(Surface)
    text += "\n* The node order of PTN was checked(%dEA)\n"%(changed)
    errsolid = []
    return nps, Surface, text, errsolid

def DistanceFromLineToNode2D(N0, nodes=[], xy=12, onlydist=0):
    x = int(xy/10)
    y = int(xy%10)

    N1=nodes[0]
    N2=nodes[1]
    N=[-1, 0, 0, 0]
    if round(N2[x]-N1[x], 8) and round(N2[y]-N1[y], 8): 
        A =  (N2[y]-N1[y])/(N2[x]-N1[x]) 
        cx = A / (A*A+1) *(N0[x]/A + N0[y] + A*N1[x] - N1[y])
        cy = A * (cx - N1[x]) + N1[y] 
        N[x] = cx
        N[y] = cy
        distance = math.sqrt((cx-N0[x])**2 + (cy-N0[y])**2)

    elif round(N2[x]-N1[x], 8) and round(N2[y]-N1[y], 8) == 0: 
        distance = abs(N0[y] - N1[y])
        N[x] = N0[x]
        N[y] = N1[y]

    else: 
        distance = abs(N0[x] - N1[x])
        
        N[x] = N1[x]
        N[y] = N0[y]
    if onlydist ==1: 
        return distance
    else: 
        return distance, N 
