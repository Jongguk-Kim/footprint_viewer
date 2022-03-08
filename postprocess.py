import numpy as np 
from basicfunctions import Sorting, Area, findout_offset, meshgrid_in_Quadrilateral,\
     Delete_Close_Points, Angle_Between_Vectors
from os.path import isfile 
from os import getcwd 
import glob, math  
from layout import NODE 
from pattern import PATTERN
import MeshViewer as cute  

import time 

def timer(func): 
    def wrapper(*args, **kwargs): 
        start = time.time()
        rv = func(*args, **kwargs)
        total = time.time() - start
        print ("### Time to generate mesh grid: %.2f"%(total))
        return rv 
    return wrapper

def Edge_Boundary(surf): 
    
    edges = []
    sfedge = []
    for sf in surf: 
        edges.append([sf[1], sf[2], 1, sf[0], 0, sf[5][0], sf[5][1] ])
        edges.append([sf[2], sf[3], 2, sf[0], 0, sf[5][1], sf[5][2]])
        sfedge.append([sf[1], sf[2], 1, sf[0], 0])
        sfedge.append([sf[2], sf[3], 2, sf[0], 0])
        if sf[4] > 10**7: 
            edges.append([sf[3], sf[4], 3, sf[0], 0, sf[5][2], sf[5][3]])
            edges.append([sf[4], sf[1], 4, sf[0], 0, sf[5][3], sf[5][0]])
            sfedge.append([sf[3], sf[4], 3, sf[0], 0])
            sfedge.append([sf[4], sf[1], 4, sf[0], 0])
        else: 
            edges.append([sf[3], sf[1], 3, sf[0], 0, sf[5][2], sf[5][0]])
            sfedge.append([sf[3], sf[1], 3, sf[0], 0])

    sfedge = np.array(sfedge)

    bde=[]
    for i, ed in enumerate(sfedge): 
        if edges[i][4] >0: 
            continue 
        ix1 = np.where(sfedge[:,:2]==ed[0])[0]
        ix2 = np.where(sfedge[:,:2]==ed[1])[0]
        ix = np.intersect1d(ix1, ix2)
        if len(ix) ==1: 
            edges[i][4] = 1
            bde.append(edges[i])
        else:
            N = len(ix)
            for x in ix:
                edges[x][4] = N 
    return bde 


class FOOTPRINT(): 
    def __init__(self) :

        ## file information 
        # self.file = ""  
        self.ptn =None ;     self.mesh = None;   self.smart = None 
        self.isPTN = False
        self.sfric = False 
        self.workingfile = None 
        self.ISLM = False 
        self.ISLM_cali  = False 
        self.tireSize = False 
        self.tirePattern = False 


        ## from footprint data 
        self.surface_values=None 

        ## initial values 
        self.px =[];         self.py =[];        self.pv =[]
        self.initial_shift = 0 
        self.vmin = 50000
        self.vmax = 500000
        self.density = 20 
        self.group = 'pcr'
        self.fitting = 6 
        self.camber = 0 

        ## current values 
        self.xs=[]; self.ys=[]; self.vs =[]
        self.centerPress=[]

        ## values on board 

        # self.shift = 0 
        self.lateralShift = 0.0 
        self.angle = 0.0

        self.totalArea = 0.0
        self.actualArea = 0.0 
        self.roundness = 0.0
        self.lengths=[]
        self.widths=[]
        self.GWShape = 0.0
        self.localGWShape = 0
        self.roundnessShape = 0
        self.squarenessShape = 0
        self.FPClines =[]
        self.boundary=[]
        self.areaPoints=[]
        self.edge_groove=[]

        self.edges =None 

        self.centerPress=None 

        self.ISLM_totalArea = 0.0
        self.ISLM_actualArea = 0.0 
        self.ISLM_roundness = 0.0
        self.ISLM_lengths=list(range(20))
        self.ISLM_widths=list(range(20))
        self.ISLM_GWShape = 0.0
        self.ISLM_localGWShape = 0
        self.ISLM_roundnessShape = 0
        self.ISLM_squarenessShape = 0
        self.ISLM_FPC=None 
        self.ISLM_FPC_detail = None 
        self.ISLM_boundary=None 
        self.ISLM_boundary_Init=None 
        self.ISLM_contourX =[] 
        self.ISLM_contourY =[] 
        self.ISLM_contourV =[] 
        self.ISLM_centerPress=None 

        self.ISLM_cali_totalArea = 0.0
        self.ISLM_cali_actualArea = 0.0 
        self.ISLM_cali_roundness = 0.0
        self.ISLM_cali_lengths=list(range(20))
        self.ISLM_cali_widths=list(range(20))
        self.ISLM_cali_GWShape = 0.0
        self.ISLM_cali_localGWShape = 0
        self.ISLM_cali_roundnessShape = 0
        self.ISLM_cali_squarenessShape = 0

        self.ISLM_caliPressure = None 
        self.ISLM_caliFPC = None 
        self.ISLM_caliFPC_detail = None 
        self.ISLM_caliboundary=None 
        self.ISLM_caliboundary_Init=None 

        self.ISLM_cali_contourX =[] 
        self.ISLM_cali_contourY =[] 
        self.ISLM_cali_contourV =[] 

    def inputBoardValues(self, vmin=0,\
         vmax=0, density=0, group="",\
              file="", angle=0, shift=0): 
        self.vmin = vmin
        self.vmax = vmax 
        self.density = density  
        self.group = group 
        self.workingfile = file 
        self.angle = 0
        self.lateralShift = shift 
    def inputFitting(self, fitting): 
        self.fitting = fitting 

    def inputISLM_FPC(self, fname): 
        # print (" ISLM FPC Input")
        with open(fname) as F: 
            lines = F.readlines()
        for line in lines: 
            if 'ActualContactArea' in line: 
                wd = line.split("=")[1].strip()
                self.ISLM_actualArea = float(wd)/10000
                print ('* ActualContactArea', self.ISLM_actualArea)
            if 'TotalContactArea' in line: 
                wd = line.split("=")[1].strip()
                self.ISLM_totalArea = float(wd)/10000
                print ('* TotalContactArea', self.ISLM_totalArea)
            if 'Roundness(%)=' in line: 
                wd = line.split("=")[1].strip()
                self.ISLM_roundness = float(wd)/100
                print ('* Roundness', self.ISLM_roundness)
            if 'DetailedContactLength(mm)' in line: 
                wd = line.split("=")[1].strip()
                ws = wd.split("/")
                # foot.lengths[3], foot.lengths[5], foot.lengths[10], foot.lengths[15], foot.lengths[17]
                self.ISLM_lengths[3] = float(ws[0].strip())
                self.ISLM_lengths[5] = float(ws[1].strip())
                self.ISLM_lengths[10] = float(ws[2].strip())
                self.ISLM_lengths[15] = float(ws[3].strip())
                self.ISLM_lengths[17] = float(ws[4].strip())
                print ("* LENGTHS", self.ISLM_lengths)
            if 'DetailedContactWidth(mm)' in line: 
                wd = line.split("=")[1].strip()
                ws = wd.split("/")
                self.ISLM_widths[3] = float(ws[0].strip())
                self.ISLM_widths[5] = float(ws[1].strip())
                self.ISLM_widths[10] = float(ws[2].strip())
                self.ISLM_widths[15] = float(ws[3].strip())
                self.ISLM_widths[17] = float(ws[4].strip())
                print ("* WIDTHS", self.ISLM_widths)
            if 'Gullwing Global Shape Factor' in line: 
                wd = line.split("=")[1].strip()
                self.ISLM_GWShape = float(wd)
                print ("* Gullwing Global Shape Factor", self.ISLM_GWShape)
            if 'Gullwing Local Shape Factor' in line: 
                wd = line.split("=")[1].strip()
                self.ISLM_localGWShape = float(wd)
                print ("* Gullwing Local Shape Factor", self.ISLM_GWShape)
            if 'Roundness Shape Factor' in line: 
                wd = line.split("=")[1].strip()
                self.ISLM_roundnessShape = float(wd)
                print ("* Roundness Shape Factor", self.ISLM_GWShape)
            if 'Squareness Shape Factor' in line: 
                wd = line.split("=")[1].strip()
                self.ISLM_squarenessShape = float(wd)
                print ("* Squareness Shape Factor", self.ISLM_GWShape)

    def inputISLM_caliFPC(self, fname): 
        # print (" ISLM FPC Input")
        
        with open(fname) as F: 
            lines = F.readlines()
        for line in lines: 
            if 'ActualContactArea' in line: 
                wd = line.split("=")[1].strip()
                self.ISLM_cali_actualArea = float(wd)/10000
                print ('* ActualContactArea', self.ISLM_cali_actualArea)
            if 'TotalContactArea' in line: 
                wd = line.split("=")[1].strip()
                self.ISLM_cali_totalArea = float(wd)/10000
                print ('* TotalContactArea', self.ISLM_cali_totalArea)
            if 'Roundness(%)=' in line: 
                wd = line.split("=")[1].strip()
                self.ISLM_cali_roundness = float(wd)/100
                print ('* Roundness', self.ISLM_cali_roundness)
            if 'DetailedContactLength(mm)' in line: 
                wd = line.split("=")[1].strip()
                ws = wd.split("/")
                # foot.lengths[3], foot.lengths[5], foot.lengths[10], foot.lengths[15], foot.lengths[17]
                self.ISLM_cali_lengths[3] = float(ws[0].strip())
                self.ISLM_cali_lengths[5] = float(ws[1].strip())
                self.ISLM_cali_lengths[10] = float(ws[2].strip())
                self.ISLM_cali_lengths[15] = float(ws[3].strip())
                self.ISLM_cali_lengths[17] = float(ws[4].strip())
                print ("* LENGTHS", self.ISLM_cali_lengths)
            if 'DetailedContactWidth(mm)' in line: 
                wd = line.split("=")[1].strip()
                ws = wd.split("/")
                self.ISLM_cali_widths[3] = float(ws[0].strip())
                self.ISLM_cali_widths[5] = float(ws[1].strip())
                self.ISLM_cali_widths[10] = float(ws[2].strip())
                self.ISLM_cali_widths[15] = float(ws[3].strip())
                self.ISLM_cali_widths[17] = float(ws[4].strip())
                print ("* WIDTHS", self.ISLM_cali_widths)
            if 'Gullwing Global Shape Factor' in line: 
                wd = line.split("=")[1].strip()
                self.ISLM_cali_GWShape = float(wd)
                print ("* Gullwing Global Shape Factor", self.ISLM_cali_GWShape)
            if 'Gullwing Local Shape Factor' in line: 
                wd = line.split("=")[1].strip()
                self.ISLM_cali_localGWShape = float(wd)
                print ("* Gullwing Local Shape Factor", self.ISLM_cali_GWShape)
            if 'Roundness Shape Factor' in line: 
                wd = line.split("=")[1].strip()
                self.ISLM_cali_roundnessShape = float(wd)
                print ("* Roundness Shape Factor", self.ISLM_cali_GWShape)
            if 'Squareness Shape Factor' in line: 
                wd = line.split("=")[1].strip()
                self.ISLM_cali_squarenessShape = float(wd)
                print ("* Squareness Shape Factor", self.ISLM_cali_GWShape)


    def inputInitalFootprint(self, px, py, pv, actualArea): 
        self.px = px 
        self.py = py 
        self.pv = pv
        self.actualArea = actualArea 

    def inputCenterPressure(self, lst): 
        self.centerPress=lst 

    def inputFPC(self, areaPts, basic, adv, lines): 
        self.totalArea = basic[3] 
        self.actualArea = basic[2] 
        self.roundness = basic[4] 
        self.lengths = basic[0]
        self.widths = basic[1]
        self.GWShape = adv[0]
        self.localGWShape = adv[1]
        self.roundnessShape = adv[2] 
        self.squarenessShape = adv[3] 
        self.FPClines = lines 
        self.areaPoints=areaPts 
    def inputBoundary(self, xs, ys): 
        self.boundary = [xs, ys]
   
   
def readDAT(footfile): 
    with open(footfile) as DAT: 
        lines = DAT.readlines()
    shift = 0; offset=0 
    cmd = None 
    nodes=[]
    els =[]
    for line in lines: 
        if "*" in line: 
            if "*NODE" in line : 
                cmd = "ND"
                
            elif "*ELEMENT" in line : 
                cmd = 'EL'
                
            elif  "*LATERAL SHIFT" in line: 
                cmd = 'LS'

            elif "*OFFSET" in line: 
                cmd = 'OF'

            else: 
                cmd = None 
        else: 
            data = line.split(",")
            for d in data: 
                d = d.strip()

            if cmd == 'ND': 
                nodes.append([int(float(data[0])), float(data[1]), float(data[2]), float(data[3]), float(data[4]), float(data[5])])
            if cmd == "EL": 
                els.append( [int(data[0]), int(data[1]), int(data[2]), int(data[3]), int(data[4])] )
            if cmd == 'LS': 
                shift = float(data[0])
            if cmd == 'OF' : 
                offset = int(data[0])

    nodes = np.array(nodes)
  
    surfaces =[]
    for el in els: 
        ix = np.where(nodes[:,0]==el[1])[0]
        if len(ix):  n1 = nodes[ix[0]]
        else: continue 
        ix = np.where(nodes[:,0]==el[2])[0]
        if len(ix):  n2 = nodes[ix[0]]
        else: continue 
        ix = np.where(nodes[:,0]==el[3])[0]
        if len(ix):  n3 = nodes[ix[0]]
        else: continue 
        ix = np.where(nodes[:,0]==el[4])[0]
        if len(ix):  n4 = nodes[ix[0]]
        else: continue 

        surfaces.append([el[0], el[1], el[2], el[3], el[4], [n1, n2, n3, n4]])

    edge_Boundary = Edge_Boundary(surfaces)
    return nodes, surfaces, offset, shift,   edge_Boundary

@timer 
def interpolation_footprint_pressure(nodes=None, surfaces=None, vmin=0, vmax=100000, rotating=0.0, cmap='rainbow', displim=0.15,\
     dotNum=5, trd=10**7, savefile="", dpi=300, post=False,  Yshift=0.0, ptn=None, Onlyribshape=False, offset=0, layoutCenterNode=None): 

    minZ = np.min(nodes[:,3])
    idx = np.where(nodes[:,0] > 10**7)[0]
    tdnodes = nodes[idx]
    totalContactForce = np.sum(tdnodes[:,5])

    if offset == 0 and post == False : 
        node_ids = nodes[:,0]
        ix = np.where(node_ids>10**7)[0]
        td_ids = node_ids[ix]
        td_Offset = findout_offset(td_ids, step=10**4, shift=-10**7)
        print (" Tread Pitch Offset=%d"%(td_Offset))
    elif offset == 0 and post : 
        td_Offset = 10000
    else: 
        td_Offset = offset

    if not isfile('topsurf.tmp'): 
        print ("## Generating topsurf.tmp ## ")
        ptns =[]
        if not ptn: 
            jobdir = getcwd()
            ptns = glob.glob(jobdir+"/*.ptn")
        if len(ptns): 
            PTN = PATTERN(ptns[0], start_number=10**7)
        else: 
            ptn = None 
                
    topsurf = 'topsurf.tmp'
    if isfile(topsurf): 
        with open(topsurf) as TF: 
            lines = TF.readlines()
        tops =[]
        for line in lines: 
            wd = line.split(",")
            tops.append([int(wd[0].strip()), int(wd[1].strip()), int(wd[2].strip()), int(wd[3].strip()), int(wd[4].strip())])
        tops = np.array(tops)

    ht = 2.0e-3
    mht = ht*2.0
    gap = 0.3e-3

    cnt = 0 
    ActualArea = 0 

    values_onSurface=[]
    
    for sf in surfaces: 
        if isfile(topsurf):   
            if (sf[5][0][3] <= minZ + ht or sf[5][1][3] <= minZ + ht or sf[5][2][3] <= minZ + ht or sf[5][3][3] <= minZ + ht) and \
                (sf[5][0][3] <= minZ + mht and sf[5][1][3] <= minZ + mht and sf[5][2][3] <= minZ + mht and sf[5][3][3] <= minZ + mht) : 

                ix1 = np.where(tops[:, 1:] - trd == int(sf[5][0][0] - trd)%td_Offset )[0]
                ix2 = np.where(tops[:, 1:] - trd == int(sf[5][1][0] - trd)%td_Offset )[0]
                ix3 = np.where(tops[:, 1:] - trd == int(sf[5][2][0] - trd)%td_Offset )[0]
                ix4 = np.where(tops[:, 1:] - trd == int(sf[5][3][0] - trd)%td_Offset )[0]
                ix = np.intersect1d(ix1, ix2)
                ix = np.intersect1d(ix,  ix3)
                ix = np.intersect1d(ix,  ix4)
                if len(ix):
                    contacting = True 
                else: 
                    gap1 = abs(sf[5][0][3] - sf[5][1][3])
                    gap2 = abs(sf[5][1][3] - sf[5][2][3])
                    gap3 = abs(sf[5][2][3] - sf[5][3][3])
                    gap4 = abs(sf[5][3][3] - sf[5][1][3])

                    if gap1<gap and gap2<gap and gap3<gap and gap4<gap: 
                        contacting = True 
                    else: 
                        contacting = False 
            else: 
                contacting = False 
        else: 
            contacting = True  
        if contacting : 
            ys = [sf[5][0][1], sf[5][1][1], sf[5][2][1], sf[5][3][1]]
            xs = [sf[5][0][2], sf[5][1][2], sf[5][2][2], sf[5][3][2]]
            vs = [sf[5][0][4], sf[5][1][4], sf[5][2][4], sf[5][3][4]]
            values_onSurface.append([xs, ys, vs])
            ActualArea += element_contactArea(xs, ys, vs, vmin)
            if cnt ==0: 
                px, py, pv = meshgrid_in_Quadrilateral(xs, ys, num=dotNum, vs=vs)
            else: 
                ipx, ipy, ipv = meshgrid_in_Quadrilateral(xs, ys, num=dotNum, vs=vs)
                px = np.append(px, ipx, axis=0)
                py = np.append(py, ipy, axis=0)
                pv = np.append(pv, ipv, axis=0)

            cnt += 1

    print ("  Surfaces=%d(%d), Contact Ht=%f"%(cnt, len(surfaces), minZ*1000))

    if rotating: 
        print (" Points are rotated by the angle %.3f"%(rotating))
        angle = math.radians(rotating)
        for x, y in zip(px, py): 
            x = math.cos(angle)*x - math.sin(angle)*y 
            y = math.sin(angle)*x + math.cos(angle)*y 

    ###############################################################
    ## Centering in X direction 
    ###############################################################
    # print ("YSHIFT", Yshift)
    Centering = True 
    if Centering: 
        if not Yshift:
            shift = 0.0 
            if  isfile('shift.tmp'): 
                    with open('shift.tmp') as SH: 
                        line = SH.readlines()
                    shift = float(line[0].strip() )
        
        if not isinstance(layoutCenterNode, type(None)):
            ixt = np.where(nodes[:,0] > 10**7)[0]
            ixt1 = np.where(nodes[:,0] < 10**7+ 10000)[0]
            ixt = np.intersect1d(ixt, ixt1)
            ixb = np.where(nodes[:,0] < 10000)[0]
            if len(ixt) > len(ixb): 
                shift = Yshift 
            else: 
                ixs = np.where(nodes[:,0] % td_Offset == layoutCenterNode[0])[0]
                xmin = 10.0
                dfnode = None 
                # print ("   Contact Ht=%.2f"%(minZ*1000))
                for ix in ixs: 
                    if abs(nodes[ix][1])<xmin and nodes[ix][3] < minZ + 0.0001: 
                        xmin = abs(nodes[ix][1])
                        dfnode = nodes[ix]
                        # print ("   Deformed", dfnode[0]-10**7, dfnode[2]*1000, dfnode[3]*1000)
                if not isinstance(dfnode, type(None)): 
                    shift = layoutCenterNode[2] - dfnode[2] 
                else: 
                    shift = Yshift
        else: 
            shift = Yshift 
    else: 
        shift = 0 
    px = px.flatten(); py.flatten(); pv.flatten()
    print ("* Footprint Lateral shift =%.2fmm"%(shift*1000) )
    px += shift 
    # px.flatten(); py.flatten(); pv.flatten()
    return px.flatten(), py.flatten(), pv.flatten(), ActualArea, shift, values_onSurface 

def element_contactArea(xs, ys, vs, vmin): 
    counting = 0 
    for m in range(4): 
        if vs[m] > vmin: 
            counting += 1 
    if counting >0  : 
        if xs[2] == xs[3] and ys[2] == ys[3] : 
            ax = xs[:-1]
            ay = ys[:-1]
            shape = 3 
            if vs[3] > vmin:   counting -= 1 
        else: 
            ax = xs
            ay = ys
            shape = 4 
        return abs(Area(ax, ay)) * counting / shape 
    else: 
        return 0.0 

def calculatingActualArea(surfs, vmin): 
    actualArea = 0 

    for sf in surfs: 
        actualArea += element_contactArea(sf[0], sf[1], sf[2], vmin)
 
    return actualArea

def footPrintWith_belt(px, py, pv, vmin=5E4, cmap='rainbow', displim=0.15, savefile="", dpi=300, nodes=None, outer=None,\
     belt=None, bead=None, carcass=None, initnode=None, shift=0): 
    size =0.3
    fig, ax = plt.subplots()
    plt.axis('equal')
    plt.xlabel("Lateral Position(m)", size=11)
    plt.ylabel("Longitudinal Position(m)", size=11)
    plt.xticks(size=10)
    plt.yticks(size=10)
    plt.title("Contact Pressure(Pa) Distribution", size=14)
    
    points = plt.scatter(px, py, c=pv, s=size, cmap=cmap, vmin=vmin, vmax=vmin*10, edgecolors=None, linewidths=0.0 )
    
    #######################################################################
    ## saving profile 
    #######################################################################
    
    fp = open(savefile+"-outerprofile.inp", 'w')
    fp.write("*NODE\n")
    for n in initnode:
        fp.write("%10d, %.7f, %.7f, %.7f\n"%(n[0], n[3], n[2], n[1]))
    iniMinZ = np.min(initnode[:, 3])
    dfmMinZ = np.max(nodes[:,3])
    zmove = iniMinZ + dfmMinZ
    ofn = 5000
    for n in nodes:
        z = math.sqrt((n[3]**2+n[1]**2))
        fp.write("%10d, %.7f, %.7f, %.7f\n"%(n[0]+ofn, z+zmove, n[2], 0.0))
    
    fp.write("*ELEMENT, TYPE=MGAX1\n")
    # cnt = 0 
    for ed in outer: 
        # cnt += 1
        # print ("out %6d :  %6d, %6d, %6d"%(cnt, ed[4], ed[0], ed[1]))
        fp.write("%10d, %10d, %10d\n"%(ed[4], ed[0], ed[1]))
        fp.write("%10d, %10d, %10d\n"%(ed[4]+ofn, ed[0]+ofn, ed[1]+ofn))
    # cnt = 0 
    for ed in belt:
        # cnt += 1
        # print ("bt %6d :  %6d, %6d, %6d"%(cnt, ed[4], ed[0], ed[1]))
        fp.write("%10d, %10d, %10d\n"%(ed[4], ed[0], ed[1]))
        fp.write("%10d, %10d, %10d\n"%(ed[4]+ofn, ed[0]+ofn, ed[1]+ofn))
    # cnt = 0 
    for ed in bead:
        # cnt += 1
        # print ("bd %6d :  %6d, %6d, %6d"%(cnt, ed[4], ed[0], ed[1]))
        fp.write("%10d, %10d, %10d\n"%(ed[4], ed[0], ed[1]))
        fp.write("%10d, %10d, %10d\n"%(ed[4]+ofn, ed[0]+ofn, ed[1]+ofn))
    # cnt = 0 
    for ed in carcass:
        # cnt += 1
        # print ("cc %6d :  %6d, %6d, %6d"%(cnt, ed[4], ed[0], ed[1]))
        fp.write("%10d, %10d, %10d\n"%(ed[4], ed[0], ed[1]))
        fp.write("%10d, %10d, %10d\n"%(ed[4]+ofn, ed[0]+ofn, ed[1]+ofn))
    fp.write("*ELSET, ELSET=initial\n")
    cnt = 0 
    for ed in outer: 
        cnt += 1 
        if cnt%10 : fp.write("%10d, "%(ed[4]))
        else:  fp.write("%10d\n"%(ed[4]))
    if cnt%10: 
        fp.write("\n")
    fp.write("*ELSET, ELSET=deformed\n")
    cnt = 0 
    for ed in outer: 
        cnt += 1 
        if cnt%10 : fp.write("%10d, "%(ed[4]+ofn))
        else:  fp.write("%10d\n"%(ed[4]+ofn))
    if cnt%10: 
        fp.write("\n")
    fp.write("*ELSET, ELSET=BT_initial\n")
    cnt = 0 
    for ed in belt: 
        cnt += 1 
        if cnt%10 : fp.write("%10d, "%(ed[4]))
        else:  fp.write("%10d\n"%(ed[4]))
    if cnt%10: 
        fp.write("\n")
    fp.write("*ELSET, ELSET=BT_deformed\n")
    cnt = 0 
    for ed in belt: 
        cnt += 1 
        if cnt%10 : fp.write("%10d, "%(ed[4]+ofn))
        else:  fp.write("%10d\n"%(ed[4]+ofn))
    if cnt%10: 
        fp.write("\n")

    fp.write("*ELSET, ELSET=BD_initial\n")
    cnt = 0 
    for ed in bead: 
        cnt += 1 
        if cnt%10 : fp.write("%10d, "%(ed[4]))
        else:  fp.write("%10d\n"%(ed[4]))
    if cnt%10: 
        fp.write("\n")
    fp.write("*ELSET, ELSET=BD_deformed\n")
    cnt = 0 
    for ed in bead: 
        cnt += 1 
        if cnt%10 : fp.write("%10d, "%(ed[4]+ofn))
        else:  fp.write("%10d\n"%(ed[4]+ofn))
    if cnt%10: 
        fp.write("\n")

    fp.write("*ELSET, ELSET=CC_initial\n")
    cnt = 0 
    for ed in bead: 
        cnt += 1 
        if cnt%10 : fp.write("%10d, "%(ed[4]))
        else:  fp.write("%10d\n"%(ed[4]))
    if cnt%10: 
        fp.write("\n")
    fp.write("*ELSET, ELSET=CC_deformed\n")
    cnt = 0 
    for ed in bead: 
        cnt += 1 
        if cnt%10 : fp.write("%10d, "%(ed[4]+ofn))
        else:  fp.write("%10d\n"%(ed[4]+ofn))
    if cnt%10: 
        fp.write("\n")

    fp.close()
    ###########################################################

    dmax = np.max(py)
    
    lw = 0.2
    
    xs=[]; ys=[]
    
    for ed in outer: 
        ix1=np.where(nodes[:, 0] == ed[0])[0][0]
        ix2=np.where(nodes[:, 0] == ed[1])[0][0]
        xs.append([nodes[ix1][2]+shift, nodes[ix2][2]+shift])
        ys.append([ -math.sqrt (nodes[ix1][1]**2+nodes[ix1][3]**2), -math.sqrt (nodes[ix2][1]**2+nodes[ix2][3]**2)])

    xs = np.array(xs); ys = np.array(ys)
    lmin = np.min(ys)
    lmax = np.max(ys)
    lift = dmax - lmin + 0.05
    ys += lift 
    for tpx, tpy in zip(xs, ys): 
        plt.plot(tpx, tpy, 'black', lw=lw)

    nymax = np.max(ys)

    layout_Y = []
    xs=[]; ys=[]
    for ed in outer: 
        ix1=np.where(initnode[:, 0] == ed[0])[0][0]
        ix2=np.where(initnode[:, 0] == ed[1])[0][0]
        xs.append([  initnode[ix1][2]+shift,  initnode[ix2][2]+shift])
        ys.append([ -initnode[ix1][3], -initnode[ix2][3]])

    xs = np.array(xs); ys = np.array(ys)
    lmax = np.max(ys)
    lift_init = nymax - lmax
    ys += lift_init 
    for tpx, tpy in zip(xs, ys): 
        plt.plot(tpx, tpy, 'lightgray', lw=lw)

    bmin=0; bmax=0 
    neg=[]; pos=[]
    
    for ed in belt: 
        ix1=np.where(nodes[:, 0] == ed[0])[0][0]
        ix2=np.where(nodes[:, 0] == ed[1])[0][0]
        xs=[nodes[ix1][2]+shift, nodes[ix2][2]+shift]
        ys=[ -math.sqrt (nodes[ix1][1]**2+nodes[ix1][3]**2) +lift,-math.sqrt (nodes[ix2][1]**2+nodes[ix2][3]**2) +lift]
        plt.plot(xs, ys, 'red', lw=lw)
        if bmin > xs[0]:  
            neg = [xs[0], ys[0]]
            bmin = xs[0]
        if bmin > xs[1]:  
            neg = [xs[1], ys[1]]
            bmin = xs[1]
        if bmax < xs[0]:  
            pos = [xs[0], ys[0]]
            bmax = xs[0]
        if bmax < xs[1]:  
            pos = [xs[1], ys[1]]
            bmax = xs[1]

    displim *=2
    
    plt.plot([neg[0], neg[0]], [neg[1], 0], 'b:', lw=lw)
    plt.plot([pos[0], pos[0]], [pos[1], 0], 'b:', lw=lw)
    BeltMaxWidth = bmax - bmin 
    plt.text(0, nymax+0.01, "Max BT Width(deformed): %.1fmm"%(BeltMaxWidth*1000), fontsize=6, color='black', ha='center')

    ix1 = np.where(py> -0.005)[0]
    ix2 = np.where(py<  0.005)[0]
    ix = np.intersect1d(ix1, ix2)
    mid_py = py[ix];            mid_px = px[ix]
    mid_mn = np.min(mid_px);    mid_mx = np.max(mid_px)
    Max_ContactWidth = mid_mx-mid_mn 
    dmin = np.min(py)

    plt.text(0, dmin-0.04,   "Contact Width: %.1fmm"%(Max_ContactWidth*1000), fontsize=6, color='black', ha='center')

    plt.ylim(dmin-0.08, nymax*1.2)
    plt.xlim(-displim*2, displim*2)
    plt.text(-displim/2-displim*0.1, nymax*1.1, 'Original Image', fontsize=10, color='r')
    # print (" display Y lim max", nymax*1.2)
    # divider = make_axes_locatable(ax)
    # cax = divider.append_axes("right", size="1%", pad=0.1)
    # cbar = plt.colorbar(points, cax)
    cbar = plt.colorbar(points)
    cbar.ax.tick_params(labelsize=5)

    plt.savefig(savefile+"-Footshape_Original_Belt.png", dpi=dpi)
    plt.clf()

def contPress(px, py, pv, vmin=50000,  displim=0.15, savefile=None, dpi=300):
    yRange=0.001
    ix1 = np.where(py>-yRange)[0]
    ix2 = np.where(py< yRange)[0]
    ix = np.intersect1d(ix1, ix2)
    dx=px[ix]; dv=pv[ix]

    pos = -displim
    step = 1e-3 
    cp=[]
    
    while pos <=displim: 
        ix1 = np.where(dx>=pos     )[0]
        ix2 = np.where(dx<=pos+step)[0]
        ix = np.intersect1d(ix1, ix2) 
        if len(ix): 
            tv = dv[ix]
            vmax = np.max(tv)/10**6  ## unit : MPa
            cp.append([pos, vmax])
        else: 
            cp.append([pos, 0.0])
        pos+= step
    cp = np.array(cp)

    if savefile: 
        plt.plot(cp[:,0], cp[:,1], c='blue', lw=1.0)
        pMax = np.max(cp[:,1])
        ix = np.where(cp[:,1]==pMax)[0][0]
        maxValue = cp[ix]
        txt = "Max: %.3fMPa at %.3f(m)"%(maxValue[1], maxValue[0])

        minX = 0.0 
        txtPosition = [minX, maxValue[1]*1.05]
        if maxValue[0] < 0:     txtAlign = 'left'
        else: txtAlign = 'right'
        plt.annotate(txt, xy=(maxValue[0], maxValue[1]), xytext=(txtPosition[0], txtPosition[1]), \
            fontsize=10, arrowprops=dict(facecolor='black', shrink=0.05, width=1.0, headwidth=5), ha=txtAlign)

        plt.title("Contact Pressure Distribution", size=14)
        plt.ylabel("Contact Pressure(MPa)", size=10)
        plt.xlabel("Lateral Position(m)", size=10)
        plt.xticks(size=9)
        plt.yticks(size=9)

        plt.grid(True)
        plt.rc('grid', linestyle="--", color='black', linewidth=0.2)
        plt.xlim(-displim, displim)
        plt.ylim(0, maxValue[1]*1.1)
        plt.savefig(savefile+"-CpressGraph.png",dpi=dpi)
        plt.clf()

    fp = open("CpressAlongCenter.txt", "w")
    for p in cp: 
        fp.write("'%3.7f\t%3.7f\t%3.7f\n"%(0.0, p[0], p[1]*10**6))
    fp.close()
    return cp 

def Cropping_footprint_dots(ptx=None, pty=None, ptv=None, xrange=None, yrange=None): 

    if xrange: 
        ix1 = np.where(ptx>=xrange[0])[0]
        ix2 = np.where(ptx<=xrange[1])[0]
        ix  = np.intersect1d(ix1, ix2) 
    if yrange: 
        iy1 = np.where(pty>=yrange[0])[0]
        iy2 = np.where(pty<=yrange[1])[0]
        iy  = np.intersect1d(iy1, iy2) 
    if xrange and yrange: 
        idx = np.intersect1d(ix, iy)
    elif xrange: 
        idx = ix 
    elif yrange: 
        idx = iy 
    else: 
        return ptx, pty, ptv 
    return ptx[idx], pty[idx], ptv[idx]

def contactLength(PointsUp, PointsDown, maxes, position=50, half=False): 
    
    maxwidth = maxes[0]-maxes[1]
    maxlength = maxes[2]-maxes[3]
    halfWidth = maxwidth/2.0 

    mid = (maxes[0]+maxes[1])/2.0 
    cloeset = maxes[0] 

    uppoint  =0.0

    positions=[position]
    for lpos in positions: 
        specificposition   = -halfWidth + lpos *  maxwidth / 100 
        sclose = 10**7
        spoint1 = 0.0
        for pt in PointsUp:
            if abs(mid - pt[0]) < cloeset : 
                cloeset = abs(mid-pt[0])
                uppoint = pt[1]
            if abs(specificposition-pt[0]) < sclose: 
                sclose = abs(specificposition-pt[0])
                spoint1 = pt[1]

        cloeset = maxes[0] 
        sclose = 10**7
        downpoint  =0.0
        spoint2 = 0.0
        for pt in PointsDown:
            if abs(mid - pt[0]) < cloeset : 
                cloeset = abs(mid-pt[0])
                downpoint = pt[1]
            if abs(specificposition-pt[0]) < sclose: 
                sclose = abs(specificposition-pt[0])
                spoint2 = pt[1]

        centerlength = uppoint - downpoint 
        slength = spoint1 - spoint2
        lx = specificposition
        ly1 = spoint1
        ly2 = spoint2
        if half == 'up': 
            return spoint1
        elif half == 'down': 
            return -spoint2
        else: 
            return slength 

def contactWidth(PointsPos, PointsNeg, maxes, position=50): 
    maxwidth = maxes[0]-maxes[1]
    maxlength = maxes[2]-maxes[3]

    mid = (maxes[0]+maxes[1])/2.0 
    cloeset = maxes[0] 

    pointsX2 = np.append( PointsPos[:, 0], PointsNeg[:, 0])
    pointsY2 = np.append( PointsPos[:, 1], PointsNeg[:, 1])

    uppoint  =0.0
    positions=[position]
    for wpos in positions:
        mid = (maxes[2]+maxes[3])/2.0 
        specificposition = maxes[2] - maxlength * wpos / 100.0
        sclose = maxes[2]
        spoint1 = 0.0
        for pt in PointsPos:
            if abs(specificposition-pt[1]) < sclose: 
                sclose = abs(specificposition-pt[1])
                spoint1 = pt[0]

        sclose = maxes[2]
        spoint2 = 0.0
        for pt in PointsNeg:
            if abs(specificposition-pt[1]) < sclose: 
                sclose = abs(specificposition-pt[1])
                spoint2 = pt[0]

        return spoint1 - spoint2
    
def FPC(px, py, pv, cp=None, savefile="", displim=0.15, fitting=6,\
     cmap='rainbow', vmin=50000, vmax=100000, dpi=300, ActualArea=0.0,\
          ProductLine='PCR', contactForce=0, comparing=False): 

    yRange=0.001
    ix1 = np.where(py>-yRange)[0]
    ix2 = np.where(py< yRange)[0]
    ix = np.intersect1d(ix1, ix2)
    dx=px[ix]; dv=pv[ix]

    xMinCnt = np.min(dx)
    xMaxCnt = np.max(dx)
    Center_ContactWidth = xMaxCnt - xMinCnt 

    xMin = np.min(px)
    xMax = np.max(px)
    Max_ContactWidth = xMax - xMin 

    ix = np.where(py>0)[0]
    upX = px[ix]; upY = py[ix]; upP = pv[ix]
    ix = np.where(py<0)[0]
    downX = px[ix]; downY = py[ix]; downP = pv[ix]

    maxY = np.max(upY)
    minY = np.min(downY)
    Max_ContactLength = maxY - minY 

    xRange = 0.5e-3   
    step = 0.5e-3
    while xRange < 30E-3: 
        ix1 = np.where(px>-xRange)[0]
        ix2 = np.where(px< xRange)[0]
        ix = np.intersect1d(ix1, ix2)
        if len(ix) > 1000: 
            cLine = py[ix]
            yMax = np.max(cLine)
            yMin = np.min(cLine)
            Center_ContactLength = yMax - yMin 
            break 
        xRange += step 

    fittingmargin = 0.001 
    fittingmargin_gap=fittingmargin/100.0
    point_gap = 3.01e-3
    
    fpcX, fpcY, fpcV = Cropping_footprint_dots(ptx=px, pty=py, ptv=pv, xrange=[xMinCnt*1.2, xMaxCnt*1.2])
    if comparing: return fpcX, fpcY

    PointsUp, PointsDown, PointsPos, PointsNeg, maxes =\
         SearchPoints(fpcX, fpcY, fitting=fitting, fittingmargin=fittingmargin, dist=point_gap, savefig=False, delete='max')

    #############################################
    gapLimit = 0.0025
    upMax = np.max (PointsUp[:,1])
    downMin = np.min(PointsDown[:,1])
    if upMax > maxY + gapLimit :   
        ixs = np.where(fpcY<maxY - abs(upMax - maxY)/2)[0]
        ixx = np.where(fpcY>minY + abs(downMin - minY)/2)[0]
        ix = np.intersect1d(ixs, ixx)
        fpcX = fpcX[ix]
        fpcY = fpcY[ix]
        fpcV = fpcV[ix]

        PointsUp, PointsDown, PointsPos, PointsNeg, maxes =\
        SearchPoints(fpcX, fpcY, fitting=fitting, fittingmargin=fittingmargin, dist=point_gap, savefig=False, delete='max')

    ##############################################

    pointsX2 = np.append( PointsUp[:, 0], PointsDown[:, 0])
    pointsY2 = np.append( PointsUp[:, 1], PointsDown[:, 1])
            
    if abs(xMaxCnt) < abs(xMinCnt): 
        maxes[0] = xMaxCnt 
        maxes[1] = -xMaxCnt 
        halfWidth = abs(xMaxCnt)
    else: 
        maxes[0] = -xMinCnt 
        maxes[1] = xMinCnt 
        halfWidth = abs(xMinCnt)


    maxes[2] = maxY
    maxes[3] = minY 

    maxwidth = maxes[0]-maxes[1]
    maxlength = maxes[2]-maxes[3]

    print ("* Max Length=%.1f, Width=%.1f"%(maxlength*1000, maxwidth*1000))
    
    mid = (maxes[0]+maxes[1])/2.0 
    cloeset = maxes[0] 

    pointsX1 = np.append( PointsUp[:, 0], PointsDown[:, 0])
    pointsY1 = np.append( PointsUp[:, 1], PointsDown[:, 1])

    pointsX2 = np.append( PointsPos[:, 0], PointsNeg[:, 0])
    pointsY2 = np.append( PointsPos[:, 1], PointsNeg[:, 1])

    uppoint  =0.0

    # positions=[15, 25, 50, 75, 85, 5, 10, 20, 30, 35, 40, 45, 55, 60, 65, 70, 80, 90, 95, 0, 100]
    positions = range(0, 101, 5)
    Lengths=[]
    lines = []
    poscounting = 0 
    maxLength = 0 
    cals = []
    for lpos in positions: 
        poscounting += 1 
        specificposition   = -halfWidth + lpos *  maxwidth / 100 
        sclose = 10**7
        spoint1 = 0.0
        for pt in PointsUp:
            if abs(mid - pt[0]) < cloeset : 
                cloeset = abs(mid-pt[0])
                uppoint = pt[1]
            if abs(specificposition-pt[0]) < sclose: 
                sclose = abs(specificposition-pt[0])
                spoint1 = pt[1]

        cloeset = maxes[0] 
        sclose = 10**7
        downpoint  =0.0
        spoint2 = 0.0
        for pt in PointsDown:
            if abs(mid - pt[0]) < cloeset : 
                cloeset = abs(mid-pt[0])
                downpoint = pt[1]
            if abs(specificposition-pt[0]) < sclose: 
                sclose = abs(specificposition-pt[0])
                spoint2 = pt[1]

        centerlength = uppoint - downpoint 
        slength = spoint1 - spoint2
        lx = specificposition
        ly1 = spoint1
        ly2 = spoint2
        # if not lpos % 5 : 
        # if poscounting <=5: 
        if True: 
            lengthX=[]; lengthY=[]
            lengthX.append(lx); lengthY.append(ly1)
            lengthX.append(lx); lengthY.append(ly2)
            # plt.plot(lengthX, lengthY, color='darkgreen', lw=0.5)
            # plt.text(lx, ly2-0.01, "L%.1f"%(slength*1000), size=5)

            lines.append([lengthX, lengthY, [lx, ly2-0.01, "L%.1f"%(slength*1000)] ])
            
            Lengths.append(slength*1000)
        if maxLength < slength: 
            maxLength = slength 
        # print ("    POS %.1f%%, length=%.2f"%(lpos, slength*1000))
        cals.append([lpos, slength*1000])
        
        ## if contour == True: 
        mid = (maxes[2]+maxes[3])/2.0 
        specificposition = maxes[2] - maxlength * 0.5
        sclose = maxes[2]
        spoint1 = 0.0
        for pt in PointsPos:
            if abs(specificposition-pt[1]) < sclose: 
                sclose = abs(specificposition-pt[1])
                spoint1 = pt[0]

        sclose = maxes[2]
        spoint2 = 0.0
        for pt in PointsNeg:
            if abs(specificposition-pt[1]) < sclose: 
                sclose = abs(specificposition-pt[1])
                spoint2 = pt[0]


        centerwidth = spoint1 - spoint2
    # txt = ''
    # for i, cal in enumerate(cals): 
    #     if (i+1) % 10 ==0: txt += '\n'
    #     txt += "%.1f%%, %.2fmm, "%(cal[0], cal[1])
    # print (txt)



    cals =[]
    print (" >> Max Length=%.1f"%(maxLength*1000))
    maxlength = maxLength 
    ### width at specific position 
    Widths=[]
    posPositions=[]; negPositions=[]
    poscounting = 0 
    maxWidth = 0 
    for wpos in positions:
        poscounting += 1
        mid = (maxes[2]+maxes[3])/2.0 
        specificposition = maxes[2] - maxlength * wpos / 100.0
        sclose = maxes[2]
        spoint1 = 0.0
        for pt in PointsPos:
            if abs(specificposition-pt[1]) < sclose: 
                sclose = abs(specificposition-pt[1])
                spoint1 = pt[0]

        sclose = maxes[2]
        spoint2 = 0.0
        for pt in PointsNeg:
            if abs(specificposition-pt[1]) < sclose: 
                sclose = abs(specificposition-pt[1])
                spoint2 = pt[0]

        swidth = spoint1 - spoint2

        posPositions.append([spoint1, specificposition])
        negPositions.append([spoint2, specificposition])
        
        cals.append([wpos, swidth*1000])
    # print (" CHECK Width Information", len(posPositions))
    print ("** Check contact width on shoulder Lug if any")
    cnt = 0     
    idx = 0
    for pw, nw in zip(posPositions, negPositions): 
        if cnt ==0 : 
            d1 = abs(posPositions[cnt+1][idx]- pw[idx])
            d2 = abs(posPositions[cnt+2][idx]- pw[idx])
            d3 = abs(posPositions[cnt+3][idx]- pw[idx])

            if d2 < d1 and d2 < d3:  
                posPositions[cnt+2][idx] = (posPositions[cnt+1][idx]+posPositions[cnt+3][idx])/2
                print ("  Width=%6.2f, Ht=%6.2f, d=%6.1f,%6.1f,%6.1f"%(posPositions[cnt+2][0]*1000, posPositions[cnt+2][1]*1000, d1*1000, d2*1000, d3*1000))
                print ("  > width =%.2f"%( (posPositions[cnt+1][idx]+posPositions[cnt+3][idx])*500.0))

            d1 = abs(negPositions[cnt+1][idx]- pw[idx])
            d2 = abs(negPositions[cnt+2][idx]- pw[idx])
            d3 = abs(negPositions[cnt+3][idx]- pw[idx])

            if d2 < d1 and d2 < d3:  
                negPositions[cnt+2][idx] = (negPositions[cnt+1][idx]+negPositions[cnt+3][idx])/2
                print ("  Width=%6.2f, Ht=%6.2f, d=%6.1f,%6.1f,%6.1f"%(negPositions[cnt+2][0]*1000, negPositions[cnt+2][1]*1000, d1*1000, d2*1000, d3*1000))
                print ("  > width =%.2f"%( (negPositions[cnt+1][idx]+negPositions[cnt+3][idx])*500.0))

            cnt += 1
            continue 
        elif cnt >= len(posPositions)-1: 
            d1 = abs(posPositions[cnt-1][idx]- pw[idx])
            d2 = abs(posPositions[cnt-2][idx]- pw[idx])
            d3 = abs(posPositions[cnt-3][idx]- pw[idx])

            if d2 < d1 and d2 < d3:  
                posPositions[cnt-2][idx] = (posPositions[cnt-1][idx]+posPositions[cnt-3][idx])/2
                print ("  Width=%6.2f, Ht=%6.2f, d=%6.1f,%6.1f,%6.1f"%(posPositions[cnt-2][0]*1000, posPositions[cnt-2][1]*1000, d1*1000, d2*1000, d3*1000))
                print ("  > width =%.2f"%( (posPositions[cnt-1][idx]+posPositions[cnt-3][idx])*500.0))

            d1 = abs(negPositions[cnt-1][idx]- pw[idx])
            d2 = abs(negPositions[cnt-2][idx]- pw[idx])
            d3 = abs(negPositions[cnt-3][idx]- pw[idx])

            if d2 < d1 and d2 < d3:  
                negPositions[cnt-2][idx] = (negPositions[cnt-1][idx]+negPositions[cnt-3][idx])/2
                print ("  Width=%6.2f, Ht=%6.2f, d=%6.1f,%6.1f,%6.1f"%(negPositions[cnt-2][0]*1000, negPositions[cnt-2][1]*1000, d1*1000, d2*1000, d3*1000))
                print ("  > width =%.2f"%( (negPositions[cnt-1][idx]+negPositions[cnt-3][idx])*500.0))

            break 

        d1 = abs(posPositions[cnt-1][idx]- pw[idx])
        d2 = abs(posPositions[cnt+1][idx]- pw[idx])
        d3 = abs(posPositions[cnt+1][idx]- posPositions[cnt-1][idx])
        
        if d3 < d1 and d3 < d2 : 
            posPositions[cnt][idx] = (posPositions[cnt+1][idx]+ posPositions[cnt-1][idx])/2.0
            print ("  Width=%6.2f, Ht=%6.2f, d=%6.1f,%6.1f,%6.1f"%(pw[0]*1000, pw[1]*1000, d1*1000, d2*1000, d3*1000))
            print ("  > width =%.2f"%((posPositions[cnt+1][idx]+ posPositions[cnt-1][idx])*500.0))

        d1 = abs(negPositions[cnt-1][idx]- nw[idx])
        d2 = abs(negPositions[cnt+1][idx]- nw[idx])
        d3 = abs(negPositions[cnt+1][idx]- negPositions[cnt-1][idx])
        
        if d3 < d1 and d3 < d2 : 
            negPositions[cnt][idx] = (negPositions[cnt+1][idx]+ negPositions[cnt-1][idx])/2.0
            print ("  Width=%6.2f, Ht=%6.2f, d=%6.1f,%6.1f,%6.1f"%(nw[0]*1000, nw[1]*1000, d1*1000, d2*1000, d3*1000))
            print ("  > width =%.2f"%((negPositions[cnt+1][idx]+ negPositions[cnt-1][idx])*500.0))

        cnt += 1

    maxWidth = 0  
    
    for pw, nw in zip(posPositions, negPositions):
        lengthX=[]; lengthY=[]
        lengthX.append(pw[0]); lengthY.append(pw[1])
        lengthX.append(nw[0]); lengthY.append(nw[1])
        sw = pw[0]-nw[0]
        lines.append([lengthX, lengthY, [pw[0]+0.01, pw[1], "W%.1f"%(sw*1000)] ])
        if maxWidth < sw: 
            maxWidth = sw 
        
        Widths.append(sw*1000)


    print (" >> Max Width=%.1f"%(maxWidth*1000))
    Max_ContactWidth = maxWidth
    updown_points_Calculation_margin = 0.03
    downnode = NODE()
    cnt =1
    for nd in PointsDown:
        # print ("%.2f < %.2f < %.2f"%((maxes[1] + maxwidth * 0.06)*1000, nd[0]*1000, ( maxes[0] - maxwidth * 0.06)*1000))
        # if nd[0] >  maxes[1] + maxwidth * updown_points_Calculation_margin and nd[0] <  maxes[0] - maxwidth * updown_points_Calculation_margin:
            downnode.Add([cnt, 0.0, nd[0], nd[1]])
            cnt+=1
    downnode.Sort(item=2, reverse=False)

    upnode = NODE()
    for nd in PointsUp:
        # if nd[0] >  maxes[1] + maxwidth * updown_points_Calculation_margin and nd[0] <  maxes[0] - maxwidth * updown_points_Calculation_margin:
            upnode.Add([cnt, 0.0, nd[0], nd[1]])
            cnt+=1
    upnode.Sort(item=2, reverse=True)

    


    posnode = NODE()
    cnt = 0 
    for pn in posPositions:
        posnode.Add([cnt, 0.0, pn[0], pn[1]])
        cnt+=1

    posnode.Sort(item=3, reverse=False)
    negnodeUp = NODE()
    negnodeDown = NODE()
    cnt = 100 
    for pn in negPositions:
        if pn[1]>0: 
            negnodeUp.Add([cnt, 0.0, pn[0], pn[1]])
        else: 
            negnodeDown.Add([cnt, 0.0, pn[0], pn[1]])
        cnt+=1
    negnodeUp.Sort(item=3, reverse=True)
    negnodeDown.Sort(item=3, reverse=True)


    areanode = NODE()
    areanode.Combine(downnode)
    areanode.Combine(posnode)
    areanode.Node = Delete_Close_Points(areanode.Node, point_gap, backward=False )
    areanode.Combine(upnode)
    areanode.Node = Delete_Close_Points(areanode.Node, point_gap, backward=True )
    areanode.Combine(negnodeUp)
    areanode.Node = Delete_Close_Points(areanode.Node, point_gap, backward=False )
    areanode.Combine(negnodeDown)
    areanode.Node = Delete_Close_Points(areanode.Node, point_gap, backward=True )

    # npx=[]; npy=[]

    i = 1 
    PL=0 
    while i < len(areanode.Node)-1:
        L1 = math.sqrt((areanode.Node[i][2] - areanode.Node[i-1][2])**2  + (areanode.Node[i][3] - areanode.Node[i-1][3])**2)
        L2 = math.sqrt((areanode.Node[i][2] - areanode.Node[i+1][2])**2  + (areanode.Node[i][3] - areanode.Node[i+1][3])**2)
        
        if i>2 and L1 > PL * 2 and L1 > point_gap*2 and L2 > point_gap*2: 
            # print ("** del", areanode.Node[i])
            # npx.append(areanode.Node[i][2]); npy.append(areanode.Node[i][3]); 
            del(areanode.Node[i]) 
            # continue 
        PL = L1 
        i += 1

    
    
    upx = []; upy=[]
    dwx = []; dwy=[]
    for nd in areanode.Node: 
        if nd[3] >= 0: 
            upx.append(nd[2])
            upy.append(nd[3])
        elif nd[3] < 0: 
            dwx.append(nd[2])
            dwy.append(nd[3])
    deletedx=[]; deletedy=[]
    unsortedNode=[]
    
    i = 1
    while i < len(upx): 
        cangle = Angle_Between_Vectors([0, 0, 1, 0], [0, 0, upx[i], upy[i]])
        
        if i > 2: 
            if  i <= len(upx)-2: 
                if abs(unsortedNode[-1][1] -  upx[i]) > abs(unsortedNode[-1][1] - upx[i+1]) and abs(unsortedNode[-1][1] -  upx[i]) > 0.01: 
                    # print("%.3f, %.3f"%( upx[i], upy[i]))
                    deletedx.append(upx[i]); deletedy.append(upy[i])
                    upx = np.delete(upx, i)
                    upy = np.delete(upy, i)
                    continue

            if cangle <pangle or abs(cangle - pangle) > 0.25  : 
                # print("%.3f, %.3f"%( upx[i], upy[i]))
                # print ("** Angle %.2f, pangle=%.2f"%(math.degrees(cangle), math.degrees(pangle)) )
                deletedx.append(upx[i]); deletedy.append(upy[i])
                upx = np.delete(upx, i)
                upy = np.delete(upy, i)
                # print ("             deleted")
                continue
        unsortedNode.append([cangle, upx[i], upy[i]])
        pangle = cangle  
        

        i += 1 

    print(" math degrees ", math.degrees(0.25))
    # print ("************************* ")
    i = 1
    while i < len(dwx): 
        cangle = Angle_Between_Vectors([0, 0, -1, 0], [0, 0, dwx[i], dwy[i]])
        # print ("* Angle %.2f"%(math.degrees(cangle)))
        if i > 2: 
            if math.degrees(cangle) < 30.0 and math.degrees(pangle) > 150: 
                pangle = 0 
            if cangle < pangle or abs(cangle - pangle) > 0.25: 
                deletedx.append(dwx[i]); deletedy.append(dwy[i])
                dwx = np.delete(dwx, i)
                dwy = np.delete(dwy, i)
                # print ("             deleted")
                continue
        pangle = cangle  
        unsortedNode.append([cangle+math.pi, dwx[i], dwy[i]])

        i += 1 

    sortedNode = Sorting(unsortedNode)

    nlist = []
    totalareapointx=[]
    totalareapointy=[]
    for n in sortedNode: 
        # print (" Angle %f"%(n[0]))
        totalareapointx.append(n[1])
        totalareapointy.append(n[2])

    totalarea = Area(totalareapointx, totalareapointy)
    print ("** Total Area = %.2f"%(totalarea*10000))
    
    Roundness = (totalarea)/(maxlength*Max_ContactWidth)
    contactRatio = ActualArea/totalarea*100
    if ActualArea : 
        AvgContactPress = contactForce/ActualArea
    else: 
        AvgContactPress = 0 
    

    # ActualArea 

    positions=[75, 77.5, 80, 82.5, 85, 87.5, 90, 92.5, 95, 97.5]
    FootPrintRatio=[]
    print ("** Contact Ratio from 75%~97.5%")
    for pos in positions: 
        length1 = contactLength(PointsUp, PointsDown, maxes, position=pos)
        pos_opp = 100 -pos
        length2 = contactLength(PointsUp, PointsDown, maxes, position=pos_opp)
        length = (length1+length2)/2
        FootPrintRatio.append(length/maxlength)
        print ("    position %d%%, Length=%.2f Footprint ratio=%.2f"%(pos, length*1000, length/maxlength))

    del(FootPrintRatio[-1])
    Gullwing_ShapeFactor = 1 - max(FootPrintRatio) 
    print ("** GullWing Shape Factor = %.3f"%(Gullwing_ShapeFactor))

    ## Gullwing Local Shape Factor : ave. slope of Shoulder Rib Leading Edge 
    length_Up=[]

    positions=[75, 77.5, 80, 82.5, 85, 87.5, 90, 92.5, 95.0]
    for pos in positions: 
        length_Up.append(contactLength(PointsUp, PointsDown, maxes, position=pos, half='up'))
        # print (pos, contactLength(PointsUp, PointsDown, maxes, position=pos, half='up'))
    delW = Center_ContactWidth  * 0.025
    slopes=[]
    for i, lg in enumerate(length_Up): 
        if i ==0: 
            plg = lg 
            continue 
        slopes.append((-lg+plg)/delW)
        plg = lg 
    # print ("slopes", slopes)
    # print (" del w", delW, Center_ContactWidth)

    Gullwing_localShapeFactor = np.average(np.array(slopes))
    print ("** GullWing Local Shape Factor = %.3f"%(Gullwing_localShapeFactor))
    ## suppose 
    ## Center_ContactWidth 


    ## Roundness shape factor 
    ## 2nd-order Polynominal fit for center 50% of footprint 
    ## Roundness = Coef of 2nd-order term 
    fpcX, fpcY, fpcV = Cropping_footprint_dots(ptx=px, pty=py, ptv=pv, xrange=[xMinCnt*0.5, xMaxCnt*0.5])
    fitting = 2 
    PointsUp, PointsDown, PointsPos, PointsNeg, _, coefs =\
         SearchPoints(fpcX, fpcY, fitting=fitting, fittingmargin=fittingmargin, coef=True, savefig=False)

    #############################################
    maxY = np.max(fpcY)
    minY = np.min(fpcY)
    upMax = np.max (PointsUp[:,1])
    downMin = np.min(PointsDown[:,1])
    if upMax > maxY + gapLimit :   
        ixs = np.where(fpcY<maxY - abs(upMax - maxY)/2)[0]
        ixx = np.where(fpcY>minY + abs(downMin - minY)/2)[0]
        ix = np.intersect1d(ixs, ixx)
        fpcX = fpcX[ix]
        fpcY = fpcY[ix]
        fpcV = fpcV[ix]

        PointsUp, PointsDown, PointsPos, PointsNeg, _, coefs =\
        SearchPoints(fpcX, fpcY, fitting=fitting, fittingmargin=fittingmargin, coef=True, savefig=False)

    ##############################################

    RoundnessShapeFactor = coefs[0][0]   
    print ("** Roundness Shape Factor=%.3f"%(RoundnessShapeFactor))

    Squareness_ShapeFactor = maxlength / Max_ContactWidth
    print ("** Squareness Shape Factor=%.3f"%(Squareness_ShapeFactor))
    Coef_variation_cPress = None 
    if not isinstance(cp, type(None)): 
        idx = np.where(cp[:,1]>0)[0]
        pres = cp[idx]
        avgPres = np.average(pres[:,1])
        stdPres = np.std(pres[:,1])
        Coef_variation_cPress = stdPres / avgPres 
        print ('** Coef. Of Variation of CPress=%.3f'%(Coef_variation_cPress))
        print ("   - Avg=%.2f, Std=%.2f"%(avgPres, stdPres ))

    areapoints=[totalareapointx, totalareapointy]
    basicFPC = [Lengths, Widths, ActualArea, totalarea, Roundness, contactRatio, AvgContactPress]
    advFPC = [Gullwing_ShapeFactor, Gullwing_localShapeFactor, RoundnessShapeFactor, Squareness_ShapeFactor, Coef_variation_cPress]

    return areapoints, basicFPC, advFPC, lines 
    # areapoints=[totalareapointx, totalareapointy]
    # basicFPC = [Lengths, Widths, ActualArea, totalarea, Roundness, contactRatio, AvgContactPress]
    # advFPC = [Gullwing_ShapeFactor, Gullwing_localShapeFactor, RoundnessShapeFactor, Squareness_ShapeFactor, Coef_variation_cPress]

   

def Filtering_preFilter(pointsX, pointsY, mht=5e-3, image=False): 
    
    initX = pointsX; initY = pointsY

    if pointsY[1] < 0: 
        negative = True 
        pointsY *= -1 
    else: negative = False 
    
    multi = 3
    gap = 3.0E-3

    N = len(pointsX)
    yMax = np.max(pointsY); yMin = np.min(pointsY)

    mid = (yMax + yMin) / 2.0 
    current = yMin
    while current < mid: 
        ix1 = np.where(pointsY>=current)[0]
        ix2 = np.where(pointsY< current + gap) [0]
        ix3 = np.where(pointsY< current + gap*2) [0]
        ix = np.intersect1d(ix1, ix2) 
        jx = np.intersect1d(ix1, ix3) 
        if len(ix) ==1 or len(ix) * multi < len(jx) or current < yMin + gap   : 
            j = 0 
            while j < len(pointsX): 
                if pointsY[j] < current + gap and pointsY[j] > current : 
                    pointsX = np.delete(pointsX, j)
                    pointsY = np.delete(pointsY, j)
                    continue 
                j += 1 
        current += gap  


    current = yMax
    while current > mid+mid/3.0: 
        ix1 = np.where(pointsY<=current)[0]
        ix2 = np.where(pointsY > current - gap) [0]
        ix3 = np.where(pointsY > current - gap*2) [0]
        ix = np.intersect1d(ix1, ix2) 
        jx = np.intersect1d(ix1, ix3) 
        if len(ix) ==1 or len(ix) * multi < len(jx)  : 
            j = 0 
            while j < len(pointsX): 
                if pointsY[j] > current - gap and pointsY[j] < current : 
                    pointsX = np.delete(pointsX, j)
                    pointsY = np.delete(pointsY, j)
                    continue 
                j += 1 

        current -= gap  

    if negative: 
        pointsY *= -1 

    # print (" ")

    if True: 
        if negative: imageName  = "fitting_dots-prefilter_N%d.png"%(len(initX))
        else: imageName= "fitting_dots-prefilter_Ps%d.png"%(len(initX))
        import matplotlib.pyplot as plt 
        plt.clf()
        plt.scatter(initX, initY, label="points for fitting", edgecolors=None, linewidths=0.0, color='black', s=30)
        plt.scatter(pointsX, pointsY, label="Fitted points", edgecolors=None, linewidths=0.0, color='red', s=15)
        # print ("Final Fitting coefficient", A)
        plt.legend(loc=4)
        plt.xlim(-0.08, 0.08)
        plt.ylim(0.0, 0.08)
        plt.axis('equal')
        plt.savefig(imageName)
        plt.clf()

    return pointsX, pointsY 


def Filtering_Up_Down_points(pointsX, pointsY, mht=5e-3, image=False): 
    # pointsX, pointsY = Filtering_preFilter(pointsX, pointsY, mht=mht, image=image)

    i= 2 
    mht = 5e-3
    while i < len(pointsX)-2: 
        dy1 = 100.0; dy2 = 100.0
        if abs(pointsX[i] - pointsX[i-1]) :#<= dist*2: 
            dy1 = abs(pointsY[i] - pointsY[i-1])
        if abs(pointsX[i+1] - pointsX[i]) :#<= dist*2: 
            dy2 = abs(pointsY[i+1] - pointsY[i])
        # pslp =abs((pointsY[i-1]-pointsY[i-2])/(pointsX[i-1]-pointsX[i-2]))
        cslp =abs((pointsY[i]-pointsY[i-1])/(pointsX[i]-pointsX[i-1]))
        nslp =abs((pointsY[i+1]-pointsY[i])/(pointsX[i+1]-pointsX[i-1]))

        if abs((nslp-cslp)/cslp) > 1.0 or dy1 <100 or dy2 <100: 
            if dy1 <100 and dy2 <100: 
                
                if dy1 > mht and dy2 > mht : 
                    pointsX = np.delete(pointsX, i)
                    pointsY = np.delete(pointsY, i)
                    continue 
            elif dy1 ==100: 
                if dy2 > mht: 
                    pointsX = np.delete(pointsX, i+1)
                    pointsY = np.delete(pointsY, i+1)
                    i += 1 
                    continue 
            elif dy2 ==100: 
                if dy1 > mht: 
                    pointsX = np.delete(pointsX, i-1)
                    pointsY = np.delete(pointsY, i-1)
                    continue 

        i += 1 
    return pointsX, pointsY 

def searchbounarypoints(px=[], py=[], pxy=[], pos="updown", dist=5.0E-3, ix=0, iy=1):
    try:
        if pos =="updown": 
            px = np.array(px); py=np.array(py)
        else:
            pt = np.array(py); py=np.array(px)
            px = pt 
    except: 
        pxy = np.array(pxy)
        if pos =="updown": 
            px = pxy[:,ix]
            py = pxy[:,iy]
        else:
            px = pxy[:,iy]
            py = pxy[:,ix]


    HighValulePoint=[]
    LowValuePoint = []
    vmax = px.max()
    vmin = px.min()
    len = abs(vmax-vmin)

    # pN = px.size

    N = int(len/dist)+1
    
    position = vmin
    for i in range(N): 
        cnt = 0 
        ymax =0; ymin = 0
        id1 = np.where(px>=position)[0]
        id2 = np.where(px<position+dist)[0]
        id = np.intersect1d(id1, id2)
        tpy = py[id]
        ymax = tpy.max()
        ymin = tpy.min()
        cnt = len(tpy)

        if cnt > 3: 
            if pos =="updown": 
                HighValulePoint.append([position, ymax])
                LowValuePoint.append([position, ymin])
            else: 
                HighValulePoint.append([ymax, position])
                LowValuePoint.append([ymin, position])
                
        position += dist 

    return HighValulePoint, LowValuePoint, N


def SearchPoints(px=[], py=[], dist=5.0E-03, fitting=6, fittingmargin=0.001, \
    coef=False, savefig=False, delete='max', nofitting=False ): 
    # print ("no for search", len(px))
    print ("** Curve Fitting Order=%d"%(fitting))
    case=0 ## case 0 is more faster 

    marginoferr =fittingmargin  # 1mm 
    margin_gap = fittingmargin / 1000.0 # 0.001mm
    removedotstoabsolute = 1E-03  #3E-3 
    # debug =1
    xmax = px.max()
    xmin = px.min()
    ymax = py.max()
    ymin = py.min()

    pN = len(px)
    LengthPointUp=[]
    LengthPointDown=[]
    WidthPositionPos=[]
    WidthPositionNeg=[]

    ################################################################
    if nofitting : dist = 1.0e-3
    xmid = (xmax + xmin) / 2.0
    xlen = abs(xmax-xmid)
    N = int(xlen/dist) + 1
    uN = N 
    # if debug == 1: print (N, "x mid ", xmid)
    position = xmid
    subN  = 3
    dd = dist/subN 
    for i in range(N): 
        cnt = 0 
        yx =0; ym = 0
        if case== 0:
            id1 = np.where(px >= position)[0]
            id2 = np.where(px < position + dist)[0]
            ids = np.intersect1d(id1, id2)
            # print (len(id1), len(id2), ids)
            tpy = py[ids]; tpx = px[ids]
            if len(ids)>0: 
                ix = np.where(tpy>0)[0];    uy = tpy[ix];  subX = tpx[ix]
                tp = position 
                subY = []
                for k in range(int(subN)): 
                    id1 = np.where(subX >= tp)[0]
                    id2 = np.where(subX < tp +dd)[0]
                    idn = np.intersect1d(id1, id2)
                    
                    if len(idn): 
                        mx = np.max(uy[idn])
                        subY.append(mx)
                    tp += dd 
                if len(subY): 
                    avgY = np.average(np.array(subY))
                    LengthPointUp.append([position+dist*0.5, avgY])
                elif len(uy): 
                    avgY = np.max(uy)
                    LengthPointUp.append([position+dist*0.5, avgY])

                ix = np.where(tpy<0)[0];    dy = tpy[ix];  subX = tpx[ix]
                tp = position 
                subY = []
                for k in range(int(subN)): 
                    id1 = np.where(subX >= tp)[0]
                    id2 = np.where(subX < tp +dd)[0]
                    idn = np.intersect1d(id1, id2)
                    if len(idn): 
                        subY.append(np.min(dy[idn]))
                    tp += dd 
                if len(subY): 
                    avgY = np.average(np.array(subY))
                    LengthPointDown.append([position+dist*0.5, avgY])
                elif len(dy): 
                    avgY = np.min(dy)
                    LengthPointDown.append([position+dist*0.5, avgY])
            position += dist 

    xlen = abs(xmid-xmin)
    N = int(xlen/dist) + 1
    dN = N
    # if debug == 1: print (N)
    position = xmid
    for i in range(N): 
        cnt = 0 
        yx =0; ym = 0#; yxx = 0; ymx = 0
        if case ==0: 
            id1 = np.where(px < position)[0]
            id2 = np.where(px >= position - dist)[0]
            ids = np.intersect1d(id1, id2)
            tpy = py[ids]; tpx = px[ids]
            if len(ids)>0:
                ix = np.where(tpy > 0)[0];    uy = tpy[ix];  subX = tpx[ix]
                tp = position 
                subY = []
                for k in range(int(subN)): 
                    id1 = np.where(subX < tp)[0]
                    id2 = np.where(subX >= tp -dd)[0]
                    idn = np.intersect1d(id1, id2)
                    if len(idn): 
                        subY.append(np.max(uy[idn]))
                    tp -= dd 
                if len(subY): 
                    avgY = np.average(np.array(subY))
                    LengthPointUp.append([position-dist*0.5, avgY])
                elif len(dy): 
                    avgY = np.max(dy)
                    LengthPointUp.append([position-dist*0.5, avgY])

                ix = np.where(tpy<0)[0];    dy = tpy[ix];  subX = tpx[ix]
                tp = position 
                subY = []
                for k in range(int(subN)): 
                    id1 = np.where(subX < tp)[0]
                    id2 = np.where(subX >= tp -dd)[0]
                    idn = np.intersect1d(id1, id2)
                    if len(idn): 
                        subY.append(np.min(dy[idn]))
                    tp -= dd 
                if len(subY): 
                    avgY = np.average(np.array(subY))
                    LengthPointDown.append([position-dist*0.5, avgY])
                elif len(dy): 
                    avgY = np.min(dy)
                    LengthPointDown.append([position-dist*0.5, avgY])
            position -= dist
    
    LengthPointUp = Sorting(LengthPointUp)
    p1 = np.array(LengthPointUp)
    pointsX = p1[:, 0]
    pointsY = p1[:, 1]

    if nofitting : 
        rawX = np.array(pointsX); rawY = np.array(pointsY)
        p2 = np.array(LengthPointDown)
        pointsX = p2[:, 0]
        pointsY = p2[:, 1]
        rawX = np.append(rawX, pointsX); rawY = np.append(rawY, pointsY)
        return rawX, rawY

    ################################################################
    ymid = (ymax + ymin) / 2.0
    ylen = abs(ymax-ymid)
    N = int(ylen/dist) + 1
    # if debug == 1: print (N, "y mid ", ymid, "Points=", pN, "dist", dist)
    position = ymid
    for i in range(N): 
        cnt = 0 
        xx =0.0; xm = 0.0#; xxy = 0; xmy = 0
        # if case ==0: 
        id1 = np.where(py >= position)[0]
        id2 = np.where(py < position + dist)[0]
        ids = np.intersect1d(id1, id2)
        tpx = px[ids]
        if len(ids): 
            try: 
                maxarg = np.argmax(tpx)
                minarg = np.argmin(tpx)
                xx = tpx[maxarg]
                xm = tpx[minarg]
                WidthPositionPos.append([xx, position+0.5*dist])
                WidthPositionNeg.append([xm, position+0.5*dist])
            except: 
                pass 
        position += dist 

    ylen = abs(ymid-ymin)
    N = int(ylen/dist) + 1
    # if debug == 1: print (N)
    position = ymid
    for i in range(N): 
        cnt = 0 
        xx =0; xm = 0
        id1 = np.where(py < position)[0]
        id2 = np.where(py > position - dist)[0]
        ids = np.intersect1d(id1, id2)
        tpx = px[ids]
        if len(ids): 
            try: 
                maxarg = np.argmax(tpx)
                minarg = np.argmin(tpx)
                xx = tpx[maxarg]
                xm = tpx[minarg]
                WidthPositionPos.append([xx, position-0.5*dist])
                WidthPositionNeg.append([xm, position-0.5*dist])
            except: 
                pass 
        position -= dist 
    # print ("W Pos", WidthPositionPos)
    # print ("W Neg", WidthPositionNeg)
    ordern = fitting
    WidthPositionPos = Sorting(WidthPositionPos)
    WidthPositionNeg = Sorting(WidthPositionNeg)

    ###################################################
    print ("** Curving Fitting Up points")
    
    orgptx = pointsX; orgpty = pointsY
    pointsX, pointsY = Filtering_Up_Down_points(pointsX, pointsY, mht=5e-3, image=False)
    A, err, err_value = curvefitting(pointsX, pointsY, order=ordern)

    # print("err", err)
    # print("max", err_value)
    ptx = pointsX
    pty = pointsY
    ermax = err.max()

    cnt = 0
    perr = 0.0
   
    N = len(pty)
    
    # print (" No of pooints - %d"%(N))
    # print(" Margin Gap", margin_gap)
    # print ("*******************************************")
    delpx =[]; delpy=[]
    # erM = 0.01 
    while abs(perr - err_value) > margin_gap:  
        
        if delete=='max': ermaxindx =  np.argmax(err)
        else: ermaxindx =  np.argmin(err)
        delpx.append(ptx[ermaxindx]); delpy.append(pty[ermaxindx])
        tx = np.delete(ptx, ermaxindx)
        ty = np.delete(pty, ermaxindx)

        perr = err_value
        A, err, err_value = curvefitting(tx, ty, order=ordern)
        ptx = tx
        pty = ty
        cnt += 1
        
        # if  np.max(err) > erM or np.min(err) < -erM: continue 

        if err_value > perr: break
        if err_value / len(err) <  marginoferr: break
        if cnt > N/2: break

    ## verification 
    nd =[]
    for x, y in zip(delpx, delpy): 
        nd.append([x, y])

    nd = sorted(nd, key=lambda x: x[0])
    cn = len(delpx)
    g = 3.01e-3; gy=2.0e-3
    cnt = 0 
    for i, n in enumerate(nd): 
        if i <= 1 or i >= cn -2: continue 
        dxp = nd[i-1][0]; dyp = nd[i-1][1]
        dx0 = nd[i  ][0]; dy0 = nd[i  ][1]
        dxn = nd[i+1][0]; dyn = nd[i+1][1]
        if abs(dx0-dxp) <= g and abs(dx0-dxn) <= g and abs(dy0-dyp) < gy and abs(dy0-dyn) < gy: 
            ptx = np.append(ptx, dx0)
            pty = np.append(pty, dy0)
            cnt += 1 
    if cnt : 
        A, err, err_value = curvefitting(ptx, pty, order=ordern)


    print ("* Fitting Coefficient: ", A)
    # if not isfile("fitting_dots.png"): savefig = True  
    # else: savefig = False 
    savefig  = False    
    if savefig:
        
        import matplotlib.pyplot as plt 
        plt.clf()
        plt.scatter( p1[:, 0],  p1[:, 1], label="points for fitting", edgecolors=None, linewidths=0.0, color='black', s=50)
        plt.scatter(pointsX, pointsY, label="after filtered", edgecolors=None, linewidths=0.0, color='blue', s=20)
        plt.scatter(delpx, delpy, label="deleted points", edgecolors=None, linewidths=0.0, color='white', s=10)
        
        plt.scatter(ptx, pty, label="Fitted points", edgecolors=None, linewidths=0.0, color='red', s=5)

        pitting_points =  True  
        if pitting_points: 
            s = np.min(px)
            widthmax = np.max(px)
            d = 0.01
            ptx = []; pty=[]
            while s < widthmax: 
                yt = 0.0
                for i in range(ordern+1):
                    yt += A[i]*pow(s, float(i))
                ptx.append(s)
                pty.append(yt)
                s += d 

            plt.plot(ptx, pty, label="calculated points",  color='red', lw=1)

        # print ("Final Fitting coefficient", A)
        plt.legend(loc=4)
        plt.xlim(-0.3, 0.3)
        plt.axis('equal')
        plt.savefig("fitting_dots%d.png"%(len(p1)))
        plt.clf()
    # print (f"Up - Fitting Error {round(err_value, 6)} (iter-{cnt-1}, {len(ptx)})") 
    print (" Up - Fitting Error %.6f, (iteration %d, current points=%d)"%(err_value, cnt-1, len(ptx)))
    # print (" Coefficients")
    # print (A)
    if coef: 
        coefs =[A]
    
    # import matplotlib.pyplot as plt 
    # plt.scatter(ptx, pty) 

    p1 =[]
    position = xmid
    # for k in range(uN-1): 
    for k in range(uN-1):
        yt = 0.0
        for i in range(ordern+1):
            yt += A[i]*pow(position, float(i))
        p1.append([position, yt]) 
        position += dist
    position = xmid
    for k in range(dN-1): 
        yt = 0.0
        for i in range(ordern+1):
            yt += A[i]*pow(position, float(i))
        p1.append([position, yt]) 
        position -= dist

    ###################################################
    print ("** Curving Fitting Down points")
    p2 = np.array(LengthPointDown)
    pointsX = p2[:, 0]
    pointsY = p2[:, 1]
    orgptx = pointsX; orgpty = pointsY
    pointsX, pointsY = Filtering_Up_Down_points(pointsX, pointsY, mht=5e-3)
    
    A, err, err_value = curvefitting(pointsX, pointsY, order=ordern)
    ptx = pointsX
    pty = pointsY
    ermax = err.min()

    cnt = 0
    perr = 0
    N = len(pty)
    delpx =[]; delpy=[]
    while abs(perr - err_value) > margin_gap:  
        tx = ptx
        ty = pty 
        if delete=='max': erminindx =  np.argmin(err)
        else: erminindx =  np.argmax(err)
        delpx.append(tx[erminindx]); delpy.append(ty[erminindx])
        tx = np.delete(tx, erminindx)
        ty = np.delete(ty, erminindx)

        perr = err_value
        A, err, err_value = curvefitting(tx, ty, order=ordern)

        ptx = tx
        pty = ty
        cnt += 1
        if err_value > perr: break
        if err_value / len(err) <  marginoferr: break
        if cnt > N/2: break

     ## verification 
    nd =[]
    for x, y in zip(delpx, delpy): 
        nd.append([x, y])

    nd = sorted(nd, key=lambda x: x[0])
    cn = len(delpx)
    g = 3.01e-3; gy = 2.0e-3
    cnt = 0 
    for i, n in enumerate(nd): 
        if i <= 1 or i >= cn -2: continue 
        dxp = nd[i-1][0]; dyp = nd[i-1][1]
        dx0 = nd[i  ][0]; dy0 = nd[i  ][1]
        dxn = nd[i+1][0]; dyn = nd[i+1][1]
        if abs(dx0-dxp) <= g and abs(dx0-dxn) <= g and abs(dy0-dyp) < gy and abs(dy0-dyn) < gy: 
            ptx = np.append(ptx, dx0)
            pty = np.append(pty, dy0)
            cnt += 1 
    if cnt : 
        A, err, err_value = curvefitting(ptx, pty, order=ordern)


    # print (f"Down - Fitting Error {round(err_value, 6)} (iter-{cnt-1}, {len(ptx)})") 
    print ("* Fitting Coefficient: ", A)
    print (" Down - Fitting Error %.6f, (iteration %d, %d)"%(err_value, cnt-1, len(ptx)))  
    # print(" ERR VALUE", err_value, delpx[-1], delpy[-1])
    # print (" Coefficients")
    # if not isfile("fitting_dots_down.png"): savefig = True  
    # else: savefig = False 
    
    if savefig:
        
        import matplotlib.pyplot as plt 
        plt.clf()
        plt.scatter( p2[:, 0],  p2[:, 1], label="points for fitting", edgecolors=None, linewidths=0.0, color='black', s=50)
        plt.scatter(pointsX, pointsY, label="after filtered", edgecolors=None, linewidths=0.0, color='blue', s=20)
        plt.scatter(delpx, delpy, label="deleted points", edgecolors=None, linewidths=0.0, color='white', s=10)
        
        plt.scatter(ptx, pty, label="Fitted points", edgecolors=None, linewidths=0.0, color='red', s=5)
        pitting_points =  True  
        if pitting_points: 
            s = np.min(px)
            widthmax = np.max(px)
            d = 0.01
            ptx = []; pty=[]
            while s < widthmax: 
                yt = 0.0
                for i in range(ordern+1):
                    yt += A[i]*pow(s, float(i))
                ptx.append(s)
                pty.append(yt)
                s += d 

            plt.plot(ptx, pty, label="calculated points",  color='red', lw=1)
        # print ("Final Fitting coefficient", A)
        plt.legend(loc=4)
        plt.xlim(-0.3, 0.3)
        plt.axis('equal')
        plt.savefig("fitting_dots_down%d.png"%(len(p2)))
        plt.clf()


    # print (A)
    if coef: 
        coefs.append(A)
    # plt.scatter(ptx, pty) 
    p2 =[]
    position = xmid
    for k in range(uN-1): 
        yt = 0.0
        for i in range(ordern+1):
            yt += A[i]*pow(position, float(i))
        p2.append([position, yt]) 
        # print ("%f, %f"%(position, yt))
        position += dist
    position = xmid
    for k in range(dN-1): 
        yt = 0.0
        for i in range(ordern+1):
            yt += A[i]*pow(position, float(i))
        p2.append([position, yt]) 
        # print ("%f, %f"%(position, yt))
        position -= dist

    ymax = max(p1[1])
    ymin = min(p2[1]) 
    if coef: 
        return np.array(p1), np.array(p2), np.array(WidthPositionPos), np.array(WidthPositionNeg), [xmax, xmin, ymax, ymin], coefs
    else: 
        return np.array(p1), np.array(p2), np.array(WidthPositionPos), np.array(WidthPositionNeg), [xmax, xmin, ymax, ymin]

def curvefitting(nx=[], ny=[], order=1):

    S = int(order) + 1
    N = len(nx)

    sx = []; sxy =[]
    for i in range(int(S*2)):
        tsx =0.0
        tsy =0.0 
        for j in range(N): 
            tsx += pow(nx[j], float(i))
            if (i<int(S)): tsy += pow(nx[j], float(i)) * ny[j]
        sx.append(tsx)
        if (i<int(S)): sxy.append(tsy)
    matrix = np.zeros((int(S), int(S)))
    for i in range(int(S)):
        for j in range(int(S)):
            matrix[i][j] = sx[i+j]
    try:
        RM = np.linalg.inv(matrix)
        A = np.matmul(RM, sxy)
    except:
        # print(matrix)
        return np.zeros(S), np.zeros(N), 0.0

    err = []
    evalue = 0.0
    for i in range(N): 
        y =0 
        for j in range(int(S)):
            y += A[j] * pow(nx[i], float(j))
        evalue += abs(ny[i]-y)
        err.append(-ny[i]+y)
        
    err = np.array(err)

    return A, err, evalue 



def extract_profile_crown(fname): 

    node, element, elset, surface, tie, xy, rims = cute.Mesh2DInformation(fname)

    crown = ['CTR', 'SUT','UTR', 'TRW']

    ctb = element.Elset("CTB") 
    if not len(ctb.Element) : 
        ctb = element.Elset("CTR") 
        crown = ['SUT','UTR', 'TRW']
    for st in crown: 
        eset = element.Elset(st)
        if len(eset.Element): 
            ctb.Combine(eset)
    
    edges = ctb.OuterEdge(node)

    bt = element.ElsetToEdge("BT2")

    if len(bt.Edge): 
        for ed in bt.Edge: 
            edges.Add(ed)

    bt1 = element.ElsetToEdge("BT1")
    if len(bt1.Edge): 
        for ed in bt1.Edge: 
            edges.Add(ed) 
    bt3 = element.ElsetToEdge("BT3")

    if len(bt3.Edge): 
        for ed in bt3.Edge: 
            edges.Add(ed) 

    return [np.array(node.Node), edges.Edge]
    

class inDoorFootprint(): 
    def __init__(self, jobFile=None) :

        with open(jobFile) as F: 
            lines = F.readlines()
        start = False 
        xs=[]; ys=[]; vs=[]
        for i, line in enumerate(lines): 
            if not start: 
                if "SIZE" in line.upper(): print(line.strip())
                if "PATTERN" in line.upper() : print(line.strip())
                if "INFL" in line.upper() : print(line.strip())
                if "TEST" in line.upper() : print(line.strip())
                if "SPEC. No." in line.upper(): print(line.strip())
                

            if "Camera Raw Data" in line:
                start = True 
                continue 
            if start : 
                dl = line.split(",")
                
                tl = []
                cal = 0 
                for j, d in enumerate(dl):
                    v = float(d.strip())
                    if j > 0:  
                        if pv >= 10  and v >=10 : 
                            ys.append(float(i))
                            xs.append(float(j))
                            vs.append(v)
                        pv = v 
                    else: 
                        pv = v 

        self.X = np.array(xs)
        self.Y = np.array(ys)
        self.P = np.array(vs)

        ## scale :
        #  신 장비의 경우 해상도 2020 * 2020 이며, 이미지의 크기는 500 mm * 500 mm
        #  구 장비의 경우 해상도 480 * 640 이며, 이미지의 크기는 292 mm * 390 mm 
        if len(lines) > 1000: 
            xdots=2022; ydots=2022; xwidth=0.5; ywidth=0.5 
        else: 
            xdots=480; ydots=640; xwidth=0.292; ywidth=0.390
        self.scaling( xdots=xdots, ydots=ydots, xWidth=xwidth, ywidth=ywidth )

    def scaling(self, xdots=2022, ydots=2022, xWidth=0.5, ywidth=0.5): 
        self.X /= xdots  / xWidth 
        self.Y /= ydots  / ywidth

    def positioning(self, value=40, add=0.01): 

        jx = np.where(self.P > value)[0]

        tx40 = self.X[jx] 
        mi = np.min(tx40); mj = np.max(tx40) 

        cx = mi; step = 1.0E-3

        start = -100.0; end=100.0
        minNumber = 100
        while cx < mj: 

            ix1 = np.where(tx40>cx)[0]
            ix2 = np.where(tx40<cx+step)[0]
            ix = np.intersect1d(ix1, ix2)
            if len(ix) > minNumber and start ==-100.0:  
                start = cx 
                break 
            cx += step 

        cx = mj 
        while cx > mi: 
            ix1 = np.where(tx40<cx)[0]
            ix2 = np.where(tx40>cx-step)[0]
            ix = np.intersect1d(ix1, ix2)

            if len(ix) > minNumber and end ==100.0:  
                end = cx 
                break 
            cx -= step 
        
        mxd = (start+end)/2.0 



        ix1 = np.where(self.X>=start-add)[0]
        ix2 = np.where(self.X<=end+add)[0]
        ix = np.intersect1d(ix1, ix2)
        self.X = self.X[ix] - mxd 
        self.Y = self.Y[ix]
        self.P = self.P[ix]


        jx = np.where(self.P > value)[0]
        ty40 = self.Y[jx] 
        mi = np.min(ty40); mj = np.max(ty40)
        cy = mi; step = 1.0E-3

        start = -100.0; end=100.0
        # minNumber = 100
        while cy < mj: 

            ix1 = np.where(ty40>cy)[0]
            ix2 = np.where(ty40<cy+step)[0]
            ix = np.intersect1d(ix1, ix2)
            if len(ix) > minNumber and start ==-100.0:  
                start = cy 
                break 
            cy += step 

        cy = mj 
        while cy > mi: 
            ix1 = np.where(ty40<cy)[0]
            ix2 = np.where(ty40>cy-step)[0]
            ix = np.intersect1d(ix1, ix2)

            if len(ix) > minNumber and end ==100.0:  
                end = cy 
                break 
            cy -= step 
        

        myd = (start+end)/2.0 
        
        ix1 = np.where(self.Y>=start-add)[0]
        ix2 = np.where(self.Y<=end+add)[0]
        ix = np.intersect1d(ix1, ix2)

        self.Y = self.Y[ix] - myd 
        self.X = self.X[ix]
        self.P = self.P[ix]



    def imageBoundary(self, gap=5): 
        
        N = len(self.P)
        gP = np.average(self.P)*2/3
        xs =[]; ys=[]; ps=[]
        for i, _ in enumerate(self.P): 
            if i < 100 or i > N -100: continue 
            if self.P[i] > gP: continue 
            if (self.P[i-1] - self.P[i]) > gap and (self.P[i-2] - self.P[i]) > gap: 
                    xs.append(self.X[i])
                    ys.append(self.Y[i])
                    ps.append(self.P[i])
                    continue 
            if (self.P[i+1] - self.P[i]) > gap and (self.P[i+2] - self.P[i]) > gap: 
                    xs.append(self.X[i])
                    ys.append(self.Y[i])
                    ps.append(self.P[i])
                    continue 
            if (self.P[i-1] - self.P[i]) > gap*2: 
                    xs.append(self.X[i])
                    ys.append(self.Y[i])
                    ps.append(self.P[i])
                    continue 

        self.X = np.array(xs)
        self.Y = np.array(ys)
        self.P = np.array(ps)




