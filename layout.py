import numpy as np 
import math 
from basicfunctions import isNumber, Angle_3Points, NodeDistance

import time 
def timer(func): 
    def wrapper(*args, **kwargs): 
        start = time.time()
        rv = func(*args, **kwargs)
        total = time.time() - start
        print ("### Time to generate mesh grid: %.2f"%(total))
        return rv 
    return wrapper

class LAYOUT: 
    def __init__(self, meshfile): 
        self.Node = NODE()
        self.Element = ELEMENT()
        self.Elset = ELSET()
        self.Surface = SURFACE()
        self.carcassGa = 0.0
        self.moldOD=0.0
        self.GD=0.0
        self.UTG=0.0
        self.edge_outer=[]
        try: 
            self.Mesh2DInformation(meshfile)
        except Exception as EX: 
            print(EX)
            print (" NOT MESH FILE!!")


    def Mesh2DInformation(self, InpFileName):
        with open(InpFileName) as INP:
            lines = INP.readlines()
        skipComment = 0 
        for line in lines:
            if line[0] == '\n': continue 
            if '**' in line:
                if 'The List of Rebar Element sets' in line: 
                    skipComment =1 
                if skipComment ==0: 
                    if "C01" in line : 
                        w = line.split(",")
                        self.carcassGa = float(w[6].strip())
                    if 'CAVITY_OD' in line and not 'CODE' in line and not 'NODE' in line: 
                        txt = line.split(":")[1].strip()
                        self.moldOD = float(txt)/1000.0
                    if 'GROOVE DEPTH' in line : 
                        txt = line.split(":")[1].strip()
                        self.GD = float(txt)/1000.0
                continue 
            if '*' in line:
                # print(line.strip())
                word = list(line.split(','))
                if '*HEADING' in word[0].upper(): 
                    spt = 'HD'
                elif '*NODE' in word[0].upper(): 
                    spt = 'ND'
                elif '*ELEMENT' in word[0].upper():
                    # EL = list(word[1].split('='))
                    # EL = EL[1].strip()
                    if 'MGAX' in line:
                        spt = 'M1'
                    elif 'CGAX3' in line :
                        spt = 'C3'
                    elif 'CGAX4' in line :
                        spt = 'C4'
                    else:
                        spt = 'NN'
                elif '*SURFACE' in word[0].upper(): 
                    spt = 'SF'
                    name = word[2].split('=')[1].strip()

                    self.Surface.AddName(name)
                #                    print ('Name=', name, 'was stored', Surface.Surface)
                elif '*TIE' in word[0].upper(): 
                    spt = 'TI'
                    name = word[1].split('=')[1].strip()
                elif '*ELSET' in word[0].upper(): 
                    spt = 'ES'
                    name = word[1].split('=')[1].strip()
                    if name != "BETWEEN_BELTS" and name != "BD1" and name != "BetweenBelts":
                        self.Elset.AddName(name)

                elif 'NSET' in word[0].upper(): 
                    spt = 'NS'
                    name = word[1].split('=')[1].strip()

                else:
                    spt = ''
            else:
                word = list(line.strip().split(','))
                if spt == 'HD':
                    pass
                if spt == 'ND':
                    self.Node.Add([int(word[0]), float(word[3]), float(word[2]), float(word[1])])
                if spt == 'M1':
                    # Element   [EL No,                  N1,          N2,  N3, N4,'elset Name', N,  Area/length, CenterX, CenterY]
                    self.Element.Add([int(word[0]), int(word[1]), int(word[2]), 0, 0, '', 2])#, math.sqrt(math.pow(C1[1] - C2[1], 2) + math.pow(C1[2] - C2[2], 2)), (C1[1] + C2[1]) / 2.0, (C1[2] + C2[2]) / 2.0])
                if spt == 'C3':
                    self.Element.Add([int(word[0]), int(word[1]), int(word[2]), int(word[3]), 0, '', 3])#, A[0], A[1], A[2]])
                if spt == 'C4':
                    self.Element.Add([int(word[0]), int(word[1]), int(word[2]), int(word[3]), int(word[4]), '', 4])#, A[0], A[1], A[2]])
                if spt == 'NS':
                    pass
                if spt == 'ES':
                    if name.upper() != "BETWEEN_BELTS" and name.upper() != "BD1":
                        # if 'BTT' in name: print (line.strip())
                        for i in range(len(word)):
                            if isNumber(word[i]) == True:
                                self.Elset.AddNumber(int(word[i]), name)
                if spt == 'SF':
                    pass

                else:
                    pass

        if not self.moldOD: 
            hmax =0  
            for n in self.Node.Node: 
                if hmax < n[3]: 
                    hmax = n[3] 
            self.moldOD = hmax * 2 

        for i in range(len(self.Elset.Elset)):
            for j in range(1, len(self.Elset.Elset[i])):
                for k in range(len(self.Element.Element)):
                    if self.Elset.Elset[i][j] == self.Element.Element[k][0]:
                        self.Element.Element[k][5] = self.Elset.Elset[i][0]
                        break

    def OuterEdge(self): 
        self.edge_outer = self.Element.OuterEdge(self.Node)

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

        npnode = np.array(self.Node)
        idx = np.where(npnode[:,0] == n)[0]
        if len(idx) ==0: 
            print ("Cannot Find Node (%d)"%(n))
            NullList = [0, 0.0, 0.0, 0.0]
            return NullList
        else: 
            return self.Node[idx[0]]
        
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
    def Combine(self, node):
        N = len(node.Node)
        for i in range(N): 
            self.Add(node.Node[i])
    
class SURFACE:
    def __init__(self):
        self.Surface=[]
    def Add(self, surf):
        self.Surface.append(surf)

    def AddName(self, name): 
        pre = 0 
        for surf in self.Surface:
            if surf[0] == name: 
                pre =1
                break 
        if pre == 0: 
            self.Surface.append([name])
    def AddSurface(self, name, number, face): 
        pre = 0 
        for i, surf in enumerate(self.Surface):
            if surf[0] == name: 
                self.Surface[i].append([number, face])
                pre =1 
                break 
        if pre ==0: 
            self.AddName(name)
            self.AddSurface(name, number, face)

class TIE:
    def __init__(self):
        self.Tie = []
    def Add(self, name, slave, master):
        self.Tie.append([name, slave, master])
class EDGE:
    def __init__(self):
        self.Edge = []
    def __del__(Self): 
        # print ("EDGE IS DELETED")
        pass

    def Help(self):
        print ("*********************************************************************************")
        print ("EDGE : Node1, Node2, Elset_Name, FacdID, Element_No, D")
        print (" D : -1= Edge, 0 = Free Edge, -2 = not Free Edge, 1= outer edges, Above 1(2~) = Tie No")
        print ("***************************************************************************")

    def Add(self, edge):
        self.Edge.append(edge)

    def Print(self, tail=0, head=0, **kwargs):
        Print_list(self.Edge, tail=tail, head=head, **kwargs)
        
    def Combine(self, iEdge):
        N = len(iEdge.Edge)
        for i in range(N): 
            self.Add(iEdge.Edge[i])

    def Sort(self, reverse=False):
        edges=[]
        e1 =[]; e2=[]
        for i, e in enumerate(self.Edge):
            edges.append([e[0], e[1], i])
            e1.append(e[0])
            e2.append(e[1])
            # print("%d, %d, %d"%(e[0], e[1], e[4]))

        edges = np.array(edges)
        e1 = np.array(e1)
        e2 = np.array(e2)

        
        if reverse == False:        
            e = np.setdiff1d(e1, e2) ## it remains the 1st node at the first element 
            col = 0
            ncol =1 
        else:                       
            e = np.setdiff1d(e2, e1)  ## it remains the 2nd node at the last element 
            col = 1
            ncol = 0 
        

        ## if it's is closed, there will be nothing 

        srted=[] 
        idx = np.where(edges[:,col] == e) [0]
        while len(idx)>0: 
            i = edges[ idx[0] ][2]
            # print ("sorting", self.Edge[i])
            srted.append(self.Edge[i])

            idx = np.where(edges[:,col] == self.Edge[i][ncol])[0]
        self.Edge = srted 
        return srted 
    def Image(self, Node, file='EDGE', marker ="", dpi=100):
        MembWidth = 0.5
        color = 'black'

        fig, ax = plt.subplots()
        ax.axis('equal')
        ax.axis('on')
        N = len(self.Edge)
        MinX = 100000.0
        MaxX = -100000.0
        MinY = 100000.0
        MaxY = -100000.0
        for i in range(N):
            N1 = Node.NodeByID(self.Edge[i][0])
            N2 = Node.NodeByID(self.Edge[i][1])
            x1 = N1[2]
            y1 = N1[3]
            x2 = N2[2]
            y2 = N2[3]
            plt.plot([x1, x2], [y1, y2], color, lw=MembWidth, marker=marker)
            if MinX > x1:
                MinX = x1
            if MaxX < x1:
                MaxX = x1
            if MinY > y1:
                MinY = y1
            if MaxY < y1:
                MaxY = y1
            if MinX > x2:
                MinX = x2
            if MaxX < x2:
                MaxX = x2
            if MinY > y2:
                MinY = y2
            if MaxY < y2:
                MaxY = y2
        plt.xlim(MinX - 0.01, MaxX + 0.01)
        plt.ylim(MinY - 0.01, MaxY + 0.01)
        plt.savefig(file, dpi=dpi)
class ELSET:
    def __init__(self):
        self.Elset = []

    def Print(self):
        print ("*************************************************************")
        print ("** [[Elset1, E11, E12, ..], {Elset2, E21, E22, ...], ...]")
        print ("*************************************************************")

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
        exist = 0
        for i in range(len(self.Elset)):
            if self.Elset[i][0] == name:
                exist = 1
                self.Elset[i].append(n)
                break
        if exist == 0:
            self.Elset.append([name, n])
class ELEMENT:
    def __init__(self):
        self.Element = []
    def Add(self, e):
        self.Element.append(e)
    def Delete(self, n):
        N = len(self.Element)
        for i in range(N):
            if self.Element[i][0] == n:
                del (self.Element[i])
                break
    def DeleteDuplicate(self, id=1):
        npe = []
        for el in self.Element: 
            npe.append([el[0], el[2], el[3], el[4], el[6]])
        npe = np.array(npe)

        npe0 = npe[:,0]
        npe1 = np.unique(npe0)
        if len(npe0) != len(npe1): 
            i = 0 
            while i < len(self.Element): 
                idx = np.where(npe[:,0]==self.Element[i][0])[0]
                if len(idx) > 1: 
                    for ix in idx: 
                        if ix != i: 
                            del(self.Element[ix])
                            npe = np.delete(npe, ix, axis=0)
                            print (" >> Element %d was deleted."%(self.Element[i][0]))
                            i -= 1 
                            break 
                i += 1 

    def Nodes(self, **args):

        Node = NODE()
        for key, value in args.items():
            if key == 'Node' or key == 'node':
                Node = value


        I = len(self.Element)
        NL = []

        nds = []
        for el in self.Element: 
            nds.append(el[1]); nds.append(el[2]); 
            if el[3]>0: nds.append(el[3])
            if el[4]>0: nds.append(el[4])
        nds = np.array(nds, dtype=np.int32)
        NL = np.unique(nds)

        if len(Node.Node) == 0:  return NL
        
        npnd = np.array(Node.Node)
        NC = NODE() 
        for nd in NL: 
            ix = np.where(npnd[:,0] == nd) [0]
            if len(ix) == 1: 
                NC.Add([ int(npnd[ix[0]][0]), npnd[ix[0]][1], npnd[ix[0]][2], npnd[ix[0]][3]] )
        return NC 

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
    def Combine(self, element):
        N=len(element.Element)
        for i in range(N): 
            self.Add(element.Element[i])
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
        OEdges = OuterEdge(FEdges, Node, self)
        return OEdges
        ## other method
        npn = np.array(Node.Node)
        free = self.FreeEdge()
        outer = EDGE()
        my = 10**7
        start = 0 
        for i, e in enumerate(free.Edge): 
            ix = np.where(npn[:,0] == e[0])[0][0]
            if npn[ix][2]> 0 and npn[ix][3]<my: 
                my = npn[ix][3]
                start = i 
        
        outer.Add(free.Edge[start])
        end = free.Edge[start][0]
        nxt = free.Edge[start]
        cnt = 0
        N=5
        while nxt[1] != end: 
            cnt +=1
            if cnt > 1000: 
                print(" Error to find outer edges")
                break 
            nt = self.NextEdge(nxt, free)
            if len(nt) ==1: 
                nxt = nt[0]
            elif len(nt) ==2: 
                fst = 0
                ne = nt[fst]
                ns = ne[0]
                for i in range(N): 
                    ne = self.NextEdge(ne, free)
                    if len(ne) ==2 : 
                        if ne[0][1] == ns or ne[1][1] == ns: 
                            fst = 1
                        break 
                    else: 
                        ne = ne[0]
                        if ne[1] == ns: 
                            fst = 1
                            break 
                nxt = nt[fst]
            outer.Add(nxt)
        return outer 
    def TieEdge(self, Node):
        FreeEdge = self.FreeEdge()
        OuterEdge(FreeEdge, Node, self)  # Don't Delete this line.
        TieNum = 1
        i = 0;        iTemp = 0;        j = 0
        connectedEdge = []
        TEdge = EDGE()
        while i < len(FreeEdge.Edge):
            if FreeEdge.Edge[i][5] < 1:
                TieNum += 1
                nodeStart = FreeEdge.Edge[i][0]
                FreeEdge.Edge[i][5] = TieNum
                TEdge.Add(FreeEdge.Edge[i])  # marked as TIE edge with No.
                iTemp = i
                while FreeEdge.Edge[iTemp][1] != nodeStart:
                    j += 1
                    if j > 100:
                        break  # in case infinite loop
                    connectedEdge = NextEdge(FreeEdge, iTemp)  # find next edge
                    if len(connectedEdge) == 1:  # in case of being found just 1 edge
                        iTemp = connectedEdge[0]
                    elif len(connectedEdge) == 2:  # when other tie is connected (2 ties are connected)
                        if FreeEdge.Edge[connectedEdge[0]][1] == nodeStart:
                            iTemp = connectedEdge[0]
                        elif FreeEdge.Edge[connectedEdge[1]][1] == nodeStart:
                            iTemp = connectedEdge[1]
                        else:
                            if FreeEdge.Edge[connectedEdge[0]][5] < 1 and FreeEdge.Edge[connectedEdge[1]][5] < 1:
                                iTemp = FindTieLoop(nodeStart, connectedEdge, FreeEdge)
                            elif FreeEdge.Edge[connectedEdge[0]][5] < 1:
                                iTemp = connectedEdge[0]
                            elif FreeEdge.Edge[connectedEdge[1]][5] < 1:
                                iTemp = connectedEdge[1]
                            else:
                                print ('[INPUT] {' + str(FreeEdge.Edge[connectedEdge[0]]) + ',' + str(FreeEdge.Edge[connectedEdge[1]]) + ' (0) TIE Conection InCompletion')
                                break
                    else:
                        print ('[INPUT] 2 or more Ties are Connected.')
                        break
                    # After finding next TIE Edge ################################
                    FreeEdge.Edge[iTemp][5] = TieNum
                    TEdge.Add(FreeEdge.Edge[iTemp])
                del connectedEdge
                connectedEdge = []
            i += 1
        return TEdge
    def MasterSlaveEdge(self, Node, Op = 0, **args):

        npn = np.array(Node.Node)
        for key, value in args.items():
            if key == 'op':
                Op = int(value)
        
        TieEdge = self.TieEdge(Node)
        iNum = 2
        mlength = 0 
        ErrRatio = 0.01
        
        NumTie = 0 
        N = len(TieEdge.Edge)
        for i in range(N):
            if TieEdge.Edge[i][5] > NumTie:
                NumTie = TieEdge.Edge[i][5]
            
        iMaster = []
        while iNum <=NumTie:
            MaxLength = 0.0
            SumLength = 0.0
            k = 0
            Save = 0
            while k < N:
                if TieEdge.Edge[k][5] == iNum:
                    # N1 = Node.NodeByID(TieEdge.Edge[k][0])
                    # N2 = Node.NodeByID(TieEdge.Edge[k][1])
                    ix = np.where(npn[:,0] == TieEdge.Edge[k][0])[0][0]; N1 = Node[ix]
                    ix = np.where(npn[:,0] == TieEdge.Edge[k][1])[0][0]; N2 = Node[ix]
                    Length = math.sqrt((N1[2]-N2[2])* (N1[2]-N2[2]) + (N1[3]-N2[3])*(N1[3]-N2[3]))
                    SumLength += Length
                    if Length > MaxLength:
                        MaxLength = Length
                        Save = k
                k += 1
            SumLength -= MaxLength
            if SumLength > MaxLength * (1+ErrRatio) or SumLength < MaxLength * (1-ErrRatio):
                print ('ERROR::PRE::TIE CREATION INCOMPLETE ON', TieEdge.Edge[Save][3])
            iMaster.append(Save)
            iNum += 1
        
        MasterEdge=EDGE()
        SlaveEdge =EDGE()
        M = len(iMaster)
        for i in range(N):
            exist = 0
            for j in range(M):
                if i == iMaster[j]:
                    MasterEdge.Add(TieEdge.Edge[i])
                    exist =1
                    break
            if exist == 0:
                SlaveEdge.Add(TieEdge.Edge[i])
        
        ## Op == 0 return MasterEdge and SlaveEdge
        ## Op == 1 return only Master Edge
        ## Op == 2 return Only Slave Edge
        if Op == 0:
            return MasterEdge, SlaveEdge
        elif Op == 1:
            return MasterEdge
        else:
            return SlaveEdge
    def Print(self, **kwargs): 
        Print_list(self.Element, **kwargs)
    def NextEdge(self, edge, edges, rev=0): 
        ed =[]
        if rev == 0: 
            for e in edges.Edge: 
                if e[0] == edge[1]: 
                    ed.append(e)
        else: 
            for e in edges.Edge: 
                if e[1] == edge[0]: 
                    ed.append(e)

        return ed 

def FreeEdge(edge):
    FEdge = EDGE()
    edges = np.array(edge.Edge)
    for i, ed in enumerate(edges): 
        ix1 = np.where(edges[:,1] == ed[0])[0]
        ix2 = np.where(edges[:,0] == ed[1])[0]
        ix = np.intersect1d(ix1, ix2)
        if len(ix) ==0: 
            edge.Edge[i][5] = 0
            FEdge.Add(edge.Edge[i])
        else: 
            edge.Edge[i][5] = -2
    return FEdge
def OuterEdge(FreeEdge, Node, Element):
    N = len(FreeEdge.Edge)
    MinY = 9.9E20
    cNodes = [0]
    npn = np.array(Node.Node)
    
    for i in range(N):
        ix = np.where(npn[:,0]==FreeEdge.Edge[i][0])[0][0]; N1 = Node.Node[ix]
        ix = np.where(npn[:,0]==FreeEdge.Edge[i][1])[0][0]; N2 = Node.Node[ix]

        if N1[3] < MinY:
            MinY = N1[3]
            cNodes[0] = N1[0]
        if N2[3] < MinY:
            MinY = N2[3]
            cNodes[0] = N2[0]
    if cNodes[0] == 0:
        zMin = np.min(npn[:, 3])
        ix = np.where(npn[:,3] == zMin)[0]

        for x in ix: 
            if npn[:,2] < 0: 
                cNodes = npn[x]
                break 
        # cNodes[0] = Node.NodeIDByCoordinate('z', 0.0, closest=1)

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

def readBodyLayout(meshfile=None): 
    # print ("Mesh file", meshfile)
    layout = LAYOUT(meshfile)
    # print(layout.Node.Node)
    layout.OuterEdge()

    beltEdge =[]   # print ("EDGE : Node1, Node2, Elset_Name, FacdID, Element_No, D")
    carcass =[]
    BDR = ELEMENT()
    BDL = ELEMENT()
    npn = np.array(layout.Node.Node)

    for el in layout.Element.Element: 
        if 'BT1' in el[5] or 'BT2' in el[5] or 'BT3' in el[5]: 
            beltEdge.append([el[1], el[2], el[5], 0, el[0], 0])
        if 'C01' in el[5]: 
            carcass.append([el[1], el[2], el[5], 0, el[0], 0])
        if 'BEAD_L' in el[5] or 'BEAD_R' in el[5]:
            ix = np.where(npn[:, 0]==el[1])[0][0]
            if npn[ix][2]>0: 
                BDR.Add(el)
            else:
                BDL.Add(el)

    edge_BD = BDR.OuterEdge(layout.Node)
    edge_BDL = BDL.OuterEdge(layout.Node)
    # edge_BD.Image(layout.Node, file="RBead.png")
    # edge_BDL.Image(layout.Node, file="LBead.png")
    # print (len(edge_BD.Edge));     print (len(edge_BDL.Edge))
    edge_BD.Combine(edge_BDL)
    # print (len(edge_BD.Edge))
    # layout.edge_outer.Image(layout.Node, file='outer.png')
    # print (layout.edge_outer)
    fp = open('outer.tmp', 'w')
    for ed in layout.edge_outer.Edge: 
        fp.write("%d, %d, 0, 0, %d, 0\n"%(ed[0], ed[1], ed[4]))
    fp.close()
    fp = open('belt.tmp', 'w')
    for ed in beltEdge: 
        fp.write("%d, %d, 0, 0, %d, 0\n"%(ed[0], ed[1], ed[4]))
    fp.close()
    fp = open('bead.tmp', 'w')
    for ed in edge_BD.Edge: 
        fp.write("%d, %d, 0, 0, %d, 0\n"%(ed[0], ed[1], ed[4]))
    for ed in edge_BDL.Edge: 
        fp.write("%d, %d, 0, 0, %d, 0\n"%(ed[0], ed[1], ed[4]))
    fp.close()
    fp = open('carcass.tmp', 'w')
    for ed in carcass: 
        fp.write("%d, %d, 0, 0, %d, 0\n"%(ed[0], ed[1], ed[4]))
    fp.close()
    try: 
        topsurf_edges = GrooveDetectionFromEdge(layout.edge_outer, layout.Node)
        no_tread = 10**7
        ffname = 'topsurf.tmp'
        tfp = open(ffname, 'w')
        for ed in topsurf_edges.Edge:
            if ed[7] == 0: 
                tfp.write("%10d, %10d, %10d, %10d, %10d\n"%(ed[4]+ no_tread, ed[0]+no_tread, ed[1]+no_tread, ed[1]+no_tread, ed[0]+no_tread))
        tfp.close()
    
    except: 
        print ("** NO Groove FOUND. Mesh was generated from axi.")
        pass 
        
    return layout.edge_outer.Edge, beltEdge, edge_BD.Edge, carcass 

def LayoutMesh_From_axi(axi="", limit=10000, output="" ):
    print (" AXI FILE : %s"%(axi))
    with open(axi) as I: 
        lines = I.readlines()
    fp = open(output, 'w')
    fp.write("**************************************\n")
    fp.write("** TIRE MESH from AXI \n")
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
        else: 

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

def GrooveDetectionFromEdge(oEdge, node, OnlyTread=1, TreadNumber=10000000, **args):
    for key, value in args.items():
        if key == 'onlytread' or key == 'tread':
            OnlyTread = int(value)
        if key == 'treadno' or key == 'Treadno' or key == 'Treadstartno' or key == 'TreadStartNo':
            TreadNumber = int(value)
        # if key == 'offset':
            # Offset = int(value)
        # if key == 'step':
            # Step = int(value)

    N = len(oEdge.Edge)
    LengthOfEdge = len(oEdge.Edge[0])

    
    # print 'oEDGE - ', oEdge.Edge[5], len(oEdge.Edge)
    # print '        ', node.NodeByID(oEdge.Edge[5][0]), node.NodeByID(oEdge.Edge[5][1])
    # Reverse = 0
     # Original Edge : [Node 1, Node 2, Elset Name, Face ID, Tie No]
    # Edge from INP(Mesh file) : Add Groove/Tread Edge Mark
    # Data to Add : Edge Length, Groove Mark for Node 1, Groove Mark for Node 2, Tread Mark]
    # Groove Mark : 0 - node not in Groove, 1 - node in Groove
    # Tread Mark : 0 - Edge not on Tread, 2 - Edge on Tread

    TreadElset = ['CTB', 'SUT', 'CTR', 'UTR']#, 'TRW']
    TN = len(TreadElset)
    cEdge = EDGE()
    counting = 0
    for i in range(N):
        TreadID = 0
        for j in range(TN):
            if oEdge.Edge[i][2] == TreadElset[j]:
                TreadID = 2
                cEdge.Add(oEdge.Edge[i])
                N1 = node.NodeByID(oEdge.Edge[i][0])
                N2 = node.NodeByID(oEdge.Edge[i][1])
                length = NodeDistance(N1, N2, xy=0)
                cEdge.Edge[counting][-1] = length
                cEdge.Edge[counting].append(0)
                cEdge.Edge[counting].append(0)
                cEdge.Edge[counting].append(TreadID)
                counting += 1
                break


    CriticalAngle = 45.0   ### if there is an error to detect grooves, check this angle... 
    N = len(cEdge.Edge)
    cEdge.Edge[0][6] = 0
    PA = 0.0  # Previous Angle
    PN = 0  # Previous Node Groove ID : 0 - Not in the groove, 1 - in the groove
    for i in range(1, N):
        N1 = node.NodeByID(cEdge.Edge[i][0])
        N2 = node.NodeByID(cEdge.Edge[i][1])
        if N2[2] - N1[2] != 0: 
            A = math.degrees(math.atan((N2[3] - N1[3]) / (N2[2] - N1[2])))
        else: 
            A = 90.0
            
        if N2[2] < N1[2]:
            A = -A
        cEdge.Edge[i][3] = round(A, 2)

        if i == 1:
            if A > CriticalAngle:
                cEdge.Edge[0][6] = 1  ## To see if TBR shoulder side part
                PN = cEdge.Edge[0][7] = 0
                PA = A
                continue
            else:
                PN = cEdge.Edge[i][7] = 0
                PA = A
                continue
        else:
            cEdge.Edge[i][6] = PN
            if PN == 0:
                if PA > CriticalAngle:
                    if A > CriticalAngle:
                        cEdge.Edge[i][6] = 1
                        PN = cEdge.Edge[i][7] = 0
                        PA = A
                        continue
                    else:  # if  (a < CriticalAngle and a > -CriticalAngle)
                        PN = cEdge.Edge[i][7] = 0
                        PA = A
                        continue
                else:
                    if A < CriticalAngle and A > -CriticalAngle:
                        PN = cEdge.Edge[i][7] = 0
                        PA = A
                        continue
                    elif A < -CriticalAngle and abs(PA - A) > CriticalAngle:
                        PN = cEdge.Edge[i][7] = 1
                        PA = A
                        continue
                    else:  # a > CriticalAngle:
                        PN = cEdge.Edge[i][7] = 1
                        PA = A
                        continue
            else:
                if PA > CriticalAngle:
                    if A < CriticalAngle and A > -CriticalAngle:
                        cEdge.Edge[i][6] = 1
                        PN = cEdge.Edge[i][7] = 0
                        PA = A
                        continue
                    else:
                        PN = cEdge.Edge[i][7] = 1
                        PA = A
                        continue
                else:
                    PN = cEdge.Edge[i][7] = 1
                    PA = A
                    continue

                    
    for i in range(N-1, 0, -1):
        if (cEdge.Edge[i][6] == 0 and cEdge.Edge[i][7] == 0):
            break
        elif cEdge.Edge[i][6] == 1 and cEdge.Edge[i][7] == 1:
            cEdge.Edge[i][6] = 1
            cEdge.Edge[i][7] = 0
        elif cEdge.Edge[i][6] == 0 and cEdge.Edge[i][7] == 1:    
            cEdge.Edge[i][6] = 1
            cEdge.Edge[i][7] = 0
            break
    
    for i in range(N): 
        N1 = node.NodeByID(cEdge.Edge[i][0])
        N2 = node.NodeByID(cEdge.Edge[i][1])
        DIST = NodeDistance(N1,N2, xy=23)
        cEdge.Edge[i][5]=DIST
        N1 = node.NodeByID(cEdge.Edge[i][0])
        N2 = node.NodeByID(cEdge.Edge[i][1])
        N3 = [9999, N1[1], N1[2]+0.1, N1[3]]

        Angle = Angle_3Points(xs=[N1[2], N2[2], N3[2]], ys=[N1[3], N2[3], N3[3]], center_index=1 )
        # Angle = Angle_3nodes(N3, N2, N1, xy=23)
        cEdge.Edge[i][3]=Angle*180.0/3.141592

    # sys.exit()
    
    return cEdge

def getNode_2dMesh(meshfile): 
    with open(meshfile) as DFM: 
        lines = DFM.readlines()
    
    cmd = None 
    NODE=[]
    for line in lines: 
        if "**" in line : continue 
        if "*" in line :
            if "*NODE" in line: 
                cmd = 'ND'
            else: 
                cmd = None 
        else: 
            if cmd == 'ND': 
                wd = line.split(",")
                NODE.append([int(wd[0].strip()), float(wd[3].strip()), float(wd[2].strip()), float(wd[1].strip())])

    return np.array(NODE)

@timer 
def makeDatFile_ODB(npn, pv, npp, filename, direction='X+'): 
    fp = open(filename, 'w')

    if "X" in direction: 
        dx = 2; dy=3;  dz = 1 
    elif 'Z' in direction: 
        dx = 2; dy = 1; dz = 3

    fp.write("*NODE\n")
    for i, n in enumerate(npn): 
        fp.write("%10d, %12.6f, %12.6f, %12.6f, %12.6f, 0.0\n"%(n[0], n[dy], n[2], n[dz], pv[i]))

    if len(npp) > 100000: 
        # fp.write("** %s\n"%(npp[0]))
        tpp =[]
        for i, n in enumerate(npn): 
            ix=np.where(npp[:,2] == n[0])[0]
            for x in ix: 
                tpp.append(npp[x])
    
        npp = np.array(tpp)

        # fp.write("** %s\n"%(npp[0]))

    # fp.write("** NO. of element %d\n"%(len(npp)))
    fp.write("*ELEMENT\n")
    for pp in npp: 
        ix1 = np.where(npp[:,0]==pp[0])[0]
        ix = np.where(npn[:,0] == npp[ix1[0]][2])[0]
        if not len(ix): 
            continue 

        if len(npp[ix1]) > 1: 
            itx = npp[ix1][1]
        else: 
            continue 

        itx_unique = np.unique(itx)
        if len(itx_unique) <3: 
            continue 
        for tx in itx_unique: 
            ix2 = np.where(npp[:,1]==tx)[0]
            ix = np.intersect1d(ix1, ix2)
            sfs = npp[ix]
            sf = np.unique(sfs, axis=0)
            if len(sf) ==3 or len(sf) ==4: 
                try: 
                    xs =[]; ys=[]; ids =[]
                    idx = np.where(npn[:,0] == sf[0][2])[0][0]; 
                    xs.append(npn[idx][dx]); ys.append(npn[idx][dy]); ids.append(npn[idx][0])
                    idx = np.where(npn[:,0] == sf[1][2])[0][0]; 
                    xs.append(npn[idx][dx]); ys.append(npn[idx][dy]); ids.append(npn[idx][0])
                    idx = np.where(npn[:,0] == sf[2][2])[0][0]
                    xs.append(npn[idx][dx]); ys.append(npn[idx][dy]); ids.append(npn[idx][0])
                    if len(sf) == 4: 
                        idx = np.where(npn[:,0] == sf[3][2])[0][0]
                        xs.append(npn[idx][dx]); ys.append(npn[idx][dy]); ids.append(npn[idx][0])
                except: 
                    # print ("No node: ", sf[0][2], sf[1][2], sf[2][2], sf[3][2])
                    continue 
                _, _, srted= points_ConvexHull(xs, ys)
                if len(srted) != len(xs): 
                    # print (sf[0][0], len(srted), int(ids[srted[0]]), int(ids[srted[1]]), int(ids[srted[2]]), ",", ids)
                    for i in range(len(xs)): 
                        for j in range(len(srted)): 
                            if ids[i] == ids[srted[j]]: 
                                break 
                        else: 
                            srted.append(i)
                            break 
                    
                if len(sf) == 4:
                    fp.write("%10d, %10d, %10d, %10d, %10d\n"%(sf[0][0], ids[srted[0]], ids[srted[1]], ids[srted[2]], ids[srted[3]]))
                else : 
                    fp.write("%10d, %10d, %10d, %10d, %10d\n"%(sf[0][0], ids[srted[0]], ids[srted[1]], ids[srted[2]], ids[srted[2]]))
    
    fp.close()

    
def points_ConvexHull (xs=None, ys=None): 
    xs = np.array(xs)
    ys = np.array(ys)

    ## start point : left bottom point 
    miny = np.min(ys)
    ix = np.where(ys==miny)[0]
    mx = xs[ix]
    arx = np.argmin(mx)
    sx = xs[ix[arx]]; sy = ys[ix[arx]] ## 
    sortedX =[sx]; sortedY = [sy]; mch=[ix[arx]]
    #################################################
    px = sx; py = sy 
    pre_angle = 0.0
    stop_angle = 6.2657320147 # math.pi*(2 - 1.0/180.0)
    dpi = np.pi * 2

    
    mchIdx=-1
    I = len(xs)
    for i in range(I): 
        current_minAngle = 10000.0
        for j in range(I): 
            if xs[j] == px and ys[j] == py : continue 
            
            cos = (xs[j]-px) / math.sqrt( (xs[j]-px)**2 + (ys[j]-py)**2) 
            if cos >1.0: cos = 1.0 
            elif cos < -1.0: cos = -1.0 
            angle = math.acos(cos)

            if py > ys[j]: 
                angle = dpi - angle 

            if current_minAngle > angle and angle >= pre_angle: 
                current_minAngle = angle 
                tx = xs[j]; ty=ys[j]; mchIdx=j 
        if sx == tx and sy == ty: 
            break 
        if angle > stop_angle: 
            break 

        sortedX.append(tx)
        sortedY.append(ty)
        px = tx; py=ty 
        pre_angle = current_minAngle
        mch.append(mchIdx)

    i = 0 
    while i < len(sortedX): 
        j = i + 1
        while j < len(sortedX): 
            if sortedX[i] == sortedX[j] and sortedY[i] == sortedY[j]:  
                del(sortedX[j])
                del(sortedY[j])
                continue 
            j += 1
        i += 1
    return sortedX, sortedY, mch

def CalculateAngleFrom3Node(N1, N2, Center, XY=13, **args):
    for key, value in args.items():
        if key == 'xy':
            XY = int(value)
    
    ix = int(XY/10)
    iy = int(XY)%10  
    
    V1x = N1[ix] - Center[ix]
    V1y = N1[iy] - Center[iy]
    V2x = N2[ix] - Center[ix]
    V2y = N2[iy] - Center[iy]
    
    L1 = math.sqrt(V1x*V1x + V1y*V1y)
    L2 = math.sqrt(V2x*V2x + V2y*V2y)
    
    # print 'L1', L1, ',V1:', V1x, ',', V1y
    # print "L2", L2, ',V2:', V2x, ',', V2y
    try:
        return math.acos(round((V1x*V2x+V1y*V2y)/(L1*L2), 10))
    except: 
        
        print ("ERROR. During Calculating Angle between Nodes - a in acos(a) is over 1.0 or L1*L2 = 0")
        print ("  Length : ", L1,',', L2)
        print ("  Node 1 : ", N1)
        print ("  Node 2 : ", N2)
        print ("  Center : ", Center)
        return 0