
import math 

import numpy as np 



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

def Sorting(np_lst, item=0, reverse=False):
    tmpNode=[]
    for nd in np_lst:    tmpNode.append(nd)  

    try: arr = np_lst[:, item]
    except: 
        np_lst = np.array(np_lst)
        arr = np_lst[:, item]

    if reverse == False: args = np.argsort(arr)
    else:                args = np.argsort(arr)[::-1]
    
    sortedlist = []
    for i, arg in enumerate(args):
        sortedlist.append(tmpNode[arg])
    return  np.array(sortedlist)

def Area(xs=[], ys=[]): 
    x =[]; y=[]
    for px, py in zip(xs, ys):
        x.append(px)
        y.append(py)
    x.append(xs[0]); y.append(ys[0])

    A = [0.0, 0.0, 0.0]
    n = len(x)-1
    for i in range(n):
        s = x[i] * y[i + 1] - x[i + 1] * y[i]
        A[0] += s
        A[1] += (x[i] + x[i + 1]) * s
        A[2] += (y[i] + y[i + 1]) * s
    A[0] = A[0] / 2.0
    return A[0]

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

def Delete_Close_Points(lst, point_gap, backward=False): 
    if backward: 
        i = 0 
        while i < len(lst): 
            j = i +1
            fd = False  
            while j < len(lst): 
                L = math.sqrt((lst[i][2] -lst[j][2])**2  + (lst[i][3] - lst[j][3])**2)
                if L <point_gap:
                    fd = True 
                    break 
                j +=1
            if fd: 
                del(lst[i])
                continue 
            i += 1 
    else: 
        i = 0 
        while i < len(lst): 
            j = i +1
            while j < len(lst): 
                L = math.sqrt((lst[i][2] -lst[j][2])**2  + (lst[i][3] - lst[j][3])**2)
                if L <point_gap:
                    del(lst[j]) 
                    continue 
                j +=1
            i += 1 
    return lst 

def NormalVector_plane(n1, n2, n3): 
    a = [0, n1[1]-n2[1], n1[2]-n2[2], n1[3]-n2[3]]
    b = [0, n3[1]-n2[1], n3[2]-n2[2], n3[3]-n2[3]]

    II = a[2]*b[3] - a[3]*b[2]
    JJ = a[3]*b[1] - a[1]*b[3]
    KK = a[1]*b[2] - a[2]*b[1]

    return [0, II, JJ, KK]

def NormalVector(n1, n2):
    x1 = n1[1]; y1 = n1[2]; z1 = n1[3]
    x2 = n2[1]; y2 = n2[2]; z2 = n2[3]

    vx = y1*z2 - z1*y2 
    vy = z1*x2 - x1*z2 
    vz = x1*y2 - y1*x2 
    norm = [0, vx, vy, vz]
    return norm 
 

def Angle_Between_Vectors(va, vb):
    la = math.sqrt(va[1]*va[1] + va[2]*va[2] + va[3]*va[3])
    lb = math.sqrt(vb[1]*vb[1] + vb[2]*vb[2] + vb[3]*vb[3])
    cos = round((va[1]*vb[1]+va[2]*vb[2]+va[3]*vb[3]) / la/lb , 8)
    return round(math.acos(cos), 10)    
   
def rotatePoints(px, py, angle=0): 
    radian_angle = math.radians(angle)
    xs = np.cos(radian_angle) * px - np.sin(radian_angle) * py
    ys = np.sin(radian_angle) * px + np.cos(radian_angle) * py 
    return xs, ys 

def isNumber(s):
  try:
    float(s)
    return True
  except ValueError:
    return False

def Angle_3Points(xs=None, ys=None, center_index=1 ,**args):
    centerX = xs[center_index]; centerY = ys[center_index]
    idx = center_index + 1 
    if idx >=3: idx -= 3 
    V1x = xs[idx]; V1y = ys[idx]
    idx = idx + 1 
    if idx >=3: idx -= 3 
    V2x = xs[idx]; V2y = ys[idx]
    
    L1 = math.sqrt(V1x*V1x + V1y*V1y)
    L2 = math.sqrt(V2x*V2x + V2y*V2y)
    
    try:
        return math.acos(round((V1x*V2x+V1y*V2y)/(L1*L2), 10))
    except: 
        print ("ERROR. During Calculating Angle between Nodes - a in acos(a) is over 1.0 or L1*L2 = 0")
        print ("  Length : ", L1,',', L2)
        print ("  Node 1 : ", N1)
        print ("  Node 2 : ", N2)
        print ("  Center : ", Center)
        return 0

def NodeDistance(N1, N2, xy=0): 
    if xy > 10: 
        x = int(xy/10); y=int(xy%10)
        return math.sqrt((N2[x] - N1[x])**2 + (N2[y] - N1[y])**2)
    else: 
        return math.sqrt((N2[1] - N1[1])**2 + (N2[2] - N1[2])**2+ (N2[3] - N1[3])**2)

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
    
    try:
        return math.acos(round((V1x*V2x+V1y*V2y)/(L1*L2), 10))
    except: 
        
        print ("ERROR. During Calculating Angle between Nodes - a in acos(a) is over 1.0 or L1*L2 = 0")
        print ("  Length : ", L1,',', L2)
        print ("  Node 1 : ", N1)
        print ("  Node 2 : ", N2)
        print ("  Center : ", Center)
        return 0



def combinationWords(ia=None, ib=None, ic=None): 
    combis= []
    if isinstance(ic, type(None)): 
        combis= []
        ti = ia.split(",")
        tj = ib.split(",")
        for n in i: 
            for t in tj: 
                combis.append([n, t])
    else: 
        
        ti = ia.split(",")
        tj = ib.split(",")
        tk = ic.split(",")

        for k in ti: 
            for l in tj: 
                for m in tk : 
                    combis.append([k, l, m])
    return combis 