import numpy as np 
import matplotlib.pyplot as plt
from PyQt5 import QtCore, QtGui, QtWidgets, QtNetwork
from os.path import isfile 
from os import getcwd
import math, glob, json

from scipy import ndimage
from scipy.ndimage.filters import convolve

from scipy import misc
import numpy as np

class cannyEdgeDetector:
    def __init__(self, imgs, sigma=1, kernel_size=5, weak_pixel=75, strong_pixel=255, lowthreshold=0.05, highthreshold=0.15):
        self.imgs = imgs
        self.imgs_final = []
        self.img_smoothed = None
        self.gradientMat = None
        self.thetaMat = None
        self.nonMaxImg = None
        self.thresholdImg = None
        self.weak_pixel = weak_pixel
        self.strong_pixel = strong_pixel
        self.sigma = sigma
        self.kernel_size = kernel_size
        self.lowThreshold = lowthreshold
        self.highThreshold = highthreshold
        return 
    
    def gaussian_kernel(self, size, sigma=1):
        size = int(size) // 2
        x, y = np.mgrid[-size:size+1, -size:size+1]
        normal = 1 / (2.0 * np.pi * sigma**2)
        g =  np.exp(-((x**2 + y**2) / (2.0*sigma**2))) * normal
        return g
    
    def sobel_filters(self, img):
        Kx = np.array([[-1, 0, 1], [-2, 0, 2], [-1, 0, 1]], np.float32)
        Ky = np.array([[1, 2, 1], [0, 0, 0], [-1, -2, -1]], np.float32)

        Ix = ndimage.filters.convolve(img, Kx)
        Iy = ndimage.filters.convolve(img, Ky)

        G = np.hypot(Ix, Iy)
        G = G / G.max() * 255
        theta = np.arctan2(Iy, Ix)
        return (G, theta)

    def non_max_suppression(self, img, D):
        M, N = img.shape
        Z = np.zeros((M,N), dtype=np.int32)
        angle = D * 180. / np.pi
        angle[angle < 0] += 180


        for i in range(1,M-1):
            for j in range(1,N-1):
                try:
                    q = 255
                    r = 255

                   #angle 0
                    if (0 <= angle[i,j] < 22.5) or (157.5 <= angle[i,j] <= 180):
                        q = img[i, j+1]
                        r = img[i, j-1]
                    #angle 45
                    elif (22.5 <= angle[i,j] < 67.5):
                        q = img[i+1, j-1]
                        r = img[i-1, j+1]
                    #angle 90
                    elif (67.5 <= angle[i,j] < 112.5):
                        q = img[i+1, j]
                        r = img[i-1, j]
                    #angle 135
                    elif (112.5 <= angle[i,j] < 157.5):
                        q = img[i-1, j-1]
                        r = img[i+1, j+1]

                    if (img[i,j] >= q) and (img[i,j] >= r):
                        Z[i,j] = img[i,j]
                    else:
                        Z[i,j] = 0


                except IndexError as e:
                    pass

        return Z

    def threshold(self, img):

        highThreshold = img.max() * self.highThreshold;
        lowThreshold = highThreshold * self.lowThreshold;

        M, N = img.shape
        res = np.zeros((M,N), dtype=np.int32)

        weak = np.int32(self.weak_pixel)
        strong = np.int32(self.strong_pixel)

        strong_i, strong_j = np.where(img >= highThreshold)
        zeros_i, zeros_j = np.where(img < lowThreshold)

        weak_i, weak_j = np.where((img <= highThreshold) & (img >= lowThreshold))

        res[strong_i, strong_j] = strong
        res[weak_i, weak_j] = weak

        return (res)

    def hysteresis(self, img):

        M, N = img.shape
        weak = self.weak_pixel
        strong = self.strong_pixel

        for i in range(1, M-1):
            for j in range(1, N-1):
                if (img[i,j] == weak):
                    try:
                        if ((img[i+1, j-1] == strong) or (img[i+1, j] == strong) or (img[i+1, j+1] == strong)
                            or (img[i, j-1] == strong) or (img[i, j+1] == strong)
                            or (img[i-1, j-1] == strong) or (img[i-1, j] == strong) or (img[i-1, j+1] == strong)):
                            img[i, j] = strong
                        else:
                            img[i, j] = 0
                    except IndexError as e:
                        pass

        return img
    
    def detect(self):
        imgs_final = []
        for i, img in enumerate(self.imgs):    
            self.img_smoothed = convolve(img, self.gaussian_kernel(self.kernel_size, self.sigma))
            self.gradientMat, self.thetaMat = self.sobel_filters(self.img_smoothed)
            self.nonMaxImg = self.non_max_suppression(self.gradientMat, self.thetaMat)
            self.thresholdImg = self.threshold(self.nonMaxImg)
            img_final = self.hysteresis(self.thresholdImg)
            self.imgs_final.append(img_final)

        return self.imgs_final


def main(): 
    fname = "19MA120425KI_13_1_60_Raw.csv"
    fname = "21MS060641KI_1_1_100_Raw.csv"
    # jobFile, _ = QtWidgets.QFileDialog.getOpenFileName(None, "Select footprint Raw", "", "File Open(*.csv)")

    with open(fname) as F: 
        lines = F.readlines()
    start = False 
    data =[]
    row = 0 
    xs=[]; ys=[]; vs=[]
    for i, line in enumerate(lines): 
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
                        ys.append(float(i)/10)
                        xs.append(float(j)/10)
                        vs.append(v)
                    pv = v 
                else: 
                    pv = v 
            
    print (" reading done")
    xs = np.array(xs)
    ys = np.array(ys)
    vs = np.array(vs)
    for i in range(10, 15, 5): ### 20 is good .. so far 
        ix = np.where(vs>float(i))[0]
        px = xs[ix]; py=ys[ix]; pv=vs[ix]
        print (" starting to print (over %d) : %d , max Value=%f"%(i, len(px), np.max(pv)))
        plt.scatter(px, py, c=pv, cmap='rainbow', vmin=float(i)+10, vmax=np.max(pv)*0.8, s=0.1, linewidths=0)   
        plt.axis('equal') 
        # plt.show()
        # plt.xlim(750, 2000)
        # plt.ylim(800, 2200)
        plt.savefig("value_%d.png"%(i), dpi=300)
        img = cannyEdgeDetector("value_%d.png"%(i))
        plt.clf()

def EdgeDection(image): 
    img = cannyEdgeDetector(image)

def load_ISLM_Original(fname): 
    fp = open(fname, 'r')
    lines = fp.readlines()
    fp.close()
    xs=[]; ys=[]; vs=[]
    for line in lines: 
        dt = line.split(',')
        ys.append(float(dt[0].strip()))
        xs.append(float(dt[1].strip()))
        vs.append(float(dt[2].strip()))
    vmin = 100000
    xs = np.array(xs); ys=np.array(ys); vs = np.array(vs)
    idx = np.where(vs>=vmin)[0]
    xs = xs[idx]; ys= ys[idx]; vs=vs[idx]
    plt.scatter(xs,ys, c=vs, cmap='jet', vmin=vmin, vmax=vmin*10, s=5, linewidths=0.0)
    plt.axis('equal')
    plt.show()
   
def contactsurfacearea(node, area): 
    with open(node) as N: 
        lines = N.readlines()
    nodes = []
    cmd = None 
    for line in lines: 
        if "*" in line: 
            if "*NODE" in line: 
                cmd = True 
            else: 
                cmd = None 
        else: 
            if cmd : 
                words = line.split(",")
                nid = int(words[0].strip())
                nx = float(words[1].strip())
                ny = float(words[2].strip())
                nz = float(words[3].strip())
                nv = float(words[4].strip())
                nodes.append([nid, nx, ny, nz, nv])
    
    # with open(area) as N: 
    #     lines = N.readlines()

    nodes = np.array(nodes)
    
    plt.scatter(nodes[:,2], nodes[:,1], c=nodes[:,4])
    plt.axis('equal')
    plt.show()
        


if __name__ == "__main__": 

    # main()
    fname = "RND-3003473VT00036-0-D101-0001-Contour_Original.dat"
    # load_ISLM_Original(fname)


    nodefile = "D:\\01_ISLM_Scripts\\temporary\\Foot\\RND-3003763VT00009-0-D101-0002-postfoot.dat"
    areafile = "D:\\01_ISLM_Scripts\\temporary\\Foot\\RND-1031373VT00003-0-D101-0001-ContactSurfaceArea_Original.txt"

    contactsurfacearea(nodefile, areafile)
    