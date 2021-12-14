import numpy as np 
import math 
import matplotlib.pyplot as plt
# from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from time import time 

from PyQt5.QtCore import QRegExp
from PyQt5.QtGui import QColor, QTextCharFormat, QFont, QSyntaxHighlighter
from numpy.core.overrides import array_function_dispatch
from numpy.lib.arraysetops import isin



def userformat(color, style=''):
    """Return a QTextCharFormat with the given attributes.
    """
    _color = QColor()
    _color.setNamedColor(color)

    _format = QTextCharFormat()
    _format.setForeground(_color)
    if 'bold' in style:
        _format.setFontWeight(QFont.Bold)
    if 'italic' in style:
        _format.setFontItalic(True)

    return _format


# Syntax styles that can be shared by all languages
STYLES = {
    'keyword': userformat('blue'),
    'operator': userformat('red'),
    'brace': userformat('darkGray'),
    'defclass': userformat('black', 'bold'),
    'string': userformat('magenta'),
    'string2': userformat('darkMagenta'),
    'comment': userformat('darkGreen', 'italic'),
    'self': userformat('black', 'italic'),
    'numbers': userformat('brown'),
}


class PythonHighlighter (QSyntaxHighlighter):
    """Syntax highlighter for the Python language.
    """
    # Python keywords
    keywords = [
        'and', 'assert', 'break', 'class', 'continue', 'def',
        'del', 'elif', 'else', 'except', 'exec', 'finally',
        'for', 'from', 'global', 'if', 'import', 'in',
        'is', 'lambda', 'not', 'or', 'pass', 'print',
        'raise', 'return', 'try', 'while', 'yield',
        'None', 'True', 'False',
    ]

    # Python operators
    operators = [
        '=',
        # Comparison
        '==', '!=', '<', '<=', '>', '>=',
        # Arithmetic
        '\+', '-', '\*', '/', '//', '\%', '\*\*',
        # In-place
        '\+=', '-=', '\*=', '/=', '\%=',
        # Bitwise
        '\^', '\|', '\&', '\~', '>>', '<<',
    ]

    # Python braces
    braces = [
        '\{', '\}', '\(', '\)', '\[', '\]',
    ]
    def __init__(self, document):
        QSyntaxHighlighter.__init__(self, document)

        # Multi-line strings (expression, flag, style)
        # FIXME: The triple-quotes in these two lines will mess up the
        # syntax highlighting from this point onward
        self.tri_single = (QRegExp("'''"), 1, STYLES['string2'])
        self.tri_double = (QRegExp('"""'), 2, STYLES['string2'])

        rules = []

        # Keyword, operator, and brace rules
        rules += [(r'\b%s\b' % w, 0, STYLES['keyword'])
            for w in PythonHighlighter.keywords]
        rules += [(r'%s' % o, 0, STYLES['operator'])
            for o in PythonHighlighter.operators]
        rules += [(r'%s' % b, 0, STYLES['brace'])
            for b in PythonHighlighter.braces]

        # All other rules
        rules += [
            # 'self'
            (r'\bself\b', 0, STYLES['self']),

            # Double-quoted string, possibly containing escape sequences
            (r'"[^"\\]*(\\.[^"\\]*)*"', 0, STYLES['string']),
            # Single-quoted string, possibly containing escape sequences
            (r"'[^'\\]*(\\.[^'\\]*)*'", 0, STYLES['string']),

            # 'def' followed by an identifier
            (r'\bdef\b\s*(\w+)', 1, STYLES['defclass']),
            # 'class' followed by an identifier
            (r'\bclass\b\s*(\w+)', 1, STYLES['defclass']),

            # From '#' until a newline
            (r'#[^\n]*', 0, STYLES['comment']),

            # Numeric literals
            (r'\b[+-]?[0-9]+[lL]?\b', 0, STYLES['numbers']),
            (r'\b[+-]?0[xX][0-9A-Fa-f]+[lL]?\b', 0, STYLES['numbers']),
            (r'\b[+-]?[0-9]+(?:\.[0-9]+)?(?:[eE][+-]?[0-9]+)?\b', 0, STYLES['numbers']),
        ]

        # Build a QRegExp for each pattern
        self.rules = [(QRegExp(pat), index, fmt)
            for (pat, index, fmt) in rules]


    def highlightBlock(self, text):
        """Apply syntax highlighting to the given block of text.
        """
        # Do other syntax formatting
        for expression, nth, userformat in self.rules:
            index = expression.indexIn(text, 0)

            while index >= 0:
                # We actually want the index of the nth match
                index = expression.pos(nth)
                length = expression.cap(nth).length()
                self.setFormat(index, length, userformat)
                index = expression.indexIn(text, index + length)

        self.setCurrentBlockState(0)

        # Do multi-line strings
        in_multiline = self.match_multiline(text, *self.tri_single)
        if not in_multiline:
            in_multiline = self.match_multiline(text, *self.tri_double)


    def match_multiline(self, text, delimiter, in_state, style):
        """Do highlighting of multi-line strings. ``delimiter`` should be a
        ``QRegExp`` for triple-single-quotes or triple-double-quotes, and
        ``in_state`` should be a unique integer to represent the corresponding
        state changes when inside those strings. Returns True if we're still
        inside a multi-line string when this function is finished.
        """
        # If inside triple-single quotes, start at 0
        if self.previousBlockState() == in_state:
            start = 0
            add = 0
        # Otherwise, look for the delimiter on this line
        else:
            start = delimiter.indexIn(text)
            # Move past this match
            add = delimiter.matchedLength()

        # As long as there's a delimiter match on this line...
        while start >= 0:
            # Look for the ending delimiter
            end = delimiter.indexIn(text, start + add)
            # Ending delimiter on this line?
            if end >= add:
                length = end - start + add + delimiter.matchedLength()
                self.setCurrentBlockState(0)
            # No; multi-line string
            else:
                self.setCurrentBlockState(in_state)
                length = text.length() - start + add
            # Apply formatting
            self.setFormat(start, length, style)
            # Look for the next match
            start = delimiter.indexIn(text, start + length)

        # Return True if still inside a multi-line string, False otherwise
        if self.currentBlockState() == in_state:
            return True
        else:
            return False


class myCanvas(FigureCanvas):
    def __init__(self, parent=None, *args, **kwargs):
        self.figure = plt.figure()
        FigureCanvas.__init__(self, self.figure)
        self.setParent(parent)
        self.canvas = FigureCanvas(self.figure)
        self.toolbar = NavigationToolbar(self.canvas, self)

        self.xs=[]; self.ys=[]
        self.mclick=0; self.clicked = 0 
        self.dots=[];   self.circle=[];     self.chars=[]; self.lines=[]
        self.achars=[]; self.cline=[]; self.llen=[] 
        self.distance=0 
        self.fpcLine=[]
        self.fpcArea=[]
        self.fpcTexts=[]

        self.fontsize = 10

        ##############################################################
        self.norm = plt.Normalize(1,4)
        self.cmap = plt.cm.RdYlGn
        self.annot = plt.annotate("", xy=(0,0), xytext=(20,20),textcoords="offset points",
                    bbox=dict(boxstyle="round", fc="w"),
                    arrowprops=dict(arrowstyle="->"))
        self.annot.set_visible(False)
        # self.ax = self.figure.add_subplot(111)
        self.c = np.random.randint(1,5,size=15)

        self.shift = 0 
        self.scale = 1.0 

        self.valuepointx=None
        self.valuepointy=None
        self.valuepoint=None
        
        ##############################################################

    def hover(self, event): 
        vis = self.annot.get_visible() 
        if event.inaxes == self.ax: 
            cont, ind = self.sc.contains(event) 
            if cont: 
                self.annot.set_visible(False)
                pos = self.sc.get_offsets()[ind["ind"][0]]
                # Pres = (pos[1] + self.shift) / self.scale * 10E-7
                gap = 0.00001
                
                ix1=np.where(self.valuepointx>pos[0]-gap)[0]
                ix2=np.where(self.valuepointx<pos[0]+gap)[0]
                ix = np.intersect1d(ix1, ix2)
                while not len(ix): 
                    gap += 0.00002
                    ix1=np.where(self.valuepointx>pos[0]-gap)[0]
                    ix2=np.where(self.valuepointx<pos[0]+gap)[0]
                    ix = np.intersect1d(ix1, ix2)

                if len(ix): 
                    Pres = self.valuepoint[ix[0]] * 10E-7
                    self.annot = plt.annotate("Press=%.2fMpa"%(Pres), xy=[pos[0]+0.01, pos[1]+0.01], xytext=[pos[0]+0.01, pos[1]+0.01],textcoords="offset points",
                        bbox=dict(boxstyle="round", fc="w"),
                        arrowprops=dict(arrowstyle="->", alpha=0.4))
                    self.annot.set_visible(True)
                    self.figure.canvas.draw_idle()
            else: 
                if vis: 
                    self.annot.set_visible(False)
                    self.figure.canvas.draw_idle()


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

    def onReleased(self, event): 

        if event.button == 2: 
            self.mclick += 1

            if self.mclick == 4: 
                self.xs=[]
                self.ys=[]
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
                self.xs.append(event.xdata)
                self.ys.append(event.ydata)
                d, = plt.plot(event.xdata, event.ydata, 'o', color='gray')
                self.dots.append(d)
                

                x1 = self.xs[0]; x2=self.xs[1]; x3=self.xs[2]
                y1 = self.ys[0]; y2=self.ys[1]; y3=self.ys[2]
                A = x1*(y2-y3) - y1 *(x2-x3) + x2*y3 - x3*y2
                B = (x1*x1 + y1*y1)*(y3-y2) +(x2**2 + y2**2)*(y1-y3) + (x3**2+y3**2)*(y2-y1)
                C = (x1**2 + y1**2)*(x2-x3)+(x2**2+y2**2)*(x3-x1) + (x3*x3 + y3*y3)*(x1-x2)
                D = (x1*x1 + y1*y1)*(x3*y2-x2*y3)+(x2*x2+y2*y2)*(x1*y3-x3*y1)+(x3*x3+y3*y3)*(x2*y1-x1*y2)

                cx = -B/A/2.0
                cy = -C/A/2.0
                R = math.sqrt(B*B + C*C - 4*A*D) / 2/abs(A)

                self.xs.append(event.xdata)
                self.ys.append(event.ydata)
                d, = plt.plot(cx, cy, 'o', color='red')
                self.dots.append(d)

                ch = plt.text((self.xs[0]+self.xs[1])/2.0, (self.ys[0]+self.ys[1])/2.0, "R="+str(round(R*1000, 2)), size=self.fontsize, color='black')
                self.chars.append(ch)

                crcl = plt.Circle((cx, cy), R, color='gray', fill=False)
                self.ax.add_artist(crcl)
                self.circle.append(crcl)

            else:
                self.xs.append(event.xdata)
                self.ys.append(event.ydata)
                d, = plt.plot(event.xdata, event.ydata, 'o', color='gray')
                self.dots.append(d)

            current_xlim=self.ax.get_xlim()
            current_ylim=self.ax.get_ylim()
            plt.xlim(current_xlim[0], current_xlim[1])
            plt.ylim(current_ylim[0], current_ylim[1])
            self.figure.canvas.draw_idle()

        elif event.button == 1: 
            self.clicked =0  
            self.mclick = 0
            self.xs=[];           self.ys=[]

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
            self.dots=[];      self.lines=[];        self.chars=[];        self.achars=[]
            self.cline=[];     self.llen=[]

            current_xlim=self.ax.get_xlim()
            current_ylim=self.ax.get_ylim()
            plt.xlim(current_xlim[0], current_xlim[1])
            plt.ylim(current_ylim[0], current_ylim[1])
            self.figure.canvas.draw_idle()
        elif event.button ==3:
            self.clicked += 1
            self.xs.append(event.xdata)
            self.ys.append(event.ydata)
            d, = plt.plot(event.xdata, event.ydata, 'o', color='red')
            self.dots.append(d)
            N = len(self.xs)-1
            if N> 0: 
                self.distance = round( math.sqrt((self.xs[N]-self.xs[N-1])**2 + (self.ys[N]-self.ys[N-1])**2 ) *1000, 2)
                ch = plt.text((self.xs[N]+self.xs[N-1])/2.0, (self.ys[N]+self.ys[N-1])/2.0, str(self.distance), size=self.fontsize)
                self.chars.append(ch)

                ln, = plt.plot([self.xs[N-1], self.xs[N]],[self.ys[N-1], self.ys[N]], color='orange')
                self.lines.append(ln)

                if self.clicked > 2: 

                    sx = 0; sy=0
                    for x, y in zip(self.xs, self.ys): 
                        sx += x
                        sy += y
                    cx = sx/(float(N)+1)
                    cy = sy/(float(N)+1)

                    area = self.Area(ix=self.xs, iy=self.ys)
                    for achar in self.achars: 
                        achar.set_visible(False)
                    self.figure.canvas.draw_idle()
                    ach= plt.text(cx, cy, "Area="+ str(round(area*1_000_0, 2)), color='gray', size=self.fontsize)
                    self.achars.append(ach)

                    for char in self.llen:
                        char.set_visible(False)

                    for line in self.cline: 
                        line.remove()
                    self.cline=[]
                    
                    self.distance = round( math.sqrt((self.xs[N]-self.xs[0])**2 + (self.ys[N]-self.ys[0])**2 ) *1000, 2)
                    ch = plt.text((self.xs[N]+self.xs[0])/2.0, (self.ys[N]+self.ys[0])/2.0, str(self.distance), color='gray', size=self.fontsize)
                    self.llen.append(ch)

                    ln, = plt.plot([self.xs[0], self.xs[N]],[self.ys[0], self.ys[N]], color='gray', linestyle="--" )
                    self.cline.append(ln)
            current_xlim=self.ax.get_xlim()
            current_ylim=self.ax.get_ylim()
            plt.xlim(current_xlim[0], current_xlim[1])
            plt.ylim(current_ylim[0], current_ylim[1])
            self.figure.canvas.draw_idle()

    def plotComparing(self, pts, legends=None, items=None, size=1, grv=None, pressure=None, colors=None, marks=None, sizes=None): 
        self.figure.clear()
        self.ax = self.figure.add_subplot(111)
        self.ax.axis('equal')
        size *= 5

        if pressure: 
            distance = 0.1
            tpress=[]
            yrange = 0.001
            for ps in pressure:
                ix1 = np.where(ps[1] > -yrange)[0]
                ix2 = np.where(ps[1] < yrange)[0]
                ix = np.intersect1d(ix1, ix2) 
                tpress.append([ps[0][ix], ps[1][ix], ps[2][ix] ])

            position = 10; maxp = 0 
            for pt, ps in zip(pts, tpress): 
                mn = np.min(pt[1])
                if position > mn: 
                    position = mn 
                mx = np.max(ps[2])
                if maxp < mx: 
                    maxp = mx  

            height_graph = 0.2 
            scale = height_graph / maxp 
            position -= height_graph
            
            self.scale = scale 
            

        EA = len(pts)
        
        vsx=[]; vsy=[]; vsv=[]
        if isinstance(grv, type(None)): 
            cnt = 0 
            for i, pt in zip(items, pts): 
                size = sizes[cnt]
                self.ax.scatter(pt[0], pt[1], c=colors[cnt], s=size*3,   marker=marks[cnt], label=legends[cnt])
                if pressure: 
                    vx, vy, vv = self.add_pressure(pressure[cnt], yrange=yrange, distance=distance, below=True, position=position, size=size, color=colors[cnt], scale=scale, mark=marks[cnt], EA=EA )
                    vsx.append(vx); vsy.append(vy); vsv.append(vv)
                cnt += 1 
        else: 
            cnt = 0 
            for i, pt, gv in zip(items, pts, grv): 
                size = sizes[cnt]
                if pt : 
                    self.ax.scatter(pt[0], pt[1], c=colors[cnt], s=size*3,   marker=marks[cnt], label=legends[cnt])
                    if len(gv): 
                        mx = np.max(pt[0])
                        for g in gv: ## g =[n1, n2, face, en, 0, N1, N2]  [sf[1], sf[2], 3, sf[0], 0, sf[5][2], sf[5][0]]
                            xs =[]; ys =[]
                            xs.append(g[5][2]); xs.append(g[6][2])
                            ys.append(g[5][1]); ys.append(g[6][1])
                            lx = abs(xs[1]-xs[0]); ly = abs(ys[1]-ys[0])
                            if ly > lx: 
                                if abs(xs[0]) < mx : 
                                    self.ax.plot(xs, ys, c=colors[cnt], linewidth=size*0.1, ls=':' )

                    if pressure: 
                        vx, vy, vv = self.add_pressure(pressure[cnt], yrange=yrange, distance=distance, below=True, position=position, size=size, color=colors[cnt], scale=scale, mark=marks[cnt], EA=EA)
                        vsx.append(vx); vsy.append(vy); vsv.append(vv)

                cnt += 1 

        plt.legend(fontsize=8, loc=1)

        if pressure: 
            valueX=[]; valueY=[]; value=[]
            for x, y, v in zip(vsx, vsy, vsv): 
                valueX += x; valueY += y; value += v 
            self.sc = self.ax.scatter(valueX, valueY, s=size*3, edgecolors=None, c=None, linewidths=0.0, marker=None)
            self.figure.canvas.mpl_connect("motion_notify_event", self.hover)
            self.valuepointx = np.array(valueX)
            self.valuepointy = np.array(valueY) 
            self.valuepoint = np.array(value )

        # lim = 0.25
        # self.ax.scatter([lim, lim, -lim, -lim], [-lim, lim, lim, -lim], edgecolors=None, linewidths=0.0, c='gray', s=0.01)
        self.figure.canvas.mpl_connect('button_release_event', self.onReleased)
        self.figure.tight_layout()
        self.figure.canvas.draw()


    def add_pressure(self, pressure, yrange=0.001, distance=0.02, below=True, position=0, size=1, color='black', scale=1, mark='*', EA=1): 
        xs = pressure[0]; ys=pressure[1]; pv = pressure[2] 

        yrange = 0.001 
        ix1 = np.where(ys>-yrange)[0]; ix2 = np.where(ys<yrange)[0]
        ix = np.intersect1d(ix1, ix2)

        xpress = xs[ix]; press = pv[ix]
        indexs = np.argsort(xpress)
        axs = xpress[indexs]
        avs = press[indexs]
        avs *= scale 

        if not below: 
            shift =   position + distance 
        else: 
            shift = - position + distance 

        self.shift = shift
        # self.ax.scatter(axs, avs-shift, c=color, s=size, edgecolors=None, linewidths=0.0, marker=mark)
        cnt = 0 
        gap = 0.001

        minValue = np.min(avs) 
        cnt = 0 
        sumdist = 0 
        j = 0 
        for lx, ly in zip(axs, avs):
            if ly > minValue*1.2: 
                tDist = abs(lx - axs[j-1])
                if avs[j-1] < minValue*1.1 or tDist > 0.005: continue 
                sumdist += tDist
                cnt +=1 
            if cnt == 10: 
                break 
            j += 1 
        avgDist = sumdist / 10 
        if EA > 1: 
            if gap < avgDist : 
                gap = avgDist * 1.5
        else: 
            gap =  0 #avgDist * 0.3
        # print (" Pressure Avg. Dist=%.2E Points Gap=%.2E"%(avgDist*1000, gap*1000))
        
        
        cnt = 0 
        points=[]
        ptx=[]; pty=[]; ptv=[]
        Pvs = press[indexs]
        gapCount = 0 
        for lx, ly, lv in zip(axs, avs, Pvs):
            if cnt ==0: 
                prex=lx 
                points.append([lx, -shift])
            if lx - prex > gap and points[-1][1] != -shift: 
                points.append([prex, -shift])
                points.append([lx, -shift])
                gapCount += 1 
            ptx.append(lx); pty.append(ly-shift); ptv.append(lv)
            points.append([lx, ly-shift])
            prex = lx 
            cnt += 1

        if EA > 1: 
            whilecount = 0 
            while gapCount > 20: 
                gap +=0.0005
                cnt = 0 
                points=[]
                ptx=[]; pty=[]; ptv=[]
                Pvs = press[indexs]
                gapCount = 0 
                for lx, ly, lv in zip(axs, avs, Pvs):
                    if cnt ==0: 
                        prex=lx 
                        points.append([lx, -shift])
                    if lx - prex > gap and points[-1][1] != -shift: 
                        points.append([prex, -shift])
                        points.append([lx, -shift])
                        gapCount += 1 
                    ptx.append(lx); pty.append(ly-shift); ptv.append(lv)
                    points.append([lx, ly-shift])
                    prex = lx 
                    cnt += 1
                
                whilecount += 1 
                if whilecount > 10: 
                    break 
            

        points.append([axs[-1], -shift])
        polygon = plt.Polygon(np.array(points), linewidth=size, edgecolor=color, facecolor='none', closed=False)
        self.ax.add_patch(polygon)

        return ptx, pty, ptv 

    def Plotting(self, xs=None, ys=None, pv=None, **args) :
        ## points = plt.scatter(px, py, c=pv, s=size, cmap=cmap, vmin=vmin, vmax=vmin*10, edgecolors=None, linewidths=0.0 )
        vmin = 50000; vmax = 500000
        size = 0.3 ; cmap = 'rainbow'
        adding = None 
        grid = False 
        files = False
        filter = True ## not islm data 
        contour = False 
        profile = False 
        lateralShift=0
        for key, value in args.items(): 
            if key == 'vmin': vmin = value 
            if key == 'vmax': vmax = value 
            if key == 'size': size = value 
            if key == 'cmap': cmap = value 
            if key == 'adding': adding = value 
            if key == 'grid' : grid = value 
            if key == 'legends': files = value 
            if key == 'filter': filter = value 
            if key == 'contour': contour = value 
            if key == 'profile': profile = value 
            if key == 'lateralShift': lateralShift=value 

        self.figure.clear()
        plt.clf()
        if not grid: 
            self.ax = self.figure.add_subplot(111)
            plt.axis('equal')   

            if not profile and len(xs): 
                if not contour : 
                    self.ax.scatter(xs, ys, c=pv, s=size, cmap=cmap, vmin=vmin, vmax=vmax, edgecolors=None, linewidths=0.0)
                else: 
                    self.ax.scatter(xs, ys, c=pv, s=size*4, cmap=cmap, vmin=vmin, vmax=vmax, edgecolors=None, linewidths=0.0)

                ymin = np.min(ys)
                distance = 0.02

                if adding: 
                    pts = adding[0]
                    lines = adding[1]

                    if len(pts): self.ax.scatter(pts[0], pts[1], s=size, c='black', edgecolors=None, linewidths=0.0)
                    if len(lines): 
                        for line in lines: 
                            self.ax.plot(line[0], line[1], linewidth = 0.5)

                yrange = 0.001 
                ix1 = np.where(ys>-yrange)[0]; ix2 = np.where(ys<yrange)[0]
                ix = np.intersect1d(ix1, ix2)

                xpress = xs[ix]; press = pv[ix]
                indexs = np.argsort(xpress)
                axs = xpress[indexs]
                avs = press[indexs]

                height_graph = 0.15 
                maxPress = np.max(avs) 
                self.scale =  height_graph / maxPress 
                avs *= self.scale
                

                if not profile: 
                    yprofile = np.max(avs)
                else: 
                    yprofile = np.min(avs)
                shift =yprofile - ymin + distance 
                avs -= shift 
                self.shift = shift  

                points=[]
                valueX=[]; valueY=[]; value=[]
                values = press[indexs]
                for ax, av, vv in zip(axs, avs, values): 
                    if len(points): 
                        points.append([ax, -shift])
                    points.append([ax, av])
                    points.append([ax, -shift])
                    valueX.append(ax); valueY.append(av); value.append(vv)

                polygon = plt.Polygon(np.array(points), linewidth=size*2, edgecolor='gray', facecolor='none', closed=False)
                self.ax.add_patch(polygon)
                self.sc = self.ax.scatter(valueX, valueY, s=size*3, edgecolors=None, c=None, linewidths=0.0, marker=None)
                self.figure.canvas.mpl_connect("motion_notify_event", self.hover)

                self.valuepointx = np.array(valueX)
                self.valuepointy = np.array(valueY) 
                self.valuepoint = np.array(value )
                
            elif len(xs): 
                npn, edges = profile
                y_pf =[]
                for ed in edges: 
                    ix = np.where(npn[:,0] == ed[0])[0][0]; n1 = npn[ix]
                    ix = np.where(npn[:,0] == ed[1])[0][0]; n2 = npn[ix]
                    X=[n1[2]+lateralShift, n2[2]+lateralShift]
                    Y=[n1[3], n2[3]]
                    line, = self.ax.plot(X, Y, color="black", linewidth=1.0)

                    y_pf.append(n1[3]); y_pf.append(n2[3])
                ypf = np.array(y_pf)
                distance = 0.0
                ymin = np.max(ypf)

                yrange = 0.001 
                ix1 = np.where(ys>-yrange)[0]; ix2 = np.where(ys<yrange)[0]
                ix = np.intersect1d(ix1, ix2)

                xpress = xs[ix]; press = pv[ix]
                indexs = np.argsort(xpress)
                axs = xpress[indexs]
                avs = press[indexs]

                height_graph = 0.2 
                maxPress = np.max(avs) 
                avs *= height_graph / maxPress 

                if not profile: 
                    yprofile = np.max(avs)
                else: 
                    yprofile = np.min(avs)
                shift =yprofile - ymin + distance 
                avs -= shift 
                # self.ax.scatter(axs, avs-shift, c='blue',s=size*10, edgecolors=None, linewidths=0.0)
                points=[]
                valueX=[]; valueY=[]; value=[]
                values = press[indexs]

                for ax, av, vv in zip(axs, avs, values): 
                    if len(points): 
                        points.append([ax, 0])
                    points.append([ax, av])
                    points.append([ax, 0])
                    valueX.append(ax); valueY.append(av); value.append(vv)

                polygon = plt.Polygon(np.array(points), linewidth=size*2, edgecolor='gray', facecolor='none', closed=False)
                self.ax.add_patch(polygon)
                self.sc = self.ax.scatter(valueX, valueY, s=size*3, edgecolors=None, c=None, linewidths=0.0, marker=None)
                self.figure.canvas.mpl_connect("motion_notify_event", self.hover)

                self.valuepointx = np.array(valueX)
                self.valuepointy = np.array(valueY) 
                self.valuepoint = np.array(value )

            else: 
                self.clearWindow()
        else: 
            # print ("PLOTTING MULTIPLE PLOTS")
            N = len(xs)
            subs =[]
            if N == 1: 
                subs = [self.figure.add_subplot(111)]
                PLOTS = 1 
            elif N <= 2: 
                subs = [self.figure.add_subplot(1,2,1), 
                        self.figure.add_subplot(1,2,2) 
                    ]
                PLOTS= 2
            elif N <= 4: 
                PLOTS = 4
                for i in range(4): 
                    subs.append(self.figure.add_subplot(2,2,i+1))
            elif N <=6: 
                PLOTS = 6
                for i in range(6): 
                    subs.append(self.figure.add_subplot(2,3,i+1))
            elif N <=9: 
                PLOTS = 9
                for i in range(9): 
                    subs.append(self.figure.add_subplot(3,3,i+1))
            elif N <=16: 
                PLOTS = 16
                for i in range(16): 
                    subs.append(self.figure.add_subplot(4,4,i+1))
            elif N <=20: 
                PLOTS = 20 
                for i in range(20): 
                    subs.append(self.figure.add_subplot(4,5,i+1))

            # t1 = time()
            lim = 0 
            xx=[]; yy=[]; vv=[]
            for tx, ty, tv in zip(xs, ys, pv): 
                value = np.min (tx) 
                if lim < abs(value): lim = abs(value)
                value = np.max (tx) 
                if lim < abs(value): lim = abs(value)
                value = np.min (ty) 
                if lim < abs(value): lim = abs(value)
                value = np.max (ty) 
                if lim < abs(value): lim = abs(value)
            #     if N > 4 and filter : 
            #         cnt = 0 
            #         txx=[]; tyy=[]; tvv=[]
            #         for x, y, v in zip(tx, ty, tv): 
            #             if cnt % N: 
            #                 txx.append(x); tyy.append(y); tvv.append(v)
            #             cnt += 1
            #         xx.append(np.array(txx)); yy.append(np.array(tyy)); vv.append(np.array(tvv))
            # if N > 4: 
            #     xs = xx; ys = yy; pv = vv 
            # # print (" TIME = %.2fsec"%(time()-t1))
            # print("PLOTS, XS", PLOTS, len(xs))

            lim *= 1.1
            cnt = 0 
            for cnt in range(PLOTS): 
                
                if cnt < len(xs): 
                    subs[cnt].axis('equal')
                    subs[cnt].axis('off')
                    subs[cnt].scatter([lim, lim, -lim, -lim], [-lim, lim, lim, -lim], edgecolors=None, linewidths=0.0, c='gray', s=0.01)
                    if not contour : 
                        subs[cnt].scatter(xs[cnt],  ys[cnt], c=pv[cnt], s=size, cmap=cmap, vmin=vmin, vmax=vmax, edgecolors=None, linewidths=0.0)
                    else: 
                        subs[cnt].scatter(xs[cnt],  ys[cnt], c=pv[cnt], s=size*4, cmap=cmap, vmin=vmin, vmax=vmax, edgecolors=None, linewidths=0.0)
                    subs[cnt].text(-lim, lim, files[cnt], fontsize=6 )
                else: 
                    subs[cnt].axis('off')
                    subs[cnt].scatter([lim, lim, -lim, -lim], [-lim, lim, lim, -lim], edgecolors=None, linewidths=0.0, c='gray', s=0.01)

        self.figure.canvas.mpl_connect('button_release_event', self.onReleased)
        self.figure.tight_layout()
        
        self.figure.canvas.draw()

    def addPoints(self, points=None, lines=None, lw=1, size=0.3): 
        if isinstance(points, type(None)): 
            return 
        try: 
            self.removePoints(action=False)
        except: 
            pass 

        xs = points[0]
        ys = points[1]

        self.fpcArea, = plt.plot(xs,ys, color='black', linewidth=lw, marker="o", markersize=size)
        if not isinstance(lines, type(None)): 
            cnt = 0 
            for line in lines: 
                ln, = plt.plot(line[0], line[1], color='blue', lw=0.5)
                self.fpcLine.append(ln)
                if cnt %2:  
                    ch = plt.text(line[2][0], line[2][1], line[2][2], size=8)
                    self.fpcTexts.append(ch)
                cnt += 1

        current_xlim=self.ax.get_xlim()
        current_ylim=self.ax.get_ylim()
        plt.xlim(current_xlim[0], current_xlim[1])
        plt.ylim(current_ylim[0], current_ylim[1])

        self.figure.canvas.draw_idle()

    def removePoints(self, action=True): 
        self.fpcArea.remove()
        for ln in self.fpcLine: 
            try: 
                ln.remove()
            except: 
                continue 
        for txt in self.fpcTexts: 
            txt.set_visible(False)
        if action : 
            current_xlim=self.ax.get_xlim()
            current_ylim=self.ax.get_ylim()
            plt.xlim(current_xlim[0], current_xlim[1])
            plt.ylim(current_ylim[0], current_ylim[1])
            self.figure.canvas.draw_idle()

    def clearWindow(self): 
        plt.clf()
        self.figure.clear()
        self.ax = self.figure.add_subplot(111)
        self.ax.axis('equal')
        lim = 0.25
        self.ax.scatter([lim, lim, -lim, -lim], [-lim, lim, lim, -lim], edgecolors=None, linewidths=0.0, c='gray', s=0.01)
        self.figure.canvas.mpl_connect('button_release_event', self.onReleased)
        self.figure.tight_layout()
        self.figure.canvas.draw()


