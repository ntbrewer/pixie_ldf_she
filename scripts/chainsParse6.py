inf=open('chains/MayCalchains.txt','r')
lines=inf.readlines()

nl=[]
for i in range(0,len(lines)):
    if lines[i]=='\n':
        nl.append(i)

blocks=[]
for i in range(0,len(nl)-1):
    blocks.append(lines[nl[i]+1:nl[i+1]])

ACB=[]
for i in blocks:
    if any('I' in j for j in i) and not any('F' in j for j in i):
        ACB.append(i)

import numpy as np
from datetime import datetime
from pyqtgraph.Qt import QtCore, QtGui
import pyqtgraph.opengl as gl
import numpy as np
import pyqtgraph as pqt 

app = QtGui.QApplication([])
w = gl.GLViewWidget()
w.opts['distance'] = 20000
w.show()
w.setWindowTitle('pyqtgraph example: GLScatterPlotItem')

#g = gl.GLGridItem()
size = QtGui.QVector3D(10000,10000,1)
g = gl.GLGridItem(size=size)
g.translate(4000,0,0)
g.setSpacing(1000,1000,1)
w.addItem(g)

def get_block_i(b):
    darr = np.array([]).reshape(0,2)
#    sarr = np.array([]).reshape(0,2) //string arr
    for i in b[1:]:
        t = i.split(' ')
        t = [x for x in t if x!='']
        d2 = eval(t[2])
        d3 = eval(t[1])
        #d0 = t[0]
        #d1 = t[3]
        #d4 = datetime.timestamp(datetime.strptime(b[0].split('\t')[0],"%a %b %d %H:%M:%S %Y"))
        darr = np.r_[darr,[[d2,d3]]]
    return(darr.T)

for i in range(0,10000):
    dd = get_block_i(ACB[i])
    zd = dd[1,1:]
    yd = np.log(dd[0,1:]-dd[0,0:-1])
    ll = len(dd[0])-1
    xd = np.empty(ll)
    xd.fill(dd[1,1])
    xd[0]=1
    wd = np.empty(ll)
    od = np.empty(ll)
    wd.fill(3000000)
    od.fill(4000)
    arr = np.array([xd,yd,zd]).T
    sp1 = gl.GLScatterPlotItem(pos=arr)
    w.addItem(sp1)


