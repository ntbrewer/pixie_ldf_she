#/bin/python3
inf=open('chains.22Feb.txt','r')
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

SFB=[]
ACSFB=[]
for i in blocks:
    if any('F' in j for j in i) and any('A' in j for j in i):
        ACSFB.append(i)
    elif any('F' in j for j in i):
        SFB.append(i)

import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from datetime import datetime

def get_block(b,darr):
    for i in range(0,5):
        darr[i].append([])
#    pxl=[eval(b[0].split(' ')[7]),eval(b[0].split(' ')[10])]
    for i in b[1:]:
        t = i.split(' ')
        t = [x for x in t if x!='']
        darr[2][-1].append(eval(t[2]))
        darr[3][-1].append(eval(t[1]))
        darr[0][-1].append(t[0])
        darr[1][-1].append(t[3])
#        darr[0][-1].append(pxl[0])
#        darr[1][-1].append(pxl[1])
        darr[4][-1].append(datetime.timestamp(datetime.strptime(b[0].split('\t')[0],"%a %b %d %H:%M:%S %Y")))       

#darr[0] = 'X'
#darr[1] = 'Y'
#darr[2] = 'dt'
#darr[3] = 'E'
#darr[4] = 'wt'


darr=[[[]],[[]],[[]],[[]],[[]]]
for b in ACSFB:
    get_block(b,darr)

fig = plt.figure()
ax = fig.add_subplot(111)
for i in range(1,len(darr[2])):
    ax.scatter(darr[3][i],darr[2][i])

plt.show()

#ax = fig.add_subplot(111, projection='3d')
#ax.scatter(darr[2],darr[1],darr[0],c=darr[3],depthshade=False)
#ax.plot_wireframe(darr[2],darr[1],darr[0])

#plt.show()

                
#fig = plt.figure()
#ax = fig.add_subplot(111, projection='3d')
#for i in range(0,len(xd)):
#    ax.scatter(zd[i],yd[i],xd[i],c=ed[i],s=50,depthshade=False)
#    ax.plot_wireframe(zd[i],yd[i],xd[i])

#plt.show()
