#!/bin/python3
inf=open('chains/MayCalchains.txt','r')
lines=inf.readlines()

nl=[]
for i in range(0,len(lines)):
    if lines[i]=='\n':
        nl.append(i)

blocks=[]
for i in range(0,len(nl)-1):
    blocks.append(lines[nl[i]+1:nl[i+1]])

#ACB=[]
#for i in blocks:
#    if any('I' in j for j in i) and not any('F' in j for j in i):
#        ACB.append(i)

ACB=[]
for i in blocks:
    if any('I' in j for j in i) and list('A' in j for j in i).count(True)>2:
        ACB.append(i)

#SFB=[]
#ACSFB=[]
#for i in blocks:
#    if any('F' in j for j in i) and any('A' in j for j in i):
#        ACSFB.append(i)
#    elif any('F' in j for j in i):
#        SFB.append(i)

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

def get_block_i(b):
    darr = np.array([]).reshape(0,2)
#    sarr = np.array([]).reshape(0,2) //string arr
    for i in range(1,len(b[1:])):
        t = b[i].split(' ')
        t = [x for x in t if x!='']
        d2 = eval(t[2])
        d3 = eval(t[1])
        if d3 > 10000 and i == 2:
           return(None)
        else:
           darr = np.r_[darr,[[d2,d3]]]
        #d0 = t[0]
        #d1 = t[3]
        #d4 = datetime.timestamp(datetime.strptime(b[0].split('\t')[0],"%a %b %d %H:%M:%S %Y"))
        #darr = np.r_[darr,[[d2,d3]]]
    return(darr.T)

#darr[0] = 'I/A/F'
#darr[1] = 'MBVE'
#darr[2] = 'dt'
#darr[3] = 'E'
#darr[4] = 'wt'


#darr=[[[]],[[]],[[]],[[]],[[]]]
#for b in ACSFB:
#    get_block(b,darr)

#darr=[[[]],[[]],[[]],[[]],[[]]]
#for b in ACB:
#    get_block(b,darr)

def cm(E):
    if 9>E>6:
        Pr=format(int(255*np.sqrt(1-(E-7.5)**2/2.25)),'x')
    else:
        Pr='44'
    if 10>E>7:
        Pg=format(int(255*np.sqrt(1-(E-8.5)**2/2.25)),'x')
    else:
        Pg='44'
    if 11>E>8:
        Pb=format(int(255*np.sqrt(1-(E-9.5)**2/2.25)),'x')
    else:
        Pb='44'
    if len(Pr)<2:
       Pr='0'+Pr
    if len(Pg)<2:
       Pg='0'+Pg
    if len(Pb)<2:
       Pb='0'+Pb
    return('#'+Pr+Pg+Pb)

fig = plt.figure()
#ax = fig.add_subplot(111, projection='3d')
ax1 = fig.add_subplot(221)
ax2 = fig.add_subplot(223)
ax3 = fig.add_subplot(122)
plt.ion()
alp=0.1
for b in ACB:
    dd = get_block_i(b)
    if dd == None:
        continue
    zd = dd[1,1:]
    yd = np.log10(dd[0,1:]-dd[0,0:-1])
    ll = len(dd[0])-1
    xd = np.empty(ll)
    xd.fill(dd[1,1])
    #xd[0]=1
    cmap = cm(dd[1,1]/1000)
    #wd = np.empty(ll)
    #od = np.empty(ll)
    #wd.fill(-5)
    #od.fill(4000)
    try:
       ax1.scatter(yd[0],xd[0],c=cmap,alpha=alp,edgecolors='none')
       ax3.scatter(xd[1:],zd[1:],c=cmap,alpha=alp,edgecolors='none')
       ax2.scatter(yd[1:],zd[1:],c=cmap,alpha=alp,edgecolors='none')
    except ValueError:
        print(b,xd,dd)
        break

ax1.set_ylim(7000,10000)
ax2.set_ylim(4000,10000)
ax3.set_ylim(4000,10000)
ax3.set_xlim(7000,10000)

"""
    try:
       ax.scatter(xd,yd,zd,c=cmap,marker='o',alpha=alp)
    #ax.plot(xd[1:],yd[1:],zd[1:],color='b')
    #ax.plot(xd[:],yd[:],zd[:],color='b')
       ax.scatter(xd,yd,1000,c=cmap,alpha=alp)
    #ax.plot(xd,yd,1000,color='0.5')
       ax.scatter(xd[1:],wd[1:],zd[1:],c=cmap,alpha=alp)
       ax.scatter(od,yd,zd,c=cmap,alpha=alp)
    #ax.plot(od,yd,zd,color='0.5')
    #plt.show()
    #input('ok?')
    except ValueError:
        print(b,xd,dd)
        break
"""
#ax.set_xscale('symlog')
#ax.set_yscale('symlog')
#ax.view_init(elev=60,azim=20)
#ax.set_zlim(1000,10000)
#ax.set_xlim(4000,12000)
#ax.set_ylim(-5,10)
plt.show()
#ax.set_zscale('symlog')
"""                
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
for i in range(1,round(len(darr[3])/100)):
    xd=[(darr[3][i][1])]
    yd=[(1)]
    zd=[(darr[2][i][1])]
    wd=[3000000]
    od=[4000]
    if len(darr[3][i])>2:
        for j in range(2,len(darr[3][i])):
            xd.append((darr[3][i][1]))
            yd.append((darr[3][i][j]))
            zd.append((darr[2][i][j]))
            wd.append(3000000)
            od.append(1000)
    ax.scatter(xd,zd,yd,c='k',marker='o')
    ax.plot(xd,zd,yd,color='b')
    ax.scatter(xd,zd,-1000,c='r')
    ax.scatter(xd,wd,yd,c='r')
    ax.scatter(od,zd,yd,c='r')
--
"""
#plt.show()
