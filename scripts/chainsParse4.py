#/bin/python3
inf=open('chains.II.Feb16.txt','r')
lines=inf.readlines()

nl=[]
for i in range(0,len(lines)):
    if lines[i]=='\n':
        nl.append(i)

blocks=[]
for i in range(0,len(nl)-1):
    blocks.append(lines[nl[i]+1:nl[i+1]])
"""
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
"""
SFB=[]
for i in blocks:
    if len(i) >3 and any('F ' in j for j in i):
        SFB.append([i[0]])
        SFB[-1].append(i[-3])
        SFB[-1].append(i[-2])
        SFB[-1].append(i[-1])

SFB2=[]
for i in SFB:
    if not any('A ' in j for j in i):
        SFB2.append(i)
print(len(SFB),len(SFB2))
import numpy as np
#from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
#from datetime import datetime

def get_block(b,darr):

    #print(b[0].split(' ')[7])
    #pxl=[eval(b[0].split(' ')[7]),eval(b[0].split(' ')[10])]
    for i in range(0,3):
        darr[i].append([])
    for idx in range(2,4):
        t = b[idx].split(' ')
        t = [x for x in t if x!='']
        pt = b[idx-1].split(' ')
        pt = [x for x in pt if x!='']
#        print(len(t),t[2])
#        if len(t)>4:
        try:
            dt1=eval(t[3])
        except:
            dt1=eval(t[2])
        try:
            dt0=eval(pt[3])
        except:
            dt0=eval(pt[2])
#        else:
#            dt1=eval(t[2])
#        if len(pt)>4:
#            dt0=eval(pt[3])
#        else:
#            dt0=eval(pt[2])
        darr[0][-1].append(dt1-dt0)
        darr[1][-1].append(eval(t[1]))
        darr[2][-1].append(t[0])
     #   darr[0][-1].append(pxl[0])
     #   darr[1][-1].append(pxl[1])
     #   darr[4][-1].append(datetime.timestamp(datetime.strptime(b[0].split('\t')[0],"%a %b %d %H:%M:%S %Y")))       

#darr[0] = 'dt'
#darr[2] = 'E'
#darr[3] = 'I/F'



#darr=[[[]],[[]],[[]],[[]],[[]]]
darr=[[[]],[[]],[[]],[[]]]
for b in SFB2:
    get_block(b,darr)

#print(len(darr[2]))
fig = plt.figure()
ax = fig.add_subplot(311,xlabel='dt1 (ms)', ylabel='dt0 (ms)')
ax1= fig.add_subplot(312, xlabel='dt1 (ms)' ,ylabel='E (keV)')
ax2= fig.add_subplot(313, xlabel='dt0 (ms)' ,ylabel='E (keV)')
for i in range(1,len(darr[1])-1):
    ax.scatter(darr[0][i][1],darr[0][i][0])
    ax1.scatter(darr[0][i][1],darr[1][i][1])
    ax2.scatter(darr[0][i][0],darr[1][i][1])

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
