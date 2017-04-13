
inf=open('chains.txt','r')
lines=inf.readlines()

nl=[]
for i in range(0,len(lines)):
    if lines[i]=='\n':
        nl.append(i)

blocks=[]
for i in range(0,len(nl)-1):
    blocks.append(lines[nl[i]+1:nl[i+1]])
vb=[]
for i in blocks:
    if any('F ' in j for j in i):
        vb.append(i)

ouf=open('SFw240PuCa_tsort.txt','w')
for i in vb:
    ouf.writelines(i)

ouf.close()

svb=[]
for i in range(0,len(vb)):
    if len(vb[i][2].split('-'))<2:
        b=vb[i][-1].split('      ')
        svb.append([eval(b[1][0:7]),i])
svb.sort()

ouf=open('SFw240PuCa_SFEsort.txt','w')
for i in svb:
    ouf.writelines(vb[i[1]])

ouf.close()

vb=[]
for i in blocks:
    if any('in ' in j for j in i):
        vb.append(i)


ouf=open('withIN240PuCa.txt','w')
for i in vb:
    ouf.writelines(i)

ouf.close()
