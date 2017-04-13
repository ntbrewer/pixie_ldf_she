#/bin/python3
import numpy as np
file=open('gui_cal_7fit.txt','r')
formatlist=file.readlines()
file.close()
for i in range(0,len(formatlist)):
    formatlist[i]=eval(formatlist[i])

file=open('Config.Jun_15.xml','w')
file.write('----' + '\n')
file.close()

maplist=(64,63,62,61,60,59,58,57,56,55,54,53,52,51,50,49,48,
47,46,45,44,43,42,41,40,39,38,37,36,35,34,33,16,15,14,13,12,
11,10,9,8,7,6,5,4,3,2,1,32,31,30,29,28,27,26,25,24,23,22,21,
20,19,18,17,128,127,126,125,124,123,122,121,120,119,118,117,
116,115,114,113,112,111,110,109,108,107,106,105,104,103,102,
101,100,99,98,97,85,95,94,93,92,91,90,89,88,87,86,96,84,83,82
,81,80,79,78,77,76,75,74,73,72,71,70,69,68,67,66,65)

for i in formatlist:
    i0=int(i[0])-3100
    if i0 < 47:
       file=open('Config.Jun_15.xml','a')
       if np.mod(i0,16) ==0:
            file.write('</Module> ' +'\n' + '<Module number="' +repr(round(i0/16))+ '">' +  '\n')
       file.write('<Channel number= " ' + repr(np.mod(i0,16)) + ' " location= " ' + repr(48-i0) + ' " subtype="dssd_front" type="dssd_front"> ' + '\n')
       file.write('<Calibration max="32000" model="linear">' + repr(i[2]) + ' ' + repr(i[1])+ ' ' + '</Calibration> </Channel>' + '\n')
       if i0==46:
           file.write('<Channel number= " ' + repr(np.mod(i0+1,16)) + ' " location= " ' + repr(48-i0-1) + ' " subtype="dssd_front" type="dssd_front"> ' + '\n')
           file.write('<Calibration max="32000" model="linear">' + repr(i[2]) + ' ' + repr(i[1]) + ' ' + '</Calibration> </Channel>' + '\n')
       file.close()
       

    if i0 > 47 and i0 < 176:
       file=open('Config.Jun_15.xml','a')
       if np.mod(i0,16) ==0:
            file.write('</Module> ' +'\n' + '<Module number="' +repr(round((i0)/16))+ '">' +  '\n')
       file.write('<Channel number= " ' + repr(np.mod((i0),16)) + ' " location= " ' + repr(maplist[i0-48]) + ' " subtype="dssd_back" type="dssd_back"> ' + '\n')
       file.write('<Calibration max="32000" model="linear">' + repr(i[2]) + ' ' + repr(i[1]) + ' ' + '</Calibration> </Channel>' + '\n')
       file.close()

    if i0 > 176:
       file=open('Config.Jun_15.xml','a')
       if np.mod(i0-1,16) ==0:
            file.write('</Module> ' +'\n' + '<Module number="' +repr(round((i0)/16))+ '">' +  '\n')
       file.write('<Channel number= " ' + repr(np.mod((i0),16)) + ' " location= " ' + repr(i0-177) + ' " subtype="si" type="si"> ' + '\n')
       file.write('<Calibration max="32000" model="linear">' + repr(i[2]) + ' ' + repr(i[1]) + ' ' + '</Calibration> </Channel>' + '\n')
       file.close() 

file=open('Config.Jun_15.xml','a')
file.write('</Module> ' +'\n')
file.close()
