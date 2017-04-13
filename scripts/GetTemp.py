
import urllib
import time




def get_temp:
    html=urllib.request.urlopen('http://inflnr.jinr.ru/mon/draw.php?tbl=280125DE030000F3&cnt=5000')
    rhtml=html.read()
    raw=rhtml.split(b'my_time')
    tdat=raw[1].split(b'(')[1].split(b')')[0].decode()
    raw=rhtml.split(b'my_data')
    xdat=raw[1].split(b'(')[1].split(b')')[0].decode()
    ouf=open('TempReport2015.dat','a')
    ouf.write(time.ctime() + '\n')
    ouf.write('xdat=' + xdat + '\n')
    ouf.write('tdat=' + tdat + '\n')
    ouf.close()








