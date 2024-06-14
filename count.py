import os
import subprocess
zmax=12000
tmax=3600
#log100_10.txt

filelist=os.listdir()
filelist.sort()

for file in filelist:
    cnt=0
    if '.py' in file:
        continue
    if 'log' in file:
        fp=open(file,'r')
        for line in fp.readlines():
            if 'fill_holes' in line:
                cnt=cnt+1
        dz=0
        dt=0
        if len(file)==13:
            dz=int(file[3:6])
            dt=int(file[7:9])
        else:
            dz=int(file[3:5])
            dt=int(file[6:8])
        print(file[0:10]+': '+str(dz/dt+.00001)[0:4]+' : '+str(cnt/(zmax/dz)/(tmax/dt))[0:5]+' : '+str(dz*dt*cnt/100000)[0:5])

