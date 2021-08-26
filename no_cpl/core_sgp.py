import os 
import shutil
import subprocess
import datetime

sgp_dir = '/stor/soteria/hydro/shared/lasso_for_tyler/'
custom_k = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,22,24,\
            26,28,30,32,34,36,38,40,45,50]
dx = 250
nx = 401

def strptime(dtstr):
    return datetime.datetime(int(dtstr[0:4]),int(dtstr[5:7]),int(dtstr[8:10]),\
           int(dtstr[11:13]),int(dtstr[14:16]))
i = 0
for folder in os.listdir(sgp_dir):
    # get folder start date, time end date, and set dx and nx
    start_time = datetime.datetime(2022,9,9,9,9)
    end_time   = datetime.datetime(2000,2,2,2,2)
    for file in os.listdir(sgp_dir+'/'+folder):
        if len(file)<25:
            continue
        time = strptime(file[-16:])
        if time < start_time:
            start_time = time
        elif time > end_time:
            end_time = time

    # set customk to 50
    fp = open('setup_parallel.py','r')
    lines = fp.readlines()
    lines[35] = "stdate  = '"+start_time.strftime('%Y-%m-%dT%H:%M:%S')+"'\n"
    lines[36] = "enddate = '"+end_time.strftime('%Y-%m-%dT%H:%M:%S')+"'\n"
    lines[33] = 'nx      = '+str(nx)+'\n'
    lines[34] = 'dx      = '+str(dx)+'\n'
    lines[31] = 'k       = '+str(50)+'\n'
    lines[24] = "sfc_dir = '"+sgp_dir+folder+"/'\n"
    lines[30] = "dirname = 'sgp_"+str(start_time.year)+str(start_time.month)+str(start_time.day)+"'\n"
    lines[43] = "kcust   = "+str(custom_k)
    fp.close()
    fp = open('sgprun_'+str(i)+'.py','w')
    fp.writelines(lines)
    fp.close()
    i = i+1

