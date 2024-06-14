import numpy as np
import netCDF4 as nc
import datetime
import os
import pickle
clubb_dir= '/home/tsw35/tyche/clubb/'
tdir = clubb_dir+'tune_feb23/'

##### PULL DATA #####
peak_val = {}
mean_val = {}
h_c_lwpm = {}
mean_rtm = {} # mean first 6 km from 12-7
h_c_rtm = {} # ratio of mean first 6 km from 12-7
rtp2 = {} # mean rtp above 1km
maxzl = {} # 95th percentile
maxzh = {} # 95th percentile
urmu = {} # mean ur>0
urmx = {} # max ur
crs = {}
cc1s = {}
incs = {}
ccc = {}
lwp = {}
urfc= {}
thlp2 = {}
wp2 = {}
met={} #circ height method


try:
    for day in os.listdir(tdir):
        peak_val[day]=[]
        mean_val[day]=[]
        h_c_lwpm[day]=[]
        mean_rtm[day]=[]
        h_c_rtm[day]=[]
        crs[day]=[]
        incs[day]=[] #this is actually k number of columns
        cc1s[day]=[]
        ccc[day]=[]
        lwp[day]=[]
        rtp2[day]=[]
        thlp2[day]=[]
        maxzl[day]=[]
        maxzh[day]=[]
        urmu[day]=[]
        urmx[day]=[]
        urfc[day]=[]
        wp2[day]=[]
        met[day]=[]

        print('*',end='',flush=True)
        for file in os.listdir(tdir+day):
            if 'k1' in file:
                try:
                    fs=nc.Dataset(tdir+day+'/'+file+'/k_1/c_1/output/arm_sfc.nc','r')
                    fm=nc.Dataset(tdir+day+'/'+file+'/k_1/c_1/output/arm_zm.nc','r')
                    ft=nc.Dataset(tdir+day+'/'+file+'/k_1/c_1/output/arm_zt.nc','r')
                except:
                    print('ERROR: '+str(day)+':'+file)
                    continue 
                crs[day].append(0)
                cc1s[day].append(0)
                incs[day].append(1)

                lwp[day].append(fs['lwp'][:,0,0,0])
                peak_val[day].append(np.max(lwp[day][-1]))
                mean_val[day].append(np.mean(lwp[day][-1]))
                h_c_lwpm[day].append(float('nan'))
                mean_rtm[day].append(np.mean(ft['rtm'][420:840,0:150,0,0]))
                h_c_rtm[day].append(float('nan'))
                rtp2[day].append(np.mean(fmcp['rtp2'][:,25:150,0,0]))
                thlp2[day].append(np.mean(fmcp['thlp2'][:,25:150,0,0]))
                wp2[day].append(np.mean(fmcp['wp2'][:,0:150,0,0]))
                
                maxzl[day].append(float('nan'))
                maxzh[day].append(float('nan'))
                urmu[day].append(float('nan'))
                urmx[day].append(float('nan'))
                urfc[day].append(float('nan'))
                ccc[day].append(float('nan'))
                met[day].append('NAN')
            else:
                try:
                    fscp=nc.Dataset(tdir+day+'/'+file+'/k_2/agg_outsfc.nc','r')
                    fsch=nc.Dataset(tdir+day+'/'+file+'/k_2/c_2/output/arm_sfc.nc','r')
                    fscc=nc.Dataset(tdir+day+'/'+file+'/k_2/c_1/output/arm_sfc.nc','r')
                    fmcp=nc.Dataset(tdir+day+'/'+file+'/k_2/agg_outzm.nc','r')
                    ftcp=nc.Dataset(tdir+day+'/'+file+'/k_2/agg_outzt.nc','r')
                    ftch=nc.Dataset(tdir+day+'/'+file+'/k_2/c_2/output/arm_zt.nc','r')
                    ftcc=nc.Dataset(tdir+day+'/'+file+'/k_2/c_1/output/arm_zt.nc','r')
                    fccp=nc.Dataset(tdir+day+'/'+file+'/k_2/clusters.nc','r')
                except:
                    print('ERROR: '+str(day)+':'+file)
                    continue

                lwp[day].append(fscp['lwp'][:,0,0,0])
                peak_val[day].append(np.max(lwp[day][-1]))
                mean_val[day].append(np.mean(lwp[day][-1]))
                h_c_lwpm[day].append(np.mean(fsch['lwp'][:,0,0,0])/np.mean(fscc['lwp'][:,0,0,0]))
                mean_rtm[day].append(np.mean(ftcp['rtm'][420:840,0:150,0,0]))
                h_c_rtm[day].append(np.mean(ftch['rtm'][420:840,0:150,0,0])/np.mean(ftcc['rtm'][420:840,0:150,0,0]))
            
                rtp2[day].append(np.mean(fmcp['rtp2'][:,25:150,0,0]))
                thlp2[day].append(np.mean(fmcp['thlp2'][:,25:150,0,0]))
                wp2[day].append(np.mean(fmcp['wp2'][:,0:150,0,0]))
                maxzl[day].append(np.percentile(fccp['z_circl'][:],95))
                maxzh[day].append(np.percentile(fccp['z_circh'][:],95))
                ur = fccp['u_r'][:]
                urmu[day].append(np.mean(ur[ur>0]))
                urmx[day].append(np.percentile(ur,98))
                urfc[day].append(np.sum(ur>0)/len(ur))

                ns=file.replace('feb23_cr','g').replace('ccc','g').replace('cc','g').replace('/','').replace('mm','g').replace('i','g').replace('k2','').split('g')
                crs[day].append(float(ns[1]))
                cc1s[day].append(float(ns[2]))
                incs[day].append(2)
                ccc[day].append(0)
                met[day].append(ns[3])
            print('.',end='',flush=True)
        peak_val[day]=np.array(peak_val[day])
        mean_val[day]=np.array(mean_val[day])
        h_c_lwpm[day]=np.array(h_c_lwpm[day])
        mean_rtm[day]=np.array(mean_rtm[day])
        h_c_rtm[day]=np.array(h_c_rtm[day])
        crs[day]=np.array(crs[day])
        incs[day]=np.array(incs[day])
        cc1s[day]=np.array(cc1s[day])
        ccc[day]=np.array(ccc[day])
        lwp[day]=np.array(lwp[day])
        rtp2[day]=np.array(rtp2[day])
        maxzl[day]=np.array(maxzl[day])
        maxzh[day]=np.array(maxzh[day])
        urmu[day]=np.array(urmu[day])
        urmx[day]=np.array(urmx[day])
        urfc[day]=np.array(urfc[day])
        thlp2[day]=np.array(thlp2[day])
        wp2[day]=np.array(wp2[day])
        met[day]=np.array(met[day])

except Exception as e:
    print(e)
    print(day)
    print(file)
    #raise e

fp=open('tunefeb23.p','wb')
pickle.dump(mean_val,fp)
pickle.dump(peak_val,fp)
pickle.dump(h_c_lwpm,fp)
pickle.dump(mean_rtm,fp)
pickle.dump(h_c_rtm,fp)
pickle.dump(crs,fp)
pickle.dump(incs,fp)
pickle.dump(cc1s,fp)
pickle.dump(lwp,fp)
pickle.dump(rtp2,fp)
pickle.dump(maxzl,fp)
pickle.dump(maxzh,fp)
pickle.dump(urmu,fp)
pickle.dump(urmx,fp)
pickle.dump(urfc,fp)
pickle.dump(ccc,fp)
pickle.dump(thlp2,fp)
pickle.dump(wp2,fp)
pickle.dump(met,fp)
fp.close()
