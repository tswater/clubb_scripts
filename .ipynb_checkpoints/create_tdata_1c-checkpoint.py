import numpy as np
import netCDF4 as nc
import datetime
import os
import pickle

clubb_dir = '/home/tsw35/tyche/clubb/'

def cts_to_string(cts):
    out=[]
    for ct in cts:
        sp=str(ct).split('0')
        print(ct)
        if(len(sp)<4):
            if 'e-04' in str(ct):
                out.append('cr20ct0.000'+str(ct)[0]+'cf0.8')
            elif 'e-05' in str(ct):
                out.append('cr20ct'+str(ct)+'cf0.8')
        else:
            out.append('cr20ct0.000'+sp[4]+'cf0.8')
    return out

ct_n =np.linspace(.00001,.0005,50,dtype='f4')
ct_s= cts_to_string(ct_n) #filename strings

daylist2=os.listdir(clubb_dir+'circ_tune/')
daylist2.sort()

tests=['dpv_lwp','dmv_lwp','dpt_lwp','dt1_lwp','dmv_tke','mdt_lwp']
tdata={}
tdata['ct_pct']={}
tdata['tests']={}
tdata['vars']={}
tdata['days']=[]

fp1c_zt=nc.Dataset('/home/tsw35/tyche/clubb/circ_tune/20150606/cr20ct0.00031cf0.8/k_2/agg_outzt.nc','r')

tdata['altitude']=fp1c_zt['altitude'][:]
for test in tests:
    tdata['ct_pct'][test]=[]
    tdata['tests'][test]=[]
nf=0
for folder in daylist2:
    print(folder,end='...',flush=True)
    if nf == 9:
        print()
        nf=-1
    nf=nf+1
    tdata['days'].append(folder)
    
    fp2c_zt=nc.Dataset(clubb_dir+'sgp_1c/sgp_'+folder+'/k_1/c_1/output/arm_zt.nc','r')
    fp2c_zm=nc.Dataset(clubb_dir+'sgp_1c/sgp_'+folder+'/k_1/c_1/output/arm_zm.nc','r')
    fp2c_sfc=nc.Dataset(clubb_dir+'sgp_1c/sgp_'+folder+'/k_1/c_1/output/arm_sfc.nc','r')
    
    lwp_2c=fp2c_sfc['lwp'][:,0,0,0]
    tke_2c=np.sum(fp2c_zm['wp2'][:,:,0,0]+fp2c_zm['up2'][:,:,0,0]+fp2c_zm['vp2'][:,:,0,0],axis=1)
    
    #### TESTS ####
    dpv_lwp=np.zeros((len(ct_s),))
    dmv_lwp=np.zeros((len(ct_s),))
    dpt_lwp=np.zeros((len(ct_s),))
    dt1_lwp=np.zeros((len(ct_s),))
    #dpv_blh=np.zeros((len(ct_s),))
    #dmv_blh=np.zeros((len(ct_s),))
    #dpt_blh=np.zeros((len(ct_s),))
    dtke   =np.zeros((len(ct_s),))
    mdt_lwp=np.zeros((len(ct_s),)) # mean difference in tile LWP
    for i in range(len(ct_s)):
        f=ct_s[i]
        fpzt=nc.Dataset(clubb_dir+'circ_tune/'+folder+'/'+f+'/k_2/agg_outzt.nc','r')
        fpzm=nc.Dataset(clubb_dir+'circ_tune/'+folder+'/'+f+'/k_2/agg_outzm.nc','r')
        fpsfc=nc.Dataset(clubb_dir+'circ_tune/'+folder+'/'+f+'/k_2/agg_outsfc.nc','r')
        fpsfc_c1 = nc.Dataset(clubb_dir+'circ_tune/'+folder+'/'+f+'/k_2/c_1/output/arm_sfc.nc','r')
        fpsfc_c2 = nc.Dataset(clubb_dir+'circ_tune/'+folder+'/'+f+'/k_2/c_2/output/arm_sfc.nc','r')
        fpclst   = nc.Dataset(clubb_dir+'circ_tune/'+folder+'/'+f+'/k_2/clusters.nc','r')
        frac=fpclst['frac'][:]
        
        lwp=fpsfc['lwp'][:,0,0,0]
        dpv_lwp[i]=np.max(lwp)-np.max(lwp_2c)
        dmv_lwp[i]=np.mean(lwp)-np.mean(lwp_2c)
        dpt_lwp[i]=np.where(lwp==np.max(lwp))[0][0]-np.where(lwp_2c==np.max(lwp_2c))[0][0]
        dt1_lwp[i]=np.where(lwp>.1*np.max(lwp))[0][0]-np.where(lwp_2c>.1*np.max(lwp_2c))[0][0]
        tke_cpl   =np.sum(fpzm['wp2'][:,:,0,0]+fpzm['up2'][:,:,0,0]+fpzm['vp2'][:,:,0,0],axis=1)
        try:
            dtke[i]   =np.mean(tke_cpl-tke_2c)
        except:
            dtke[i]   =0
        lwp_c1=fpsfc_c1['lwp'][:,0,0,0]
        lwp_c2=fpsfc_c2['lwp'][:,0,0,0]
        dlwp_cpl  =np.abs(lwp_c2-lwp_c1)/((lwp_c2*frac[1]+lwp_c1*frac[0])/2)
        try:
            mdt_lwp[i]=np.mean(dlwp_cpl)
        except:
            mdt_lwp[i]=0
    tdata['tests']['dpv_lwp'].append(dpv_lwp)
    tdata['tests']['dmv_lwp'].append(dmv_lwp)
    tdata['tests']['dpt_lwp'].append(dpt_lwp)
    tdata['tests']['dt1_lwp'].append(dt1_lwp)
    tdata['tests']['dmv_tke'].append(dtke)
    tdata['tests']['mdt_lwp'].append(mdt_lwp)

tdata['ntests']={}
for test in tests:
    tdata['ntests'][test]=[]
nf=0
i=0
print()
print('FIRST STEP DONE')
print()

for folder in daylist2:
    print(folder,end='...',flush=True)
    if nf == 9:
        print()
        nf=-1
    nf=nf+1
    
    fp2c_zt=nc.Dataset(clubb_dir+'sgp_1c/sgp_'+folder+'/k_1/c_1/output/arm_zt.nc','r')
    fp2c_zm=nc.Dataset(clubb_dir+'sgp_1c/sgp_'+folder+'/k_1/c_1/output/arm_zm.nc','r')
    fp2c_sfc=nc.Dataset(clubb_dir+'sgp_1c/sgp_'+folder+'/k_1/c_1/output/arm_sfc.nc','r')
    
    lwp_2c=fp2c_sfc['lwp'][:,0,0,0]
    
    tke_2c    =np.sum(fp2c_zm['wp2'][:,:,0,0]+fp2c_zm['up2'][:,:,0,0]+fp2c_zm['vp2'][:,:,0,0],axis=1)
    
    tdata['ntests']['dpv_lwp'].append(tdata['tests']['dpv_lwp'][i]/np.max(lwp_2c))
    tdata['ntests']['dmv_lwp'].append(tdata['tests']['dmv_lwp'][i]/np.mean(lwp_2c))
    tdata['ntests']['dpt_lwp'].append(tdata['tests']['dpt_lwp'][i]/100)
    tdata['ntests']['dt1_lwp'].append(tdata['tests']['dt1_lwp'][i]/100)
    tdata['ntests']['dmv_tke'].append(tdata['tests']['dmv_tke'][i]/np.mean(tke_2c))
    tdata['ntests']['mdt_lwp'].append(tdata['tests']['mdt_lwp'][i])
    
    i=i+1

print()
print('STEP 2 DONE')
print()

ctpcts=[5,10,50,100]
for pct in ctpcts:
    print(pct)
    name='ct_pct'+str(pct)
    tdata[name]={}
    for test in tests:
        tdata[name][test]=[]
        for i in range(len(daylist2)):
            wheres=np.where(tdata['ntests'][test][i]>=(pct/100))[0]
            if len(wheres)==0:
                tdata[name][test].append(float('NaN'))
            else:
                tdata[name][test].append(wheres[0])

pickle.dump(tdata,open("tdata_base1c.p","wb"))
