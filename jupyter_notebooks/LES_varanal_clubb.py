# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.15.2
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# %%
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc
import datetime
import os
import seaborn as sns
from matplotlib.animation import FuncAnimation
from IPython.display import HTML
sns.set_theme()

# %%
days=    [20160625,20160716,20160719,20160720,20170609,
          20170626,20170627,20170629,20170705,20170709,
          20170712,20170716,20170717,20170719,20170720,
          20170728,20170826,20170922,20170923,20170924,
          20180522,20180530,20180618,20180619,20180704,
          20180705,20180523,20180707,20180709,20180710,
          20180711,20180712,20180809,20180811,20180916,
          20180917,20190707,20190709,20190714,20190804,
          20190805]


# %%
def vpt(T,r,p):
    theta=T*(100000/(p))**(2/7)
    thv=theta*(1+.61*r)
    return thv


# %%

# %%
# TIMES
# LES is minutes since 2016-01-01 with 10 min dt
# Varanal is seconds since 2012-5-1 with hourly dt
# clubb time is wrong; 0 is 5am local (10am utc) and dt is 1 min
def index_time(date,clb_st=5,les_st=7,var_st=datetime.datetime(2012,5,1,0,0)):
    clb_i=(date.hour-clb_st)*60+date.minute
    if clb_i<0:
        clb_i=clb_i+60*24
    les_i=int((date.hour-les_st)*6+date.minute/10)
    var_i=int((date-var_st).total_seconds()/60/60)+5
    return clb_i,les_i,var_i


# %%

# %%
les_1c ='/home/tsw35/tyche/data/LES_1C/'
clubb_hf='/home/tsw35/tyche/clubb/sgp_cpl_d4/'
clubb_2c='/home/tsw35/tyche/clubb/sgp_2c_d/'
clubb_1c='/home/tsw35/tyche/clubb/sgp_1c_d/'
clubb_hf2='/home/tsw35/tyche/clubb/sgp_2c_hf/'
clubb_nf2='/home/tsw35/tyche/clubb/sgp_2c_nf/'
clubb_f2='/home/tsw35/tyche/clubb/sgp_nocpl_i/'

varanal='/home/tsw35/soteria/clubb/data/sgp60varanarap_2012-2019.nc'

# %%
fpv=nc.Dataset(varanal,'r')
var_p=fpv['lev'][:]

# %%

# %%
i=0
j=0
plt.figure(figsize=(18,6))
altl=np.linspace(0,5000,166)
fp1=nc.Dataset(clubb_1c+'sgp_'+str(days[0])+'/k_1/agg_outzt.nc','r')
altc=fp1['altitude'][0:166]
for day in days:
    if j==6:
        plt.figure(figsize=(18,6))
        j=0
    dys=str(day)
    fp1=nc.Dataset(clubb_1c+'sgp_'+str(day)+'/k_1/c_1/output/arm_zt.nc','r')
    fp2=nc.Dataset(clubb_2c+'sgp_'+str(day)+'/k_2/agg_outzt.nc','r')
    fp3=nc.Dataset(clubb_hf+'sgp_'+str(day)+'/k_2/agg_outzt.nc','r')
    fplm=nc.Dataset(les_1c+'trimfr2_'+str(day)+'_01.nc','r')
    fplt=nc.Dataset(les_1c+'trimfr2_'+str(day)+'_00.nc','r')
    alt2=fp1['altitude'][0:166]
    plt.subplot(1,6,j+1)
    hr=20
    st_time= datetime.datetime(int(dys[0:4]),int(dys[4:6]),int(dys[6:8]),hr)
    clb_t,les_t,var_t=index_time(st_time)
    #plt.plot(var_p)
    altv=np.flip(np.interp(np.flip(var_p*100),np.flip(fp2['p_in_Pa'][clb_t,0:166,0,0]),np.flip(altc)))
    plt.plot(fpv['u'][var_t,:-20],altv[0:-20],'b-')
    plt.plot(fpv['v'][var_t,:-20],altv[0:-20],'r-')
    plt.plot(fplm['u'][les_t,0:166],altl,'b--')
    plt.plot(fplm['v'][les_t,0:166],altl,'r--')
    j=j+1
    print('.',end='',flush=True)

# %%
i=0
j=0
plt.figure(figsize=(18,6))
altl=np.linspace(0,5000,166)
fp1=nc.Dataset(clubb_1c+'sgp_'+str(days[0])+'/k_1/agg_outzt.nc','r')
altc=fp1['altitude'][0:166]
for day in days:
    if j==6:
        plt.figure(figsize=(18,6))
        j=0
    dys=str(day)
    fp1=nc.Dataset(clubb_1c+'sgp_'+str(day)+'/k_1/c_1/output/arm_zt.nc','r')
    fp2=nc.Dataset(clubb_2c+'sgp_'+str(day)+'/k_2/agg_outzt.nc','r')
    fp3=nc.Dataset(clubb_hf+'sgp_'+str(day)+'/k_2/agg_outzt.nc','r')
    fplm=nc.Dataset(les_1c+'trimfr2_'+str(day)+'_01.nc','r')
    fplt=nc.Dataset(les_1c+'trimfr2_'+str(day)+'_00.nc','r')
    alt2=fp1['altitude'][0:166]
    plt.subplot(1,6,j+1)
    hr=18
    st_time= datetime.datetime(int(dys[0:4]),int(dys[4:6]),int(dys[6:8]),hr)
    clb_t,les_t,var_t=index_time(st_time)
    
    altv=np.flip(np.interp(np.flip(var_p*100),np.flip(fp2['p_in_Pa'][clb_t,0:166,0,0]),np.flip(altc)))
    vpt_v=vpt(fpv['T'][var_t,:],fpv['q'][var_t,:]/1000,var_p*100)
    
    plt.title(day)

    
    vpt_1  =fp1['thvm'][clb_t,0:166,0,0]
    vpt_nf =fp2['thvm'][clb_t,0:166,0,0]
    vpt_hf =fp3['thvm'][clb_t,0:166,0,0]
    vpt_l  =fplm['thv'][les_t,0:166]
    
    plt.plot(vpt_1-vpt_1[0],altc,':')
    plt.plot(vpt_nf-vpt_nf[0],altc,':')
    plt.plot(vpt_hf-vpt_hf[0],altc,'--')
    plt.plot(vpt_l-vpt_l[0],altl)
    #plt.plot(fplt['thv'][les_t,0:166],altl)
    plt.plot(vpt_v-vpt_v[5],altv,'o-',mfc='none',linewidth=1)
    if j ==0:
        plt.ylabel('Elevation (m)')
    plt.legend(['1C','2C','CPL','L','VAR'])
    #plt.xlim(np.min(fp2['thvm'][clb_t,0:166,0,0])-2.5,fp1['thvm'][clb_t,64,0,0]+.5)
    plt.ylim(0,6000)
    plt.xlim(-1,15)
    plt.xlabel('$Delta\ \theta_v$ (K)')

    '''
    plt.plot(fp1['T_in_K'][clb_t,0:166,0,0],altc)
    plt.plot(fp2['T_in_K'][clb_t,0:166,0,0],altc)
    #plt.plot(fplm['th'][les_t,0:166],altl)
    #plt.plot(fplt['th'][les_t,0:166],altl)
    plt.plot(fpv['T'][var_t,:],altv,'o-',mfc='none',linewidth=1)
    plt.legend(['C1','C2','L1','L2','VAR'])
    plt.ylim(0,6000)
    plt.xlim(250,305)
    '''
    
    
    i=i+1
    j=j+1

# %%
np.max(alt2)

# %%
####
#### DIRECT COMPARISON ####
####

i=0
j=0
plt.figure(figsize=(18,6))
altl=np.linspace(0,5000,166)
fp1=nc.Dataset(clubb_1c+'sgp_'+str(days[0])+'/k_1/agg_outzt.nc','r')
altc=fp1['altitude'][0:166]
for day in days:
    if j==6:
        plt.figure(figsize=(18,6))
        j=0
    dys=str(day)
    fp1=nc.Dataset(clubb_1c+'sgp_'+str(day)+'/k_1/c_1/output/arm_zt.nc','r')
    fp2=nc.Dataset(clubb_2c+'sgp_'+str(day)+'/k_1/c_1/output/arm_zt.nc','r')
    fp3=nc.Dataset(clubb_hf+'sgp_'+str(day)+'/k_1/c_1/output/arm_zt.nc','r')
    fplm=nc.Dataset(les_1c+'trimfr2_'+str(day)+'_01.nc','r')
    fplt=nc.Dataset(les_1c+'trimfr2_'+str(day)+'_00.nc','r')
    alt2=fp1['altitude'][0:166]
    plt.subplot(1,6,j+1)
    hr=15
    st_time= datetime.datetime(int(dys[0:4]),int(dys[4:6]),int(dys[6:8]),hr)
    clb_t,les_t,var_t=index_time(st_time)
    
    altv=np.flip(np.interp(np.flip(var_p*100),np.flip(fp2['p_in_Pa'][clb_t,0:166,0,0]),np.flip(altc)))
    vpt_v=vpt(fpv['T'][var_t,:],fpv['q'][var_t,:]/1000,var_p*100)
    
    vpt_vc=np.interp(altc,altv,vpt_v)
    vpt_vl=np.interp(altl,altv,vpt_v)
    
    vpt_cnf=fp2['thvm'][clb_t,0:166,0,0]
    vpt_cf =fp1['thvm'][clb_t,0:166,0,0]
    vpt_chf=fp3['thvm'][clb_t,0:166,0,0]
    vpt_l  =fplm['thv'][les_t,0:166]
    
    
    plt.plot(vpt_vc-vpt_cf,altc,':')
    plt.plot(vpt_vc-vpt_cnf,altc,':')
    plt.plot(vpt_vc-vpt_chf,altc,'--')
    plt.plot(vpt_vl-vpt_l,altl)
    plt.plot([0,0],[0,6000],'k--',alpha=.5)
    
    plt.legend(['CF','CNF','CHF','L','0'])
    plt.title(day)
    plt.xlim(-5,5)
    plt.ylim(0,6000)
    
    i=i+1
    j=j+1
    

# %%
clb_t

# %%
j=0
plt.figure(figsize=(18,6))
for day in days:
    
    if j==6:
        plt.figure(figsize=(18,6))
        j=0
    plt.subplot(1,6,j+1)
    
    
    fpf1=nc.Dataset(clubb_f2+'sgp_'+str(day)+'/k_2/c_1/output/arm_zt.nc','r')
    fpnf1=nc.Dataset(clubb_nf2+'sgp_'+str(day)+'/k_2/c_1/output/arm_zt.nc','r')
    fphf1=nc.Dataset(clubb_hf2+'sgp_'+str(day)+'/k_2/c_1/output/arm_zt.nc','r')
    fpf2=nc.Dataset(clubb_f2+'sgp_'+str(day)+'/k_2/c_2/output/arm_zt.nc','r')
    fpnf2=nc.Dataset(clubb_nf2+'sgp_'+str(day)+'/k_2/c_2/output/arm_zt.nc','r')
    fphf2=nc.Dataset(clubb_hf2+'sgp_'+str(day)+'/k_2/c_2/output/arm_zt.nc','r')
    
    var='thvm'
    
    vptf1=fpf1[var][800,0:166,0,0]
    vptn1=fpnf1[var][800,0:166,0,0]
    vpth1=fphf1[var][800,0:166,0,0]
    vptf2=fpf2[var][800,0:166,0,0]
    vptn2=fpnf2[var][800,0:166,0,0]
    vpth2=fphf2[var][800,0:166,0,0]
    
    plt.plot(vptf1,altc,'b:')
    plt.plot(vptf2,altc,'r:')
    plt.plot(vptn1,altc,'b--')
    plt.plot(vptn2,altc,'r--')
    plt.plot(vpth1,altc,'b-')
    plt.plot(vpth2,altc,'r-')
    
    #plt.xlim(np.min(vpth1)-1,np.min(vpth1)+10)
    plt.ylim(0,4000)
    
    plt.legend(['F_C','F_H','NF_C','NF_H','HF_C','HF_H'])
    plt.title(day)
    j=j+1

# %%
j=0
plt.figure(figsize=(18,6))
a_n=[]
a_f=[]
a_h=[]
hr=16
time=(hr-5)*60
for day in days:
    
    if j==6:
        plt.figure(figsize=(18,6))
        j=0
    plt.subplot(1,6,j+1)
    
    
    fpf1=nc.Dataset(clubb_f2+'sgp_'+str(day)+'/k_2/c_1/output/arm_zt.nc','r')
    fpnf1=nc.Dataset(clubb_nf2+'sgp_'+str(day)+'/k_2/c_1/output/arm_zt.nc','r')
    fphf1=nc.Dataset(clubb_hf2+'sgp_'+str(day)+'/k_2/c_1/output/arm_zt.nc','r')
    fpf2=nc.Dataset(clubb_f2+'sgp_'+str(day)+'/k_2/c_2/output/arm_zt.nc','r')
    fpnf2=nc.Dataset(clubb_nf2+'sgp_'+str(day)+'/k_2/c_2/output/arm_zt.nc','r')
    fphf2=nc.Dataset(clubb_hf2+'sgp_'+str(day)+'/k_2/c_2/output/arm_zt.nc','r')
    
    var='p_in_Pa'
    
    vptf1=fpf1[var][time,0:166,0,0]
    vptn1=fpnf1[var][time,0:166,0,0]
    vpth1=fphf1[var][time,0:166,0,0]
    vptf2=fpf2[var][time,0:166,0,0]
    vptn2=fpnf2[var][time,0:166,0,0]
    vpth2=fphf2[var][time,0:166,0,0]
    
    a_f.append(altc[np.argmax((vptf2-vptf1)[1:]>0)])
    a_n.append(altc[np.argmax((vptn2-vptn1)[1:]>0)])
    a_h.append(altc[np.argmax((vpth2-vpth1)[1:]>0)])
    
    plt.plot(vptf2-vptf1,altc,'k:')
    plt.plot(vptn2-vptn1,altc,'k--')
    plt.plot(vpth2-vpth1,altc,'k-')
    
    #plt.xlim(np.min(vpth1)-1,np.min(vpth1)+10)
    plt.ylim(0,4000)
    
    plt.legend(['F_C','NF_C','HF_C'])
    plt.title(day)
    j=j+1

# %%

# %%
plt.scatter(np.linspace(1,len(a_f),len(a_f)),a_f)
plt.scatter(np.linspace(1,len(a_f),len(a_f)),a_n)
plt.scatter(np.linspace(1,len(a_f),len(a_f)),a_h)

# %%
print(np.mean(a_f))
print(np.mean(a_n))
print(np.mean(a_h))

# %%
plt.plot(vpt_v,altv)
plt.xlim(300,330)

# %%
plt.plot(fp1['p_in_Pa'][clb_t,0:166,0,0],altc)

# %%
print(float(fpv['time'][var_t]))

# %%
plt.plot(vpt_v,altv)
plt.xlim(300,400)

# %%
plt.plot(vpt_v,altv,'o-')
plt.xlim(350,400)

# %%
(fpv['time'][1])

# %%
plt.plot(fp1['T_in_K'][0,0:166,0,0],altc)
plt.plot(fp1['T_in_K'][clb_t,0:166,0,0],altc)

# %%

plt.imshow(np.transpose(fpv['T_adv_v'][var_t-10:var_t+24*12-10,1:20]+fpv['T_adv_h'][var_t-10:var_t+24*12-10,1:20]),origin='lower',cmap='terrain',extent=(0,12,0,5.5))
plt.colorbar()

# %%
plt.imshow(np.transpose(fpv['dTdt'][var_t-10:var_t+24*12-10,1:20]),origin='lower',cmap='terrain',extent=(0,12,0,5.5))
plt.colorbar()

# %%
p2=fp2['p_in_Pa'][clb_t,0:166,0,0]

# %%
altv=np.interp(np.flip(var_p)*100,np.flip(p2),np.flip(altc))

# %%
var_p[-1]

# %%
dts=datetime.datetime(2016,1,1,0,0)
a=[]
for t in fplm['time'][:]:
    a.append(dts+datetime.timedelta(seconds=t*60))

# %%
plt.plot(a)

# %%
plt.plot(fpv['sw_dn_srf'][0+5:0+20])

# %%
plt.plot(altv)

# %%
plt.scatter(fpv['dTdt'][48+20:48*100:24,5]+,fpv['T_adv_h'][48+20:48*100:24,5],alpha=.25)

# %%
var_t=var_t-4
print(fpv['T'][var_t+1,10]-fpv['T'][var_t,10])
print(fpv['dTdt'][var_t,10])

# %%
print(fp2['time'][:].shape)
print(fplm['time'][:].shape)

# %%
plt.hist(fpv['dqdt'][48+20:48*100,5],bins=100,density=True)
plt.xlim(-1,1)
print(np.std(fpv['dqdt'][48+20:48*100,5]))

# %%
plt.hist(fpv['q_adv_h'][48+20:48*100,5],bins=100,density=True)
plt.xlim(-1,1)
print(np.std(fpv['q_adv_h'][48+20:48*100,5]))
print(np.std(fpv['dqdt'][48+20:48*100,5]))
print(np.mean(fpv['q_adv_h'][48+20:48*100,5]))
print(np.mean(fpv['dqdt'][48+20:48*100,5]))

# %%
a=np.array([0,1,2,3,4,5,6,7,8,9,10])
b=np.interp(np.linspace(0,10,5),np.linspace(0,10,11),a)
print(b)

# %%
np.linspace(0,10,5)

# %%
