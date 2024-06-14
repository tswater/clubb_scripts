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
#     display_name: Python 3
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
plt.rcParams.update({'figure.max_open_warning': 0})

# %%
clubb_dir= '/home/tsw35/tyche/clubb/'

# %%

# %%
plt.figure()
fscp_d2 = nc.Dataset(clubb_dir+'test_cpl/k_2/agg_outsfc.nc','r')
fscp_d = nc.Dataset(clubb_dir+'test_cpl2/k_2/agg_outsfc.nc','r')
fsn = nc.Dataset(clubb_dir+'test_ncpl/k_2/agg_outsfc.nc','r')
plt.plot(fscp_d2['lwp'][:,0,0,0])
plt.plot(fscp_d['lwp'][:,0,0,0])
plt.plot(fsn['lwp'][:,0,0,0])
plt.legend(['D2','D','2C'])

# %%
fscp_d2h = nc.Dataset(clubb_dir+'test_cpl/k_2/c_2/output/arm_sfc.nc','r')
fscp_d2c = nc.Dataset(clubb_dir+'test_cpl/k_2/c_1/output/arm_sfc.nc','r')
fscp_dh = nc.Dataset(clubb_dir+'test_cpl2/k_2/c_2/output/arm_sfc.nc','r')
fscp_dc = nc.Dataset(clubb_dir+'test_cpl2/k_2/c_1/output/arm_sfc.nc','r')
fsnh = nc.Dataset(clubb_dir+'test_ncpl/k_2/c_2/output/arm_sfc.nc','r')
fsnc = nc.Dataset(clubb_dir+'test_ncpl/k_2/c_1/output/arm_sfc.nc','r')

# %%
plt.plot(fscp_d2h['lwp'][:,0,0,0]-fscp_d2c['lwp'][:,0,0,0])
plt.plot(fscp_dh['lwp'][:,0,0,0]-fscp_dc['lwp'][:,0,0,0])
plt.plot(fsnh['lwp'][:,0,0,0]-fsnc['lwp'][:,0,0,0])
plt.legend(['D2','D','2C'])

# %%
fccp_d2=nc.Dataset(clubb_dir+'test_cpl/k_2/clusters.nc','r')
fccp_d=nc.Dataset(clubb_dir+'test_cpl2/k_2/clusters.nc','r')


# %%
def dz2dz(dz_):
    out1=[]
    out2=[]
    for i in range(len(dz_)):
        if dz_[i]==0:
            out1.append(0)
            out2.append(0)
        sdz=str(dz_[i])
        out1.append(float(sdz[0:2]))
        out2.append(float(sdz[3:5]))
    out1=np.array(out1)
    out2=np.array(out2)
    return out1,out2


# %%
plt.figure(figsize=(15,5))
plt.subplot(1,4,1)
plt.plot(fccp_d2['u_r'][:])
plt.plot(fccp_d['u_r'][:])
plt.subplot(1,4,2)
plt.plot(fccp_d2['z_circh'][:])
plt.plot(fccp_d['z_circh'][:])
plt.subplot(1,4,3)
plt.plot(fccp_d2['z_circh'][:]-fccp_d2['z_circl'][:])
plt.plot(fccp_d['z_circh'][:]-fccp_d['z_circl'][:])
plt.subplot(1,4,4)
dz_c2a,dz_c2b=dz2dz(fccp_d2['dz_circ'][:])
dz_ca,dz_cb=dz2dz(fccp_d['dz_circ'][:])
plt.plot(dz_c2a*40,'*')
plt.plot(dz_ca*40,'*')
plt.plot(dz_cb*40,'*')

# %%
(A*dz*300+1*dz*W*290*300)/(A*dz+W*dz*300)

# %%
L=10000
300-10/L*300

# %%

# %%
fmnh = nc.Dataset(clubb_dir+'test_ncpl/k_2/c_2/output/arm_zm.nc','r')
fmnc = nc.Dataset(clubb_dir+'test_ncpl/k_2/c_1/output/arm_zm.nc','r')
fmn  = nc.Dataset(clubb_dir+'test_ncpl/k_2/agg_outzm.nc','r')

# %%
fmch = nc.Dataset(clubb_dir+'test_cpl2/k_2/c_2/output/arm_zm.nc','r')
fmcc = nc.Dataset(clubb_dir+'test_cpl2/k_2/c_1/output/arm_zm.nc','r')
fmc  = nc.Dataset(clubb_dir+'test_cpl2/k_2/agg_outzm.nc','r')

# %%
data=fmc['rtp2'][:,25:100,0,0]-fmn['rtp2'][:,25:100,0,0]
plt.imshow(data.T,origin='lower',cmap='coolwarm',extent=(7,20,0,6))

# %%
data2 = fmc['rtp2'][:,25:100,0,0]
plt.plot(np.mean(data2,axis=1))
data2 = fmn['rtp2'][:,25:100,0,0]  
plt.plot(np.mean(data2,axis=1))

# %%
dt1=np.mean(fmc['rtp2'][:,25:150,0,0],axis=1)
dt2=np.mean(fmn['rtp2'][:,25:150,0,0],axis=1)
dt2=dt2/np.max(dt1)
dt1=dt1/np.max(dt1)
print(np.mean(dt1-dt2))
plt.plot(dt1)
plt.plot(dt2)

# %%
plt.plot(fmn['rtp2'][600,:,0,0],'g')
plt.plot(fmc['rtp2'][600,:,0,0],'r-')
#plt.plot(fmnc['rtp2'][600,:,0,0],'b-')


# %%

# %%
for k in fccp_d2.variables:
    print(k)

# %%
tst=str(fccp_d['dz_circ'][100])

# %%
print(tst)

# %%
fsn.close()
fsnh.close()
fsnc.close()

# %%
tst[3:4]

# %%
import pickle
fp=open('tunedp.p','rb')
mean_val=pickle.load(fp)
peak_val=pickle.load(fp)
h_c_lwpm=pickle.load(fp)
mean_rtm=pickle.load(fp)
h_c_rtm=pickle.load(fp)
crs=pickle.load(fp)
incs=pickle.load(fp)
cc1s=pickle.load(fp)
lwp =pickle.load(fp)
rtp2 = pickle.load(fp)
maxzl = pickle.load(fp)
maxzh = pickle.load(fp)
urmu = pickle.load(fp)
urmx = pickle.load(fp)
urfc = pickle.load(fp)
ccc0 = pickle.load(fp)
thlp2 = pickle.load(fp)
wp2 = pickle.load(fp)

fp.close()

# %%
days=list(mean_val.keys())
days.sort()
ccc={}
dn=0
for day in days:
    ccc[day]=[]
    for i in range(len(crs[day])):
        if incs[day][i]==1:
            ccc[day].append(0)
            dn=1
        else:
            ccc[day].append(ccc0[day][i-dn])


# %%

# %%
def plot_test(param,test,name,c2=0,day=0):
    param=np.copy(param)
    #param=param[cc1s[day]==3]
    #test=test[cc1s[day]==3]
    gu=np.unique(param)
    gu.sort()
    test_line=[]
    for g in gu:
        test_line.append(np.nanmean(test[param==g]))
    try:
        c2 = c2+1
        plt.scatter(param,test,alpha=.25)
    except:
        c2=np.array(c2)
        c2=c2[cc1s[day]==3]
        color=plt.cm.coolwarm((np.array(c2)-.5)/(1.5-.5))
        plt.scatter(param,test,c=color,alpha=.25)
    plt.plot(gu,test_line,'b-')
    xaxh=np.nanmax(gu[gu>0])-np.nanmin(gu[gu>0])
    plt.xlim(np.nanmin(gu[gu>0])-xaxh*.1,np.nanmax(gu[gu>0])+xaxh*.3)
    yaxh=np.nanmax(test[test>0])-np.nanmin(test[test>0])
    ymin=np.nanmin(test[test>0])-.1*yaxh
    ymax=np.nanmax(test[test>0])+.3*yaxh
    plt.ylim(ymin,ymax)
    plt.title(name)


# %%
days=list(mean_val.keys())
days.sort()
import warnings
warnings.filterwarnings("ignore")
for day in days:
    if len(crs[day])<110:
        continue
    plt.figure(figsize=(24,12))
    
    # CRS
    plt.subplot(3,8,1)
    plot_test(crs[day],wp2[day],'WP2')
    plt.ylabel('CR')
    plt.subplot(3,8,2)
    plot_test(crs[day],mean_val[day],'Mean LWP')
    plt.subplot(3,8,3)
    plot_test(crs[day],h_c_lwpm[day],'H/C Ratio LWP')
    plt.subplot(3,8,4)
    plot_test(crs[day],rtp2[day],'RTP2')
    plt.subplot(3,8,5)
    plot_test(crs[day],thlp2[day],'THLP2')
    plt.subplot(3,8,6)
    plot_test(crs[day],peak_val[day],'PEAK LWP')
    plt.subplot(3,8,7)
    plot_test(crs[day],urmx[day],'MAX UR')
    
    plt.subplot(3,8,8)
    colors=['red','white','gainsboro','lightgrey','darkgrey','dimgrey','black']
    gs=np.unique(crs[day])
    gs.sort() 
    leglist=[]
    for idx in range(len(gs)):
        g=gs[idx]
        c=colors[idx]
        m_LWP=np.mean(lwp[day][crs[day]==g,:],axis=0)
        plt.plot(m_LWP,color=c)
        leglist.append(str(g))
    plt.xlim(120,1000)
    plt.legend(leglist,loc='upper left',labelspacing=.25,borderaxespad=0.1)
    plt.title('Mean LWP in Time')
    
    # CC1S
    print(day)
    plt.subplot(3,8,9)
    plot_test(cc1s[day],wp2[day],'WP2')
    plt.ylabel('CC1')
    plt.subplot(3,8,10)
    plot_test(cc1s[day],mean_val[day],'Mean LWP')
    plt.subplot(3,8,11)
    plot_test(cc1s[day],h_c_lwpm[day],'H/C Ratio LWP')
    plt.subplot(3,8,12)
    plot_test(cc1s[day],rtp2[day],'RTP2')
    plt.subplot(3,8,13)
    plot_test(cc1s[day],thlp2[day],'THLP2')
    plt.subplot(3,8,14)
    plot_test(cc1s[day],peak_val[day],'PEAK LWP')
    plt.subplot(3,8,15)
    plot_test(cc1s[day],urmx[day],'MAX UR')
    
    plt.subplot(3,8,16)
    colors=['red','lightgrey','darkgrey','black']
    gs=np.unique(cc1s[day])
    gs.sort() 
    leglist=[]
    for idx in range(len(gs)):
        g=gs[idx]
        c=colors[idx]
        m_LWP=np.mean(lwp[day][cc1s[day]==g,:],axis=0)
        plt.plot(m_LWP,color=c)
        leglist.append(str(g))
    plt.xlim(120,1000)
    plt.legend(leglist,loc='upper left',labelspacing=.25,borderaxespad=0.1)
    plt.title('Mean LWP in Time')
    
    # CC1S
    plt.subplot(3,8,17)
    plot_test(ccc[day],wp2[day],'WP2')
    plt.ylabel('CCC')
    plt.subplot(3,8,18)
    plot_test(ccc[day],mean_val[day],'Mean LWP')
    plt.xlabel('Parameter Value')
    plt.subplot(3,8,19)
    plot_test(ccc[day],h_c_lwpm[day],'H/C Ratio LWP')
    plt.subplot(3,8,20)
    plot_test(ccc[day],rtp2[day],'RTP2')
    plt.xlabel('Parameter Value')
    plt.subplot(3,8,21)
    plot_test(ccc[day],thlp2[day],'THLP2')
    plt.subplot(3,8,22)
    plot_test(ccc[day],peak_val[day],'PEAK LWP')
    plt.xlabel('Parameter Value')
    plt.subplot(3,8,23)
    plot_test(ccc[day],urmx[day],'MAX UR')
    
    plt.subplot(3,8,24)
    colors=['red','white','gainsboro','lightgrey','darkgrey','dimgrey','black','black']
    gs=np.unique(ccc[day])
    gs.sort() 
    leglist=[]
    for idx in range(len(gs)):
        g=gs[idx]
        c=colors[idx]
        m_LWP=np.mean(lwp[day][ccc[day]==g,:],axis=0)
        plt.plot(m_LWP,color=c)
        leglist.append(str(g))
    plt.xlim(120,1000)
    plt.legend(leglist,loc='upper left',labelspacing=.25,borderaxespad=0.1)
    plt.title('Mean LWP in Time')
    
    plt.suptitle(day)
    plt.subplots_adjust(wspace=.2,hspace=.35)

# %%
cc1s[day][cc1s[day]==3]

# %%
for day in days:
    try:
        test1='/home/tsw35/tyche/clubb/circ_tune2d/'+str(day)+'/dp_cr3k2cc3ccc1'
        test2='/home/tsw35/tyche/clubb/circ_tune2d/'+str(day)+'/dp_cr3k2cc3ccc1.5'
        #fth=nc.Dataset(test+'/k_2/c_2/output/arm_zt.nc','r')
        #ftc=nc.Dataset(test+'/k_2/c_1/output/arm_zt.nc','r')
        fc1 =nc.Dataset(test1+'/k_2/clusters.nc','r')
        fc2 = nc.Dataset(test2+'/k_2/clusters.nc','r')
    except:
        continue
    plt.figure()
    dt1=fc1['dvpt'][:]
    dt2=fc2['dvpt'][:]
    ur=fc1['u_r'][:]
    ur2=fc2['u_r'][:]
    '''
    ur0=fc['u_r0'][:]
    dtt=dtt[ur0>0]
    cnt=len(dtt)/len(ur0)
    ur0=ur0[ur0>0]
    a=np.mean(ur0/dtt)
    plt.scatter(dtt,[0]*len(dtt))
    dtt=np.linspace(0,2)
    ur1=a*(dtt)**(.5)
    ur2=a*(dtt)**(.75)
    ur3=a*(dtt)**(1)
    ur4=a*(dtt)**(1.25)
    ur5=a*(dtt)**(1.5)
    ur6=a*(dtt)**(2)
    #dtt[dtt<=0]=float('nan')
    
    plt.plot(dtt,ur1)
    plt.plot(dtt,ur2)
    plt.plot(dtt,ur3)
    plt.plot(dtt,ur4)
    plt.plot(dtt,ur5)
    plt.plot(dtt,ur6)
    plt.title(str(day)+': '+str(cnt))
    plt.legend(['r','.5','.75','1','1.25','1.5','2'])
    #plt.plot(np.linspace(0,200,len(fth['thlm'][:,0,0,0])),dtt*10)
    '''
    plt.subplot(1,2,1)
    plt.plot(dt1)
    plt.plot(dt2)
    plt.subplot(1,2,2)
    plt.plot(ur)
    plt.plot(ur2)

# %%
fc1

# %%
a=np.array([1,2,3,4,5])
print(a[3])
print(a[4-3:4])

# %%
fc

# %%
plt.imshow((fth['thlm'][:,0:75,0,0]-ftc['thlm'][:,0:75,0,0]).T,cmap='coolwarm',origin='lower',extent=(0,22,0,25))
plt.colorbar()

# %%
test={}
test['hello']=[1,2,3]
def changeit(dick_):
    dick_['hello'].append(4)
changeit(test)
print(test)

# %%
a='um'
b=a[:-1]+'p2'

# %%
b

# %%
days=    [20160625,20160716,20160719,20160720,20170609,
          20170626,20170627,20170629,20170705,20170709,
          20170712,20170716,20170717,20170719,20170720,
          20170728,20170826,20170922,20170923,20170924,
          20180522,20180530,20180618,20180619,20180704,
          20180705,20180523,20180707,20180709,20180710,
          20180711,20180712,20180809,20180811,20180916,
          20180917,20190707,20190709,20190714]

# %%
fplm.variables

# %%
dlwpc=[]
dlwpc2=[]
dlwpl=[]
drtpc=[]
drtpl=[]
dtkec=[]
dtkel=[]

# %%

altl=np.linspace(0,5000,166)
les_1c ='/home/tsw35/tyche/data/LES_1C/'

for day in days:
    plt.figure(figsize=(18,5))
    dir2c = clubb_dir+'sgp_2c_d/sgp_'+str(day)+'/'
    dir1c = clubb_dir+'sgp_1c_d/sgp_'+str(day)+'/'
    dircp = clubb_dir+'sgp_cpl_d/sgp_'+str(day)+'/'
    
    fmcp = nc.Dataset(dircp+'k_2/agg_outzm.nc','r')
    fmcpc = nc.Dataset(dircp+'k_2/c_1/output/arm_zm.nc','r')
    fmcph = nc.Dataset(dircp+'k_2/c_2/output/arm_zm.nc','r')
    fm2c = nc.Dataset(dir2c+'k_2/agg_outzm.nc','r')
    fm1c = nc.Dataset(dir1c+'k_1/c_1/output/arm_zm.nc','r')
    fplm=nc.Dataset(les_1c+'trimfr2_'+str(day)+'_01.nc','r')
    fplt=nc.Dataset(les_1c+'trimfr2_'+str(day)+'_00.nc','r')
    
    alt=fmcp['altitude'][0:150]
    
    var1='rtp2'#'thlp2'
    var2='qv2'#'thl2'
        
    for i in range(1,7):
        lest=6*6+(i-1)*6
        clbt=8*60+(i-1)*60
        
        plt.subplot(1,6,i)
        '''
        plt.plot(fmcp[var1][clbt,0:150,0,0],alt,'b-')
        plt.plot(fm1c[var1][clbt,0:150,0,0],alt,'b--')
        plt.plot(fplt[var2][lest,0:166],altl,'g-')
        plt.plot(fplm[var2][lest,0:166],altl,'g--')
        '''
        
        plt.plot(np.mean(fmcp[var1][clbt-5:clbt+5,0:150,0,0],axis=0),alt,'b-')
        plt.plot(np.mean(fm1c[var1][clbt-5:clbt+5,0:150,0,0],axis=0),alt,'b--')
        #plt.plot(np.mean(fm2c[var1][clbt-5:clbt+5,0:150,0,0],axis=0),alt,'grey')
        plt.plot(fplt[var2][lest,0:166],altl,'g-')
        plt.plot(fplm[var2][lest,0:166],altl,'g--')
        
        '''
        plt.plot(fmcp[var1][clbt,0:150,0,0]-fm1c[var1][clbt,0:150,0,0],alt,'b-')
        plt.plot(fplt[var2][lest,0:166]-fplm[var2][lest,0:166],altl,'g-')
        '''
        
        '''
        plt.plot(np.mean(fmcp[var1][clbt-5:clbt+5,0:150,0,0]-fm1c[var1][clbt-5:clbt+5,0:150,0,0],axis=0),alt,'b-')
        plt.plot(fplt[var2][lest,0:166]-fplm[var2][lest,0:166],altl,'g-')
        '''
        
        plt.title(str(day)+': '+str(clbt/60+5))
        
        i=i+1
    break
        

# %%
for day in days:
    plt.figure()
    dircp = clubb_dir+'sgp_cpl_d/sgp_'+str(day)+'/'
    fcc=nc.Dataset(dircp+'k_2/clusters.nc','r')
    plt.plot()


# %%
def moving_average(a, n=30) :
    ret = np.cumsum(a, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    return ret[n - 1:] / n


# %%
days2=[20170627, 20190709, 20180704, 20170728, 20170922, 20170717]

# %%
days_=os.listdir(clubb_dir+'sgp_2c_d/')
days_.sort()
daysa=[]
for day in days_:
    if '20190804' in day:
        break
    daysa.append(day[4:])

# %%
j=0
for day in [20150801,20160625,20180707]:
    i = 5
    lest=6*6+(i-1)*6
    clbt=8*60+(i-1)*60
    dir2c = clubb_dir+'mar15fit_2c/sgp_'+str(day)+'/'
    dir1c = clubb_dir+'mar15fit_1c/sgp_'+str(day)+'/'
    dircp = clubb_dir+'mar17fit_cpl/sgp_'+str(day)+'/'
    #dircp2 = clubb_dir+'sgp_cpl_dhi/sgp_'+str(day)+'/'
    try:
        fmcp = nc.Dataset(dircp+'k_2/agg_outzm.nc','r')
        #fmcpc = nc.Dataset(dircp+'k_2/c_1/output/arm_zm.nc','r')
        #fmcph = nc.Dataset(dircp+'k_2/c_2/output/arm_zm.nc','r')
        #fm2c = nc.Dataset(dir2c+'k_2/agg_outzm.nc','r')
        fm1c = nc.Dataset(dir1c+'k_1/c_1/output/arm_zm.nc','r')
        fplm=nc.Dataset(les_1c+'trimfr2_'+str(day)+'_01.nc','r')
        fplt=nc.Dataset(les_1c+'trimfr2_'+str(day)+'_00.nc','r')
        fscp=nc.Dataset(dircp+'k_2/agg_outsfc.nc','r')
        #fscp2=nc.Dataset(dircp+'k_2/agg_outsfc.nc','r')
        fs1c=nc.Dataset(dir1c+'k_1/c_1/output/arm_sfc.nc','r')
    except:
        continue
    plt.figure(figsize=(6,4),dpi=300)
    
    lwp2=moving_average(fscp['lwp'][:,0,0,0])
    lwp1=moving_average(fs1c['lwp'][:,0,0,0])
    lwpht=fplt['lwp'][:]
    lwphm=fplm['lwp'][:]
    plt.plot(np.linspace(5,22,lwp1.shape[0]),lwp2,'b-')
    plt.plot(np.linspace(5,22,lwp1.shape[0]),lwp1,'b--')
    plt.plot(np.linspace(7,22,lwpht.shape[0]),lwpht,'g-')
    plt.plot(np.linspace(7,22,lwphm.shape[0]),lwphm,'g--')
    plt.ylabel(r"$LWP (kg/m^{-2})$",fontsize=20,weight='bold')
    plt.yticks(fontsize=18)
    #plt.title(day)
    plt.xticks(np.arange(5,25,5),fontsize=18)
    if j == 0:
        plt.title(r"LWP",fontsize=30,weight='bold')
        #plt.yticks(np.arange(0,.0011,.0003),fontsize=18,rotation=45)
    #else:
    #    plt.yticks(np.arange(0,.014,.004),fontsize=18,rotation=45)
    if j == 2:
        plt.xlabel(r"$Time\ (hr)$",fontsize=24,weight='bold')
    
    j=j+1

# %%
j=0
for day in [20150801,20160625,20180707]:
    i = 5
    lest=6*6+(i-1)*6
    clbt=8*60+(i-1)*60
    dir2c = clubb_dir+'tall_2c/sgp_'+str(day)+'/'
    dir1c = clubb_dir+'tall_1c/sgp_'+str(day)+'/'
    dircp = clubb_dir+'tall2_cpl/sgp_'+str(day)+'/'
    #dircp2 = clubb_dir+'sgp_cpl_dhi/sgp_'+str(day)+'/'
    try:
        fmcp = nc.Dataset(dircp+'k_2/agg_outzm.nc','r')
        #fmcpc = nc.Dataset(dircp+'k_2/c_1/output/arm_zm.nc','r')
        #fmcph = nc.Dataset(dircp+'k_2/c_2/output/arm_zm.nc','r')
        #fm2c = nc.Dataset(dir2c+'k_2/agg_outzm.nc','r')
        fm1c = nc.Dataset(dir1c+'k_1/c_1/output/arm_zm.nc','r')
        fplm=nc.Dataset(les_1c+'trimfr2_'+str(day)+'_01.nc','r')
        fplt=nc.Dataset(les_1c+'trimfr2_'+str(day)+'_00.nc','r')
        fscp=nc.Dataset(dircp+'k_2/agg_outsfc.nc','r')
        #fscp2=nc.Dataset(dircp+'k_2/agg_outsfc.nc','r')
        fs1c=nc.Dataset(dir1c+'k_1/c_1/output/arm_sfc.nc','r')
    except:
        continue
    plt.figure(figsize=(4,4),dpi=300)
    
    var1='thlp2'#'rtp2'#'thlp2'
    var2='thl2'#'qv2'#'thl2'
    plt.plot(np.mean(fmcp[var1][clbt-5:clbt+5,0:125,0,0],axis=0),alt[0:125],'b-')
    plt.plot(np.mean(fm1c[var1][clbt-5:clbt+5,0:125,0,0],axis=0),alt[0:125],'b--')
    #plt.plot(np.mean(fm2c[var1][clbt-5:clbt+5,0:150,0,0],axis=0),alt,'grey')
    plt.plot(fplt[var2][lest,0:166],altl,'g-')
    plt.plot(fplm[var2][lest,0:166],altl,'g--')
    #ax.ticklabel_format(scilimits=(0,0))
    #formatter.set_powerlimits((0, 0))
    #plt.yticks(np.arange(0,5,1.5),fontsize=18)
    #plt.xticks(np.arange(0,1.5,.3),fontsize=18)
    if j == 0:
        plt.title(r"$\theta'^2$",fontsize=30,weight='bold')
        plt.xlim(-.08,1.5)
    #if j == 1:
        #plt.xlim(-.045,.75)
    if j == 2:
        plt.xlabel(r"$\theta'^2\ (K^2)$",fontsize=24,weight='bold')
        plt.xlim(-.06,1)
    
    j=j+1

# %%
j=0
for day in [20150801,20160625,20180707]:
    i = 5
    lest=6*6+(i-1)*6
    clbt=8*60+(i-1)*60
    dir2c = clubb_dir+'tall_2c/sgp_'+str(day)+'/'
    dir1c = clubb_dir+'tall_1c/sgp_'+str(day)+'/'
    dircp = clubb_dir+'tall2_cpl/sgp_'+str(day)+'/'
    #dircp2 = clubb_dir+'sgp_cpl_dhi/sgp_'+str(day)+'/'
    try:
        fmcp = nc.Dataset(dircp+'k_2/agg_outzm.nc','r')
        #fmcpc = nc.Dataset(dircp+'k_2/c_1/output/arm_zm.nc','r')
        #fmcph = nc.Dataset(dircp+'k_2/c_2/output/arm_zm.nc','r')
        #fm2c = nc.Dataset(dir2c+'k_2/agg_outzm.nc','r')
        fm1c = nc.Dataset(dir1c+'k_1/c_1/output/arm_zm.nc','r')
        fplm=nc.Dataset(les_1c+'trimfr2_'+str(day)+'_01.nc','r')
        fplt=nc.Dataset(les_1c+'trimfr2_'+str(day)+'_00.nc','r')
        fscp=nc.Dataset(dircp+'k_2/agg_outsfc.nc','r')
        #fscp2=nc.Dataset(dircp+'k_2/agg_outsfc.nc','r')
        fs1c=nc.Dataset(dir1c+'k_1/c_1/output/arm_sfc.nc','r')
    except:
        continue
    fig,ax = plt.subplots(figsize=(4.5,4),dpi=300)
    var1='rtp2'#'thlp2'
    var2='qv2'#'thl2'
    plt.plot(np.mean(fmcp[var1][clbt-5:clbt+5,0:125,0,0],axis=0),alt[0:125],'b-')
    plt.plot(np.mean(fm1c[var1][clbt-5:clbt+5,0:125,0,0],axis=0),alt[0:125],'b--')
    #plt.plot(np.mean(fm2c[var1][clbt-5:clbt+5,0:150,0,0],axis=0),alt,'grey')
    plt.plot(fplt[var2][lest,0:166],altl,'g-')
    plt.plot(fplm[var2][lest,0:166],altl,'g--')
    ax.ticklabel_format(scilimits=(0,0))
    #formatter.set_powerlimits((0, 0))
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    if j == 0:
        plt.title(r"$q'^2$",fontsize=30,weight='bold')
        plt.legend(['2 COLUMN','1 COLUMN','LES HET','LES HMG'],fontsize=18)
        plt.ylabel('2017-06-09\n Elevation (km)',fontsize=20,weight='bold')
    if j == 1:
        plt.ylabel('2017-07-09\n Elevation (km)',fontsize=20,weight='bold')
    if j == 2:
        plt.xlabel(r"$q'^2\ (g/g)^2$",fontsize=24,weight='bold')
        plt.ylabel('2018-07-07\n Elevation (km)',fontsize=20,weight='bold')
        plt.xticks(np.arange(0,.000015,.000003),fontsize=18)
    j=j+1

# %%

altl=np.linspace(0,5000,166)
les_1c ='/home/tsw35/tyche/data/LES_1C/'
dlwpc=[]
dlwpc2=[]
dlwpl=[]
drtpc=[]
drtpl=[]
dtkec=[]
dtkel=[]
j= 6
for day in daysa:
    print('.',end='')
    if j == 6:
        #plt.figure(figsize=(15,7))
        plt.figure(figsize=(18,5))
        j=0
    i = 5
    lest=6*6+(i-1)*6
    clbt=8*60+(i-1)*60
    dir2c = clubb_dir+'tall_2c/sgp_'+str(day)+'/'
    dir1c = clubb_dir+'tall_1c/sgp_'+str(day)+'/'
    dircp = clubb_dir+'tall_cpl/sgp_'+str(day)+'/'
    #dircp2 = clubb_dir+'sgp_cpl_dhi/sgp_'+str(day)+'/'
    try:
        fmcp = nc.Dataset(dircp+'k_2/agg_outzm.nc','r')
        fmcpc = nc.Dataset(dircp+'k_2/c_1/output/arm_zm.nc','r')
        fmcph = nc.Dataset(dircp+'k_2/c_2/output/arm_zm.nc','r')
        fm2c = nc.Dataset(dir2c+'k_2/agg_outzm.nc','r')
        fm1c = nc.Dataset(dir1c+'k_1/c_1/output/arm_zm.nc','r')
        fplm=nc.Dataset(les_1c+'trimfr2_'+str(day)+'_01.nc','r')
        fplt=nc.Dataset(les_1c+'trimfr2_'+str(day)+'_00.nc','r')
        fscp=nc.Dataset(dircp+'k_2/agg_outsfc.nc','r')
        #fscp2=nc.Dataset(dircp+'k_2/agg_outsfc.nc','r')
        fs1c=nc.Dataset(dir1c+'k_1/c_1/output/arm_sfc.nc','r')
    except:
        continue
    
    alt=fmcp['altitude'][0:150]
    
    var1='thlp2'#'rtp2'#'thlp2'
    var2='thl2'#'qv2'#'thl2'
    
    tkecp=np.sum(fmcp['wp2'][:,:,0,0],axis=1)*.5
    tke1c=np.sum(fm1c['wp2'][:,:,0,0],axis=1)*.5
    
    dtkec.append(np.mean(tkecp-tke1c))
    dtkel.append(np.mean(fplt['w2'][:]-fplm['w2'][:]))
    
    dlwpc.append(np.mean(fscp['lwp'][:,0,0,0]-fs1c['lwp'][:,0,0,0]))
    #dlwpc2.append(np.mean(fscp2['lwp'][:,0,0,0]-fs1c['lwp'][:,0,0,0]))
    dlwpl.append(np.mean(fplt['lwp'][:]-fplm['lwp'][:]))
    
    drtpc.append(np.mean(fmcp[var1][clbt-5:clbt+5,25:125,0,0]-fm1c[var1][clbt-5:clbt+5,25:125,0,0]))
    drtpl.append(np.mean(fplt[var2][lest,33:166]-fplm[var2][lest,33:166]))
    
    '''
    plt.subplot(2,3,j+1)
    plt.plot(np.linspace(0,1,tkecp.shape[0]),tkecp,'b-')
    plt.plot(np.linspace(0,1,tkecp.shape[0]),tke1c,'b--')
    plt.plot(np.linspace(0,1,fplt['tke'][:].shape[0]),np.sum(fplt['w2'][:],axis=1),'g-')
    plt.plot(np.linspace(0,1,fplt['tke'][:].shape[0]),np.sum(fplm['w2'][:],axis=1),'g--')
    plt.legend(['CLUBB CIRC','CLUBB 1C','LES HET','LES HMG'])
    
    '''
    '''
    ax=plt.subplot(1,6,j+1)
    plt.plot(np.mean(fmcp[var1][clbt-5:clbt+5,0:150,0,0],axis=0),alt,'b-')
    plt.plot(np.mean(fm1c[var1][clbt-5:clbt+5,0:150,0,0],axis=0),alt,'b--')
    #plt.plot(np.mean(fm2c[var1][clbt-5:clbt+5,0:150,0,0],axis=0),alt,'grey')
    plt.plot(fplt[var2][lest,0:166],altl,'g-')
    plt.plot(fplm[var2][lest,0:166],altl,'g--')
    ax.ticklabel_format(scilimits=(0,0))
    #formatter.set_powerlimits((0, 0))
    plt.legend(['CLUBB CIRC','CLUBB 1C','LES HET','LES HMG'])
    '''
    '''
    plt.plot(np.mean(fmcp[var1][clbt-5:clbt+5,0:150,0,0]-fm1c[var1][clbt-5:clbt+5,0:150,0,0],axis=0),alt,'b-')
    plt.plot(fplt[var2][lest,0:166]-fplm[var2][lest,0:166],altl,'g-')
    plt.legend(['CLUBB','LES'])
    '''
    if j == 0:
        plt.ylabel('Elevation (m)')
    plt.title(str(day))
    plt.xlabel('$\sigma^2_{\theta}$ (K)')
    j=j+1

# %%
fplt

# %%
plt.plot([-1,1],[-1,1],'grey',alpha=.5)
plt.scatter(drtpc,drtpl)
plt.xlim(-.0000001,.0000015)
plt.ylim(-.0000001,.0000015)
print(metrics.r2_score(drtpc,drtpl))
print(stats.pearsonr(drtpc,drtpl))
plt.figure()
plt.plot([-1,1],[-1,1],'grey',alpha=.5)
plt.scatter(dtkec,dtkel)
plt.xlim(-.05,.6)
plt.ylim(-.05,.6)
print(metrics.r2_score(dtkec,dtkel))
print(stats.pearsonr(dtkec,dtkel))
plt.figure()
plt.plot([-1,1],[-1,1],'grey',alpha=.5)
plt.scatter(dlwpc,dlwpl)
plt.xlim(-.005,.02)
plt.ylim(-.005,.02)
print(metrics.r2_score(dlwpc,dlwpl))
print(stats.pearsonr(dlwpc,dlwpl))

# %%
plt.hist(np.array(dlwpc2)-np.array(dlwpc))

# %%
for day in daysa:
    if j == 6:
        #plt.figure(figsize=(15,7))
        plt.figure(figsize=(18,5))
        j=0
    i = 5
    lest=6*6+(i-1)*6
    clbt=8*60+(i-1)*60
    dir2c = clubb_dir+'sgp_2c_d/sgp_'+str(day)+'/'
    dir1c = clubb_dir+'sgp_1c_d/sgp_'+str(day)+'/'
    dircp = clubb_dir+'sgp_cpl_d4/sgp_'+str(day)+'/'
    try:
        fc = nc.Dataset(dir2c+'k_2/clusters.nc','r')
    except:
        continue
    break

# %%
fc

# %%

altl=np.linspace(0,5000,166)
les_1c ='/home/tsw35/tyche/data/LES_1C/'

j= 6
for day in daysa:
    if j == 6:
        #plt.figure(figsize=(15,7))
        plt.figure(figsize=(18,5))
        j=0
    i = 5
    lest=6*6+(i-1)*6
    clbt=8*60+(i-1)*60
    dir2c = clubb_dir+'sgp_2c_d/sgp_'+str(day)+'/'
    dir1c = clubb_dir+'sgp_1c_d/sgp_'+str(day)+'/'
    dircp = clubb_dir+'sgp_cpl_d4/sgp_'+str(day)+'/'
    try:
        ftcpc = nc.Dataset(dir2c+'k_2/c_1/output/arm_zt.nc','r')
        ftcph = nc.Dataset(dir2c+'k_2/c_2/output/arm_zt.nc','r')
    except:
        continue
    
    alt=ftcpc['altitude'][0:150]

    plt.subplot(1,6,j+1)
    datac=np.mean(ftcpc['thvm'][clbt-5:clbt+5,0:150,0,0],axis=0)
    plt.plot(datac,alt,'b-')
    plt.plot(np.mean(ftcph['thvm'][clbt-5:clbt+5,0:150,0,0],axis=0),alt,'r-')
    
    plt.xlim(np.min(datac)-1,np.min(datac)+8)
    
    plt.title(str(day)+': '+str(clbt/60+5))
    j=j+1
    if j>50:
        break

# %%
'''
def vspeed(thvh,thvc,cur=.1,a=1,l=10000):
    ur=cur*l**(1/2)*9.8**(1/2)*((thvh-thvc)/300)**a
    t0=np.mean(thvh[:,3:5],axis=1)
    cutt=.2
    zur=ur.copy()
    for t in range(ur.shape[0]):
        dE=0
        first=True
        for j in range(1,ur.shape[1]):
            if (thvh[t,j]-t0)<cutt:
                zur[t,j]=zur[t,j-1]+zur[t,j]
            elif thvh[t,j]>thvc[t,j]:
                zur[t,j]=zur[t,j-1]
            elif thvh[t,j]<thvc[t,j]:
                
'''        


# %%
def vspeed(thvh,thvc,W,A,dz,r,P,Tk,cur=2.5,a=1,l=20000):
    ur=cur*l**(1/2)*9.8**(1/2)*((thvh-thvc)/300)**a
    KE=np.zeros(ur.shape)
    KEh=np.zeros(ur.shape)
    BE=np.zeros(ur.shape)
    w=np.zeros(ur.shape)
    
    # first level
    msk=KE[:,0]==0
    KEh[:,0]=.5*W*np.mean(ur[:,0:5],axis=1)**3*dz
    KE[:,0]=KEh[:,0]
    w[:,0]=(2*KEh[:,0]/A)**(1/3)
    dT=[]
    for x in range(1,ur.shape[1]):
        if x<=5:
            KEh[msk,x]=.5*W*np.mean(ur[msk,0:5],axis=1)**3*dz
        else:
            KEh[msk,x]=.5*W*ur[msk,x]**3*dz
        if x>4:
            tenv=(thvh[msk,x-1]+thvh[msk,x])/2
            #rh=RH(r[msk,x],P[msk,x],Tk[msk,x-1])
            #ad=9.808*(1+(2501000*r[msk,x]/rh)/(287*Tk[msk,x-1]))/(1003.5+(2501000**2*r[msk,x]/rh)/(461.5*Tk[msk,x-1]**2))
            #ad[rh<1]=.0098
            #tprc=(Tk[msk,x-1]-ad*dz)*(1+.61*r[msk,x-1])*(100000/P[msk,x])**.286
            tprc=thvh[msk,x-1]
            try:
                dT.append(np.nanmean(tprc-tenv))
            except:
                pass
            BE[msk,x]=A*9.8*(tprc-tenv)/(2*tenv)*w[msk,x-1]*dz
        KE[msk,x]=BE[msk,x]+KEh[msk,x]+KE[msk,x-1]
        msk=msk&(KE[:,x]>0)
        w[msk,x]= (2*KE[msk,x]/A)**(1/3)
        ur[~msk,x]=0
    #print(np.nanmean(np.abs(dT)))
    return w, KE, KEh, BE,ur,dT


# %%
w=2
lpse=7/2000
d=np.sqrt(w**2*300/(9.8*lpse))
print(d*lpse)
print(W*dz/A)

# %%
print(list(range(-1,-1,-1)))

# %%
for x in range(5,-1,-1):
    print(x)

# %%
.5*(12725434)*(.01)**2


# %%
def dTt(thvh,thvc,W,A,dz,r,P,Tk,Tk2,cur=1,a=1,l=20000,decay=2500):
    vpts=np.zeros(thvh.shape)
    ur=cur*l**(1/2)*9.8**(1/2)*((thvh-thvc)/300)**a
    out=np.zeros(thvh.shape)
    tprcp=np.max(Tk[:],axis=1)
    tprc=np.max(Tk[:],axis=1)
    urt=np.cumsum(ur,axis=1)
    rr=r[:,0].copy()
    vpts[:,0]=tprc*(1+.61*rr)*(100000/P[:,0])**.286
    #mc=np.zeros(rr.shape,dtype=bool)
    for x in range(1,thvh.shape[1]):
        tenv=(thvh[:,x]+thvh[:,x-1])/2#(thvh[:,x-1]/(100000/P[:,x-1])**.286+thvh[:,x]/(100000/P[:,x])**.286)
        #rr=r[:,x]
        rh=RH(rr,P[:,x],tprcp)
        #mc=((tprc>thvh[:,x])&(rh>1))|mc
        rr2=rr.copy()
        rr2[rh>1]=rr2[rh>1]/rh[rh>1]
        ad=9.808*(1+(2501000*rr2)/(287*tprcp))/(1003.5+(2501000**2*rr2)/(461.5*tprcp**2))
        ad[rh<1]=.0098
        msk=ur[:,x]>0
        tprc=(tprcp-ad*dz)
        tprc[msk]=((tprcp[msk]-ad[msk]*dz)*urt[msk,x-1]+ur[msk,x]*Tk[msk,x])/urt[msk,x]
        tprc=(tprc+Tk[:,x]*(dz/decay))/(dz/decay+1)
        rr[msk]=(rr[msk]*urt[msk,x-1]+ur[msk,x]*r[msk,x])/urt[msk,x]
        rr=(rr+r[:,x]*(dz/decay))/(dz/decay+1)
        T1=tprcp*(1+.61*rr2)*(100000/P[:,x])**.286
        #*(1+.61*r2[:])*(100000/P[:,x])**.286)
        tprcp=tprc.copy()#/(1+.61*r2[:])/(100000/P[:,x])**.286
        tprc=tprc*(1+.61*rr2)*(100000/P[:,x])**.286
        #b=(1+.61*r[:,x-1])*(100000/P[:,x])**.286
        #out[:,x]=b*Tk[:,x-1]*(np.log(dt+thvh[:,x-1])-np.log(thvh[:,x-1]))*dz/dt-\
        #         ad*b*dz**2/dt**2*(dt-thvh[:,x-1]*(np.log(dt+thvh[:,x-1])-np.log(thvh[:,x-1])))-dz
        
        vpts[:,x]=tprc
        out[:,x]=(tprc-thvh[:,x])/(tenv)*9.8+(T1-thvh[:,x-1])*9.8/tenv
        #out[mc,x]=0
    return out,vpts


# %%
def vspeed(thvh,thvc,W,A,dz,r,P,Tk,Tk2,cur=1,a=1,l=20000):
    vpts=np.zeros(thvh.shape)
    ur=cur*l**(1/2)*9.8**(1/2)*((thvh-thvc)/300)**a
    print(np.mean(ur[500:700,5:10]))
    BE=np.zeros(thvh.shape)
    tprcp=Tk[:,0]
    tprc=Tk[:,0]
    vpts[:,0]=tprc
    urt=np.cumsum(ur,axis=1)
    tm=60
    #A=np.zeros(thvh.shape)
    ket=.5*ur**3*W*dz/A*tm
    rr=r[:,0].copy()
    w=np.zeros(thvh.shape)
    w[:,0]=np.sqrt(2*ket[:,0]/tm)
    #mc=np.zeros(rr.shape,dtype=bool)
    for x in range(1,thvh.shape[1]):
        tenv=(thvh[:,x]+thvh[:,x-1])/2#(thvh[:,x-1]/(100000/P[:,x-1])**.286+thvh[:,x]/(100000/P[:,x])**.286)
        #rr=r[:,x]
        rh=RH(rr,P[:,x],tprcp)
        #mc=((tprc>thvh[:,x])&(rh>1))|mc
        rr2=rr.copy()
        rr2[rh>1]=rr2[rh>1]/rh[rh>1]
        ad=9.808*(1+(2501000*rr2/rh)/(287*tprcp))/(1003.5+(2501000**2*rr2/rh)/(461.5*tprcp**2))
        ad[rh<1]=.0098
        msk=ur[:,x]>0
        tprc=(tprcp-ad*dz)
        tprc[msk]=((tprcp[msk]-ad[msk]*dz)*urt[msk,x-1]+ur[msk,x]*Tk[msk,x])/urt[msk,x]
        rr[msk]=(rr[msk]*urt[msk,x-1]+ur[msk,x]*r[msk,x])/urt[msk,x]
        T1=tprcp*(1+.61*rr2)*(100000/P[:,x])**.286
        #*(1+.61*r2[:])*(100000/P[:,x])**.286)
        tprcp=tprc.copy()#/(1+.61*r2[:])/(100000/P[:,x])**.286
        tprc=tprc*(1+.61*rr2)*(100000/P[:,x])**.286
        #b=(1+.61*r[:,x-1])*(100000/P[:,x])**.286
        #out[:,x]=b*Tk[:,x-1]*(np.log(dt+thvh[:,x-1])-np.log(thvh[:,x-1]))*dz/dt-\
        #         ad*b*dz**2/dt**2*(dt-thvh[:,x-1]*(np.log(dt+thvh[:,x-1])-np.log(thvh[:,x-1])))-dz
        
        vpts[:,x]=tprc
        BE[:,x]=(tprc-thvh[:,x])/(tenv)*9.8+(T1-thvh[:,x-1])*9.8/tenv
        E=2/tm*((w[:,x-1]*tm*BE[:,x]*dz)+ket[:,x])
        e=.00001
        w[E[:]<0,x]=0
        w[E[:]>0,x]=np.sqrt(E[E>0])
        #out[mc,x]=0
    
    return w


# %%
def vspeed(thvh,thvc,W,A,dz,r,P,Tk,Tk2,cur=1,a=1,l=20000):
    vpts=np.zeros(thvh.shape)
    ur=cur*l**(1/2)*9.8**(1/2)*((thvh-thvc)/300)**a
    print(np.mean(ur[500:700,5:10]))
    BE=np.zeros(thvh.shape)
    tprcp=Tk[:,0]
    tprc=Tk[:,0]
    vpts[:,0]=tprc
    urt=np.cumsum(ur,axis=1)
    tm=60
    #A=np.zeros(thvh.shape)
    ket=.5*ur**3*W*dz/A*tm
    rr=r[:,0].copy()
    w=np.zeros(thvh.shape)
    w[:,0]=np.sqrt(2*ket[:,0]/tm)
    decay=2500
    #mc=np.zeros(rr.shape,dtype=bool)
    for x in range(1,thvh.shape[1]):
        tenv=(thvh[:,x]+thvh[:,x-1])/2#(thvh[:,x-1]/(100000/P[:,x-1])**.286+thvh[:,x]/(100000/P[:,x])**.286)
        #rr=r[:,x]
        rh=RH(rr,P[:,x],tprcp)
        #mc=((tprc>thvh[:,x])&(rh>1))|mc
        rr2=rr.copy()
        rr2[rh>1]=rr2[rh>1]/rh[rh>1]
        ad=9.808*(1+(2501000*rr2)/(287*tprcp))/(1003.5+(2501000**2*rr2)/(461.5*tprcp**2))
        ad[rh<1]=.0098
        msk=ur[:,x]>0
        tprc=(tprcp-ad*dz)
        tprc[msk]=((tprcp[msk]-ad[msk]*dz)*urt[msk,x-1]+ur[msk,x]*Tk[msk,x])/urt[msk,x]
        tprc=(tprc+Tk[:,x]*(dz/decay))/(dz/decay+1)
        rr[msk]=(rr[msk]*urt[msk,x-1]+ur[msk,x]*r[msk,x])/urt[msk,x]
        rr=(rr+r[:,x]*(dz/decay))/(dz/decay+1)
        T1=tprcp*(1+.61*rr2)*(100000/P[:,x-1])**.286
        #*(1+.61*r2[:])*(100000/P[:,x])**.286)
        tprcp=tprc.copy()#/(1+.61*r2[:])/(100000/P[:,x])**.286
        tprc=tprc*(1+.61*rr2)*(100000/P[:,x])**.286
        #b=(1+.61*r[:,x-1])*(100000/P[:,x])**.286
        #out[:,x]=b*Tk[:,x-1]*(np.log(dt+thvh[:,x-1])-np.log(thvh[:,x-1]))*dz/dt-\
        #         ad*b*dz**2/dt**2*(dt-thvh[:,x-1]*(np.log(dt+thvh[:,x-1])-np.log(thvh[:,x-1])))-dz
        
        vpts[:,x]=tprc
        BE[:,x]=(tprc-thvh[:,x])/(tenv)*9.8+(T1-thvh[:,x-1])*9.8/tenv
        E=2/tm*((w[:,x-1]*tm*BE[:,x]*dz)+ket[:,x])
        e=.00001
        w[E[:]<0,x]=0
        w[E[:]>0,x]=np.sqrt(E[E>0])
        #out[mc,x]=0
    
    return w


# %%
def dTt2(thvh,thvc,W,A,dz,r,P,Tk,Tk2,cur=1,a=1,l=20000,decay=2500):
    vpts=np.zeros(thvh.shape)
    ur=cur*l**(1/2)*9.8**(1/2)*((thvh-thvc)/300)**a
    out=np.zeros(thvh.shape)
    tprcp=Tk[:,0]
    tprc=Tk[:,0]
    urt=np.cumsum(ur,axis=1)
    rr=r[:,0].copy()
    vpts[:,0]=tprc[:]*(1+.61*rr)*(100000/P[:,0])**.286
    #mc=np.zeros(rr.shape,dtype=bool)
    for x in range(1,thvh.shape[1]):
        tenv=(thvh[:,x]+thvh[:,x-1])/2#(thvh[:,x-1]/(100000/P[:,x-1])**.286+thvh[:,x]/(100000/P[:,x])**.286)
        #rr=r[:,x]
        rh=RH(rr,P[:,x],tprcp)
        #mc=((tprc>thvh[:,x])&(rh>1))|mc
        rr2=rr.copy()
        rr2[rh>1]=rr2[rh>1]/rh[rh>1]
        ad=9.808*(1+(2501000*rr2)/(287*tprcp))/(1003.5+(2501000**2*rr2)/(461.5*tprcp**2))
        ad[rh<1]=.0098
        msk=ur[:,x]>0
        tprc=(tprcp-ad*dz)
        #tprc[msk]=((tprcp[msk]-ad[msk]*dz)*urt[msk,x-1]+ur[msk,x]*Tk[msk,x])/urt[msk,x]
        tprc=(tprc+Tk[:,x]*(dz/decay))/(dz/decay+1)
        #rr[msk]=(rr[msk]*urt[msk,x-1]+ur[msk,x]*r[msk,x])/urt[msk,x]
        rr=(rr+r[:,x]*(dz/decay))/(dz/decay+1)
        T1=tprcp*(1+.61*rr2)*(100000/P[:,x])**.286
        #*(1+.61*r2[:])*(100000/P[:,x])**.286)
        tprcp=tprc.copy()#/(1+.61*r2[:])/(100000/P[:,x])**.286
        tprc=tprc*(1+.61*rr2)*(100000/P[:,x])**.286
        #b=(1+.61*r[:,x-1])*(100000/P[:,x])**.286
        #out[:,x]=b*Tk[:,x-1]*(np.log(dt+thvh[:,x-1])-np.log(thvh[:,x-1]))*dz/dt-\
        #         ad*b*dz**2/dt**2*(dt-thvh[:,x-1]*(np.log(dt+thvh[:,x-1])-np.log(thvh[:,x-1])))-dz
        
        vpts[:,x]=tprc[:]
        out[:,x]=(tprc-thvh[:,x])/(tenv)*9.8+(T1-thvh[:,x-1])*9.8/tenv
        #out[mc,x]=0
    return out,vpts


# %%
def height(thvh,thvc,W,A,dz,r,P,Tk,Tk2,cur=1,a=1,l=20000,decay=2500):
    vpts=np.zeros(thvh.shape)
    ur=cur*l**(1/2)*9.8**(1/2)*((thvh-thvc)/300)**a
    out=np.zeros(thvh.shape)
    tprcp=Tk[:,0]
    tprc=Tk[:,0]
    urt=np.cumsum(ur,axis=1)
    rr=r[:,0].copy()
    vpts[:,0]=tprc[:]*(1+.61*rr)*(100000/P[:,0])**.286
    #mc=np.zeros(rr.shape,dtype=bool)
    for x in range(1,thvh.shape[1]):
        tenv=(thvh[:,x]+thvh[:,x-1])/2#(thvh[:,x-1]/(100000/P[:,x-1])**.286+thvh[:,x]/(100000/P[:,x])**.286)
        #rr=r[:,x]
        rh=RH(rr,P[:,x],tprcp)
        #mc=((tprc>thvh[:,x])&(rh>1))|mc
        rr2=rr.copy()
        rr2[rh>1]=rr2[rh>1]/rh[rh>1]
        ad=9.808*(1+(2501000*rr2)/(287*tprcp))/(1003.5+(2501000**2*rr2)/(461.5*tprcp**2))
        ad[rh<1]=.0098
        msk=ur[:,x]>0
        tprc=(tprcp-ad*dz)
        #tprc[msk]=((tprcp[msk]-ad[msk]*dz)*urt[msk,x-1]+ur[msk,x]*Tk[msk,x])/urt[msk,x]
        tprc=(tprc+Tk[:,x]*(dz/decay))/(dz/decay+1)
        #rr[msk]=(rr[msk]*urt[msk,x-1]+ur[msk,x]*r[msk,x])/urt[msk,x]
        rr=(rr+r[:,x]*(dz/decay))/(dz/decay+1)
        T1=tprcp*(1+.61*rr2)*(100000/P[:,x])**.286
        #*(1+.61*r2[:])*(100000/P[:,x])**.286)
        tprcp=tprc.copy()#/(1+.61*r2[:])/(100000/P[:,x])**.286
        tprc=tprc*(1+.61*rr2)*(100000/P[:,x])**.286
        #b=(1+.61*r[:,x-1])*(100000/P[:,x])**.286
        #out[:,x]=b*Tk[:,x-1]*(np.log(dt+thvh[:,x-1])-np.log(thvh[:,x-1]))*dz/dt-\
        #         ad*b*dz**2/dt**2*(dt-thvh[:,x-1]*(np.log(dt+thvh[:,x-1])-np.log(thvh[:,x-1])))-dz
        
        6
        acc=(tprc-thvh[:,x])/(tenv)*9.8+(T1-thvh[:,x-1])*9.8/tenv
        #out[mc,x]=0
    return out,vpts


# %%

altl=np.linspace(0,5000,166)
les_1c ='/home/tsw35/tyche/data/LES_1C/'
import warnings
warnings.filterwarnings("ignore")
save=True
j= 6
for day in daysa:
    print('.',end='')
    if j == 6:
        #plt.figure(figsize=(15,7))
        plt.figure(figsize=(15,7))
        j=0
    i = 5
    lest=6*6+(i-1)*6
    clbt=8*60+(i-1)*60
    dir2c = clubb_dir+'sgp_2c_d/sgp_'+str(day)+'/'
    dir1c = clubb_dir+'sgp_1c_d/sgp_'+str(day)+'/'
    dircp = clubb_dir+'sgp_2c_d/sgp_'+str(day)+'/'
    try:
        ftcpc = nc.Dataset(dircp+'k_2/c_1/output/arm_zt.nc','r')
        ftcph = nc.Dataset(dircp+'k_2/c_2/output/arm_zt.nc','r')
        fc    = nc.Dataset(dircp+'k_2/clusters.nc','r')
    except:
        continue
    alt=ftcph['altitude'][0:150]

    plt.subplot(2,3,j+1)
    #data=np.ftcpc['thvm'][:,0:150,0,0]-ftcph['thvm'][:,0:150,0,0]
    #abmax=1.5
    #im=plt.imshow(data.T,origin='lower',cmap='coolwarm',vmin=-abmax,vmax=abmax,extent=(5,22,0,6))
    thvc=ftcpc['thvm'][:,0:150,0,0]
    thvh=ftcph['thvm'][:,0:150,0,0]
    Tk=ftcph['T_in_K'][:,0:150,0,0]
    Tk2=ftcpc['T_in_K'][:,0:150,0,0]
    dz=40
    W=fc['W'][0,1]
    A=fc['frac'][1]*100000**2
    P=ftcph['p_in_Pa'][:,0:150,0,0]
    r_=ftcph['rtm'][:,0:150,0,0]
    out,vpts=dTt(thvh,thvc,W,A,dz,r_,P,Tk,Tk2,decay=3000)#np.cumsum(thvh-thvc,axis=1)
    abmax=np.nanpercentile(np.abs(out[0:780,0:50]),99)
    abmin=-abmax
    #out[out<=0]=float('nan')
    #im=plt.imshow(out[0:780,:].T,vmax=.1,vmin=-.1,origin='lower',cmap='coolwarm',extent=(5,18,0,6))
    #im=plt.imshow(out.T,vmax=abmax,origin='lower',cmap='terrain',extent=(5,18,0,6))
    #plt.colorbar()
    plt.plot(thvc[600,:],alt)
    plt.plot(thvh[600,:],alt)
    plt.plot(vpts[600,:],alt)
    plt.title(str(day)+': '+str(clbt/60+5))
    j=j+1
    if j == 6:
        if save:
            save=False
        else:
            break

# %%
#plt.plot(thvc[200,0:75],alt[0:75])
plt.plot(thvh[780,:],alt)
plt.plot(vpts[780,:],alt)
#plt.plot(vpts[780,0:75]-thvh[780,0:75],alt[0:75]) 
#plt.xlim(305,315)

# %%
np.sqrt(-1)

# %%

altl=np.linspace(0,5000,166)
les_1c ='/home/tsw35/tyche/data/LES_1C/'
import warnings
warnings.filterwarnings("ignore")

j= 6
for day in daysa:
    print('.',end='')
    if j == 6:
        #plt.figure(figsize=(15,7))
        plt.figure(figsize=(15,7))
        j=0
    i = 5
    lest=6*6+(i-1)*6
    clbt=8*60+(i-1)*60
    dir2c = clubb_dir+'sgp_2c_d/sgp_'+str(day)+'/'
    dir1c = clubb_dir+'sgp_1c_d/sgp_'+str(day)+'/'
    dircp = clubb_dir+'sgp_cpl_d4/sgp_'+str(day)+'/'
    try:
        ftcpc = nc.Dataset(dircp+'k_2/c_1/output/arm_zt.nc','r')
        ftcph = nc.Dataset(dircp+'k_2/c_2/output/arm_zt.nc','r')
        fc    = nc.Dataset(dircp+'k_2/clusters.nc','r')1
    except:
        continue
    alt=ftcph['altitude'][0:150]

    plt.subplot(2,3,j+1)
    #data=np.ftcpc['thvm'][:,0:150,0,0]-ftcph['thvm'][:,0:150,0,0]
    #abmax=1.5
    #im=plt.imshow(data.T,origin='lower',cmap='coolwarm',vmin=-abmax,vmax=abmax,extent=(5,22,0,6))
    thvc=ftcpc['thvm'][:,0:150,0,0]
    thvh=ftcph['thvm'][:,0:150,0,0]
    Tk=ftcph['T_in_K'][:,0:150,0,0]
    dz=40
    W=fc['W'][0,1]
    A=fc['frac'][1]*100000**2
    P=ftcph['p_in_Pa'][:,0:150,0,0]
    data,k,kh,b,ur=vspeed(thvh,thvc,W,A,dz,ftcph['rtm'][:,0:150,0,0],P,Tk)#np.cumsum(thvh-thvc,axis=1)
    data[data<=0]=float('nan')
    #ur[ur>=0]=float('nan')
    abmax=np.nanpercentile(np.abs(ur[0:780,:]),98)
    abmin=-abmax #0*np.nanpercentile(data[0:780,:],2)
    im=plt.imshow(ur[0:780,:].T,origin='lower',cmap='coolwarm',vmin=abmin,vmax=abmax,extent=(5,18,0,6))
    plt.colorbar()
    
    #abmax=np.percentile(np.abs(ur[:,0:100]),98)
    #plt.imshow(ur[:,0:100].T,origin='lower',vmin=-abmax,vmax=abmax,cmap='coolwarm',extent=(5,22,0,6))
    #plt.imshow(kh[500:550,0:12].T,origin='lower',cmap='coolwarm',extent=(0,50,0,50))
    #plt.colorbar()
    #plt.plot(k[600,0:100-1],alt[1:100])
    #plt.plot(kh[600,:100],alt[:100])
    #plt.plot(b[600,:100],alt[:100])
    #plt.legend(['V','X','B'])
    #plt.plot(np.linspace(5,22,len(thvh[:,0])),thvh[:,0])
    
    #dt=np.abs(b[600,:100]/kh[600,:100])
    #plt.plot(dt,alt[:100])
    #plt.xlim(np.nanpercentile(dt,10),np.nanpercentile(dt,90))
    
    plt.title(str(day)+': '+str(clbt/60+5))
    j=j+1
    #if j == 6:
    #    break

# %%
A

# %%
testx=np.array([-1000,-100,-10,-1,-.1,-.01,-.001,-.0001,0,.0001,.001,.01,.1,1,10,100,1000])/5
roots=np.zeros((3,len(testx)))
dt=60
B=-.01
for x in testx:
    a=np.roots([dt,dz,-(2*B*dz+2*x)*dt,-2*x*dz])
    print(str(x)+': '+str(a))

# %%
x=-1.9
a=np.roots([dt,dz,-(2*B*dz+2*x)*dt,-2*x*dz])
print(a)

# %%
np.roots([60,dz,0,0])

# %%
4*.2**2-2**3

# %%
4-2

# %%
plt.imshow(ur.T,origin='lower',cmap='coolwarm',vmin=np.min(np.abs(ur)),vmax=np.max(np.abs(ur)))
plt.colorbar()


# %%

altl=np.linspace(0,5000,166)
les_1c ='/home/tsw35/tyche/data/LES_1C/'
import warnings
warnings.filterwarnings("ignore")

j= 6
for day in daysa:
    if j == 6:
        #plt.figure(figsize=(15,7))
        plt.figure(figsize=(15,7))
        j=0
    i = 5
    lest=6*6+(i-1)*6
    clbt=8*60+(i-1)*60
    dir2c = clubb_dir+'sgp_2c_d/sgp_'+str(day)+'/'
    dir1c = clubb_dir+'sgp_1c_d/sgp_'+str(day)+'/'
    dircp = clubb_dir+'sgp_cpl_d4/sgp_'+str(day)+'/'
    try:
        ftcpc = nc.Dataset(dircp+'k_2/c_1/output/arm_zt.nc','r')
        ftcph = nc.Dataset(dircp+'k_2/c_2/output/arm_zt.nc','r')
        fscpc = nc.Dataset(dircp+'k_2/c_1/output/arm_sfc.nc','r')
        fscph = nc.Dataset(dircp+'k_2/c_2/output/arm_sfc.nc','r')
        #fc    = nc.Dataset(dircp+'k_2/clusters.nc','r')
    except:
        continue
    alt=ftcph['altitude'][0:150]

    plt.subplot(2,3,j+1)
    thvc=ftcpc['T_in_K'][:,0:150,0,0]
    tm=np.linspace(5,22,len(thvc[:,0]))
    #thvh=ftcph['T_in_K'][:,0:150,0,0]
    plt.plot(tm,fscpc['sh'][:,0,0,0])
    plt.plot(tm,fscph['sh'][:,0,0,0])
    
    plt.title(str(day)+': '+str(clbt/60+5))
    j=j+1

# %%
tm

# %%

# %%
ass = np.linspace(1,2,11)
print(ass)
xs=np.linspace(0,1)
for a in ass:
    #plt.plot(xs,((xs/300)**(a)-((1/300)**a*(xs)**(1))))
    plt.plot(xs,(xs/300)**(a))
    #plt.plot(xs,(1/300)**a*(xs)**(1))

# %%
ass = np.linspace(1,2,11)
print(ass)
xs=np.linspace(0,1)
for a in ass:
    plt.plot(xs,(1/300)**a*(xs)**(1))

# %%
dk=np.mean(ftcpc['T_in_K'][:,0:10,0,0]-ftcph['T_in_K'][:,0:10,0,0],axis=1)
dv=np.mean(ftcpc['thvm'][:,0:10,0,0]-ftcph['thvm'][:,0:10,0,0],axis=1)
plt.plot(dk)
plt.plot(dv)

# %%
diffs=[]
for day in daysa:
    print('.',end='')
    dir2c = clubb_dir+'sgp_2c_d/sgp_'+str(day)+'/'
    dir1c = clubb_dir+'sgp_1c_d/sgp_'+str(day)+'/'
    dircp = clubb_dir+'sgp_cpl_d4/sgp_'+str(day)+'/'
    

    ft2cc = nc.Dataset(dir2c+'k_2/c_1/output/arm_zt.nc','r')
    ft2ch = nc.Dataset(dir2c+'k_2/c_2/output/arm_zt.nc','r')
    diffs.append(np.mean(ft2ch['thvm'][:,0:25,0,0]-ft2cc['thvm'][:,0:25,0,0]))
    

# %%
plt.hist(diffs)

# %%
from sklearn import metrics
from scipy import stats


# %%
def compute_LE(LE,clst_):
    LEclst=np.zeros((LE.shape[0],2))
    for t in range(LE.shape[0]):
        LEclst[t,0]=np.mean(LE[t][clst_==0])
        LEclst[t,1]=np.mean(LE[t][clst_==1])
    return LEclst


# %%
def get_lhet(run_dir):
    fp_=open(run_dir+'/tw_run_param.txt','r')
    for line in fp_:
        if 'l_het' in line:
            return float(line.split(' ')[4])
    return 5


# %%
from scipy import stats

# %%
ln=len(v_comps.keys())
vlst=list(v_comps.keys())
cors=np.zeros((ln,ln))
for i in range(ln):
    v1=vlst[i]
    for j in range(ln):
        v2=vlst[j]
        cors[i,j]=stats.spearmanr(v_comps[v1],v_comps[v2])[0]


# %%
fig, ax = plt.subplots(figsize=(5,6))
im = ax.imshow(cors.T,cmap='bwr',vmin=-1,vmax=1)
plt.colorbar(im)
plt.grid(False)
ax.set_yticks(np.arange(len(v_comps_l)))
ax.set_yticklabels(v_comps_l)
ax.set_xticks(np.arange(len(v_comps_l)))
ax.set_xticklabels(v_comps_l)
plt.setp(ax.get_xticklabels(), rotation=45, ha="right",
     rotation_mode="anchor")
print()

# %%
lhet=[]
for day in daysa:
    
    lhet.append(np.sqrt(get_lhet(clubb_dir+'sgp_cpl_d4/sgp_'+day)))

# %%
plt.hist(lhet)

# %%

# %%
np.sqrt(32968)

# %%
plt.plot(lhet)

# %%
for i in range(ln):
    v1=vlst[i]
    plt.figure()
    plt.scatter(v_comps[v1][:],v_comps['drtp2'][:])
    plt.ylim(np.min(v_comps['drtp2'][:]),np.max(v_comps['drtp2'][:]))
    if v1=='drtp2':
        plt.xlim(np.min(v_comps[v1][:]),np.max(v_comps[v1][:]))
    plt.title(v1)

# %%
fccp.variables.keys()

# %%
plt.plot(fccp['u_r'][:])


# %%
def RH(r,P,T):
    e=r*P/(.62197+r)
    esat=.61078*np.exp((17.27*(T-273.15))/(T-273.15+237.3))*1000
    #print(str(e)+','+str(esat))
    return(e/esat)


# %%
def adprof(r,P,T0,alt,T):
    tout=[]
    tout.append(T0*(100000/P[0])**.286)
    ad=.0098
    ads=[]
    Tprev=T0
    rhs=[]
    for i in range(len(r)):
        if i ==0:
            continue
        Tprev = Tprev-ad*(alt[2]-alt[1])
        rh=RH(r[0],P[i],Tprev)
        rhs.append(rh)
        if rh>=1:
            ad=9.808*(1+(2501000*r[0]/rh)/(287*Tprev))/(1003.5+(2501000**2*r[0]/rh)/(461.5*Tprev**2))
            ads.append(ad)
            tout.append((Tprev)*(100000/P[i])**.286*(1+.61*r[0]/rh))
        else:
            ad=.0098
            tout.append((Tprev)*(100000/P[i])**.286*(1+.61*r[0]))
        #print(rh)
    tout=np.array(tout)
    #print(np.mean(rhs))
    #print(np.mean(ads))
    return tout


# %%
plt.plot(Ttc*(1+.61*rtmc)*(100000/pac)**.286-vptc)

# %%

altl=np.linspace(0,5000,166)
les_1c ='/home/tsw35/tyche/data/LES_1C/'

j= 6
for day in daysa:
    print('.',end='')
    if j == 6:
        #plt.figure(figsize=(15,7))
        plt.figure(figsize=(18,5))
        j=0
    i = 5
    lest=4*6+(i-1)*6
    clbt=6*60+(i-1)*60
    dir2c = clubb_dir+'sgp_2c_d/sgp_'+str(day)+'/'
    dir1c = clubb_dir+'sgp_1c_d/sgp_'+str(day)+'/'
    dircp = clubb_dir+'sgp_cpl_d4/sgp_'+str(day)+'/'
    try:
        ftcpc = nc.Dataset(dircp+'k_2/c_1/output/arm_zt.nc','r')
        ftcph = nc.Dataset(dircp+'k_2/c_2/output/arm_zt.nc','r')
    except:
        continue
    
    alt=ftcpc['altitude'][0:150]
    Ttc=np.mean(ftcpc['T_in_K'][clbt-5:clbt+5,0:150,0,0],axis=0)
    Tth=np.mean(ftcph['T_in_K'][clbt-5:clbt+5,0:150,0,0],axis=0)
    Tc=np.mean(ftcpc['thlm'][clbt-5:clbt+5,0:150,0,0],axis=0)
    Th=np.mean(ftcph['thlm'][clbt-5:clbt+5,0:150,0,0],axis=0)
    rtmc=np.mean(ftcpc['rtm'][clbt-5:clbt+5,0:150,0,0],axis=0)
    rtmh=np.mean(ftcph['rtm'][clbt-5:clbt+5,0:150,0,0],axis=0)
    pac=np.mean(ftcpc['p_in_Pa'][clbt-5:clbt+5,0:150,0,0],axis=0)
    pah=np.mean(ftcph['p_in_Pa'][clbt-5:clbt+5,0:150,0,0],axis=0)
    
    
    #Tca = adprof(rtmc,pac,Ttc[0],alt,Ttc)
    Tha = adprof(rtmh,pah,Tth[0],alt,Tth)
    
    vptc=Tc*(1+.61*rtmh)
    vpth=Th*(1+.61*rtmh)
    #vptca=Tca*(1+.61*rtmc[0])
    vptha=Tha

    plt.subplot(1,6,j+1)
    #plt.plot((vptc+vpth)/2,alt,'grey',alpha=.5)
    plt.plot(vptc,alt,'b-')
    #plt.plot(vptca,alt,'b--')
    plt.plot(vpth,alt,'r-')
    plt.plot(vptha,alt,'r--')
    plt.title(str(day)+': '+str(clbt/60+5))
    #plt.ylim(0,5000)
    #plt.xlim(vptc[0]-1,vptc[0]+20)
    #ce=np.cumsum((vptha-vpth)/vpth*40)
    #plt.plot(ce,alt)
    #plt.xlim(-2,2)

    
    j=j+1

# %%
.0098*(1-rtmc[0]*.85)

# %%
a=np.mean(ftcpc['thlm'][clbt-5:clbt+5,0:150,0,0],axis=0)/np.mean(ftcpc['T_in_K'][clbt-5:clbt+5,0:150,0,0],axis=0)

# %%
vptca[0]-vptc[0]

# %%
(100000/pac[0])**(.286)

# %%
j=0
for day in [20160625]:
    i = 5
    lest=6*6+(i-1)*6
    clbt=8*60+(i-1)*60
    dir2c = clubb_dir+'tall_2c/sgp_'+str(day)+'/'
    dir1c = clubb_dir+'tall_1c/sgp_'+str(day)+'/'
    dircp = clubb_dir+'tall2_cpl/sgp_'+str(day)+'/'
    #dircp2 = clubb_dir+'sgp_cpl_dhi/sgp_'+str(day)+'/'
    try:
        fmcp = nc.Dataset(dircp+'k_2/agg_outzm.nc','r')
        #fmcpc = nc.Dataset(dircp+'k_2/c_1/output/arm_zm.nc','r')
        #fmcph = nc.Dataset(dircp+'k_2/c_2/output/arm_zm.nc','r')
        #fm2c = nc.Dataset(dir2c+'k_2/agg_outzm.nc','r')
        fm1c = nc.Dataset(dir1c+'k_1/c_1/output/arm_zm.nc','r')
        #fplm=nc.Dataset(les_1c+'trimfr2_'+str(day)+'_01.nc','r')
        #fplt=nc.Dataset(les_1c+'trimfr2_'+str(day)+'_00.nc','r')
        fscp=nc.Dataset(dircp+'k_2/agg_outsfc.nc','r')
        #fscp2=nc.Dataset(dircp+'k_2/agg_outsfc.nc','r')
        fs1c=nc.Dataset(dir1c+'k_1/c_1/output/arm_sfc.nc','r')
    except:
        continue
    plt.figure(figsize=(4.5,4.5),dpi=300)
    
    lwp2=moving_average(fscp['lwp'][:,0,0,0])
    lwp1=moving_average(fs1c['lwp'][:,0,0,0])
    #lwpht=fplt['lwp'][:]
    #lwphm=fplm['lwp'][:]
    plt.plot(np.linspace(5,22,lwp1.shape[0]),lwp2,'b-',linewidth=3)
    plt.plot(np.linspace(5,22,lwp1.shape[0]),lwp1,'g--',linewidth=3)
    plt.ylabel(r"LWP $(kg/m^{-2})$",fontsize=20,weight='bold')
    plt.yticks(fontsize=18)
    plt.legend(['Circulation','No Circulation'],fontsize=18)
    #plt.title(day)
    if j == 0:
        #plt.title(r"LWP",fontsize=30,weight='bold')
        #plt.yticks(np.arange(0,.0011,.0003),fontsize=18,rotation=45)
        plt.xlabel(r"Time $(hr)$",fontsize=24,weight='bold')
        plt.xticks(np.arange(5,25,2),fontsize=18)
        plt.xlim(12,21)
        plt.ylim(-0.005,.18)
   
    
    j=j+1

# %%

# %%
