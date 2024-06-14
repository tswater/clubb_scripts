import numpy as np
import netCDF4 as nc
import os
import copy
import argparse
# ----------- #
# USER INPUTS #
# ----------- #
main_dir = '/home/tsw35/tyche/clubb/test_cpl/'
force_agg = True
dt = 60
out_var = {'zm': {'latitude':[],'longitude':[],'altitude':[],'time':[],\
                  'thlp2_ta':[],'thlp2_tp':[],'thlp2_dp1':[],'thlp2_dp2':[],\
                  'thlp2_ma':[],'thlp2_forcing':[],'thlp2':[],'wp2':[],'up2':[],\
                  'vp2':[],'Richardson_num':[],'bv_freq_sqd':[],'wprtp':[],\
                  'wpthlp':[],'rtp2':[]},\
        'zt': {'latitude':[],'longitude':[],'altitude':[],'time':[],'thlm':[],\
                  'um':[],'vm':[],'p_in_Pa':[],'rtm':[],'rvm':[],'T_in_K':[],'thvm':[],\
                  'ug':[],'vg':[],'wm':[],'cloud_cover':[],'rho':[],'rcm':[],\
                  'wp3':[],'thlp3':[],'rtp3':[],'wpthlp2':[],'wp2thlp':[],\
                  'wprtp2':[],'wp2rtp':[],'Skw_zt':[],'thlm_ma':[],'rtm_ma':[],\
                  'rtm_forcing':[],'thlm_forcing':[],'mixt_frac':[],'w_1':[],
                  'varnce_w_1':[],'varnce_w_2':[],'varnce_thl_1':[],'varnce_thl_2':[],\
                  'varnce_rt_1':[],'varnce_rt_2':[],'w_2':[],'thl_1':[],'thl_2':[],
                  'rt_1':[],'rt_2':[],'corr_rt_thl_1':[],'corr_rt_thl_2':[]},\
        'sfc':{'latitude':[],'longitude':[],'altitude':[],'time':[],'sh':[],\
               'lh':[],'lwp':[],'ustar':[],'cc':[],'z_cloud_base':[],'T_sfc':[],\
               'z_inversion':[]}}
remove_L0 = True #remove restart files

#### ARG PARSER ####
prs = argparse.ArgumentParser(description='Short sample app')
prs.add_argument('-i', action='store', dest='dirct', default=main_dir)
args = prs.parse_args()
main_dir = args.dirct


# ---------------- #
# HELPER FUNCTIONS #
# ---------------- #






# ----------- #
# SETUP STUFF #
# ----------- #
ks = []
rst_cpl = False
agg_clst = True
n_rst = 0




for kfile in os.listdir(main_dir):
    if 'k_' in kfile:
        ki =int(kfile[2:])
        ks.append(ki)
        for cfile in os.listdir(main_dir+kfile):
            if 'agg' in cfile:
                agg_clst = False
                continue
            
            elif cfile[0:2]!='c_':
                continue
            
            if 'restart' in os.listdir(main_dir+kfile+'/'+cfile):
                rst_cpl = True
                old_dir = main_dir+kfile+'/'+cfile+'/output/old/'
                for ncfile in os.listdir(old_dir):
                    sp1 = ncfile.split('_')
                    sp2 = sp1[2].split('.')
                    num = int(sp2[0])
                    if num > n_rst:
                        n_rst = num
# get initial and final time and dt
fp = open(main_dir+'k_'+str(ks[len(ks)-1])+'/c_1/input/arm_model.in','r')
for line in fp.readlines():
    if line[0:12]=='time_initial':
        lp = line.split(' ')
        lp = [x for x in lp if x !='']
        t_init = float(lp[2])
    elif line[0:10]=='time_final':
        lp = line.split(' ')
        lp = [x for x in lp if x !='']
        t_final = float(lp[2])
    elif line[0:7]=='dt_main':
        lp = line.split(' ')
        lp = [x for x in lp if x !='']
        dt = float(lp[2])

print(n_rst)
print(ks)
print(rst_cpl)

# ---------------- # 
# COLLAPSE RESTART #
# ---------------- #
print('Beginning temporal restart aggregation')

if rst_cpl:
    for ki in ks:
        print('K: '+str(ki))
        for ci in range(1,ki+1):
            print('C: '+str(ci))
            ovar=copy.deepcopy(out_var)
            old_dir = main_dir+'k_'+str(ki)+'/c_'+str(ci)+'/output/old/'
            out_dir = main_dir+'k_'+str(ki)+'/c_'+str(ci)+'/output/'
            for f in ovar.keys():
                print('f: '+str(f))
                # iterate through all timesteps
                for i in range(n_rst):
                    try:
                        fp_in = nc.Dataset(old_dir+'arm_'+f+'_'+str(i+1)+'.nc','r')
                    except:
                        continue
                    for var in ovar[f].keys():
                        dims = fp_in[var].dimensions
                        shp = fp_in[var].shape
                        if len(ovar[f][var])==0:
                            if dims[0] == 'time':
                                if len(dims)==1:
                                    ovar[f][var]=np.zeros((shp[0]*n_rst))
                                elif len(dims)==2:
                                    ovar[f][var]=np.zeros((shp[0]*n_rst,shp[1]))
                                elif len(dims)==3:
                                    ovar[f][var]=np.zeros((shp[0]*n_rst,shp[1],shp[2]))
                                elif len(dims)==4:
                                    ovar[f][var]=np.zeros((shp[0]*n_rst,shp[1],shp[2],shp[3]))
                            else:
                                ovar[f][var]=fp_in[var][:]
                        if dims[0] == 'time':
                            i1 = i*shp[0]
                            i2 = (i+1)*shp[0]
                            if len(dims)>1:
                                ovar[f][var][i1:i2,:,:,:] = fp_in[var][0:,:,:,:]
                            else:
                                ovar[f][var][i1:i2] = fp_in[var][0:]
                        else:
                            continue
                
                # create output netcdf
                fpout = nc.Dataset(out_dir+'arm_'+f+'.nc','w')
                for name,dim in fp_in.dimensions.items():
                    if name == 'time':
                        fpout.createDimension(name,n_rst*fp_in['time'].shape[0])
                    else:
                        fpout.createDimension(name,len(dim))
                for var in ovar[f].keys():
                    if ovar[f][var].shape == (1,):
                        continue
                    if var == 'time':
                        fpout.createVariable(var,fp_in[var].datatype,('time'))
                    elif var == 'altitude':
                        fpout.createVariable(var,fp_in[var].datatype,('altitude'))
                    else:
                        fpout.createVariable(var,fp_in[var].datatype,('time','altitude','latitude','longitude'))
                    fpout[var][:]=ovar[f][var][:]
 

# ------------------ #
# AGGREGATE CLUSTERS #
# ------------------ #
print('Beginning cluster aggregation')
if agg_clst or force_agg:
    times = nc.Dataset(main_dir+'k_'+str(ks[0])+'/c_1/output/arm_zm.nc','r')['time'][:]
    alts  = nc.Dataset(main_dir+'k_'+str(ks[0])+'/c_1/output/arm_zm.nc','r')['altitude'][:]
    for ki in ks:
        print('K: '+str(ki))
        ovar=copy.deepcopy(out_var)
        out_dir = main_dir+'k_'+str(ki)+'/'
        
        # read in clusters
        fp_clst = nc.Dataset(main_dir+'k_'+str(ki)+'/clusters.nc','r')
        k_masks_grid = fp_clst['cluster'][:]
        shp = fp_clst['cluster'].shape
        k_masks = np.reshape(k_masks_grid,(shp[0],shp[1]*shp[2]))
        fp_clst.close()
        
        nt = int((t_final-t_init)/3600)+1

        btsmp = [] # timestamp for surface info
        for t in range(nt):
            btsmp.append(t_init+t*3600)
        
        print(nt)

        # generate weights
        wts_0 = np.zeros((shp[0],ki))
        for j in range(ki):
            for t in range(shp[0]):
                wts_0[t,j] = np.sum(k_masks[t,:]==j)/k_masks.shape[1]
        wts_z = np.zeros((len(times),len(alts),ki))
        wts_sfc = np.zeros((len(times),ki))
        idx = 0
        delta = btsmp[1]-btsmp[0]
        for t in range(len(times)):
            if (btsmp[idx]+delta/2)<t:
                idx = idx+1
            for j in range(ki):
                wts_z[t,:,j] = wts_0[idx,j]
                wts_sfc[t,j] = wts_0[idx,j]
        
        # aggregate
        for f in ovar.keys():
            agfile = out_dir+'agg_out'+f+'.nc'
            print('f: '+str(f))
            for ci in range(1,ki+1):
                fp_in = nc.Dataset(main_dir+'k_'+str(ki)+'/c_'+str(ci)+'/output/arm_'+f+'.nc','r')
                if ci == 1:
                    for var in ovar[f].keys():
                        if var in ['latitude','longitude']:
                            continue
                        if (var == 'altitude') and (f=='sfc'):
                            continue
                        shp = fp_in[var].shape
                        if var in ['altitude','time']:
                            ovar[f][var]=fp_in[var][:]
                        else:
                            ovar[f][var]=np.zeros(shp)
                for var in ovar[f].keys():
                    if var in ['latitude','longitude','altitude','time']:
                        continue
                    if f == 'sfc':
                        ovar[f][var][:,0,0,0]=ovar[f][var][:,0,0,0]+fp_in[var][:,0,0,0]*wts_sfc[:,ci-1]
                    else:
                        ovar[f][var][:,:,0,0]=ovar[f][var][:,:,0,0]+fp_in[var][:,:,0,0]*wts_z[:,:,ci-1]
            
            # Fix Variances
            if (ki == 2) and (f=='zm'):
                for varb in ['um','vm','wm','thlm','rtm']:
                    varp=varb[:-1]+'p2'
                    fp1t=nc.Dataset(main_dir+'k_'+str(ki)+'/c_1/output/arm_zt.nc','r')
                    fp2t=nc.Dataset(main_dir+'k_'+str(ki)+'/c_2/output/arm_zt.nc','r')
                    xm1=fp1t[varb][:,:,0,0]
                    xm2=fp2t[varb][:,:,0,0]
                    f1=wts_z[:,:,0]
                    f2=wts_z[:,:,1]
                    ovar[f][varp][:,:,0,0]=ovar[f][varp][:,:,0,0]+f1*(xm1-(f1*xm1+f2*xm2))**2


            # output everything to files
            fp_out=nc.Dataset(agfile,'w')
            fp_out.createDimension('time',size=len(times))
            fp_out.createDimension('latitude',size=1)
            fp_out.createDimension('longitude',size=1)
            if f =='sfc':
                fp_out.createDimension('altitude',size=1)
            else:
                fp_out.createDimension('altitude',size=len(alts))
            dim = ('time','altitude','latitude','longitude')
            fp_out.createVariable('time','d',dimensions=('time'))
            fp_out.createVariable('altitude','d',dimensions=('altitude'))
            for var in ovar[f]:
                if (var == 'altitude') or (var=='time'):
                    fp_out[var][:]=ovar[f][var][:]
                    continue
                if var in ['latitude','longitude']:
                    continue
                fp_out.createVariable(var,'d',dimensions=dim)
                fp_out[var][:]=ovar[f][var][:]
        


