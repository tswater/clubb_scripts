
import numpy as np
import os
import datetime

sgp_dir = '/stor/soteria/hydro/shared/data/tylersclutter/surfaces201589/'
day_dir = '/stor/soteria/hydro/shared/data/tylersclutter/doz_tylertrim/'
sgp_dir2= '/stor/soteria/hydro/shared/lasso_for_tyler/'

#### READ IN SURFACE DATA ####
def read_sfc_data(var,nt,stdt,override='X'):
    nx=401
    bnm = 'jss'+var+'_bdy_02_' #base name
    varout_g = np.zeros((nt,nx,nx))
    varout_v = np.zeros((nt,nx*nx))
    if nt == 1:
        file=override
        with open(file,'r') as fp:
            varout_v = np.array([float(i) for i in fp.readlines()])
        varout_g=np.reshape(varout_v,(nx,nx))
        return varout_g,varout_v
    for t in range(nt):
        dt = stdt+datetime.timedelta(hours=t)
        tnm = dt.strftime('%Y-%m-%d-%H-%M')
        file = sfc_dir+bnm+tnm
        with open(file,'r') as fp:
            var_v = np.array([float(i) for i in fp.readlines()])
        varout_g[t,:,:] = np.reshape(var_v,(nx,nx))
        varout_v[t,:] = var_v[:]
    return varout_g,varout_v

#################


#################

def get_lengthscale(data,dirc=-1,samp=5000,cut=.25):
    #dirc: direction; -1 means no direction
    a,b = data.shape
    dx=250
    x_H=np.zeros((a,b))
    y_H=np.zeros((a,b))
    for i in range(a):
        for j in range(b):
            x_H[i,j]=i*dx
            y_H[i,j]=j*dx

    H_flat = data.flatten()
    x_flat = x_H.flatten()
    y_flat = y_H.flatten()
    idx = np.random.choice(len(H_flat),size=samp,replace=False)
    H_sb = H_flat[idx]
    x_Hsb = x_flat[idx]
    y_Hsb = y_flat[idx]
    mu=np.mean(H_sb)
    a = (H_sb[:,np.newaxis].T - mu)*(H_sb[:,np.newaxis]-mu)
    h = ((x_Hsb[:,np.newaxis] - x_Hsb.T)**2 + \
        (y_Hsb[:,np.newaxis] - y_Hsb.T)**2)**0.5
    Qf = a.flatten()
    hf = h.flatten()
    
    bins = np.linspace(0,75000,76)
    means=np.zeros((len(bins)-1,))
    for i in range(len(bins)-1):
        means[i]=np.mean(Qf[(hf>bins[i])&(hf<bins[i+1])])

    return bins[0:-1][means<=(.25*means[0])][0]

#################

test_time = '15'
lens = []
for file in os.listdir(day_dir):
    datest = file[8:12]+'-'+file[12:14]+'-'+file[14:16]+'-'+test_time
    date = datetime.datetime(int(file[8:12]),int(file[12:14]),int(file[14:16]),int(test_time))
    if file [18]=='1':
        try:
            Hg,Hv = read_sfc_data('sh',1,date,sgp_dir+'jsssh_bdy_02_'+datest+'-00')
        except:
            Hg,Hv = read_sfc_data('sh',1,date,sgp_dir2+'sgp_'+file[8:16]+\
                                  '_00/'+'jsssh_bdy_02_'+datest+'-00')
        lens.append(get_lengthscale(Hg))
        print(datest+':  '+str(lens[-1]))

print()
print('MEDIAN:'+str(np.median(lens)))
print('MEAN:  '+str(np.mean(lens)))
print('MAX:   '+str(np.max(lens)))
print('MIN:   '+str(np.min(lens)))

