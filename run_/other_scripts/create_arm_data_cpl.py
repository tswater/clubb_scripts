import xarray as xr
import datetime
import numpy as np

# ----------- #
# USER INPUTS #
# ----------- #
in_file = '/home/tsw35/soteria/clubb/data/sgo60varanarap.nc' #input data
model_file = '../../original_arm_model.in' # input file for temporal data
out_dir = '../input/case_setups/' #output location


# --------- #
# FUNCTIONS #
# --------- #

#### LOAD IN INFORMATION FROM MODEL.IN ####
def model_run_info():
 fp = open(model_file,'r')
 for line in fp:
  if line[0:3]=='day':
   day   = int(line[8:10])
  elif line[0:5]=='month':
   month = int(line[8:10])
  elif line[0:4]=='year':
   year  = int(line[8:12])
  elif line[0:12]=='time_initial':
   t_i   = float(line[15:24])
  elif line[0:10]=='time_final':
   t_f   = float(line[13:24])+3600 #(add a buffer_
  else:
   continue
 return day,month,year,t_i,t_f



#### PREPARE SURFACE DATA ####
def prepare_surface_data(fp,date,idate):

 #Write the data as requested by CLUBB
 lh = np.array(fpi['LH'])
 sh = np.array(fpi['SH'])
 tskin = np.array(fpi['T_skin'])

 #Write out the data
 fp = open(out_dir+'arm_sfc.in','w')

 #Write hearder information
 fp.write('! $Id$\n')
 fp.write('!\n')
 fp.write('! Note that T_sfc is included here, but it is not currently used. It is\n')
 fp.write('! included here to facilitate a transition to using T_sfc in the future\n')
 fp.write('! if needed.\n')
 fp.write('Time[s]    latent_ht[W\m^2]   sens_ht[W\m^2]   T_sfc[K]\n')

 dt = datetime.timedelta(hours=1)
 for it in range(len(lh)):
  timestamp = (date - idate).total_seconds()
  tmp = '%.2f %.2f %.2f %.2f\n' % (timestamp,lh[it],sh[it],tskin[it]+273.15)
  fp.write(tmp)
  date += dt

 return


#### PREPARE FORCING DATA ####
def prepare_forcing_data(fpi,date,idate,fdate):

 #Write the data as requested by CLUBB
 lev = 100*np.array(fpi['lev']) #Pa
 dTdt = np.array(fpi['dTdt'])[:,:]/3600.0 # K/hour -> K/s
 dqdt = np.array(fpi['dqdt'])[:,:]/1000.0/3600.0 #g/kg/hour -> kg/kg/s
 u = np.array(fpi['u'])[:,:]
 v = np.array(fpi['v'])[:,:]
 omega = np.array(fpi['omega'])[:,:]*100/3600.0 #mb/hr -> Pa/s

 #Write out the data
 fp = open(out_dir+'arm_forcings.in','w')
 #fp = open('tmp.txt','w')

 #Write hearder information
 fp.write('! $Id$\n')
 fp.write('! The vertical coordinate is entered in the first column as height (z[m]) or pressure (Press[Pa]).\n')
 fp.write('! Temperature can be entered as liquid water potential temperature (thlm_f[K\s]), potential temperature (thm_f[K\s]), or absolute temperature(T_f[K\s]).\n')
 fp.write('! Prescribed, time-independent vertical velocity can be entered as velocity (w[m\s]) or pressure velocity (omega[Pa\s] or omega[mb\hr])\n')
 fp.write('Press[Pa]       T_f[K\s]        rtm_f[kg\kg\s]   um_ref[m\s]   vm_ref[m\s]      um_f[m\s^2]     vm_f[m\s^2]     omega[Pa\s]          ug[m\s]         vg[m\s]\n')
 #fp.write('Press[Pa]       T_f[K\s]        rtm_f[kg\kg\s]   um_ref[m\s]   vm_ref[m\s]      um_f[m\s^2]     vm_f[m\s^2]     w[m\s]          ug[m\s]         vg[m\s]\n')
 dt = datetime.timedelta(hours=1)
 it = 0
 while date < fdate:
  print(date)
  timestamp = (date - idate).total_seconds()
  tmp = '%.2f %d\n' % (timestamp,lev.size)
  fp.write(tmp)
  for il in np.arange(lev.size):
   tmp = '%.4f %.4e %.4e %.4f %.4f %.4f %.4f %.4f %.4f %.4f\n' % (lev[il],dTdt[it,il],dqdt[it,il],
         u[it,il],v[it,il],-999.9,-999.9,omega[it,il],-999.9,-999.9)
   #      u[it,il],v[it,il],-999.9,-999.9,1.0,u[it,il],v[it,il])
   fp.write(tmp)
  date += dt
  it += 1

 return



#### PREPARE SOUNDING DATA ####
def prepare_sounding_data(fpi):

 #Write the data as requested by CLUBB
 lev = 100*np.array(fpi['lev']) #Pa
 T = np.array(fpi['T'])[0,:] # absolute temperature (T[K])
 r = np.array(fpi['q'])[0,:]/1000.0 #g/kg -> kg/kg
 u = np.array(fpi['u'])[0,:]
 v = np.array(fpi['v'])[0,:]
 omega = np.array(fpi['omega'])[0,:]*100/3600.0 #mb/hr -> Pa/s
 print(lev.size,T.size,r.size,u.size,v.size,omega.size)

 #Write out the data
 fp = open(out_dir+'arm_sounding.in','w')

 #Write hearder information
 fp.write('! The vertical coordinate is entered in the first column as height (z[m]) or pressure (Press[Pa]).\n')
 fp.write('! Temperature can be entered as liquid water potential temperature (thlm[K]), potential temperature (thm[K]), or absolute temperature(T[K]).\n')
 fp.write('! Prescribed, time-independent vertical velocity can be entered as velocity (w[m\s]) or pressure velocity (omega[Pa\s])\n')
 fp.write('! If vertical velocity is time dependent, insert missing values throughout column.\n')
 fp.write('! In all cases, missing values are denoted by -999.9\n')
 #fp('z[m]     thm[K]   rt[kg\kg]   u[m\s]   v[m\s]   w[m\s]   ug[m\s]   vg[m\s]\n')
 fp.write('Press[Pa]     T[K]   rt[kg\kg]   u[m\s]   v[m\s]   omega[Pa\s]   ug[m\s]   vg[m\s]\n')
 #fp.write('Press[Pa]     T[K]   rt[kg\kg]   u[m\s]   v[m\s]   w[m\s]   ug[m\s]   vg[m\s]\n')
 for il in np.arange(lev.size):
  tmp = '%.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f\n' % (lev[il],T[il],r[il],u[il],v[il],omega[il],-999.9,-999.9)
  fp.write(tmp)

 #Close the file
 fp.close()

 return


# ------------- #
# ACTUAL SCRIPT #
# ------------- #

#Determine Model Run Period
day,month,year,t_i,t_f=model_run_info()

tdt_i = datetime.datetime(year,month,day)+datetime.timedelta(seconds=t_i)
tdt_f = datetime.datetime(year,month,day)+datetime.timedelta(seconds=t_f)
tdt_b = datetime.datetime(year,month,day,0)

tstr_i= tdt_i.strftime('%Y-%m-%dT%H:%M:%S')
tstr_f= tdt_f.strftime('%Y-%m-%dT%H:%M:%S')

print(tstr_i)
print(tstr_f)

#Extract the data for the desired time period
fpi = xr.open_dataset(in_file)
fpi = fpi.loc[dict(time=slice(tstr_i, tstr_f))]

#Prepare the surface data
prepare_surface_data(fpi,tdt_b,tdt_i)

#Prepare the sounding data
prepare_sounding_data(fpi)

#Prepare the forcing data
prepare_forcing_data(fpi,tdt_i,tdt_b,tdt_f)
