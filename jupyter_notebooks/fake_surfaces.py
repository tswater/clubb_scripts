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
import scipy as sci
import scipy.ndimage
import re
import os
import rasterio
import seaborn as sns
import matplotlib.patheffects as pe
import matplotlib as mpl
mpl.rcParams['figure.dpi'] = 300
sns.set_theme()
plt.rcParams.update({'figure.max_open_warning': 0})

# %%

# %%
ex_1=np.zeros((5,5))
ex_1[0,:]=[5,3,2,1.2,2.7]
ex_1[1,:]=[6,5,2.2,2,3]
ex_1[2,:]=[6,3,1.8,.8,2]
ex_1[3,:]=[6.5,4.4,1.5,1,1.7]
ex_1[4,:]=[5.2,6.2,4,2,2]
#ex_1=7-ex_1
ex_1=ex_1-np.min(ex_1)
ex_1=ex_1/np.max(ex_1)

# %%
#ex_2=np.zeros((5,5))
#ex_2[0,:]=[.25,.25,.25,.25,.25]+(np.random.rand(1,5)/2-.25)
#ex_2[1,:]=[.75,.75,.25,.75,.25]+(np.random.rand(1,5)/2-.25)
#ex_2[2,:]=[.25,.25,.25,.25,.25]+(np.random.rand(1,5)/2-.25)
#ex_2[3,:]=[.75,.75,.75,.25,.25]+(np.random.rand(1,5)/2-.25)
#ex_2[4,:]=[.75,.75,.75,.75,.25]+(np.random.rand(1,5)/2-.25)

# %%

# %%
#### STEP 1 ####
plt.figure(figsize=(3,5))
plt.subplot(1,2,1)
plt.imshow(ex_1,cmap='coolwarm')
plt.xticks([-.5,.5,1.5,2.5,3.5,4.5],[])
plt.yticks([-.5,.5,1.5,2.5,3.5,4.5],[])
plt.title('Example 1',fontsize=14)

plt.subplot(1,2,2)
plt.imshow(ex_2,cmap='coolwarm')
plt.xticks([-.5,.5,1.5,2.5,3.5,4.5],[])
plt.yticks([-.5,.5,1.5,2.5,3.5,4.5],[])
plt.title('Example 2',fontsize=14)


# %%
delta_1=0
i_1=0
i_2=0
cut_1=0
cut_2=0
delta_2=0
for i in range(50,80):
    cut=np.percentile(ex_1,i)
    delta=np.mean(ex_1[ex_1>cut])-np.mean(ex_1[ex_1<=cut])
    if delta>delta_1:
        cut_1=cut
        delta_1=delta
        i_1=i
        
    cut=np.percentile(ex_2,i)
    delta=np.mean(ex_2[ex_2>cut])-np.mean(ex_2[ex_2<=cut])
    if delta>delta_2:
        cut_2=cut
        delta_2=delta
        i_2=i
p_1=np.zeros((5,5))
p_1[ex_1>cut_1]=np.mean(ex_1[ex_1>cut_1])
p_1[ex_1<=cut_1]=np.mean(ex_1[ex_1<=cut_1])

p_2=np.zeros((5,5))
p_2[ex_2>cut_2]=np.mean(ex_2[ex_2>cut_2])
p_2[ex_2<=cut_2]=np.mean(ex_2[ex_2<=cut_2])

print(str(cut_1)+': '+str(i_1)+': '+str(delta_1))
print(str(cut_2)+': '+str(i_2)+': '+str(delta_2))

# %%
#### STEP 2 ####
plt.figure(figsize=(3,5))
plt.subplot(1,2,1)
plt.imshow(p_1,cmap='coolwarm',vmin=np.min(ex_1),vmax=np.max(ex_1))
plt.xticks([-.5,.5,1.5,2.5,3.5,4.5],[])
plt.yticks([-.5,.5,1.5,2.5,3.5,4.5],[])
plt.plot([-.5,.5,.5,1.5,1.5],[.5,.5,3.5,3.5,4.5],linewidth=4,color='gold',
        path_effects=[pe.Stroke(linewidth=5, foreground='w'), pe.Normal()])
plt.title('Example 1',fontsize=14)

plt.subplot(1,2,2)
plt.imshow(p_2,cmap='coolwarm',vmin=np.min(ex_2),vmax=np.max(ex_2))
plt.xticks([-.5,.5,1.5,2.5,3.5,4.5],[])
plt.yticks([-.5,.5,1.5,2.5,3.5,4.5],[])
plt.plot([-.5,1.5,1.5,-.5],[.5,.5,1.5,1.5],linewidth=4,color='gold',
         path_effects=[pe.Stroke(linewidth=5, foreground='w'), pe.Normal()])
plt.plot([2.5,3.5,3.5,2.5,2.5],[.5,.5,1.5,1.5,.5],linewidth=4,color='gold',
         path_effects=[pe.Stroke(linewidth=5, foreground='w'), pe.Normal()])
#plt.plot([-.5,.5,.5],[3.5,3.5,4.5],linewidth=5,color='gold',
#         path_effects=[pe.Stroke(linewidth=6, foreground='w'), pe.Normal()])
plt.plot([-.5,2.5,2.5,3.5,3.5],[2.5,2.5,3.5,3.5,4.5],linewidth=4,color='gold',
         path_effects=[pe.Stroke(linewidth=5, foreground='w'), pe.Normal()])
         

plt.title('Example 2',fontsize=14)

# %%
ah_1=np.sum(p_1>cut_1)
ac_1=np.sum(p_1<=cut_1)
ah_2=np.sum(p_2>cut_2)
ac_2=np.sum(p_2<=cut_2)
L_1=6
L_2=15
hh_1=int(100*np.round(ah_1/L_1,2))
hc_1=int(100*np.round(ac_1/L_1,2))
hh_2=int(100*np.round(ah_2/L_2,2))
hc_2=int(100*np.round(ac_2/L_2,2))

r_1=np.zeros((L_1*100,hh_1+hc_1))
r_2=np.zeros((L_2*100,hh_2+hc_2))
r_1[:,0:hh_1+1]=p_1[p_1>cut_1][0]
r_1[:,hh_1+1:-1]=p_1[p_1<=cut_1][0]
r_2[:,0:hh_2+1]=p_2[p_2>cut_2][0]
r_2[:,hh_2+1:-1]=p_2[p_2<=cut_2][0]

# %%
print(r_1.shape)
print(r_2.shape)
print(hh_2)

# %%
###### STEP 3 ####
a=3
fig=plt.figure(figsize=(5.3/3*a,1*a))
plt.subplot(2,1,1)
plt.imshow(r_1.T,cmap='coolwarm',vmin=np.min(ex_1),vmax=np.max(ex_1))
plt.title('Example 1',fontsize=14)
plt.yticks([16,116,216,316,416],[])
plt.xticks(np.linspace(0,L_1*100,L_1+1),[])
plt.plot([0,L_1*100],[83,83],linewidth=4,color='gold',
         path_effects=[pe.Stroke(linewidth=5, foreground='w'), pe.Normal()])
plt.xlim(0,600)
#plt.xlim(0,1700)
#plt.ylim(0,416)

plt.subplot(2,1,2)
plt.imshow(r_2.T,cmap='coolwarm',vmin=np.min(ex_2),vmax=np.max(ex_2))
plt.title('Example 2',fontsize=14)
plt.yticks([67,167],[])
plt.xticks(np.linspace(0,L_2*100,L_2+1),[])
plt.plot([0,L_2*100],[67,67],linewidth=4,color='gold',
         path_effects=[pe.Stroke(linewidth=5, foreground='w'), pe.Normal()])
plt.xlim(0,1500)
#plt.ylim(0,416)
plt.show()


# %%
ac=np.min(p_1)
ah=np.max(p_1)
lc=[ah*.9,ah*.9,ah*.95,ah*1.3,ah*1.4,ah*1.5,ah*1.6]
lh=[ah,ah,ah,ah*1.1,ah*1.3,ah*1.5,ah*1.6]

hs=[0,1.2,1.5,2,2.3,3,4]

# %%
#### STEP 4 ####

plt.figure(figsize=(2,3.5))
ylabels=[r'$0$',r'$z_{circ}$',r'$z_{crit}$',r'$z_{lim_2}$',r'$z_{lim_1}$']
xlabels=[r'$\theta_{crit}$',r'$\theta_{lim}$']
plt.plot(lh,hs,c='indianred')
plt.plot(lc,hs,c='cornflowerblue')
plt.yticks([0,1,1.6,2,2.3],ylabels,fontsize=16)
plt.ylim(-.1,3.1)
plt.xticks([ah*1.02,ah*1.3],xlabels,fontsize=16)

# %%
ah*1.3

# %%
#### STEP 5 ####
plt.figure(figsize=(2,2))
plt.imshow(p_2,cmap='coolwarm',vmin=np.min(ex_2),vmax=np.max(ex_2))
plt.xticks([-.5,.5,1.5,2.5,3.5,4.5],[])
plt.yticks([-.5,.5,1.5,2.5,3.5,4.5],[])
plt.title(r'$u_{R_0}=3$',fontsize=14)
dx=.35
dy=.35
wd=.15
hl=.2
cf=.25

plt.arrow(0,0+cf,0,dy,width=wd,head_length=hl,facecolor='black')
plt.arrow(0,2-cf,0,-dy,width=wd,head_length=hl,facecolor='black')
plt.arrow(1,0+cf,0,dy,width=wd,head_length=hl,facecolor='black')
plt.arrow(1,2-cf,0,-dy,width=wd,head_length=hl,facecolor='black')
plt.arrow(3,0+cf,0,dy,width=wd,head_length=hl,facecolor='black')
plt.arrow(3,2-cf,0,-dy,width=wd,head_length=hl,facecolor='black')

plt.arrow(0,2+cf,0,dy,width=wd,head_length=hl,facecolor='black')
plt.arrow(1,2+cf,0,dy,width=wd,head_length=hl,facecolor='black')
plt.arrow(2,2+cf,0,dy,width=wd,head_length=hl,facecolor='black')
plt.arrow(3,3+cf,0,dy,width=wd,head_length=hl,facecolor='black')

plt.arrow(2+cf,1,dx,0,width=wd,head_length=hl,facecolor='black')
plt.arrow(2-cf,1,-dx,0,width=wd,head_length=hl,facecolor='black')
plt.arrow(4-cf,1,-dx,0,width=wd,head_length=hl,facecolor='black')
plt.arrow(4-cf,4,-dx,0,width=wd,head_length=hl,facecolor='black')
plt.arrow(3-cf,3,-dx,0,width=wd,head_length=hl,facecolor='black')

# %%
#### STEP 5 ####
plt.figure(figsize=(2,2))
plt.imshow(p_2,cmap='coolwarm',vmin=np.min(ex_2),vmax=np.max(ex_2))
plt.xticks([-.5,.5,1.5,2.5,3.5,4.5],[])
plt.yticks([-.5,.5,1.5,2.5,3.5,4.5],[])
plt.title(r'$u_{b}$=(0,4)',fontsize=14)
dx=.35
dy=.35
wd=.15
hl=.2
cf=.25

plt.arrow(2+cf,1,dx,0,width=wd,head_length=hl,facecolor='black')
plt.arrow(2-cf,1,-dx,0,width=wd,head_length=hl,facecolor='black')
plt.arrow(4-cf,1,-dx,0,width=wd,head_length=hl,facecolor='black')
plt.arrow(4-cf,4,-dx,0,width=wd,head_length=hl,facecolor='black')
plt.arrow(3-cf,3,-dx,0,width=wd,head_length=hl,facecolor='black')

# %%
#### STEP 7 ####
plt.figure(figsize=(4,1))
plt.imshow(r_2.T,cmap='coolwarm',vmin=np.min(ex_2),vmax=np.max(ex_2))
plt.title(r'$u_R=1,\ \ \ \ \ \ w_R(0)=\frac{u_R*IW*dz}{A_{hot}}$',fontsize=14,pad=10)
plt.yticks([67,167],[])
plt.xticks(np.linspace(0,L_2*100,L_2+1),[])
wd2=20
hl2=25
plt.plot([0,L_2*100],[67,67],linewidth=4,color='gold',
         path_effects=[pe.Stroke(linewidth=5, foreground='w'), pe.Normal()])
plt.xlim(0,1500)

plt.arrow(150,120-cf*100,0,-dy*100,width=wd2,head_length=hl2,facecolor='black',zorder=2)
plt.arrow(450,120-cf*100,0,-dy*100,width=wd2,head_length=hl2,facecolor='black',zorder=2)
plt.arrow(750,120-cf*100,0,-dy*100,width=wd2,head_length=hl2,facecolor='black',zorder=2)
plt.arrow(1050,120-cf*100,0,-dy*100,width=wd2,head_length=hl2,facecolor='black',zorder=2)
plt.arrow(1350,120-cf*100,0,-dy*100,width=wd2,head_length=hl2,facecolor='black',zorder=2)
plt.show()

# %%
#### STEP 7 ####
plt.figure(figsize=(4,1))
plt.imshow(r_2.T,cmap='coolwarm',vmin=np.min(ex_2),vmax=np.max(ex_2))
plt.title(r'$\lambda$',fontsize=14,pad=10)
plt.yticks([67,167],[])
plt.xticks(np.linspace(0,L_2*100,L_2+1),[])
wd2=20
hl2=25
plt.plot([0,L_2*100],[67,67],linewidth=4,color='gold',
         path_effects=[pe.Stroke(linewidth=5, foreground='w'), pe.Normal()])
plt.xlim(0,1500)

plt.arrow(150,120-cf*100,0,-dy*100,width=wd2,head_length=hl2,facecolor='black',zorder=2)
plt.arrow(450,120-cf*100,0,-dy*100,width=wd2,head_length=hl2,facecolor='black',zorder=2)
plt.arrow(750,120-cf*100,0,-dy*100,width=wd2,head_length=hl2,facecolor='black',zorder=2)
plt.arrow(1050,120-cf*100,0,-dy*100,width=wd2,head_length=hl2,facecolor='black',zorder=2)
plt.arrow(1350,120-cf*100,0,-dy*100,width=wd2,head_length=hl2,facecolor='black',zorder=2)
plt.show()

# %%

# %%
