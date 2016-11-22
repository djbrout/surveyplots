import numpy as np
import matplotlib.pyplot as plt
import array
import math
from astropy.cosmology import FlatLambdaCDM
cosmo = FlatLambdaCDM(H0=70, Om0=0.3)

fig, ax = plt.subplots(2,1,sharex=True)

pos=[]
for i in range(1,299):
    pos.append(i/100.0)


list1, idsurvey,z1, mass1,x11,c1,mu1, mu1e,mures1 = np.loadtxt('../DATA/SALT2mu_fitopt.fitres', usecols=(0,3, 7,13,20,22,37,39,41), unpack=True, dtype='string', skiprows=12)
z1 = z1.astype(float)
mu1 = mu1.astype(float)
mu1e = mu1e.astype(float)


idsurvey2,z2,mu2,mu2e=np.loadtxt('/project/rkessler/dscolnic/HST_analysis/FITOPT000+SALT2mu.FITRES',usecols=(3,7,35,37), unpack=True, dtype='string', skiprows=12)
z2 = z2.astype(float)
mu2 = mu2.astype(float)
mu2e = mu2e.astype(float)
print z2, mu2,mu2e
print idsurvey2


z1=np.append(z1,z2)
mu1=np.append(mu1,mu2)
mu1e=np.append(mu1e,mu2e)
idsurvey=np.append(idsurvey,idsurvey2)



x=cosmo.luminosity_distance(pos).value
y=5.0*(np.log10(x))+25.0

x=cosmo.luminosity_distance(z1).value
mu_syn=5.0*(np.log10(x))+25.0
line, = ax[0].plot(pos, y, lw=2)
line, = ax[1].plot(pos, y*0, lw=2)


survey=1
tempsurv=[]
tempz=[]
tempmu=[]
tempmue=[]
uniqsurv=set(idsurvey)
uniqsurv=['61','62','63','64','65','66', '5','1', '4', '101','106','100','15']
labels=['CFA1','CFA2','CFA3S','CFA3K','CFA4p1','CFA4p2','CSP','SDSS','SNLS','SNAP','CANDELS','GOODS','PS1']
xpos=[0.03,0.03,0.03,0.03,0.03,0.03,0.03,.2,.45,.6,.8,.8,.45,.35,.4]
ypos=[37,38.2,39.4,40.6,41.8,43,44.2,45.6,36.8,39,40,41,42,43,44,43,45,46,43,42,43]
fo=[14,14,14,14,14,14,14,20,20,20,20,20,20,20]
col=['Red', 'Purple','DarkRed','RosyBrown','LightSalmon','Pink','DarkKhaki','Green','DimGray','Gold','Orange','Purple','Blue']

co2=0
for uniqsu in uniqsurv:
        xx=np.where((idsurvey==uniqsu))        
        if len(xx[0]>1):
            print 'co2', len(xx[0])
            tempz=z1[xx[0]]
            tempmu=mu1[xx[0]]
            tempmu2=mu1[xx[0]]-mu_syn[xx[0]]
            tempmue=mu1e[xx[0]]
            #print 'col2', co2, len(col)
            ax[0].errorbar(tempz, tempmu, yerr=tempmue, fmt='o', color=col[co2],ecolor=col[co2])
            ax[0].text(xpos[co2],ypos[co2],labels[co2],fontdict={'fontsize':fo[co2]}, color=col[co2])
                      
        co2=co2+1
    #stop


bins = np.linspace(0, 1, 11)
bins=np.append(bins,[1.15,1.4,1.7,2.3])

digitized = np.digitize(z1, bins)
bin_means = [mu1[digitized == i].mean() for i in range(0, len(bins))]
bin_z = [z1[digitized == i].mean() for i in range(0, len(bins))]
bin_std = [np.std(mu1[digitized == i])/np.sqrt(len(mu1[digitized == i])) for i in range(0, len(bins))]

        #ax1.errorbar(bin_z, bin_means, yerr=bin_std, fmt='ko', ecolor='k')
        
        
    #plt.annotate('local max', xy=(.5, 43), xytext=(3, 1.5),
        #            arrowprops=dict(facecolor='black', shrink=0.05),
        #            )
    
        
ax[0].set_ylabel('m-M (mag)')
ax[0].set_xlim(0.01,2)
ax[0].set_xscale('log')

digitized = np.digitize(z1, bins)
bin_means = [(mu1[digitized == i]-mu_syn[digitized == i]).mean() for i in range(0, len(bins))]
bin_z = [z1[digitized == i].mean() for i in range(0, len(bins))]
bin_std = [np.std(mu1[digitized == i]-mu_syn[digitized == i])/np.sqrt(len(mu1[digitized == i])) for i in range(0, len(bins))]

#ax2.errorbar(bin_z, bin_means, yerr=bin_std, fmt='ko', ecolor='k')




idsurvey,z1, mu1, mu1e = np.loadtxt('DJ_PS1.fitres', usecols=(2, 6,44,45), unpack=True, dtype='string', skiprows=12)
z1 = z1.astype(float)
mu1 = mu1.astype(float)
mu1e = mu1e.astype(float)


idsurvey2,z2,mu2,mu2e=np.loadtxt('/project/rkessler/rhounsell/IFU_SDT_COSMO/IFU/FITOPT000+SALT2mu.FITRES',usecols=(3,6,45,46), unpack=True, dtype='string', skiprows=12)
z2 = z2.astype(float)
mu2 = mu2.astype(float)
mu2e = mu2e.astype(float)
idsurvey2[idsurvey2=='15']='61'

z1=np.append(z1,z2)
mu1=np.append(mu1,mu2)
mu1e=np.append(mu1e,mu2e)
idsurvey=np.append(idsurvey,idsurvey2)

idsurvey3,z3,nn_itype3,mu3,mu3e,zphot3 = np.loadtxt('/project/rkessler/rkessler/DES/analysis/NN_SALT2/out02_SN+hostZspec//NNSTAGE07_FINALGRID0_fitDATA+SIM3/NNFIT_DATA//FITOPT000+SALT2mu.FITRES', usecols=(3,7, 15,37,39,27), unpack=True, dtype='string', skiprows=20)
print np.sum(nn_itype3=='1'), np.sum(nn_itype3!='1')
z3=z3.astype(float)
mu3=mu3.astype(float)
mu3e=mu3e.astype(float)

print 'mu1',mu1
print 'mu2',mu2
print 'mu3',mu3



print z2, mu2,mu2e
print idsurvey3


z1=np.append(z1,z3)
mu1=np.append(mu1,mu3)
mu1e=np.append(mu1e,mu3e)
idsurvey=np.append(idsurvey,idsurvey3)



x=cosmo.luminosity_distance(pos).value
y=5.0*(np.log10(x))+25.0

x=cosmo.luminosity_distance(z1).value
mu_syn=5.0*(np.log10(x))+25.0
line, = ax[0].plot(pos, y, lw=2)
line, = ax[1].plot(pos, y, lw=2)

        

survey=1
tempsurv=[]
tempz=[]
tempmu=[]
tempmue=[]
uniqsurv=set(idsurvey)
print uniqsurv

uniqsurv=['15','10','103']
labels=['PS1','DES','WFIRST']
xpos=[.2,.45,.6,.8,.8,.45,.35,.4]
ypos=[39,40,41,42,43,44,43,45,46,43,42,43]
fo=[14,14,14,14,14,14,14,20,20,20,20,20,20,20]
col=['Green','DimGray','Gold','Orange','Purple','Blue']

co2=0
for uniqsu in uniqsurv:
            xx=np.where((idsurvey==uniqsu))
            if len(xx[0]>1):
                print 'co2', len(xx[0])
                tempz=z1[xx[0]]
                tempmu=mu1[xx[0]]
                tempmu2=mu1[xx[0]]-mu_syn[xx[0]]
                tempmue=mu1e[xx[0]]
                #print 'col2', co2, len(col)
                ax[1].errorbar(tempz, tempmu, yerr=tempmue, fmt='o', color=col[co2],ecolor=col[co2])
                ax[1].text(xpos[co2],ypos[co2],labels[co2],fontdict={'fontsize':fo[co2]}, color=col[co2])
            co2=co2+1

ax[1].set_ylabel('m-M (mag)')
ax[1].set_xlim(0.01,2)
ax[1].set_xscale('log')
ax[1].set_ylim(32,48)

plt.tight_layout()
plt.show()
                                                            
    
plt.savefig('hubble_nw.png')
    
