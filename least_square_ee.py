import numpy as np
import csv
import matplotlib.pyplot as plt
import json
import math as m
import os
import pandas as pd
import scipy.interpolate as sp
from scipy import interpolate
test = pd.read_csv('/Users/bhagyasubrayan/Desktop/Explosion Parameters/Plasticc/plasticc_test_lightcurves_01.csv')
f_0 = 27.5
dmod = 40.977
true_peak = 60499.461
#Working with Takashi's models
with open('/Users/bhagyasubrayan/Desktop/Explosion Parameters/modelnames.txt') as f:
    stripped = (line.strip() for line in f)
    lines = (line.split() for line in stripped if line)
    with open('model.csv', 'w') as out_file:
        writer = csv.writer(out_file)
        #writer.writerow(('name', 'mass'))
        writer.writerows(lines)
f = open('model.csv','r')
content = f.read()
splitcontent = content.split()
#print(type(splitcontent[9]))
n = []
mass = []
exp_ener = []
mass_loss = []
csm_edge =[]
beta=[]
p_days=[]
for i in range(len(splitcontent)):
    a2 = splitcontent[i].split(',')
    a2.remove('Msun')
    a2.remove('erg')
    if( len(a2) >5):
        a2.remove('cm')
        a2.remove('Msun/yr')
    else:
        a2.append(0)
        a2[-2] = 0
        a2[-3] = 0
    n.append(a2[0])
    mass.append(a2[1])
    exp_ener.append(a2[2])
    mass_loss.append(a2[3])
    csm_edge.append(a2[4])
    beta.append(a2[5])
    model_arr = pd.DataFrame()
    model_arr['Model Name'] = n
    model_arr['Progenitor mass'] = mass
    model_arr['Explosion Energy'] = exp_ener
    model_arr['Mass Loss'] = mass_loss
    model_arr['CSM edge'] = csm_edge
    model_arr['Beta'] = beta
print(model_arr)
group = model_arr.groupby(['Progenitor mass','Mass Loss','CSM edge','Beta']).get_group(('14','1e-2','1e15','5'))
print(group)
chi_m_arr = pd.DataFrame()
chi_m_arr['Model name'] = group['Model Name']
chi_m_arr['Explosion Energy'] = group['Explosion Energy']
chi_m_arr = chi_m_arr.reset_index(drop=True)
print(chi_m_arr)
sl_list = []
path = "/Users/bhagyasubrayan/Desktop/Explosion Parameters/Takashi's models/multicolorlightcurves/"
#Collecting magnitudes from these models in group
type(group.iloc[1,0])
len(group)
for k in range(0,len(group)):
    chib_list_k = []
    m_k = group.iloc[k,0]
    mod = pd.read_table(path + group.iloc[k,0]+'.sdss2', skiprows = 2,names = ['epoch','u','g','r','i','z','kepler'], sep='\s+')
    #print(mod)
#Filtering and making sense of the observed data
    for i in range(0,6):
        t_g = test.groupby(['object_id', 'passband']).get_group((13,i))
        #print(t_g)
        am_m_i =[]
        mag_err_i = []
        ndays_i = []
        new_i = t_g.filter(['mjd','flux','flux_err'])
        new_i = new_i.drop(new_i[new_i.mjd < (true_peak - 150) ].index)
        #print(new_i)
        new_i = new_i.drop(new_i[new_i.flux < 0].index)
        new_i = new_i.reset_index(drop=True)
        #print(new_i)
        for j in range(len(new_i)):
            a1 = -2.5*np.log10(new_i.iloc[j]['flux'])+f_0
            abs = a1 - dmod
            mag_err = np.abs(2.5*(new_i.iloc[j]['flux_err']/new_i.iloc[j]['flux'])*np.log(10))
            am_m_i.append(abs)
            mag_err_i.append(mag_err)
        new_i['Abs_mag'] = am_m_i
        new_i['Mag_error'] = mag_err_i
        print(new_i)
        new_i = new_i.drop(new_i[new_i.Mag_error > 10].index)
        new_i = new_i.reset_index(drop= True)
        print((new_i))
        tag = new_i.loc[new_i['Abs_mag'].idxmin()]
        #print(tag[0])
        for k in range(0,len(new_i)):
            d = new_i.iloc[k]['mjd'] - true_peak
            #print(d)
            ndays_i.append(d)
        #print(len(ndays))
        new_i['t_max_afterdays'] = ndays_i
        print(new_i)
        x_i = new_i['t_max_afterdays'].values
        y_obs_chi_i = new_i['Abs_mag'].values
        stdev_i = y_obs_chi_i.std()
        #print(stdev_i)
        #plt.errorbar(x_i,y_obs_chi_i, yerr = new_i['flux_err'],fmt= 'o')
        #plt.gca().invert_yaxis()
        #plt.show()
        err = new_i['Mag_error']
        plt.errorbar(x_i,y_obs_chi_i,yerr = err, fmt = 'ro')
        #plt.title('PlasTicc Light curves for a IIP in bands'+str(i))
        #plt.xlabel('days after explosion')
        #plt.ylabel('Absolute Magnitude ')
        #plt.gca().legend(('u','g','r','i','z','y'))
        pdays_i=[]
        b_arr = mod.iloc[:,[0,i+1]].drop_duplicates()
        g = b_arr.groupby('epoch').mean()
        epo = b_arr['epoch'].drop_duplicates()
        mag = g.values.squeeze()
        arr = pd.DataFrame()
        arr['Epoch'] = epo
        arr['Magnitude'] = mag
        #print(arr)
        #fil_max= arr.loc[arr['Magnitude'].idxmin()]
        #print(fil_max)
        #for k in range(0,len(arr)):
           #d1 = arr.iloc[k]['Epoch'] - fil_max[0]
           #print(d)
           #pdays_i.append(d1)
        #xnew = np.linspace(pdays_i[0],pdays_i[-1],500)
        #print(len(pdays_j))
        #print(len(mag))
        f1 = interpolate.interp1d(epo, mag, kind='cubic')
        y_pred_i = f1(x_i)
        new_i['pred_mag'] = y_pred_i
        #print(new_i)
        #plt.plot(x_i,f1(x_i),'*')
        #plt.gca().legend(('u','g','r','i','y','kepler'))
        #print((y_obs_chi))
        #print((y_pred_j))
        chi_i=[]
        c = 0
        for q in range(0,len(y_obs_chi_i)):
           c += ((y_obs_chi_i[q] - y_pred_i[q]) / stdev_i)**2
        chib_list_k.append(c)
        #print(chib_list_k)
        plt.scatter(x_i,y_pred_i,color = 'g')
        plt.legend(['Predicted','Observed'])
        plt.title('Comparison in band'+ str(i)+m_k)
        plt.gca().invert_yaxis()
        #plt.ylim(-10,-50)
        plt.show()
    sl_list.append(chib_list_k)
print(len(sl_list))
sl_array = pd.DataFrame(sl_list,columns = ('u','g','r','i','y','kepler'))
print(sl_array)
result = pd.concat([chi_m_arr,sl_array],axis=1)
print(result)
#plt.gca().invert_yaxis()
for i in range(0,6):
    plt.plot(result['Explosion Energy'],result.iloc[:,i+2],'o')
    plt.xlabel('Explosion Energy')
    plt.ylabel(u'$\u03C7^2$')
    plt.title(u'$\u03C7^2$'+'\t'+'Test')
    plt.gca().legend(('u','g','r','i','y','kepler'))
plt.show()
