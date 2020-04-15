import numpy as np
import csv
import matplotlib.pyplot as plt
import json
import math as m
import os
import pandas as pd
import scipy.interpolate as sp
from scipy import interpolate
test = pd.read_csv('/Users/bhagyasubrayan/Desktop/Plastic/plasticc_test_lightcurves_01.csv')
f_0 = 27.5
dmod = input('Enter the distance modulus:')
true_peak = input('Enter the mjd for true peak:')
model_num = input('Enter the model number:')
#Working with Takashi's models
list = os.listdir("/Users/bhagyasubrayan/Desktop/Plastic/public_data/")
lines = open('/Users/bhagyasubrayan/Desktop/Plastic/public_data/'+ list[0]).readlines()
a = open('modelnames.txt', 'w').writelines(lines[2:])
with open('modelnames.txt', 'r') as in_file:
    stripped = (line.strip() for line in in_file)
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
#print(model_arr)
group = model_arr.groupby(['Explosion Energy','Mass Loss','CSM edge','Beta']).get_group(('1e51','1e-2','1e15','5'))
chi_m_arr = pd.DataFrame()
chi_m_arr['Model name'] = group['Model Name']
chi_m_arr['P Mass'] = group['Progenitor mass']
chi_m_arr = chi_m_arr.reset_index(drop=True)
print(chi_m_arr)
sl_list = []
#Collecting magnitudes from these models in group
for i in range(0,1):
    for k in range(0,len(group)):
        mod = pd.read_table('mod'+'_'+ group.iloc[k,0]+'.sdss2.txt', names = ['epoch','u','g','r','i','z','kepler'], sep='\s+')
        t_g = test.groupby(['object_id', 'passband']).get_group((int(model_num),i))
        print(t_g)
        am_m_i =[]
        ndays_i = []
        new_i = t_g.filter(['mjd','flux','flux_err'])
        new_i = new_i.reset_index(drop=True)
        new_i = new_i.drop(new_i[new_i.mjd < (int(true_peak) - 150) ].index)
        #print(new_i)
        for j in range(len(new_i)):
            if (new_i.iloc[j]['flux'] < 0 ):
                new_i.iloc[j]['flux'] = 0.01
            a1 = -2.5*m.log((new_i.iloc[j]['flux']/ f_0),10)
            #print(type(m.log(new_i.iloc[j]['flux']/f_0)))
            abs = a1 - 5*m.log((float(dmod)/10),10)
            am_m_i.append(abs)
        new_i['Abs_mag'] = am_m_i
        new_i.reset_index(drop= True)
        #print((new_i))
        tag = new_i.loc[new_i['Abs_mag'].idxmin()]
        #print(tag[0])
        for k in range(0,len(new_i)):
            d = new_i.iloc[k]['mjd'] - tag[0]
            #print(d)
            ndays_i.append(d)
        #print(len(ndays))
        new_i['t_max_afterdays'] = ndays_i
        print(new_i)
        x_i = new_i['t_max_afterdays'].values
        y_obs_chi_i = new_i['Abs_mag'].values
        stdev_i = y_obs_chi_i.std()
        print(stdev_i)
        err = new_i['flux_err']
        #plt.errorbar(x_i,y_obs_chi_i,yerr = err, fmt = 'o')
        #plt.title('PlasTicc Light curves for a IIP in different bands')
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
        fil_max= arr.loc[arr['Magnitude'].idxmin()]
        for k in range(0,len(arr)):
           d1 = arr.iloc[k]['Epoch'] - fil_max[0]
           #print(d)
           pdays_i.append(d1)
        xnew = np.linspace(pdays_i[0],pdays_i[-1],500)
        print(y_obs_chi_i)
        print(f1(x_i))
        f1_i = interpolate.interp1d(pdays_i, mag, kind='cubic')
        y_pred_i = f1_i(x_i)
        plt.plot(x_i,y_obs_chi_i,'o',x_i,f1(x_i),'*')
        plt.gca().invert_yaxis()
    plt.show()
