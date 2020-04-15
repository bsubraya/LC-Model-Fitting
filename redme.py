import numpy as np
import csv
import matplotlib.pyplot as plt
import json
import math as m
import os
import pandas as pd
import scipy.interpolate as sp
from scipy import interpolate
from scipy.stats import chisquare
import ast
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
print(type(splitcontent[9]))
n = []
mass = []
exp_ener = []
mass_loss = []
csm_edge =[]
beta=[]
for i in range(len(splitcontent)):
    a = splitcontent[i].split(',')
    a.remove('Msun')
    a.remove('erg')
    if( len(a) >5):
        a.remove('cm')
        a.remove('Msun/yr')
    else:
        a.append(0)
        a[-2] = 0
        a[-3] = 0
    n.append(a[0])
    mass.append(a[1])
    exp_ener.append(a[2])
    mass_loss.append(a[3])
    csm_edge.append(a[4])
    beta.append(a[5])
    model_arr = pd.DataFrame()
    model_arr['Model Name'] = n
    model_arr['Progenitor mass'] = mass
    model_arr['Explosion Energy'] = exp_ener
    model_arr['Mass Loss'] = mass_loss
    model_arr['CSM edge'] = csm_edge
    model_arr['Beta'] = beta
print(model_arr)



group = model_arr.groupby(['Explosion Energy','Mass Loss','CSM edge','Beta']).get_group(('1e51','1e-2','1e15','1'))
for i in range(0,len(group)):
    mod = pd.read_table('mod'+'_'+ group.iloc[i,0]+'.sdss2.txt', names = ['epoch','u','g','r','i','z','kepler'], sep='\s+')
    a = mod.filter(['epoch','r'])
    b = a.groupby('epoch').mean()
    #print(df2.shape)
    z = mod['epoch'].drop_duplicates().values
    p = b.values.squeeze()
    arr = pd.DataFrame()
    arr['epoch'] = z
    arr['Magnitude'] = b
    #print(pred_arr)
    fil_max= arr.loc[arr['Magnitude'].idxmin()]
    #print(fil_max)
    #print(fil_max[0])
    #print(fil_max[1])
    pdays=[]
    for k in range(0,len(arr)):
        d = arr.iloc[k]['epoch'] - fil_max[0]
        #print(d)
        pdays.append(d)
    plt.plot(pdays,p,'o')
plt.gca().invert_yaxis()
plt.show()

obs = np.array([-2,-4,-9.0,-6])
pred =np.array( [-3,-4,-8,-8])
a = ((obs-pred)**2)
chi = chisquare(obs.astype(np.float64), pred.astype(np.float64))
print(chi)
a =4
print(a**2)


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
    a = splitcontent[i].split(',')
    a.remove('Msun')
    a.remove('erg')
    if( len(a) >5):
        a.remove('cm')
        a.remove('Msun/yr')
    else:
        a.append(0)
        a[-2] = 0
        a[-3] = 0
    n.append(a[0])
    mass.append(a[1])
    exp_ener.append(a[2])
    mass_loss.append(a[3])
    csm_edge.append(a[4])
    beta.append(a[5])
    model_arr = pd.DataFrame()
    model_arr['Model Name'] = n
    model_arr['Progenitor mass'] = mass
    model_arr['Explosion Energy'] = exp_ener
    model_arr['Mass Loss'] = mass_loss
    model_arr['CSM edge'] = csm_edge
    model_arr['Beta'] = beta
print(model_arr)
group = model_arr.groupby(['Explosion Energy','Mass Loss','CSM edge','Beta']).get_group(('1e51','1e-2','1e15','5'))
chi_m_arr = pd.DataFrame()
chi_m_arr['Model name'] = group['Model Name']
chi_m_arr['P Mass'] = group['Progenitor mass']
chi_m_arr = chi_m_arr.reset_index(drop=True)
print(chi_m_arr)
sl_list = []
