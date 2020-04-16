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
test = pd.read_csv('/Users/bhagyasubrayan/Desktop/Explosion Parameters/Plasticc/plasticc_test_lightcurves_01.csv')
f_0 = 27.5
dmod = 40.977
true_peak = 60499.461
#t_g = test.groupby(['object_id']).get_group((13))
#lc = t_g.filter(['mjd','flux'])
#print(type(lc))
chi_metric = []
p_metric= []
for i in range(0,1):
    t_g = test.groupby(['object_id', 'passband']).get_group((13,i))
    #print(t_g)
    am_m_i =[]
    ndays_i = []
    new_i = t_g.filter(['mjd','flux','flux_err'])
    new_i = new_i.reset_index(drop=True)
    new_i = new_i.drop(new_i[new_i.mjd < (true_peak - 150) ].index)
    plt.errorbar(new_i['mjd'],new_i['flux'],yerr=new_i['flux_err'],fmt='o')
    plt.show()
    #print(t_g)
    am_m_i =[]
    ndays = []
    #print(new_i)
    for j in range(len(new_i)):
        if (new_i.iloc[j]['flux'] < 0 ):
            new_i.iloc[j]['flux'] = 0.01
        a = -2.5*m.log(new_i.iloc[j]['flux']/ f_0)
        abs = a - 5*m.log(dmod/10)
        am_m_i.append(abs)
    new_i['Abs_mag'] = am_m_i
    new_i.reset_index(drop= True)
    #print(len(new_i))
    tag = new_i.loc[new_i['Abs_mag'].idxmin()]
    #print(tag[0])
    for k in range(0,len(new_i)):
        d = new_i.iloc[k]['mjd'] - tag[0]
        #print(d)
        ndays.append(d)
    #print(len(ndays))
    new_i['t_max_afterdays'] = ndays
    print(new_i)
    x = new_i['t_max_afterdays']
    #print(x)
    y_obs_chi = new_i['Abs_mag'].values
    #print(y_obs_chi)
    err = new_i['flux_err']
    plt.errorbar(x,y_obs_chi,yerr = err, fmt = 'o')
    plt.gca().invert_yaxis()
    plt.xlabel('days')
    plt.ylabel('Absolute Magnitude in r ')
plt.show()
list = os.listdir("/Users/bhagyasubrayan/Desktop/Plastic/public_data/multicolorlightcurves/")
list.remove('README.txt')
list.sort()
list = np.array(list)
#print(list.shape)
list_num = np.linspace(1, len(list),len(list))
print(list_num.shape)
#print(list)
#print(len(list))
for i in range(0,len(list)):
    #print(list[i])
    pred_days=[]
    lines = open('/Users/bhagyasubrayan/Desktop/Plastic/public_data/multicolorlightcurves/'+ list[i]).readlines()
    a = open('mod'+'_'+ list[i]+ '.txt', 'w').writelines(lines[2:])
    df = pd.read_table('mod'+'_'+ list[i]+'.txt', names = ['epoch','u','g','r','i','z','kepler'], sep='\s+')
    #print(df)
    df1 = df.filter(['epoch','r'])
    df2 = df1.groupby('epoch').mean()
    #print(x1)
    x1 = df['epoch'].drop_duplicates().values
    y1 = df2.values.squeeze()
    pred_arr = pd.DataFrame()
    pred_arr['epoch'] = x1
    pred_arr['Magnitude'] = y1
    #print(pred_arr)
    fil_max= pred_arr.loc[pred_arr['Magnitude'].idxmin()]
    #print(fil_max)
    #print(fil_max[0])
    #print(fil_max[1])
    for k in range(0,len(pred_arr)):
        d = pred_arr.iloc[k]['epoch'] - fil_max[0]
        #print(d)
        pred_days.append(d)
    #print(pred_days)
    xnew = np.linspace(pred_days[0],pred_days[-1],500)
    f1 = interpolate.interp1d(pred_days, y1, kind='cubic')
    #plt.plot(pred_days,y1,'o',xnew,f1(xnew),'-')
    #x_up = np.linspace(0,int(x1[-1]), 50).squeeze()
    #plt.title(list[i])
    #plt.xlim(0,150)
    y_pred_chi = f1(x)
    #print(y_pred_chi)
    plt.plot(x, y_pred_chi ,'*')
    #plt.gca().legend(('data','interpolate','from_observed'))
    for k in range(0,len(y_obs_chi)):
        y_obs_chi[k] = m.fabs(y_obs_chi[k])
        y_pred_chi[k] = m.fabs(y_pred_chi[k])
    #print((y_obs_chi[10]))
    #print((y_pred_chi[10]))
    chi = chisquare(y_obs_chi,y_pred_chi,ddof = 5)
    #print(chi)
    c = chi[0]
    p = chi[1]
    p_metric.append(p)
    chi_metric.append(c)
plt.gca().invert_yaxis()
plt.show()
chi_metric = np.array(chi_metric)
print(len(chi_metric))
chi_m_arr = pd.DataFrame()
chi_m_arr['Model Number'] = list_num
chi_m_arr['Chi_squared'] = chi_metric
chi_m_arr['p-value'] = p_metric
print(chi_metric.shape)
print(list_num.shape)
print(chi_m_arr)
tag2 = chi_m_arr.loc[chi_m_arr['Chi_squared'].idxmin()]
tag3 = tag2[0] - 1
print(tag3)
print('The best model would be:'+list[int(tag3)])
print(str(list[int(tag3)]))
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
model_arr['Reduced chi_squared'] = chi_metric
model_arr['p-value'] = p_metric
print(model_arr)
#print(num)
#for k in range(0, len(list)):
 #  list[k] = list[k].replace('yoonb','')
  # list[k] = list[k].replace('.sdss2','')
plt.plot(list_num,chi_metric,'*')
plt.xlabel('Model Number')
plt.ylabel(u'$\u03C7^2$')
plt.title(u'$\u03C7^2$'+'\t'+'Test')
plt.savefig("out.png")
plt.show()
