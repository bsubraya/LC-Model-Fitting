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
def write_file(df,filename):
    csv = df.to_csv(name,header=0,sep =' ',index = 0)
filters = ['u','g','r','i','z','y']
#t_g = test.groupby(['object_id']).get_group((13))
#lc = t_g.filter(['mjd','flux'])
#print(type(lc))
chi_metric = []
p_metric= []
for i in range(0,6):
    #t_g = test.groupby(['object_id', 'passband']).get_group((13,i))#for any band
    t_g = test.groupby(['object_id', 'passband']).get_group((13,i))
    #print(t_g)
    am_m_i =[]
    mag_err_i = []
    ndays_i = []
    app_mag = []
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
        app_mag.append(a1)
    new_i['Apparent'] = app_mag
    new_i['Abs_mag'] = am_m_i
    new_i['Mag_error'] = mag_err_i
    #print(new_i)
    new_i = new_i.drop(new_i[new_i.Mag_error > 10].index)
    new_i = new_i.reset_index(drop= True)
    #print((new_i))
    tag = new_i.loc[new_i['Abs_mag'].idxmin()]
    #print(tag[0])
    for k in range(0,len(new_i)):
        d = new_i.iloc[k]['mjd'] - true_peak
        #print(d)
        ndays_i.append(d)
    #print(len(ndays))
    new_i['t_max_afterdays'] = ndays_i
    print(new_i)
    file_for_sb = new_i[['mjd','Abs_mag','Mag_error']]
    print(file_for_sb)
    name = '/Users/bhagyasubrayan/Desktop/superbol-master/Plasticc_abs_'+filters[i]+'.txt'
    write_file(file_for_sb, name)
 
