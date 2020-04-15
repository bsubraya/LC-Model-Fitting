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
dmod = 40.977
true_peak = 60499.461
#t_g = test.groupby(['object_id']).get_group((13))
#lc = t_g.filter(['mjd','flux'])
#print(type(lc))
t_g = test.groupby(['object_id', 'passband',]).get_group((13,2))
am_m_i =[]
m_m_i = []
err_i = []
#print(t_g)
new_i = t_g.filter(['mjd','flux','flux_err'])
new_i = new_i.reset_index(drop=True)
#print(new_i.iloc[10]['mjd'])
#print(new_i)
for j in range(len(new_i)):
    if (new_i.iloc[j]['flux'] < 0 ):
        new_i.iloc[j]['flux'] = 0.01
    if (new_i.iloc[j]['mjd'] > (true_peak - 100)):
      a = -2.5*m.log(new_i.iloc[j]['flux']/ f_0)
      abs = a - 5*m.log(dmod/10)
      am_m_i.append(abs)
      m_m_i.append(new_i.iloc[j]['mjd'])
      err_i.append(new_i.iloc[j]['flux_err'])
    #print(m_m_i)
    #print(len(m_m_i))
print(len(m_m_i))
xnew = np.linspace(m_m_i[0], m_m_i[-1],15)
f1 = interpolate.interp1d(m_m_i, am_m_i, kind='cubic')
print(f1[60550])

plt.errorbar(m_m_i,am_m_i,yerr = err_i, fmt = 'o')
plt.gca().invert_yaxis()
plt.plot(m_m_i, am_m_i, 'o',xnew,f1(xnew),'+-')
plt.xlabel('mjd')
plt.ylabel('Absolute Magnitude')
plt.show()
