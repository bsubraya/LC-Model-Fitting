import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import json
import math as m
import os

test = pd.read_csv('/bsubraya/Desktop/PlasTiCc/plasticc_test_lightcurves_01.csv')
f_0 = 27.5 
am_m=[]
m_m = []
t_g = test.groupby(['object_id']).get_group((13))
lc = t_g.filter(['mjd','flux'])
for i in range(0, 6):
    t_g = test.groupby(['object_id', 'passband',]).get_group((13,i))
    #print(t_g)
    new_i = t_g.filter(['mjd','flux'])
    #print(new_i)
    for j in range(len(new_i)):
        b = lc.iloc[j,1]
        c = lc.iloc[j,0]
        if (b/f_0 >= 0):
            a = -2.5*m.log(b/ f_0)
            am_m.append(a)
            m_m.append(c)
            #print(am_m)
            #print(m_m)
    #plt.scatter(m_m,am_m)
    #plt.show()
