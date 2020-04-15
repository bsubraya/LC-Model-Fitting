import numpy as np
import csv
import matplotlib.pyplot as plt
import json
import math as m
import os
import pandas as pd
import scipy.interpolate as sp
from scipy import interpolate
list = os.listdir("/Users/bhagyasubrayan/Desktop/Plastic/public_data/multicolorlightcurves/")
list.remove('README.txt')
list.sort()
print(len(list))
for i in range(0,1):
    print(list[i])
    lines = open('/Users/bhagyasubrayan/Desktop/Plastic/public_data/multicolorlightcurves/'+ list[i]).readlines()
    a = open('mod'+'_'+ list[i]+ '.txt', 'w').writelines(lines[2:])
    df = pd.read_table('mod'+'_'+ list[i]+'.txt', names = ['epoch','u','g','r','i','z','kepler'], sep='\s+')
    df1 = df.filter(['epoch','r'])
    df2 = df1.groupby('epoch').mean()
    x = df['epoch'].drop_duplicates().values
    print(int(x[-1]))
    y = df2.values.squeeze()
    xnew = np.linspace(x[0],x[-1],500)
    f1 = interpolate.interp1d(x, y, kind='cubic')
    plt.plot(x,y,'o',xnew,f1(xnew),'*')
    x_up = np.linspace(0,int(x[-1]), 50).squeeze()
    plt.plot(x_up, f1(x_up),'-')
    plt.gca().invert_yaxis()
    plt.show()
