import numpy as np
import csv
import matplotlib.pyplot as plt
import json
import math as m
import os
import pandas as pd
import scipy.interpolate as sp
from scipy import interpolate
path = "/Users/bhagyasubrayan/Desktop/Explosion Parameters/Takashi's models/multicolorlightcurves/"
list = os.listdir(path)
list.remove('README.txt')
list.sort()
print(len(list))
for i in range(0,3):
    print(type(list[i]))
    df = pd.read_csv(path+list[i], skiprows = 2,names = ['epoch','u','g','r','i','z','kepler'], sep='\s+')
    df1 = df.filter(['epoch','r'])
    df2 = df1.groupby('epoch').mean()
    x = df['epoch'].drop_duplicates().values
    y = df2.values.squeeze()
    plt.scatter(x,y)
    plt.gca().invert_yaxis()
    plt.show()
    xnew = np.linspace(x[0],x[-1],500)
    f1 = interpolate.interp1d(x, y, kind='cubic')
    plt.plot(x,y,'o',xnew,f1(xnew),'*')
    x_up = np.linspace(0,int(x[-1]), 50).squeeze()
    plt.plot(x_up, f1(x_up),'-')
    plt.gca().invert_yaxis()
    plt.show()
fil_max= df2.loc[df2['r'].idxmin()]
fil_max
