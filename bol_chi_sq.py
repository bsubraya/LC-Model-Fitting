import numpy as np
import csv
import matplotlib.pyplot as plt
import json
import math as m
import os
import pandas as pd
import scipy.interpolate as sp
from scipy import interpolate
import os
from astropy.io import ascii
path = '/Users/bhagyasubrayan/Downloads/Melina Models/'
legend = []
files = os.listdir(path)
for i in range(0,10):
    #print(filename)
    legend.append(files[i])
    with open(path+files[i],'rb') as f:
        data = pd.read_table(f,names = ['epoch','flux','flux_err','unkn'], sep='\s+')
        #print(data)
        #plt.errorbar(data['epoch'],data['flux'],yerr = data['flux_err'],fmt= 'ro')
        #plt.figure(figsize = (8,8))
        plt.scatter(data['epoch'],data['flux'])
        plt.xlim(-10,100)
        plt.ylim(40,43)
        #plt.title(filename)
plt.legend(legend)
plt.xlabel('Days')
plt.ylabel('Flux in log_10')
plt.gcf().set_size_inches((10, 10))
#plt.savefig("/Users/bhagyasubrayan/Desktop/melina.png")
plt.show()
data
