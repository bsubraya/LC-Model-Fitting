import numpy as np
import csv
import matplotlib.pyplot as plt
import json
import math as m
import os
import pandas as pd
import scipy.interpolate as sp
from scipy import interpolate
from astropy.time import Time
from datetime import datetime
from astropy.cosmology import FlatLambdaCDM
from astropy.coordinates import SkyCoord
from astropy import units as u
import sfdmap
def event_data(filename):
    with open(filename) as f:
        file_data = json.load(f)
        data.append(file_data)
def process(filename):
    file = pd.read_csv(filename)
    df_g = file.groupby(['band']).get_group('g')
    df_r = file.groupby(['band']).get_group('r')
    return df_g.reset_index(drop =True), df_r.reset_index(drop=True)
def metadata(filename):
    for i in range(0,len(data)):
        for j in range(0,len(data[i])):
            if (data[i][j]['lasairname'] == filename):
                #print(data[i][j])
                max_date = data[i][j]['Max date']
                max_date = datetime.strptime(max_date, '%Y/%m/%d').strftime('%Y-%m-%d')
                #print(max_date)
                peak_mjd = Time(max_date, format='isot', scale='utc')
                redshift = data[i][j]['z']
                ra,dec = data[i][j]['R.A.'],data[i][j]['Declination']
    if(redshift == "n/a"):
        redshift = input('Redshift not found, please enter manually for event:'+ event_list[i])
    #print('z =',redshift)
                #print(peak_mjd.mjd)
    return float(redshift),ra,dec
def distance_mod(z):
    cosmo = FlatLambdaCDM(H0=70, Om0=0.3, Tcmb0=2.725)
    dmod = cosmo.distmod(z)
    return dmod.value
def predicted_mag_at_observed_epochs(df,days,df_ztf):
    b_arr = df.drop_duplicates()
    #print(b_arr)
    arr = df.groupby('epoch').mean()
    #print(arr)
    epo = b_arr['epoch'].drop_duplicates()
    #print(type(list(epo)[-1]))
    mag = arr.values.squeeze()
    n_arr = pd.DataFrame()
    n_arr['Epoch'] = epo
    n_arr['Magnitude'] = mag
    n_arr = n_arr.reset_index(drop = 'True')
    #print(n_arr)
    #print(list(epo)[-1])
    days = [x for x in days if x <= list(epo)[-1]]
    #print(days)
    mod_ztf = df_ztf.query("epoch in @days")
    #print(mod_ztf)
    f1 = interpolate.interp1d(epo, mag, kind='cubic')
    pred_mag = f1(days)
    #print(pred_mag)
    mod_ztf = mod_ztf.assign(Pred_mag =pred_mag)
    return mod_ztf
def model_param(model_name):
    best = model_arr.groupby(['Model Name']).get_group(model_name)
    return best
def final_obs_df(df,t,band):
    #Correcting for extinction
    #print(df)
    df['epoch'] = df.MJD - t
    if band == 'g':
        df['Ab_obs_mag'] = df.magpsf - dmod -Ag
        #print(True)
    else:
        df['Ab_obs_mag'] = df.magpsf - dmod - Ar
    df['mag_error'] = df.sigmapsf + s
    return df
def fit_all_models(df,band):
    chi = []
    for k in range(0,len(all_group)):
        m_n = all_group[k]
        #print(m_n)
        mod = pd.read_table(path + all_group[k] +'.sdss2', skiprows = 2,names = ['epoch','u','g','r','i','z','kepler'], sep='\s+')
        #print(mod)
        df_pred = mod[['epoch',band]]
        #To fix the problem of interpolation remove elements greater than the interpolation
        #range and modifying the corresponding ztf_obs_dataframe
        df_final = predicted_mag_at_observed_epochs(df_pred,list(df.epoch),df)
        chi2 = np.sum((((df_final.Ab_obs_mag - df_final.Pred_mag)/df_final.mag_error)**2))/(len(df_final)-5)
        #chi2_r = np.sum((((df_r_ztf.Ab_obs_mag - df_r_ztf.Pred_mag)/df_r_ztf.sigmapsf)**2))/(len(df_r_ztf)-5)
        #chi_g.append(chi2_g)
        #chi_r.append(chi2_r)
        #print(chi2)
        if(abs(chi2 - 1) <= 5.0):
            plt.plot(df_pred['epoch'],df_pred[band],'--',label = m_n)
        chi.append(chi2)#Value of reduced chi^^2 for each model
    return chi
event = ['ZTF19abcneik']
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
pd.options.display.max_rows = 999
print((model_arr))
#Defining path for multicolor wavelength data from Takashi's Models
path = "/Users/bhagyasubrayan/Desktop/Explosion Parameters/Takashi's models/multicolorlightcurves/"
#Collecting all the metadata of events throughout the three years
data=[]
event_data('./ztf-v2/2018_sne_table.json')
event_data('./ztf-v2/2019_sne_table.json')
event_data('./ztf-v2/2020_sne_table.json')
all_group = list(model_arr['Model Name'])
for i in range(0,len(event)):
    print(event[i])
    z,ra,dec = metadata(event[i])
    print('Redshift : ',z)
    print(ra,dec)
    dmod = distance_mod(z)
    df_g_ztf, df_r_ztf = process(event[i]+'.csv')
    dustmap = sfdmap.SFDMap("/Users/bhagyasubrayan/Desktop/sncosmo/sfddata-master")
    c = SkyCoord(ra+dec, unit=(u.hourangle, u.deg))
    ebv = dustmap.ebv(c)
    print('Ebv:',ebv)
    Ar = 2.751*ebv
    Ag = 3.793*ebv
    print('g extinction :',Ag)
    print('r extinction :',Ar)
    #Correcting dmod unceratinity while drawing from a normal distribution with mu and sigmaps
    mu, sigma = 0.005, 0.001 # mean and standard deviation
    s = np.random.normal(mu, sigma, 1)
    print('dmod uncertainity:',s)
    #plt.figure(figsize=(10,8))
    print(df_g_ztf)
    df_g_obs =  final_obs_df(df_g_ztf,df_g_ztf.MJD[0]-4.5,band = 'g')
    df_r_obs =  final_obs_df(df_r_ztf,df_r_ztf.MJD[0]-4.5,band = 'r')
    print(df_g_obs)
    fig = plt.figure(figsize=(15,6))
    #plt.title(event[i],fontsize='15')
    a = fig.add_subplot(1,2,1)
    chi_value_g = fit_all_models(df_g_obs,band= 'g')
    plt.errorbar(df_g_obs.epoch,df_g_obs.Ab_obs_mag,yerr = df_g_obs.mag_error,fmt = 'go',label = 'Observed g band',ms ='10',mec='black')
    plt.legend()
    plt.xlabel('Time since explosion(days)',fontsize='15')
    plt.title(event[i],fontsize='15')
    plt.ylabel('Absolute magnitude',fontsize='15')
    plt.rc('xtick', labelsize= 15)
    plt.rc('ytick', labelsize= 15)
    plt.gca().invert_yaxis()
    plt.ylim(-14,-18)
    plt.xlim(0,150)
    a = fig.add_subplot(1,2,2)
    chi_value_r = fit_all_models(df_r_obs,band= 'r')
    plt.errorbar(df_r_obs.epoch,df_r_obs.Ab_obs_mag,yerr = df_r_obs.mag_error,fmt = 'ro',label = 'Observed r band',ms ='10',mec='black')
    plt.legend()
    plt.xlabel('Time since explosion(days)',fontsize='15')
    plt.title(event[i],fontsize='15')
    plt.ylabel('Absolute magnitude',fontsize='15')
    plt.rc('xtick', labelsize= 15)
    plt.rc('ytick', labelsize= 15)
    plt.gca().invert_yaxis()
    plt.ylim(-14,-18)
    plt.xlim(0,150)
    plt.savefig('/Users/bhagyasubrayan/Desktop/degeneracy.jpg')
    plt.show()
