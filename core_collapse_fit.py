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
pd.options.mode.chained_assignment = 'raise'
#pd.options.display.max_rows = 999
#Ordering Takashi's Models
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
#print((model_arr))
#Defining path for multicolor wavelength data from Takashi's Models
path = "/Users/bhagyasubrayan/Desktop/Explosion Parameters/Takashi's models/multicolorlightcurves/"
#Collecting all the metadata of events throughout the three years
data=[]
def event_data(filename):
    with open(filename) as f:
        file_data = json.load(f)
        data.append(file_data)
event_data('./ztf-v2/2018_sne_table.json')
event_data('./ztf-v2/2019_sne_table.json')
event_data('./ztf-v2/2020_sne_table.json')
#Trying to give weights to bands, and fitting the other band based on this band
#First trying for g_events

event_g = ['ZTF19abqyouo','ZTF19abqyoxj','ZTF19aceckht','ZTF19acxowrr','ZTF18acyybvg','ZTF19abcneik','ZTF19aavkptg','ZTF20aanadof']
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
all_group = list(model_arr['Model Name'])
def model_param(model_name):
    best = model_arr.groupby(['Model Name']).get_group(model_name)
    return best
def plot_model(df,name,band,t):#Add upper df if upper limits
    plt.errorbar(df['epoch'],df['Ab_obs_mag'],yerr = df['mag_error'],fmt = band+'o',label = 'Observed '+band+' with texp : '+ str(t))
    mod = pd.read_table(path + name +'.sdss2', skiprows = 2,names = ['epoch','u','g','r','i','z','kepler'], sep='\s+')
    df_pred = mod[['epoch',band]]
    plt.plot(df_pred['epoch'],df_pred[band],band+'--',label = 'Best_fit : '+name)
    #plt.plot(upper['epoch'],upper['Ab_obs_mag'],band+'*',label = 'Upper Limits')
    plt.legend()
    plt.xlim(-40,200)
    plt.ylim(-20, -10)

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
        chi.append(chi2)#Value of reduced chi^^2 for each model
    return chi
def subset_df(df):
    new_df= df[['MJD','magpsf','sigmapsf','band']]
    return new_df
def best_texp(df,band,min,max):
    texp = []
    fit_best_model = []
    fit_best_chi = []
    fit_time = pd.DataFrame(columns = ['texplosion','Best_fit_model','Reduced_chi'])
    for k in range(min,max):
        print('Factor of texplosion:',k)
        chi_m_arr = pd.DataFrame()
        chi_m_arr['Model name'] = all_group
        chi_m_arr = chi_m_arr.reset_index(drop=True)
        texp.append(k)
        texp_date = df['MJD'][0] - (k)
        #print(texp_date)
        df_obs =  final_obs_df(df,texp_date,band)
        #print(df_obs)
        chi_value = fit_all_models(df_obs,band)
        #print(chi_value)
        chi_m_arr['Reduced_chi'] = chi_value
        #print('Reduced Chi-Squared array in '+band+' for texplosion '+ str(k+1)+ ' : ')
        chi_m_arr['Difference'] = list(abs(chi_m_arr.Reduced_chi - 1))
        #print(chi_m_arr)
        #if method == 'difference':
        best_model = chi_m_arr.loc[chi_m_arr['Difference'].idxmin()]
        #else:
            #best_model = chi_m_arr.loc[chi_m_arr['Reduced_chi'].idxmin()]
        fit_best_model.append(best_model[0])
        fit_best_chi.append(best_model[1])
    fit_time['texplosion'] = texp
    fit_time['Best_fit_model'] = fit_best_model
    fit_time['Reduced_chi'] = fit_best_chi
    fit_time['Difference'] = list(abs(fit_time.Reduced_chi - 1))
    #if method == 'difference':
    best_fit = fit_time.loc[fit_time['Difference'].idxmin()]
    #else:
        #best_fit = fit_time.loc[fit_time['Reduced_chi'].idxmin()]
    #print(best_fit)
    return fit_time, best_fit[0] , best_fit[1],best_fit[2]

def analysis(df,band,min,max):#Add upper if upper limits
    original_obs = subset_df(df)
    time_array, texplosion, fit_model, fit_value = best_texp(df,band,min,max)
    mod_df = final_obs_df(original_obs,original_obs.MJD[0] - texplosion,band = band)
    #mod_upper = final_obs_df(upper,original_obs.MJD[0]-texplosion)
    #print(mod_df)
    plot_model(mod_df, name = fit_model,band = band, t = texplosion)
    best_fit_param = model_param(fit_model)
    print(best_fit_param)
    print('Reduced_chi_value for '+ band+ ':', fit_value )
    return texplosion, fit_model
#With weight for g_band_only
#method = input('Choose the method (difference/minimum):')
for i in range(0,len(event_g)):
    print(event_g[i])
    z,ra,dec = metadata(event_g[i])
    print('Redshift : ',z)
    print(ra,dec)
    dmod = distance_mod(z)
    df_g_ztf, df_r_ztf = process(event_g[i]+'.csv')
    #print(df_g_ztf)
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
    #g_upper, r_upper = process(event_g[i]+'_upper.csv')
    plt.figure(figsize=(10,8))
    texplo_g, model_g = analysis(df_g_ztf,band = 'g',min= 1,max = 15)#upper = g_upper
    texplo_r,model_r = analysis(df_r_ztf,band = 'r', min = texplo_g - 3 ,max = texplo_g + 5)#upper = r_upper,
    plt.title(event_g[i]+ ' ')
    plt.ylabel('Absolute Magnitude')
    plt.xlabel('Days from texplosion')
    plt.gca().invert_yaxis()
    plt.show()
