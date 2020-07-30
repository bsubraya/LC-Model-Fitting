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
def process(df,band,source):
    if source == '0':
        #print(True)
        new_df = df.groupby('band').get_group(band)
    else:
        new_df = df.groupby(['band','source']).get_group((band,source))
    new_df = new_df.reset_index(drop =True)
    df_req = new_df[['time','magnitude','e_magnitude','band']]
    df_req.columns = ['MJD','magpsf','sigmapsf','band']
    return df_req
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
def plot_model(df,name,band,t,color):#Add upper df if upper limits
    plt.errorbar(df['epoch'],df['Ab_obs_mag'],yerr = df['mag_error'],fmt = color+'o',label = 'Observed '+band+' with texp : '+ str(t))
    mod = pd.read_table(path + name +'.sdss2', skiprows = 2,names = ['epoch','u','g','r','i','z','kepler'], sep='\s+')
    df_pred = mod[['epoch',band]]
    plt.plot(df_pred['epoch'],df_pred[band],color+'--',label = 'Best_fit : '+name)
    #plt.plot(upper['epoch'],upper['Ab_obs_mag'],band+'*',label = 'Upper Limits')
    plt.legend()
    plt.xlim(-40,200)
    plt.ylim(-20, -10)
def final_obs_df(eventname,df,t,band):
    #Correcting for extinction
    #print(df)
    df['epoch'] = df.MJD - t
    if band == 'g':
        df['Ab_obs_mag'] = df.magpsf - dmod -Ag
        #print(True)
    elif band == 'r':
        df['Ab_obs_mag'] = df.magpsf - dmod - Ar
    elif band =='i':
        df['Ab_obs_mag'] = df.magpsf - dmod - Ai
    elif band == 'z':
        df['Ab_obs_mag'] = df.magpsf - dmod - Az
    df['mag_error'] = df.sigmapsf + s
    return df
def metadata(eventname):
    with open('/Users/bhagyasubrayan/Downloads/'+eventname+'.json') as f:
        file_data = json.load(f)
    dec = file_data[eventname]['dec'][0]['value']
    ra = file_data[eventname]['ra'][0]['value']
    z = file_data[eventname]['redshift'][0]['value']
    dmod = distance_mod(float(z))
    print('Event:',eventname)
    print('Distance Modulus:',dmod)
    print('Redshift:', z )
    print('RA,Dec :',ra + ','+dec)
    Ag,Ar,Ai,Az,ebv = extinction(ra,dec)
    print('ebv:',ebv)
    mu, sigma = 0.005, 0.001 # mean and standard deviation
    s = np.random.normal(mu, sigma, 1)
    print('dmod uncertainity:',s)
    return z, ra, dec, Ag,Ar,Ai, Az,ebv,s,dmod
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
def best_texp(eventname,df,band,min,max):
    texp = []
    fit_best_model = []
    fit_best_chi = []
    fit_time = pd.DataFrame(columns = ['texplosion','Best_fit_model','Reduced_chi'])
    for k in np.arange(min,max,0.25):
        print('Factor of texplosion:',k)
        chi_m_arr = pd.DataFrame()
        chi_m_arr['Model name'] = all_group
        chi_m_arr = chi_m_arr.reset_index(drop=True)
        texp.append(k)
        texp_date = df['MJD'][0] - (k)
        #print(texp_date)
        df_obs =  final_obs_df(eventname,df,texp_date,band)
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
    #print(fit_time)
    best_fit = fit_time.loc[fit_time['Difference'].idxmin()]
    #else:
        #best_fit = fit_time.loc[fit_time['Reduced_chi'].idxmin()]
    #print(best_fit)
    return fit_time, best_fit[0] , best_fit[1],best_fit[2]
def analysis(eventname,df,band,color,min,max):#Add upper if upper limits
    original_obs = subset_df(df)
    time_array, texplosion, fit_model, fit_value = best_texp(eventname,df,band,min,max)
    mod_df = final_obs_df(eventname,original_obs,original_obs.MJD[0] - texplosion,band = band)
    #mod_upper = final_obs_df(upper,original_obs.MJD[0]-texplosion)
    #print(mod_df)
    plot_model(mod_df, name = fit_model,band = band, t = texplosion,color=color)
    best_fit_param = model_param(fit_model)
    print(best_fit_param)
    print('Reduced_chi_value for '+ band+ ':', fit_value )
    return texplosion, fit_model
def extinction(ra,dec):
    dustmap = sfdmap.SFDMap("/Users/bhagyasubrayan/Desktop/sncosmo/sfddata-master")
    c = SkyCoord(ra+dec, unit=(u.hourangle, u.deg))
    ebv = dustmap.ebv(c)
    Ar = 2.751*ebv   #1998 prescription
    Ag = 3.793*ebv
    Ai = 2.086*ebv
    Az = 1.479*ebv
    return Ag,Ar,Ai,Az,ebv
events = ['SN2012aw','SN2014cx','SN2013ab','SN2013fs']
def group(df,eventname,bands,source):
    band_df = []
    #print(source)
    for item in bands:
        d = process(df,band = item,source = source)
        print(d)
        band_df.append(d)
    if eventname == 'SN2014cx' or eventname == 'SN2013ab':
        return band_df[0],band_df[1],band_df[2]
    elif eventname == 'SN2013fs' or eventname == 'SN2012aw' :
        return band_df[0],band_df[1],band_df[2],band_df[3]

def event(eventname):
    sn_data = pd.read_csv('/Users/bhagyasubrayan/'+ eventname+'.csv')
    print(sn_data)
    plt.figure(figsize=(10,8))
    if eventname == 'SN2013fs':
        #sn_data_1 = pd.read_csv('/Users/bhagyasubrayan/'+ eventname+'.csv')
        #print(sn_data_1)
        bands = ['g','r','i','z']
        g, r, i ,zband = group(sn_data,eventname,bands,source ='2016MNRAS.459.3939V')
        texplo_g, model_g = analysis(eventname,g,band = 'g',color = 'g',min= 0,max = 5)#upper = g_upper
        texplo_r,model_r = analysis(eventname,r,band = 'r',color = 'r', min = texplo_g - 1 ,max = texplo_g + 3)
        texplo_i,model_i = analysis(eventname,i,band = 'i',color = 'y', min = texplo_g - 3 ,max = texplo_g + 3)
        texplo_z,model_z = analysis(eventname,zband,band = 'z',color = 'b', min = texplo_g - 3 ,max = texplo_g + 3)
    elif eventname == 'SN2014cx':
        #sn_data_2 = pd.read_csv('/Users/bhagyasubrayan/'+ eventname+'.csv')
        #print(sn_data_2)
        bands = ['g','r','i']
        g,r,i = group(sn_data,eventname,bands,source ='2016MNRAS.459.3939V')
        texplo_g, model_g = analysis(eventname,g,band = 'g',color = 'g',min= 0,max = 5)#upper = g_upper
        texplo_r,model_r = analysis(eventname,r,band = 'r',color = 'r', min = texplo_g - 1 ,max = texplo_g + 1)
        texplo_i,model_i = analysis(eventname,i,band = 'i',color = 'y', min = texplo_g - 1 ,max = texplo_g + 1)
    elif eventname =='SN2013ab':
        #sn_data_3 = pd.read_csv('/Users/bhagyasubrayan/'+ eventname+'.csv')
        #print(sn_data_3)
        bands = ['g','r','i']
        g,r,i = group(sn_data,eventname,bands,source = '2015MNRAS.450.2373B')
        texplo_r,model_r = analysis(eventname,r,band = 'r',color = 'r', min = 0 ,max = 5)
        texplo_g, model_g = analysis(eventname,g,band = 'g',color = 'g',min= texplo_r - 1 ,max = texplo_r + 1 )#upper = g_upper
        texplo_i,model_i = analysis(eventname,i,band = 'i',color = 'y', min = texplo_r - 1 ,max = texplo_r + 1)
    elif eventname == 'SN2012aw':
        #sn_data_4 = pd.read_csv('/Users/bhagyasubrayan/'+ eventname+'.csv')
        #print(sn_data_4)
        bands = ['g','r','i','z']
        g, r, i ,zband = group(sn_data,eventname,bands, source = '0')
        texplo_g, model_g = analysis(eventname,g,band = 'g',color = 'g',min= 0,max = 5)#upper = g_upper
        texplo_r,model_r = analysis(eventname,r,band = 'r',color = 'r', min = texplo_g - 1 ,max = texplo_g + 3)
        texplo_i,model_i = analysis(eventname,i,band = 'i',color = 'y', min = texplo_g - 3 ,max = texplo_g + 3)
        texplo_z,model_z = analysis(eventname,zband,band = 'z',color = 'b', min = texplo_g - 3 ,max = texplo_g + 3)
    #upper = r_upper,
    plt.title(eventname+ ' ')
    plt.xlabel('Days from texplosion')
    plt.gca().invert_yaxis()
    plt.show()
for i in range(0,len(events)):
    z, ra, dec, Ag,Ar,Ai, Az,ebv,s,dmod = metadata(events[i])
    event(eventname = events[i])
