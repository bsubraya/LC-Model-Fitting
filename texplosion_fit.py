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
#all_group = list(model_arr['Model Name'])
#len(all_group)
#chi_m_arr = pd.DataFrame()
#chi_m_arr['Model name'] = all_group
path = "/Users/bhagyasubrayan/Desktop/Explosion Parameters/Takashi's models/multicolorlightcurves/"
#Managing the observed ZTF data
#event_list = {'ZTF19abcneik','ZTF19aavkptg','ZTF18abrlljc','ZTF18abckutn'}
event_list = ['ZTF19abcneik','ZTF19aavkptg','ZTF18abrlljc','ZTF18abckutn','ZTF19acxowrr','ZTF19aceckht','ZTF20aanadof','ZTF19abqyoxj','ZTF18acyybvg']
#event_list = ['ZTF20aanadof','ZTF19abqyoxj','ZTF18acyybvg' ]
#event_list = ['ZTF20aanadof','ZTF19aceckht']
#'ZTF19abcneik','ZTF19aavkptg','ZTF18abrlljc'
data=[]
with open('./ztf-v2/2018_sne_table.json') as f:
    data1 = json.load(f)
    data.append(data1)
with open('./ztf-v2/2019_sne_table.json') as f:
    data2 = json.load(f)
    data.append(data2)
with open('./ztf-v2/2020_sne_table.json') as f:
    data3 = json.load(f)
    data.append(data3)
def process(filename):
    data = pd.read_csv(filename)
    df_g = data.groupby(['band']).get_group('g')
    df_r = data.groupby(['band']).get_group('r')
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
    if(redshift == "n/a"):
        redshift = input('Redshift not found, please enter manually for event:'+ event_list[i])
    #print('z =',redshift)
                #print(peak_mjd.mjd)
    return float(redshift)
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
    mod_ztf['Pred_mag'] = pred_mag
    return mod_ztf
all_group = list(model_arr['Model Name'])
def model_param(model_name):
    best = model_arr.groupby(['Model Name']).get_group(model_name)
    return best
#Collecting magnitudes from these models in group
#Creating the texplosion DataFrame
for i in range(0,len(event_list)):
    print(event_list[i])
    z = metadata(event_list[i])
    print('Redshift : ',z)
    #chi_g  = []
    #chi_r = []
    texp = []
    fit_best_g_model = []
    fit_best_r_model = []
    fit_best_chi_g = []
    fit_best_chi_r = []
    #len(all_group)
    #chi_m_arr = pd.DataFrame()
    #chi_m_arr['Model name'] = all_group
    #chi_m_arr['P Mass'] = group['Progenitor mass']
    #chi_m_arr = chi_m_arr.reset_index(drop=True)
    #print(chi_m_arr)
    #z = metadata(event_list[i])
    #print(event_list[i])
    dmod = distance_mod(z)
    #print(z)
    #print(dmod)
    fit_time = pd.DataFrame(columns = ['texplosion','Best_fit_g_model','Red_chi_g','Best_fit_r_model','Red_chi_r'])
    df_g_ztf, df_r_ztf = process(event_list[i]+'.csv')
    g_df = []
    r_df = []
    for k in range(0,80):
        print('Factor of texplosion:',k+1)
        good_mod = []
        chi_g  = []
        chi_r = []
        chi_m_arr = pd.DataFrame()
        chi_m_arr['Model name'] = all_group
        #chi_m_arr['P Mass'] = group['Progenitor mass']
        chi_m_arr = chi_m_arr.reset_index(drop=True)
        texp.append(k+1)
        texp_g = df_g_ztf['MJD'][0] - (k+1)
        texp_r = df_r_ztf['MJD'][0] - (k+1)
        df_g_ztf['epoch'] = df_g_ztf.MJD - texp_g
        df_r_ztf['epoch'] = df_r_ztf.MJD - texp_r
        #print(df_r_ztf)
        #print(df_g_ztf)
        #df_g_ztf = df_g_ztf.drop(df_g_ztf[df_g_ztf.epoch >= 150.0].index)
        #df_r_ztf = df_r_ztf.drop(df_r_ztf[df_r_ztf.epoch >= 150.0].index)
        df_g_ztf['Ab_obs_mag'] = df_g_ztf.magpsf - dmod
        df_r_ztf['Ab_obs_mag'] = df_r_ztf.magpsf - dmod
        df_g_ztf = df_g_ztf.reset_index(drop =True)
        df_r_ztf = df_r_ztf.reset_index(drop =True)
        #print(df_r_ztf)
        #print(df_g_ztf)
        for k in range(0,len(all_group)):
            m_n = all_group[k]
            #print(m_n)
            mod = pd.read_table(path + all_group[k] +'.sdss2', skiprows = 2,names = ['epoch','u','g','r','i','z','kepler'], sep='\s+')
            #print(mod)
            df_g_pred = mod[['epoch','g']]
            df_r_pred = mod[['epoch','r']]
            #list(df_r_ztf.epoch)
            #To fix the problem of interpolation remove elements greater than the interpolation
            #range and modifying the corresponding ztf_obs_dataframe
            df_g_ztf = predicted_mag_at_observed_epochs(df_g_pred,list(df_g_ztf.epoch),df_g_ztf)
            df_r_ztf = predicted_mag_at_observed_epochs(df_r_pred,list(df_r_ztf.epoch),df_r_ztf)
            chi2_g = np.sum((((df_g_ztf.Ab_obs_mag - df_g_ztf.Pred_mag)/df_g_ztf.sigmapsf)**2))/(len(df_g_ztf)-5)
            chi2_r = np.sum((((df_r_ztf.Ab_obs_mag - df_r_ztf.Pred_mag)/df_r_ztf.sigmapsf)**2))/(len(df_r_ztf)-5)
            chi_g.append(chi2_g)
            chi_r.append(chi2_r)
        g_df.append(df_g_ztf)
        r_df.append(df_r_ztf)
        chi_m_arr['reduced_chi_g'] = chi_g
        chi_m_arr['reduced_chi_r'] = chi_r
        #print(chi_m_arr)
        #chi.append(chi_m_arr)
        #plt.errorbar(df_g_ztf['epoch'],df_g_ztf['Ab_obs_mag'],yerr = df_g_ztf['sigmapsf'],fmt='go')
        #Plotting best model for the event in g band
        #print('Best Model in g-band:')
        best_model_g = chi_m_arr.loc[chi_m_arr['reduced_chi_g'].idxmin()]
        #plt.errorbar(df_g_ztf['epoch'],df_g_ztf['Ab_obs_mag'],yerr = df_g_ztf['sigmapsf'],fmt='go')
        best_mod_g = pd.read_table(path + best_model_g[0] +'.sdss2', skiprows = 2,names = ['epoch','u','g','r','i','z','kepler'], sep='\s+')
        #print(best_model_g)
        best_g = model_arr.groupby(['Model Name']).get_group(best_model_g[0])
        print(best_g)
        best_mod_val_g = best_mod_g[['epoch', 'g']]
        best_mod_val_g_r = best_mod_g[['epoch','r']]
        #plt.title(event_list[i])
        #plt.plot(best_mod_val_g['epoch'], best_mod_val_g['g'],'g--',label = best_model_g[0]+'_g_best_fit' )
        #plt.plot(best_mod_val_g_r['epoch'], best_mod_val_g_r['r'],'-',color= 'darkred',label = best_model_g[0]+'_r' )
        #plt.errorbar(df_g_ztf['epoch'],df_g_ztf['Ab_obs_mag'],yerr = df_g_ztf['sigmapsf'],fmt='go')
        #Plotting best model for the event in r band
        #print('Best Model in r-band:')
        best_model_r = chi_m_arr.loc[chi_m_arr['reduced_chi_r'].idxmin()]
        #plt.errorbar(df_r_ztf['epoch'],df_r_ztf['Ab_obs_mag'],yerr = df_r_ztf['sigmapsf'],fmt='ro')
        best_mod_r = pd.read_table(path + best_model_r[0] +'.sdss2', skiprows = 2,names = ['epoch','u','g','r','i','z','kepler'], sep='\s+')
        #print(best_model_r)
        best_r = model_arr.groupby(['Model Name']).get_group(best_model_r[0])
        print(best_r)
        best_mod_val_r = best_mod_r[['epoch', 'r']]
        best_mod_val_r_g = best_mod_r[['epoch', 'g']]
        #plt.title(event_list[i])
        #plt.plot(best_mod_val_r['epoch'], best_mod_val_r['r'],'r--',label = best_model_r[0]+'_r_best_fit')
        #plt.plot(best_mod_val_r_g['epoch'], best_mod_val_r_g['g'],'-',color = 'limegreen',label = best_model_r[0]+'_g')
        #plt.xlim(0,150)
        #plt.ylim(-20, -10)
        #plt.legend()
        #plt.gca().invert_yaxis()
        #plt.show()
        fit_best_chi_g.append(best_model_g[1])
        fit_best_chi_r.append(best_model_r[2])
        fit_best_g_model.append(best_model_g[0])
        fit_best_r_model.append(best_model_r[0])
    fit_time['texplosion'] = texp
    fit_time['Best_fit_g_model'] = fit_best_g_model
    fit_time['Red_chi_g'] = fit_best_chi_g
    fit_time['Best_fit_r_model'] = fit_best_r_model
    fit_time['Red_chi_r'] = fit_best_chi_r
    print(fit_time)
    fit_time_g = fit_time.loc[fit_time['Red_chi_g'].idxmin()]
    fit_time_r = fit_time.loc[fit_time['Red_chi_r'].idxmin()]
    print('Best fit for g and corresponding r value:',fit_time_g)
    print('Best fit for r and corresponding g value:',fit_time_r)
    plt.errorbar(g_df[fit_time_g[0]-1]['epoch'],g_df[fit_time_g[0]-1]['Ab_obs_mag'],yerr = g_df[fit_time_g[0]-1]['sigmapsf'],fmt='go')
    fit_mod_g = pd.read_table(path + fit_time_g[1] +'.sdss2', skiprows = 2,names = ['epoch','u','g','r','i','z','kepler'], sep='\s+')
    fit_val_g = fit_mod_g[['epoch', 'g']]
    plt.title(event_list[i])
    plt.plot(fit_val_g['epoch'], fit_mod_g['g'],'g--',label = fit_time_g[1]+'_g_best_fit_texp'+ str(fit_time_g[0]) )
    plt.errorbar(r_df[fit_time_r[0]-1]['epoch'],r_df[fit_time_r[0]-1]['Ab_obs_mag'],yerr = r_df[fit_time_r[0]-1]['sigmapsf'],fmt='ro')
    fit_mod_r = pd.read_table(path + fit_time_r[1] +'.sdss2', skiprows = 2,names = ['epoch','u','g','r','i','z','kepler'], sep='\s+')
    fit_val_r = fit_mod_r[['epoch', 'r']]
    #print(fit_time_r[0]-1)
    plt.plot(fit_val_r['epoch'], fit_mod_r['r'],'r--',label = fit_time_r[1]+'_r_best_fit_texp'+str(fit_time_r[0]) )
    #Plotting r band in the best fit of texplosion of g_band and the corresponding r model in the best fit of g
    plt.errorbar(r_df[fit_time_g[0]-1]['epoch'],r_df[fit_time_g[0]-1]['Ab_obs_mag'],yerr = r_df[fit_time_g[0]-1]['sigmapsf'],fmt='o',color = 'orange')
    fit_r_in_best_g = pd.read_table(path + fit_time_g[3] +'.sdss2', skiprows = 2,names = ['epoch','u','g','r','i','z','kepler'], sep='\s+')
    fit_above_val = fit_r_in_best_g[['epoch','r']]
    plt.plot(fit_above_val['epoch'], fit_above_val['r'],'--',label = fit_time_g[3]+'_r_best_fit_texp'+str(fit_time_g[0]),color = 'orange' )
    print('Best g_fit with texp'+ str(fit_time_g[0])+':',model_param(fit_time_g[1]))
    print('Best r_fit with texp'+ str(fit_time_g[0])+':',model_param(fit_time_g[3]))
    print('Best g_fit with texp'+ str(fit_time_r[0])+':',model_param(fit_time_r[1]))
    print('Best r_fit with texp'+ str(fit_time_r[0])+':',model_param(fit_time_r[3]))
    plt.legend()
    plt.xlim(0,150)
    plt.ylim(-20, -10)
    plt.gca().invert_yaxis()
    plt.show()
