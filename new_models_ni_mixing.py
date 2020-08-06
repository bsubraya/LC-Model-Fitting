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
#Path to all model data
path ='/Users/bhagyasubrayan/Desktop/Explosion Parameters/Takashi\'s models/new_models_SDSS/'
model_names = os.listdir(path)
model_names.remove('.DS_Store')
model_names = np.sort(model_names)
#len(model_names)
#To query the metadata of ZTF events
data=[]
def event_data(filename):
    with open(filename) as f:
        file_data = json.load(f)
        data.append(file_data)
event_data('./ztf-v2/2018_sne_table.json')
event_data('./ztf-v2/2019_sne_table.json')
event_data('./ztf-v2/2020_sne_table.json')
#Function to group the band data from a source paper or just from the given dataframe
def process(filename):
    file = pd.read_csv(filename)
    df_g = file.groupby(['band']).get_group('g')
    df_r = file.groupby(['band']).get_group('r')
    return df_g.reset_index(drop =True), df_r.reset_index(drop=True)
def process_known(df,band,source):
    if source == '0':
        #print(True)
        new_df = df.groupby('band').get_group(band)
    else:
        new_df = df.groupby(['band','source']).get_group((band,source))
    new_df = new_df.reset_index(drop =True)
    df_req = new_df[['time','magnitude','e_magnitude','band']]
    df_req.columns = ['MJD','magpsf','sigmapsf','band']
    return df_req
#Function to calculate distance modulus
def distance_mod(z):
    cosmo = FlatLambdaCDM(H0=70, Om0=0.3, Tcmb0=2.725)
    dmod = cosmo.distmod(z)
    return dmod.value
#Function designes for established events
def metadata_known(eventname):
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
#Corrects the dmod and ebv values for the given observed photometry
def final_obs_df_known(eventname,df,t,band):
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
#Interpolates for a given observed data for the epochs and
#gives the predicted value from the corresponding model
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
#Plotting function that plots the observed data to the model at the fit texplosion
def plot_model_event(df,name,band,t,color):#Add upper df if upper limits
    plt.errorbar(df['epoch'],df['Ab_obs_mag'],yerr = df['mag_error'],fmt = color+'o',label = 'Observed '+band+' with delay time : '+ str(t))
    mod = pd.read_table(path + name +'.sdss2', skiprows = 2,names = ['epoch','u','g','r','i','z','kepler'], sep='\s+')
    df_pred = mod[['epoch',band]]
    plt.plot(df_pred['epoch'],df_pred[band],color+'--',label = 'Best_fit : '+name)
    #plt.plot(upper['epoch'],upper['Ab_obs_mag'],band+'*',label = 'Upper Limits')
    plt.legend()
    plt.xlim(-40,250)
    plt.ylim(-20, -10)
def plot_model(df,name,band,t):#Add upper df if upper limits
    plt.errorbar(df['epoch'],df['Ab_obs_mag'],yerr = df['mag_error'],fmt = band+'o',label = 'Observed '+band+' with delay time : '+ str(t))
    mod = pd.read_table(path + name , skiprows = 2,names = ['epoch','u','g','r','i','z','kepler'], sep='\s+')
    df_pred = mod[['epoch',band]]
    plt.plot(df_pred['epoch'],df_pred[band],band+'--',label = 'Best_fit : '+name)
    #plt.plot(upper['epoch'],upper['Ab_obs_mag'],band+'*',label = 'Upper Limits')
    plt.legend(prop={"size":10})
    plt.xlim(-40,200)
    plt.ylim(-20, -10)

#Routine to return the chi^2 value for all the models
def fit_all_models(df,band):
    chi = []
    for k in range(0,len(model_names)):
        m_n = model_names[k]
        #print(m_n)
        mod = pd.read_table(path + model_names[k], skiprows = 2,names = ['epoch','u','g','r','i','z','kepler'], sep='\s+')
        #print(mod)
        df_pred = mod[['epoch',band]]
        #To fix the problem of interpolation remove elements greater than the interpolation
        #range and modifying the corresponding ztf_obs_dataframe
        df_final = predicted_mag_at_observed_epochs(df_pred,list(df.epoch),df)
        chi2 = np.sum((((df_final.Ab_obs_mag - df_final.Pred_mag)/df_final.mag_error)**2))/(len(df_final)- 7)
        #chi2_r = np.sum((((df_r_ztf.Ab_obs_mag - df_r_ztf.Pred_mag)/df_r_ztf.sigmapsf)**2))/(len(df_r_ztf)-5)
        #chi_g.append(chi2_g)
        #chi_r.append(chi2_r)
        #print(chi2)
        chi.append(chi2)#Value of reduced chi^^2 for each model
    return chi
#Creating a copy for final plotting
def subset_df(df):
    new_df= df[['MJD','magpsf','sigmapsf','band']]
    return new_df
#Routine to fit the best explosion time after fitting the best model in each band
def best_texp(df,band,min,max):
    texp = []
    fit_best_model = []
    fit_best_chi = []
    fit_time = pd.DataFrame(columns = ['texplosion','Best_fit_model','Reduced_chi'])
    for k in np.arange(min,max,1):
        print('Factor of texplosion:',k)
        chi_m_arr = pd.DataFrame()
        chi_m_arr['Model name'] = model_names
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
    print(fit_time)
    return fit_time, best_fit[0] , best_fit[1],best_fit[2]
#Final called function
def analysis(df,band,min,max):#Add upper if upper limits
    original_obs = subset_df(df)
    time_array, texplosion, fit_model, fit_value = best_texp(df,band,min,max)
    mod_df = final_obs_df(original_obs,original_obs.MJD[0] - texplosion,band = band)
    #mod_upper = final_obs_df(upper,original_obs.MJD[0]-texplosion)
    #print(mod_df)
    plot_model(mod_df, name = fit_model,band = band, t = texplosion)
    #best_fit_param = model_param(fit_model)
    #print(best_fit_param)
    print('Reduced_chi_value for '+ band+ ':', fit_value )
    return texplosion, fit_model
event_g = ['ZTF19abcneik','ZTF19aavkptg']
for i in range(0,2):
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
    texplo_g, model_g = analysis(df_g_ztf,band = 'g',min= 0,max = 8)#upper = g_upper
    texplo_r,model_r = analysis(df_r_ztf,band = 'r', min = texplo_g - 3 ,max = texplo_g + 5)#upper = r_upper,
    plt.title(event_g[i]+ ' ',fontsize = '15')
    plt.ylabel('Absolute Magnitude',fontsize = '15')
    plt.xlabel('Days from texplosion',fontsize = '15')
    plt.rc('xtick', labelsize= 15)
    plt.rc('ytick', labelsize= 15)
    plt.gca().invert_yaxis()
    plt.show()
