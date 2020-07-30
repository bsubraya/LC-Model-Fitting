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
data = []
with open('./ztf-v2/2018_sne_table.json') as f:
    data1 = json.load(f)
    data.append(data1)
with open('./ztf-v2/2019_sne_table.json') as f:
    data2 = json.load(f)
    data.append(data2)
with open('./ztf-v2/2020_sne_table.json') as f:
    data3 = json.load(f)
    data.append(data3)

sn_name = []
ztf_name= []
ztf_type = []
for i in range(0,len(data)):
    for j in range(0,len(data[i])):
        sn_name.append(data[i][j]['Supernova Name'])
        ztf_name.append(data[i][j]['lasairname'])
        ztf_type.append(data[i][j]['Type'])
event_df = pd.DataFrame(columns = ['ZTF_event','Type'])
event_df['ZTF_event'] = ztf_name
event_df['Type'] = ztf_type
ztf_type
type_2 = event_df.drop(event_df[event_df['Type'] != 'IIP'].index)
type_2 = type_2.sort_values(by='ZTF_event',ascending=True)
type_2 = type_2.reset_index(drop =True)
type_2
len(type_2['ZTF_event'])
event_list = ['ZTF19aadnwbv','ZTF19abcneik','ZTF19aavkptg','ZTF18abrlljc','ZTF18abckutn','ZTF19acxowrr','ZTF19aceckht']
def plotting(fid,df):
    if (('g' in fid) and ('r' in fid)):
        df_g = df.groupby(['band']).get_group(('g'))
        df_r = df.groupby(['band']).get_group(('r'))
        df_g = df_g.reset_index(drop =True)
        df_r = df_r.reset_index(drop =True)
        print('Both bands')
        if 'sigmapsf' in df.columns:
            plt.errorbar(df_g['MJD'],df_g['magpsf'],yerr = df_g['sigmapsf'],fmt='go')
            plt.errorbar(df_r['MJD'],df_r['magpsf'],yerr = df_r['sigmapsf'],fmt='ro')
            plt.title(type_2['ZTF_event'][j])
        else:
            plt.plot(df_g['MJD'],df_g['magpsf'],'gx')
            plt.plot(df_r['MJD'],df_r['magpsf'],'rx')
            plt.title(type_2['ZTF_event'][j])
    elif('g' not in fid):
        df_r = df.groupby(['band']).get_group(('r'))
        print('Only r band')
        if 'sigmapsf' in df.columns:
            plt.errorbar(df_r['MJD'],df_r['magpsf'],yerr = df_r['sigmapsf'],fmt='ro')
            plt.title(type_2['ZTF_event'][j])
        else:
            plt.plot(df_r['MJD'],df_r['magpsf'], 'rx')
            plt.title(type_2['ZTF_event'][j])
    elif('r' not in fid):
        df_g = df.groupby(['band']).get_group(('g'))
        print('Only g band')
        if 'sigmapsf' in df.columns:
            plt.errorbar(df_g['MJD'],df_g['magpsf'],yerr = df_g['sigmapsf'],fmt='go')
            plt.title(type_2['ZTF_event'][j])
        else:
            plt.plot(df_g['MJD'],df_g['magpsf'], 'gx')
            plt.title(type_2['ZTF_event'][j])
type_2
for j in range(0,len(type_2['ZTF_event'])):
    mjd = []
    magpsf = []
    sigmapsf = []
    fid = []
    upperlimits = []
    mjd_upper = []
    u_fid = []
    with open('./ztf-v2/'+type_2['ZTF_event'][j]+'.json') as f:
        event_details = json.load(f)
    if 'candidates' in event_details:
        print(type_2['ZTF_event'][j],type_2['Type'][j])
        for m in range(0,len(event_details['candidates'])):
            if 'sigmapsf' in event_details['candidates'][m]:
                mjd.append(event_details['candidates'][m]['mjd'])
                magpsf.append(event_details['candidates'][m]['magpsf'])
                sigmapsf.append(event_details['candidates'][m]['sigmapsf'])
                fid.append(event_details['candidates'][m]['fid'])
            if ('magpsf' in event_details['candidates'][m] and 'sigmapsf' not in event_details['candidates'][m]):
                mjd_upper.append(event_details['candidates'][m]['mjd'])
                upperlimits.append(event_details['candidates'][m]['magpsf'])
                u_fid.append(event_details['candidates'][m]['fid'])
        filter = {1:'g',2:'r'}
        fid = list(map(filter.get, fid))
        df = pd.DataFrame(columns = ['MJD','magpsf','sigmapsf','band'])
        upper_df = pd.DataFrame(columns=['MJD','magpsf','band'])
        upper_df['MJD'] = mjd_upper
        upper_df['magpsf'] = upperlimits
        u_fid = list(map(filter.get, u_fid))
        upper_df['band'] = u_fid
        df['MJD'] = mjd
        df['magpsf'] = magpsf
        df['sigmapsf'] = sigmapsf
        df['band'] = fid
        df = df.sort_values(by='MJD',ascending=True)
        df = df.reset_index(drop =True)
        upper_df = upper_df.sort_values(by='MJD',ascending=True)
        upper_df = upper_df.reset_index(drop =True)
        print(upper_df)
        print(df)
        df.to_csv(type_2['ZTF_event'][j]+'.csv',index = False)
        upper_df.to_csv(type_2['ZTF_event'][j]+'_upper.csv',index = False)
        plotting(fid,df)
        #if len(upper_df) != 0 :
            #plotting(u_fid,upper_df)
        plt.gca().invert_yaxis()
        plt.ylim(24,16)
        plt.show()
