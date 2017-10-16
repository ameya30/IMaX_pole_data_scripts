import os
import csv
import glob
import datetime
import sufi_imax_align 

import numpy as np

from   astropy.io import fits
################################# FUCNTIONS #################################

# Finds the avegare time

def avg_time(times):
    avg = 0

    for elem in times:
        avg += elem.second + 60*elem.minute + 3600*elem.hour 

    avg /= len(times)
    rez = '2009-06-12T' + str(avg/3600) + ':' + str((avg%3600)/60) + ':' + str(avg%60)

    return datetime.datetime.strptime(rez, '%Y-%m-%dT%H:%M:%S')

# Gives average time for IMaX cycle

def find_imax_time(filename):

    time_array_cam1 = []
    
    with open(filename, 'rb') as cam1:
        date_cam1 = csv.reader(cam1, delimiter = ' ', quotechar = '|')
    
        for row in date_cam1:
            date = row[1].split("'")[1]
            
            time = datetime.datetime.strptime(date,'%Y-%m-%dT%H:%M:%S.%f')
            time_array_cam1.append(time)
    
    ave_cam1 = avg_time(time_array_cam1)
    
    return ave_cam1

# Computes the difference between the acquisition times of two images

def diff_time(t1, t2):
    h1, m1, s1 = t1.hour, t1.minute, t1.second
    h2, m2, s2 = t2.hour, t2.minute, t2.second
    t1_secs = s1 + (60*m1)
    t2_secs = s2 + (60*m2)

    return (t2_secs - t1_secs)

# Selects closest (in time) sufi to imax

def close_sufi(imax_time, sufi_list):

    files = glob.glob('../Data/sufi397/*.fits.gz')
    time_array = np.zeros(len(files))
    n = 0
    for sufi_file in files:
        sufi_im = fits.open(sufi_file, ignore_missing_end = True)
        sufi_header = sufi_im[0].header
        sufi_date = sufi_header['DATE_OBS']
        sufi_time = datetime.datetime.strptime(sufi_date,'%Y-%m-%dT%H:%M:%S.%f').time()
        
        d = diff_time(imax_time, sufi_time)
        time_array[n] = abs(d)
        n = n+1
        sufi_min = files[np.where(time_array == time_array.min())[0][0]] 
        
    return sufi_min

################################# BEGIN PROGRAM ################################# 

csv_list = glob.glob('../Data/date_obs_cam1_*.csv')

with open('../Data/imax_num_closest_sufi_397.csv', 'wb') as csvfile:
    for el in csv_list:
           
        cycle_num = el.split('_')[3]
        cycle_num = cycle_num.split('.')[0]
        
        imax_time = find_imax_time(el)

        info = csv.writer(csvfile, delimiter = ' ', quotechar = '|', quoting = csv.QUOTE_MINIMAL)
        info.writerow('cycle num: ' + cycle_num + ', sufi_im: ' + close_sufi(imax_time, '../Data/sufi397'))
    #print el

    #cycle_num = el.split('_')[3]
    #cycle_num = cycle_num.split('.')[0]
    #
    #imax_time = find_imax_time(el)
    #
    #with open('../Data/cycle_num_' + str(cycle_num) + '.csv', 'wb') as csvfile:
    #    info = csv.writer(csvfile, delimiter = ' ', quotechar = '|', quoting = csv.QUOTE_MINIMAL)
    #    info.writerow('cycle num: ' + cycle_num + ', sufi_im: ' + close_sufi(imax_time, '../Data/sufi397'))
