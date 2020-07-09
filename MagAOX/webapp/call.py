#!/usr/bin/env python
# coding: utf-8

# #### Goal: Create a prototype web tool that turns a list of Simbad-resolvable astronomical object names into a properly formatted Magellan catalog entry.

# > Please consult the Las Campanas Observatory website: http://www.lco.cl/telescopes-information/magellan/instruments/observing-catalogs/observing-catalogs-format/?searchterm=catalog

# #### Imports

# In[ ]:

import sys
import numpy as np
import pandas as pd
from astroquery.simbad import Simbad
from subprocess import *
from flask import Flask
from flask import request, redirect, render_template

def create(objects,mode,offset):
    serEntered_list = objects
    InstRotOffset = offset
    InstRotOffsetMode = mode #InstRotOffsetMode is changed from 'OFF' to mode
    equinox = 2000
    epoch_pm = 2020.7

    userEntered_list = list(serEntered_list.split(","))
    # > There will not be a guide star specified for MagAOX observations, so fields 10-15 are 0's. These 6 required fields are to specify positions of the guide probes, which there will be none.

    # #### Checking if names are resolvable

    # In[ ]:


    customSimbad = Simbad()

    target_list = []

    x = 0

    for i in userEntered_list:
        star = customSimbad.query_object(userEntered_list[x])
        #print(star)
        if type(star) == type(None):
            print('\033[1m' + 'ERROR: The star [ '
                  + userEntered_list[x] + ' ] was not resolvable in Simbad. Please check the spelling of the name.')
            print('\033[0m')

        if type(star) != type(None):
            target_list.append(userEntered_list[x])

        x = x + 1

    target_list

    # #### Magellan Catalog Entry Format

    # In[ ]:


    MagCatEntry = pd.DataFrame(columns=['#',
                                        ' ',
                                        'RA',
                                        'Dec',
                                        'equinox',
                                        'RApm',
                                        'Decpm',
                                        'offset',
                                        'rot',
                                        'RA_probe1',
                                        'Dec_probe1',
                                        'equinox_probe1',
                                        'RA_probe2',
                                        'Dec_probe2',
                                        'equinox_probe2',
                                        'pm_epoch'],
                               index=[' '])

    # Appending Header
    MagCatEntry.at['0', '#'] = '###'
    MagCatEntry.at['0', ' '] = 'name'
    MagCatEntry.at['0', 'RA'] = 'hh:mm:ss.s'
    MagCatEntry.at['0', 'Dec'] = 'sdd:mm:ss'
    MagCatEntry.at['0', 'equinox'] = 'yyyy.0'
    MagCatEntry.at['0', 'RApm'] = 's.ss'
    MagCatEntry.at['0', 'Decpm'] = 's.ss'
    MagCatEntry.at['0', 'equinox'] = 'yyyy.0'
    MagCatEntry.at['0', 'offset'] = 'angle'
    MagCatEntry.at['0', 'rot'] = 'mode'
    MagCatEntry.at['0', 'RA_probe1'] = 'hh:mm:ss.s'
    MagCatEntry.at['0', 'Dec_probe1'] = 'sdd:mm:ss'
    MagCatEntry.at['0', 'equinox_probe1'] = 'yyyy.0'
    MagCatEntry.at['0', 'RA_probe2'] = 'hh:mm:ss.s'
    MagCatEntry.at['0', 'Dec_probe2'] = 'sdd:mm:ss'
    MagCatEntry.at['0', 'equinox_probe2'] = 'yyyy.0'
    MagCatEntry.at['0', 'pm_epoch'] = 'yyyy.0'

    x = 0
    for i in target_list:
        MagCatEntry = MagCatEntry.append(pd.Series(), ignore_index=True)
        x = x + 1

    MagCatEntry = MagCatEntry.drop([0])
    MagCatEntry = MagCatEntry.reset_index()
    MagCatEntry = MagCatEntry.drop(columns=['index'])

    # MagCatEntry


    # #### Custom Simbad, used to query Simbad

    # In[ ]:


    # customSimbad = Simbad()
    customSimbad.add_votable_fields('ra(s)')
    customSimbad.add_votable_fields('dec(s)')
    customSimbad.add_votable_fields('pmra')
    customSimbad.add_votable_fields('pmdec')

    # #### Fields 1 – 7

    # In[ ]:


    columns = {'MAIN_ID': [0], 'RA': [0], 'DEC': [0], 'PMRA': [0], 'PMDEC': [0]}
    targets = pd.DataFrame(data=columns)

    x = 0
    for i in target_list:
        star = customSimbad.query_object(target_list[x])
        star.keep_columns(['MAIN_ID', 'RA', 'DEC', 'PMRA', 'PMDEC'])
        star = star.to_pandas()
        targets = targets.append(star)
        x = x + 1

    targets = targets.iloc[1:]
    targets = targets.values.tolist()
    # targets


    # #### Appending the Series

    # #### Appending Rows 5, 8 – 16 (static)

    # In[ ]:


    x = 1
    for i in target_list:
        MagCatEntry.at[x, 'equinox'] = equinox
        MagCatEntry.at[x, 'offset'] = InstRotOffset
        MagCatEntry.at[x, 'rot'] = InstRotOffsetMode
        MagCatEntry.at[x, 'pm_epoch'] = epoch_pm
        x = x + 1

    # MagCatEntry


    # #### Appending Rows 1-4

    # In[ ]:


    x = 1
    y = 0
    for i in target_list:
        MagCatEntry.at[x, '#'] = '{0:03}'.format(x)
        MagCatEntry.at[x, ' '] = targets[y][0].decode('utf-8')
        MagCatEntry.at[x, 'RA'] = targets[y][1]
        MagCatEntry.at[x, 'Dec'] = targets[y][2]
        MagCatEntry.at[x, 'equinox'] = equinox
        MagCatEntry.at[x, 'RApm'] = targets[y][3]
        MagCatEntry.at[x, 'Decpm'] = targets[y][4]
        x = x + 1
        y = y + 1

    MagCatEntry = MagCatEntry.fillna(0)
    MagCatEntry

    # #### Saving as .txt

    # In[ ]:


    np.savetxt('static/catalog.txt', MagCatEntry, fmt="%-15s",
               header='               name           RA              Dec             equinox         RApm            Decpm           offset          rot             RA_probe1       Dec_probe1      equinox_probe1         RA_probe2       Dec_probe2      equinox_probe2         pm_epoch')
