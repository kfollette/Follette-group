import numpy as np, pandas as pd, scipy.stats as st, matplotlib.pyplot as plt
#get_ipython().run_line_magic('matplotlib', 'inline')
import matplotlib
matplotlib.use('Agg')
from astroquery.simbad import Simbad
from astropy.io import ascii
from astropy import constants as const
from astropy.time import Time
from astropy import units as u
from astropy.coordinates import *
from astropy.modeling.blackbody import blackbody_lambda, blackbody_nu
from astropy.table import Table
from astroplan import (AltitudeConstraint, AirmassConstraint, AtNightConstraint)
from astroplan import FixedTarget, Observer
from astroplan import is_observable, is_always_observable, months_observable
from astroplan import observability_table
from astroplan import FixedTarget, Observer
from astropy.coordinates import EarthLocation
from astroplan.plots import plot_airmass
from astroplan.plots import plot_sky
from astropy.time import Time
import scipy.optimize as optimization
from scipy.optimize import curve_fit
from scipy.stats import linregress
from sklearn.neighbors import KernelDensity
from subprocess import *
import matplotlib.image as mpimg


import sys
import numpy as np
import pandas as pd
from astroquery.simbad import Simbad
from subprocess import *
from flask import Flask
from flask import request, redirect, render_template
from call import create #change import call to "from call.py import create"
from Visibility import plot, vis
from SED import gen
from helper import *

app = Flask(__name__)
app.config['SEND_FILE_MAX_AGE_DEFAULT'] = 0 #attempts to prevent caching


@app.route('/test2.html') #if home if spressed, render home template
def home():
    return render_template("test2.html")

<<<<<<< HEAD
@app.route('/about') #if home if spressed, render home template
def about_page():
    return render_template("about.html")

@app.route('/resources') #if home if spressed, render home template
=======
@app.route('/about') #"about us" page
def about_page():
    return render_template("about.html")

@app.route('/resources') #resources page
>>>>>>> 7505f87630859e40d0fca89b5ed7d5897cc9bdd2
def resources_page():
    return render_template("resources.html")

#@app.route("/test2", methods=["GET", "POST"])
@app.route('/', methods=["GET", "POST"])
def server():

    if request.method == "POST":

        req = request.form


<<<<<<< HEAD
        objects = req.get("Objects")
        mode = req["Mode"]  #change req["mode"] to req.get("Objects")
=======
        objects = req.get("Objects") #get list of objects from 
        mode = req["Mode"]
>>>>>>> 7505f87630859e40d0fca89b5ed7d5897cc9bdd2
        offset = request.form["Offset"]
        date = req.get("Date")


<<<<<<< HEAD
        create(objects,mode,offset) #change call(objects,mode,offset) to create(objects,mode,offset)
        plot(date,objects)
        gen(objects)

        #img = mpimg.imread('static/airmass.png')
        table = txt_to_df('static/data/catalog.txt')
        
        print(table)
        
        vis(date,objects,table)
        visib = txt_to_df('static/data/visibility.txt')
=======
        non_resolvable = create(objects,mode,offset) #create catalog.txt and return list of names not resolvable by Simbad
        #remove non-resolvable names from objects
        for i in non_resolvable:
                objects = objects.replace(i,"").replace(",,", ",")
        if objects[0] == ',':
                objects = objects[1:]
        
        plot(date,objects) #create airmass.png
        gen(objects) #create a SED plot for each object (<object_name>_sed.png)

        table = txt_to_df('static/data/catalog.txt') #generate a Pandas DataFrame from the table in catalog.txt
        
        print(table)
        
        vis(date,objects,table) #create visibility.txt
        visib = txt_to_df('static/data/visibility.txt') #generate a Pandas DataFrame from the table in visibility.txt
>>>>>>> 7505f87630859e40d0fca89b5ed7d5897cc9bdd2
        
        print(visib)
        
        obj_list=list(objects.split(","))
<<<<<<< HEAD
        #return render_template("log.html", text=content, IM=img)
        return render_template("log.html", t1=table.to_html(index=False), t2=visib.to_html(index=False), obj_list=obj_list, len=len(obj_list))
=======
        return render_template("log.html", t1=table.to_html(index=False), t2=visib.to_html(index=False), obj_list=obj_list, len=len(obj_list), nr=non_resolvable, len_nr=len(non_resolvable))
>>>>>>> 7505f87630859e40d0fca89b5ed7d5897cc9bdd2

    else:
        return render_template("test2.html")

if __name__ == "__main__":
    app.run()
