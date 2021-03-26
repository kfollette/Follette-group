#!/usr/bin/env python
# coding: utf-8

# In[1]:

import accretion as a
import accretion_objects as objects
#from importlib import reload
#reload(a)
#reload(objects)


# In[2]:


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
#get_ipython().run_line_magic('matplotlib', 'inline')
import scipy.stats as st
import seaborn as sb
from astropy import constants as const
import random
import astropy.constants as const
import math
from tqdm import tqdm
import extinction as ex
import pdb
import glob
import scipy.optimize as optimization
from matplotlib.animation import FuncAnimation

pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)

import warnings
warnings.filterwarnings('ignore')

import plotly.offline as pyo
import plotly.express as px
import plotly.graph_objs as go
import plotly.tools as tls

# In[3]:


observed1 = pd.read_csv('accdb_updated.csv')
observed1['Upper Limit'] = observed1['Upper Limit'].fillna('No')
nolimit = observed1[observed1['Upper Limit']=='No']


# In[4]:


n = objects.AccretionDistribution(nolimit)
n.bootstrap()
n.UVExcessErrorProp(0, 0, 0, 0, 3, 0, 0, 1, variability=0, age_scatter=False)
df = n.create_df()


# In[5]:


#df


# In[6]:


def build_figure(observed,df):
    
    fig = go.Figure()

    fig.add_trace(go.Scatter(
        x=df["Mass (M$_\odot$)"],
        y=df["Mdot (M$_\odot$)"],
        mode='markers',
        name="Simulated",
        opacity=0.7,
        marker=dict(
        color='salmon'),
    ))


    fig.add_trace(go.Scatter(
        x=observed["Object Mass M_Solar"],
        y=observed["Accretion Rate M_solar yr-1"], 
        mode='markers',
        name='Observed',
        opacity=0.6,
        marker=dict(
        symbol='x',
        color='#44749D'),
    ))

    fig.add_trace(go.Scatter(
        x=observed["Object Mass M_Solar"],
        y=a.empiricalMdot(observed["Object Mass M_Solar"]), 
        mode='lines',
        name='Empirical Relationship',
        opacity=0.9,
        marker=dict(
        color='gray'),
    ))

    fig.update_layout(
        width=1200,
        height=675,
        
        title={
        'text': "Accretion Monte Carlo Error Propagation",
        'y':0.9,
        'x':0.46,
        'xanchor': 'center',
        'yanchor': 'top',
        },

        xaxis_title="Mass (solar masses)",
        yaxis_title="Mass Accretion Rate (solar masses/yr)",
        legend_title=None,

        font=dict(
        ),

        xaxis = dict(
            tickmode = 'array',
            tickvals = [0, 0.01, 0.05, 0.1, 0.5, 1, 1.5, 2],
            showgrid = False
        ),

        yaxis = dict(
            showexponent='all',
            tickmode = 'array',
            tickvals = [1e-13, 1e-12, 1e-11, 1e-10, 1e-9, 1e-8, 1e-7, 1e-6],
            ticktext = [1e-13, 1e-12, 1e-11, 1e-10, 1e-9, 1e-8, 1e-7, 1e-6],
            showgrid = True
        ),
    )
    fig.update_xaxes(type="log")
    fig.update_yaxes(type="log")
    
    return fig


# In[7]:


fig = build_figure(n.observed,df)


# In[8]:


def build_residuals(observed, df):
    
    residual_fig = go.Figure()
    
    observed_difference = np.log10(observed['Accretion Rate M_solar yr-1']) - np.log10(a.empiricalMdot(observed['Object Mass M_Solar']))
    simulated_difference = np.log10(df['Mdot (M$_\\odot$)']) - np.log10(a.empiricalMdot(df['Mass (M$_\\odot$)']))
    
    residual_fig.add_trace(go.Scatter(
        x=df["Mass (M$_\odot$)"],
        y=simulated_difference,
        mode='markers',
        name="Simulated",
        opacity=0.7,
        marker=dict(
        color='salmon'),
    ))


    residual_fig.add_trace(go.Scatter(
        x=observed["Object Mass M_Solar"],
        y=observed_difference, 
        mode='markers',
        name='Observed',
        opacity=0.6,
        marker=dict(
        symbol='x',
        color='#44749D'),
    ))

    residual_fig.add_trace(go.Scatter(
        x=observed["Object Mass M_Solar"],
        y=np.zeros(len(observed["Object Mass M_Solar"])), 
        mode='lines',
        name='Empirical Relationship',
        opacity=0.9,
        marker=dict(
        color='gray'),
    ))

    residual_fig.update_layout(
        width=1200,
        height=575,
        
        title={
        'text': "Residuals",
        'y':0.9,
        'x':0.46,
        'xanchor': 'center',
        'yanchor': 'top',
        },

        xaxis_title="Mass (solar masses)",
        yaxis_title="Residual (log space)",
        legend_title=None,

        font=dict(
        ),

        xaxis = dict(
            tickmode = 'array',
            tickvals = [0, 0.01, 0.05, 0.1, 0.5, 1, 1.5, 2],
            showgrid = False
        ),

        yaxis = dict(
            showexponent='all',
            tickmode = 'array',
            tickvals = [-3, -2, -1, 0, 1, 2, 3],
            showgrid = True
        ),
    )
    residual_fig.update_xaxes(type="log")
    
    return residual_fig


# In[9]:


residual_fig = build_residuals(n.observed, df)


# In[10]:


# current accretion distribution object: n


# In[11]:


import dash
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output


# In[ ]:


external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']

#create app
app = dash.Dash(__name__, external_stylesheets=external_stylesheets)

#create server
server = app.server

#Layout app using HTML
app.layout = html.Div(children=[
    html.Center(children=[
    html.H2(children='Simulated Scatter: Computational Modeling of (Sub)Stellar Accretion Rates'),

    html.Div(children='''
        Joe Palmo
    '''),
    ]),
    
    # Dropdown to choose whether user wants age scatter or not
    html.Div([
            #Dropdown menu to choose whether user wants Age Scatter
            html.Center(children=[
            html.H3("Uncertainties", style={"text-decoration": "underline"}, ),    
            ]),
        
            html.H5("Age Scatter:", style={"font-size": "20px", 'width': '30%', 'display': 'inline-block'}, ),
            html.Div([
            dcc.Dropdown(
                id='age_scatter',
                options=[{'label': i, 'value': i} for i in ['Age Scatter On', 'Age Scatter Off']],
                value='Age Scatter Off'
            )], style={'width': '68%', 'display': 'inline-block'}),
            
            #Input for variability
            html.H5("Variability (dex):", style={"font-size": "20px", 'width': '80%', 'display': 'inline-block'}, ),
            html.Div([
            dcc.Input(
                id='variability',
                type='number',
                value=0.0
            )], style={'width': '5%', 'display': 'inline-block'}),
        
            #Input for spec type uncertainty
            html.H5("Spectral Type Uncertainty (subclasses):", style={"font-size": "20px", 'width': '80%', 'display': 'inline-block'}, ),
            html.Div([
            dcc.Input(
                id='SpTy_uncertainty',
                type='number',
                value=0.0
            )], style={'width': '5%', 'display': 'inline-block'}),
        
            #Input for distance uncertainty
            html.H5("Distance Uncertainty (pc):", style={"font-size": "20px", 'width': '80%', 'display': 'inline-block'}, ),
            html.Div([
            dcc.Input(
                id='distance_uncertainty',
                type='number',
                value=0
            )], style={'width': '5%', 'display': 'inline-block'}),
        
            #Input for age uncertainty
            html.H5("Age Uncertainty (Myr):", style={"font-size": "20px", 'width': '80%', 'display': 'inline-block'}, ),
            html.Div([
            dcc.Input(
                id='age_uncertainty',
                type='number',
                value=0.0
            )], style={'width': '5%', 'display': 'inline-block'}),
        
            #Input for Av uncertainty
            html.H5("Av Uncertainty (mag):", style={"font-size": "20px", 'width': '80%', 'display': 'inline-block'}, ),
            html.Div([
            dcc.Input(
                id='Av_uncertainty',
                type='number',
                value=0.0
            )], style={'width': '5%', 'display': 'inline-block'}),
            
            #Input for bc uncertainty
            html.H5("Bolometric Correction Uncertainty (scale factor):", style={"font-size": "20px", 'width': '80%', 'display': 'inline-block'}, ),
            html.Div([
            dcc.Input(
                id='bc_uncertainty',
                type='number',
                value=0.0
            )], style={'width': '5%', 'display': 'inline-block'}),
                
            #Input for UV Excess uncertainty
            html.H5("UV Excess Uncertainty (% error (as decimal)):", style={"font-size": "20px", 'width': '80%', 'display': 'inline-block'}, ),
            html.Div([
            dcc.Input(
                id='observable_uncertainty',
                type='number',
                value=0.0
            )], style={'width': '5%', 'display': 'inline-block'}),
                
            #Input for Rin uncertainty
            html.H5("Rin Uncertainty (scale factor (of radius)):", style={"font-size": "20px",'width': '80%', 'display': 'inline-block' }, ),
            html.Div([
            dcc.Input(
                id='Rin_uncertainty',
                type='number',
                value=0.0
            )], style={'width': '5%', 'display': 'inline-block'}),
    ],
        style={'width': '40%', 'display': 'inline-block'}),
    
    
    #Eventual Option to choose different inputs
    
    #html.Div([
    
    #dcc.RadioItems(
    #id='input-method',
    #options=[
    #{'label': 'Bootstrap from observed values', 'value': 'bootstrap'},
    #{'label': 'Manual Input', 'value': 'sfr'},],
    #)
    #],style={'width': '60%', 'display': 'inline-block'}),
    
    
    
    #M vs Mdot plot
    dcc.Graph(
        id='MC',
        figure=fig
    ),
    
    #residuals plot
    html.Div(
        dcc.Graph(
        id='residual',
        figure=residual_fig
    ),
    ),
])

###### NEXT STEP #######
#Update Figure with the age_scatter input
@app.callback(
    Output('MC', 'figure'),
    Input('variability', 'value'),
    Input('age_scatter', 'value'),
    Input('SpTy_uncertainty', 'value'),
    Input('distance_uncertainty', 'value'),
    Input('age_uncertainty', 'value'),
    Input('Av_uncertainty', 'value'),
    Input('bc_uncertainty', 'value'),
    Input('observable_uncertainty', 'value'),
    Input('Rin_uncertainty', 'value'),)
def update_figure(variability, age_scatter, SpTy_uncertainty, distance_uncertainty, age_uncertainty, Av_uncertainty,
                  bc_uncertainty, observable_uncertainty, Rin_uncertainty):
    
    #n = objects.AccretionDistribution(nolimit)
    #n.bootstrap()
    
    if age_scatter == 'Age Scatter On':
        n.UVExcessErrorProp(SpTy_uncertainty, distance_uncertainty, age_uncertainty, Av_uncertainty, 3, bc_uncertainty, observable_uncertainty, 1, RinUnc=Rin_uncertainty, variability=variability, age_scatter=True)
    else:
        n.UVExcessErrorProp(SpTy_uncertainty, distance_uncertainty, age_uncertainty, Av_uncertainty, 3, bc_uncertainty, observable_uncertainty, 1, RinUnc=Rin_uncertainty, variability=variability, age_scatter=False)
        
    df = n.create_df()
    
    fig = build_figure(n.observed, df)
    
    residual_fig = build_residuals(n.observed, df)
    
    fig.update_layout(transition_duration=500)

    return fig

@app.callback(
    Output('residual', 'figure'),
    Input('variability', 'value'),
    Input('age_scatter', 'value'),
    Input('SpTy_uncertainty', 'value'),
    Input('distance_uncertainty', 'value'),
    Input('age_uncertainty', 'value'),
    Input('Av_uncertainty', 'value'),
    Input('bc_uncertainty', 'value'),
    Input('observable_uncertainty', 'value'),
    Input('Rin_uncertainty', 'value'),)
def update_residuals(variability, age_scatter, SpTy_uncertainty, distance_uncertainty, age_uncertainty, Av_uncertainty,
                  bc_uncertainty, observable_uncertainty, Rin_uncertainty):
    
    #n = objects.AccretionDistribution(nolimit)
    #n.bootstrap()
    
    if age_scatter == 'Age Scatter On':
        n.UVExcessErrorProp(SpTy_uncertainty, distance_uncertainty, age_uncertainty, Av_uncertainty, 3, bc_uncertainty, observable_uncertainty, 1, RinUnc=Rin_uncertainty, variability=variability, age_scatter=True)
    else:
        n.UVExcessErrorProp(SpTy_uncertainty, distance_uncertainty, age_uncertainty, Av_uncertainty, 3, bc_uncertainty, observable_uncertainty, 1, RinUnc=Rin_uncertainty, variability=variability, age_scatter=False)
       
    df = n.create_df()
    
    residual_fig = build_residuals(n.observed, df)
    
    residual_fig.update_layout(transition_duration=500)

    return residual_fig

if __name__ == '__main__':
    app.run_server(debug=False)
    


# In[ ]:





# In[ ]:




