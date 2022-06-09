#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 20 16:50:41 2020

@author: koealaquillage
"""

#*********************************#
# Script to observe the dispersion#
# Ellipses of my data. Taken from #
# A matplotlib tutorial at https://matplotlib.org/devdocs/gallery/statistics/confidence_ellipse.html #
# By G. Koenig, the 21/12/2020    #
#*********************************#

#*******Packages import***********#
#***Packages import**********#
import numpy as np
from scipy import io
import scipy.interpolate
import utide
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse,Rectangle
import matplotlib.transforms as transforms
#*****Useful functions************#
def import_data_Cristele(eta_data) :
    """ Function to import the data from the matlab files of Cristele.
    INPUTS :
    eta_data : Dataframe of the .mat file
    OUTPUTS :
    data_frame : Dataframe with the date, the longitude or latitude positions and
    the pressure"""
    
    # We load data from the Matlab files "
    eta=scipy.io.loadmat(eta_data)
    # Now we format the time, so that it is easier to read #
    
    df = pd.DataFrame({'year':eta['Temps']['year'][0][0][:,0],
                  'month':eta['Temps']['month'][0][0][:,0],
                   'day':eta['Temps']['day'][0][0][:,0],
                    'hour':eta['Temps']['hour'][0][0][:,0],
                    'minute':eta['Temps']['minute'][0][0][:,0],
                     'second':eta['Temps']['year'][0][0][:,0]})
    
    
    return pd.DataFrame({'Latitude' : eta['P'][0,0][0][0][0],
                          'Longitude' : eta['P'][0,0][1][0][0],
                          'Depth_surf' : eta['P'][0,0][2][:,0],
                          'Depth_1' : eta['P_Adcp'][0][0][0][0][0],
                          'Depth_2' : eta['P_Adcp'][0][0][0][0][1],
                          'Depth_3' : eta['P_Adcp'][0][0][0][0][2],
                          'Depth_4' : eta['P_Adcp'][0][0][0][0][3],
                          'Depth_5' : eta['P_Adcp'][0][0][0][0][4],
                          'Depth_6' : eta['P_Adcp'][0][0][0][0][5],
                          'Depth_7' : eta['P_Adcp'][0][0][0][0][6],
                          'Depth_8' : eta['P_Adcp'][0][0][0][0][7],
                          'Depth_9' : eta['P_Adcp'][0][0][0][0][8],
                          'Depth_10' : eta['P_Adcp'][0][0][0][0][9],
                          'Depth_11' : eta['P_Adcp'][0][0][0][0][10],
                          'Depth_12' : eta['P_Adcp'][0][0][0][0][11],
                          'Depth_13' : eta['P_Adcp'][0][0][0][0][12],
                          'Depth_14' : eta['P_Adcp'][0][0][0][0][13],
                          'Depth_15' : eta['P_Adcp'][0][0][0][0][14],
#                          'Depth_16' : eta['P_Adcp'][0][0][0][0][15],
#                          'Depth_17' : eta['P_Adcp'][0][0][0][0][16],
#                          'Depth_18' : eta['P_Adcp'][0][0][0][0][17],
#                          'Depth_19' : eta['P_Adcp'][0][0][0][0][18],
#                          'Depth_20' : eta['P_Adcp'][0][0][0][0][19],
#                          'Depth_21' : eta['P_Adcp'][0][0][0][0][20],
                          'U_1': eta['vitesse'][0][0][0][:,0],
                          'U_2': eta['vitesse'][0][0][0][:,1],
                          'U_3': eta['vitesse'][0][0][0][:,2],
                          'U_4': eta['vitesse'][0][0][0][:,3],
                          'U_5': eta['vitesse'][0][0][0][:,4],
                          'U_6': eta['vitesse'][0][0][0][:,5],
                          'U_7': eta['vitesse'][0][0][0][:,6],
                          'U_8': eta['vitesse'][0][0][0][:,7],
                          'U_9': eta['vitesse'][0][0][0][:,8],
                          'U_10': eta['vitesse'][0][0][0][:,9],
                          'U_11': eta['vitesse'][0][0][0][:,10],
                          'U_12': eta['vitesse'][0][0][0][:,11],
                          'U_13': eta['vitesse'][0][0][0][:,12],
                          'U_14': eta['vitesse'][0][0][0][:,13],
                          'U_15': eta['vitesse'][0][0][0][:,14],
#                          'U_16': eta['vitesse'][0][0][0][:,15],
#                          'U_17': eta['vitesse'][0][0][0][:,16],
#                          'U_18': eta['vitesse'][0][0][0][:,17],
#                          'U_19': eta['vitesse'][0][0][0][:,18],
#                         'U_20': eta['vitesse'][0][0][0][:,19],
#                          'U_21': eta['vitesse'][0][0][0][:,20],
                          'V_1': eta['vitesse'][0][0][1][:,0],
                          'V_2': eta['vitesse'][0][0][1][:,1],
                          'V_3': eta['vitesse'][0][0][1][:,2],
                          'V_4': eta['vitesse'][0][0][1][:,3],
                          'V_5': eta['vitesse'][0][0][1][:,4],
                          'V_6': eta['vitesse'][0][0][1][:,5],
                          'V_7': eta['vitesse'][0][0][1][:,6],
                          'V_8': eta['vitesse'][0][0][1][:,7],
                          'V_9': eta['vitesse'][0][0][1][:,8],
                          'V_10': eta['vitesse'][0][0][1][:,9],
                          'V_11': eta['vitesse'][0][0][1][:,10],
                          'V_12': eta['vitesse'][0][0][1][:,11],
                          'V_13': eta['vitesse'][0][0][1][:,12],
                          'V_14': eta['vitesse'][0][0][1][:,13],
                          'V_15': eta['vitesse'][0][0][1][:,14],
#                          'V_16': eta['vitesse'][0][0][1][:,15],
#                          'V_17': eta['vitesse'][0][0][1][:,16],
#                          'V_18': eta['vitesse'][0][0][1][:,17],
#                          'V_19': eta['vitesse'][0][0][1][:,18],
#                          'V_20': eta['vitesse'][0][0][1][:,19],
#                          'V_21': eta['vitesse'][0][0][1][:,20],
                          },index=pd.to_datetime(df))
    
       
# And a function for the ellipses
def confidence_ellipse(x, y, ax, n_std=3.0, facecolor='none', **kwargs):
    """
    Create a plot of the covariance confidence ellipse of *x* and *y*.

    Parameters
    ----------
    x, y : array-like, shape (n, )
        Input data.

    ax : matplotlib.axes.Axes
        The axes object to draw the ellipse into.

    n_std : float
        The number of standard deviations to determine the ellipse's radiuses.

    **kwargs
        Forwarded to `~matplotlib.patches.Ellipse`

    Returns
    -------
    matplotlib.patches.Ellipse
    """
    if x.size != y.size:
        raise ValueError("x and y must be the same size")

    cov = np.cov(x, y)
    pearson = cov[0, 1]/np.sqrt(cov[0, 0] * cov[1, 1])
    # Using a special case to obtain the eigenvalues of this
    # two-dimensionl dataset.
    ell_radius_x = np.sqrt(1 + pearson)
    ell_radius_y = np.sqrt(1 - pearson)
    ellipse = Ellipse((0, 0), width=ell_radius_x * 2, height=ell_radius_y * 2,
                      facecolor=facecolor, **kwargs)

    # Calculating the stdandard deviation of x from
    # the squareroot of the variance and multiplying
    # with the given number of standard deviations.
    scale_x = np.sqrt(cov[0, 0]) * n_std
    mean_x = np.mean(x)

    # calculating the stdandard deviation of y ...
    scale_y = np.sqrt(cov[1, 1]) * n_std
    mean_y = np.mean(y)

    transf = transforms.Affine2D() \
        .rotate_deg(45) \
        .scale(scale_x, scale_y) \
        .translate(mean_x, mean_y)

    ellipse.set_transform(transf + ax.transData)
    return ax.add_patch(ellipse)
#*******Now I define the paths to the data and the data I want to see
# List of station
       
# The data position
path_data = '/home/koenig.g/These/DONNEES/DonneesPropres_2016/'
            
data_file = 'Vit_GDigo.mat'

# The dictionnary to store data
dic_stat = []


# A list to store the different coefficient
list_coeff = []

# I need a list of levels and dictionnary for those levels
list_levels=['1','5','10']

# And dictionnary for those
dic_labels = {'1':'First level',
              '5':'Fifth level',
              '10':'Tenth level'}
# And a dictionnary for colors
dic_colors = {'1':'firebrick',
              '5':'lightseagreen',
              '10':'deeppink'}
# And a dictionnary for linestyle
dic_linestyle = {'1':'-',
                 '5':':',
                 '10':'-.'}

#*******Now I try and import it*********#
dic_stat = import_data_Cristele(path_data+data_file)
    
#*** And if I try the plotting
fig, ax = plt.subplots()

# First I do nice lines to identify the center
ax.axvline(c='grey', lw=1)
ax.axhline(c='grey', lw=1)

for lev in list_levels :
    # I scatter my data
    ax.scatter(dic_stat['U_'+lev][::5].values/1000.,
               dic_stat['V_'+lev][::5].values/1000.,
               s=0.15,color=dic_colors[lev],alpha=0.15)

for lev in list_levels :
    # And then the ellipse 
    confidence_ellipse(dic_stat['U_'+lev][::5].values/1000., 
                       dic_stat['V_'+lev][::5].values/1000.,
                       ax, n_std=1,edgecolor=dic_colors[lev],
                       label=dic_labels[lev])

# And the legend
ax.legend(loc='best')
# We add some labels
ax.set_ylabel(r'Meridional velocity ($m.s^{-1}$)')
ax.set_xlabel(r'Zonal velocity ($m.s^{-1}$)')
ax.set_title('Dispersion ellipse')

ax.autoscale()
plt.show()