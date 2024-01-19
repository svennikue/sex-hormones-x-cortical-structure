#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan 14 19:23:42 2024

@author: skuech

glasser parcellation sex differences

"""


from scipy.io import loadmat
import numpy as np
import pandas as pd
from brainspace.plotting import plot_hemispheres
from brainspace.datasets import load_conte69
from brainspace.utils.parcellation import map_to_labels



data_dir = "/Users/skuech/Documents/my_projects/female_gradients/data/glasser_hcp_ind.mat"
lh, rh = load_conte69()
data = loadmat(f"{data_dir}")

intensity = data['glasser_hcp_ind'][0][0][0]
mask = np.zeros(1206)
for i in range(1206):
  if intensity[i].mean() < 1:
    mask[i] = 0
  else:
    mask[i] = 1
all_1206 = pd.read_csv('hcp1200_id_sex_age.csv')
intensity = intensity[mask==1]
phenotype = all_1206[mask==1]
glasser = np.genfromtxt('glasser.csv')
glasser_l = glasser[:32492]
glasser_r = glasser[32492:]
glasser_r[glasser_r==180]=0
glasser = np.concatenate((glasser_l,glasser_r))
mask_brain = glasser != 0
plot = map_to_labels(intensity.mean(axis=0), glasser, mask=mask_brain)
plot[plot==0]=np.nan
def plot_surface_lr(lh, rh,
                    data,
                    size,
                    cmap,
                    color_range,
                    filename):
  plot_hemispheres(lh, rh, array_name = data, nan_color = (0,0,0,1),size = size,
                   cmap = cmap, color_bar = True, color_range=color_range,
                   interactive = False, zoom = 1.5, embed_nb = True, transparent_bg=True,
                   cb__labelTextProperty={'fontSize': 24}, screenshot=True, filename=filename)
  fig = plot_hemispheres(lh, rh, array_name = data, nan_color = (0,0,0,1),size = size,
                         cmap = cmap, color_bar = True, color_range=color_range,
                         cb__labelTextProperty={'fontSize': 24}, interactive = False, zoom = 1.5, embed_nb = True)
  return fig
plot_surface_lr(lh, rh, data = plot,
                   size = (1200, 200), color_range=(-1,1),
                   cmap = 'BuPu', filename = 'mean.png') 