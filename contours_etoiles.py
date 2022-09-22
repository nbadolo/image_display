#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 14 16:44:26 2022

@author: nbadolo
"""
# =============================================================================
# realisation de contours
# =============================================================================
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as opt
import pylab as plt

#%%


import numpy as np

def twoD_GaussianWithTilt(xdata_tuple, amplitude, xo, yo, sigma_x, sigma_y, theta, offset, tiltX, tiltY):
    ''' twoD_GaussianWithTilt(xdata_tuple, amplitude, xo, yo, sigma_x, sigma_y, theta, offset, tiltX, tiltY) 
    '''
    (x,  y, scale) = xdata_tuple
    xo = float(xo)
    yo = float(yo)
    tiltX = float(tiltX)
    tiltY = float(tiltY)
    a = (np.cos(theta)**2)/(2*sigma_x**2) + (np.sin(theta)**2)/(2*sigma_y**2)
    b = -(np.sin(2*theta))/(4*sigma_x**2) + (np.sin(2*theta))/(4*sigma_y**2)
    c = (np.sin(theta)**2)/(2*sigma_x**2) + (np.cos(theta)**2)/(2*sigma_y**2)
    g = offset + amplitude*np.exp( - (a*((x-xo)**2) + 2*b*(x-xo)*(y-yo) 
                            + c*((y-yo)**2))) + tiltX*x + tiltY*y
    return g.ravel()

def twoD_Gaussian(xdata_tuple, amplitude, xo, yo, sigma_x, sigma_y, theta, offset):
    '''twoD_Gaussian(xdata_tuple, amplitude, xo, yo, sigma_x, sigma_y, theta, offset)
    '''
    (x,  y) = xdata_tuple
    xo = float(xo)
    yo = float(yo)
    a = (np.cos(theta)**2)/(2*sigma_x**2) + (np.sin(theta)**2)/(2*sigma_y**2)
    b = -(np.sin(2*theta))/(4*sigma_x**2) + (np.sin(2*theta))/(4*sigma_y**2)
    c = (np.sin(theta)**2)/(2*sigma_x**2) + (np.cos(theta)**2)/(2*sigma_y**2)
    g = offset + amplitude*np.exp( - (a*((x-xo)**2) + 2*b*(x-xo)*(y-yo) 
                            + c*((y-yo)**2)))
    return g.ravel()

def simple_twoD_Gaussian(xdata_tuple, amplitude, xo, yo, sigma_x, sigma_y, offset):
    '''twoD_Gaussian(xdata_tuple, amplitude, xo, yo, sigma_x, sigma_y, offset)
    '''
    (x,  y) = xdata_tuple
    xo = float(xo)
    yo = float(yo)
    a = 1/(2*sigma_x**2)
    b = 0
    c = 1/(2*sigma_y**2)
    g = offset + amplitude*np.exp( -(a*((x-xo)**2) + c*((y-yo)**2)))
    return g.ravel()


    # source https://stackoverflow.com/questions/21566379/fitting-a-2d-gaussian-function-using-scipy-optimize-curve-fit-valueerror-and-m