""":module:: wham_make_file
   :platform: Unix, Windows
   :synopsis: reads the wham file

.. moduleauthor:: Alex Travesset <trvsst@ameslab.gov>, November 2016
"""


import numpy as np


def make_data(f_name, conv_unit = 1):
    """
    converts the data into a numpy file
    
    :param f_name: name of the file
    :param conv_unit: conversion unit
    :return: array with the data
    :rtype: numpy array
    """

    fcn = open(f_name, 'r')

    l_entry = fcn.readlines()

    data_file = []
    data = []

    for ln in l_entry:
        if float(ln) == -3.0:
            data_file.append(data)
            data = []
        else:
            data.append(float(ln))
            
    return conv_unit*np.array(data_file)
    
    
def make_histograms(data_file, delta_f):
    """ computes the function :math: '$H_i(R)$'
    
    :param data_file: data array
    :param n_size: number of bins per file
    :return: function $H_i(R)$'
    :rtype: list of numpy array
    """
    d_min = np.amin(data_file)
    d_max = np.amax(data_file)
    htol=20
    num = int((d_max-d_min)/delta_f)    
    
    bin_delt = np.linspace(d_min, d_max, num, endpoint = True)    
  
    h_fun = np.zeros([data_file.shape[0], bin_delt.shape[0]-1])
    
    for ind in range(h_fun.shape[0]):
        hist, bin_edges = np.histogram(data_file[ind], bins = bin_delt)
        h_fun[ind] = hist

	for i in range(len(h_fun[ind])):
        	if h_fun[ind][i] < htol:
			h_fun[ind][i] = 0.0

    return (h_fun, bin_edges)
