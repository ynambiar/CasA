#coding: utf-8

#OBSERVATION-WISE (spatial anaylsis)
#Written by Yamini Nambiar, 1/6/2013
from astropy.io import fits as pyfits
from numpy import *
import math

import numpy as np
file_path_prefix = '/home/yamini/Data/'



###Constructs fully qualified filename of flux image files to read
def construct_flux_filename(observation, ratio, bin, flux):
	if bin == 4 or bin == 8:
		flux_filename = file_path_prefix + str(observation) + '/repro/f00' + str(bin) + '/' + flux + '_flux.img'
	else:
		flux_filename = file_path_prefix + str(observation) + '/repro/f0' + str(bin) + '/' + flux + '_flux.img'
		return flux_filename

###Constructs fully qualified filename of original image files to read
def construct_img_filename(observation, ratio, bin, flux):
	if bin == 4 or bin == 8:
		flux_filename = file_path_prefix + str(observation) + '/repro/f00' + str(bin) + '/' + flux + '.img'
	else:
		flux_filename = file_path_prefix + str(observation) + '/repro/f0' + str(bin) + '/' + flux + '.img'
	return flux_filename


###Sums soft flux, medium flux, and hard flux into broad flux image
def create_broad(observation, ratio, bin):
	broad_flux = soft_flux_image + medium_flux_image + hard_flux_image
	return broad_flux


###This function prints the full & complete ratio filenames so that we can see it
def print_filename(ratio, bin):
	a = construct_ratio_filename(observation,ratio, bin)
	print a 	

###Reads in a single ratio image
def read_single_img(fqName):
	hdulist = pyfits.open(fqName)
	temp_img = zeros((rows_global,columns_global), float64)
	for i in range(0,temp_img.shape[0]):
		for j in range(0, temp_img.shape[1]):
			temp_img[i,j] = hdulist[0].data[i,j]
	return temp_img

###MASK Calculates mean of given pixel's flux values over all observations
def calc_flux_mean(row, col, ratio, bin):
	mtrx_flux_mean = zeros((1, len(obsid)), float64)
	for i in range (0, len(obsid)):
		#if g_broad_flux_img_list[i][row,col] > 0:
		mtrx_flux_mean[0,i] = g_broad_flux_img_list[i][row,col]
	flux_mean = mean(mtrx_flux_mean, axis=1, dtype=float64)
	return flux_mean


###MASK Creates mask values - filters through values in matrix and determines those that are meaningful for analysis
def construct_mask(total_rows, total_cols, ratio, bin):
	for row in range(0, total_rows):
		for col in range(0, total_cols):
			flux_mean = calc_flux_mean(row, col, ratio, bin)
		#	if flux_mean > mask_bin_values[bin]:		
			if flux_mean > mask_bin_values[bin] and check_band_count_sum(row, col, ratio, bin) == 1 and check_obsid_count_sum(row, col, ratio, bin) == 1:		
		#	if check_obsid_count_sum(row, col, ratio, bin)== 1 and check_band_count_sum(row, col, ratio, bin) == 1:		
				masked_rows.append(row)	
				masked_cols.append(col)
	return

###MASK Read in flux images into appropriate lists
def read_in_images_to_list(ratio, bin):
	for observation in obsid:
		#READ IN SOFT, MEDIUM, HARD FLUX IMAGES
		fq_soft_flux_filename = construct_flux_filename(observation, ratio, bin, bands[0])
		soft_flux_image = read_single_img(fq_soft_flux_filename)
		g_soft_flux_img_list.append(soft_flux_image)
		
		fq_medium_flux_filename = construct_flux_filename(observation, ratio, bin,  bands[1])
		medium_flux_image = read_single_img(fq_medium_flux_filename)
		g_medium_flux_img_list.append(medium_flux_image)

		fq_hard_flux_filename = construct_flux_filename(observation, ratio, bin, bands[2])
		hard_flux_image = read_single_img(fq_hard_flux_filename)
		g_hard_flux_img_list.append(hard_flux_image)
	
		broad_flux_image  = soft_flux_image + medium_flux_image + hard_flux_image
		g_broad_flux_img_list.append(broad_flux_image)
	
		soft = construct_img_filename(observation, ratio, bin, bands[0])
		g_soft_img_list.append(read_single_img(soft))

		medium = construct_img_filename(observation, ratio, bin, bands[1])
		g_medium_img_list.append(read_single_img(medium))
		
		hard = construct_img_filename(observation, ratio, bin, bands[2])
		g_hard_img_list.append(read_single_img(hard))	

		broad = read_single_img(soft) + read_single_img(medium) + read_single_img(hard)
		g_broad_img_list.append(broad)
		

		if ratio == ratios[0]:
			ratio_image = soft_flux_image/medium_flux_image
			g_ratio_img_list.append(ratio_image)
		else:		
			ratio_image = medium_flux_image/hard_flux_image
			g_ratio_img_list.append(ratio_image)

	return

###MASK Returns true if there is at least 1 count for a pixel across all observations 
def check_band_count_sum(row, col, ratio, bin):
	for i in range (0, len(obsid)):
		soft =	g_soft_img_list[i][row, col]
		medium = g_medium_img_list[i][row, col]
		hard = g_hard_img_list[i][row, col]
		if (soft + medium + hard) < 1:
			return 0
	return 1

###MASK Returns true if there are at least 200 counts for a pixel across all observations,
def check_obsid_count_sum(row, col, ratio, bin):	
	broad_list = []
	for i in range (0, len(obsid)):
		broad_value= g_broad_img_list[i][row,col]
		broad_list.append(broad_value)
	count_sum = sum (broad_list)
	if count_sum > 200:
		return 1
	else:
		return 0

###Calculates matrix of log (base 10) values of given lists long as the value does not equal zero
def construct_log_list(mask_list):
	log_list = []
	for i in range(0,len (mask_list)):
			if mask_list[i] != 0:
				log_list.append(log10(mask_list[i]))
			else:
				print '****MASK ERROR: VALUE 0****'
	return log_list

###Calculates mean of given list
def calc_mean(mask_list):
	mtrx = zeros((1, len(mask_list)), float64)
	for i in range(0, len(mask_list)):
		mtrx[0,i] = mask_list[i]
	mtrx_mean = mean(mtrx, axis=1, dtype=float64)
	return mtrx_mean

###Calculates standard deviation of given list
def calc_stddev(mask_list):
	mtrx_stddev = zeros ((1, len (mask_list)), float64)
	for i in range (0, len (mask_list)):
		mtrx_stddev[0,i] = mask_list[i]
	stddev = std(mtrx_stddev, axis = 1, dtype=float64) #to take the mean of each row (there’s only one)
	return stddev

###Renormalizes values of a given list (of logged values) by dividing each value by the given mean 
def construct_renorm_list(mask_list, mean):
	renorm_list = []
	for i in range (0, len (mask_list)):
		renorm_list.append(mask_list[i] - mean)
	return renorm_list
		
###Returns list of most extreme masks values 
def flag_vals(values_list, rows_list, cols_list, stddev, mean):
	for i in range(0, len (values_list)):
		if values_list[i] > (mean+3*stddev) or values_list[i] < (mean - 3*stddev):
			#flagged_vals.append(values_list[i])
			flagged_rows.append(rows_list[i])
			flagged_cols.append(cols_list[i])	
	return


###Transfers flagged values into an array corresponding to correct coordinates 
def transfer_vals_array(rows_list, cols_list):
	flagged_values_array = zeros([rows_global, columns_global], float64)
	for i in range(0, len (rows_list)):
		flagged_vals_array[rows_list[i], cols_list[i]] = 1
	return flagged_values_array

###Draws mask
def draw_mask(data, ratio, bin):
	import Image
	rescaled = (255.0 / data.max() * (data - data.min())).astype(np.uint8)
	im = Image.fromarray(rescaled)
	im.save("mask_" + ratio + str(bin) + "_observation.png")

###Draws flagged regions
def draw_flagged_regions(data, ratio, bin):
	import Image
	rescaled = (255.0 / data.max() * (data - data.min())).astype(np.uint8)
	im = Image.fromarray(rescaled)
	im.save("regions_" + ratio + str(bin) + "_observation.png")

		

			
#main
import sys 

if len(sys.argv) < 5:
	print "Four parameters are required: a bin size, a ratio, number of rows, number of columns"
	exit (1)

bin_global=sys.argv[1]
ratio_global=sys.argv[2]
rows_global = int(sys.argv[3])
columns_global = int(sys.argv[4])

mask_array = zeros([rows_global, columns_global], float64)
flagged_vals_array = zeros([rows_global, columns_global], float64)

soft_flux_image = zeros([rows_global, columns_global], float64)
medium_flux_image = zeros([rows_global, columns_global], float64)
hard_flux_image = zeros([rows_global, columns_global], float64)


obsid = [198,1952, 5196, 6745, 9117,10935, 10936, 14229]
bands = ['soft', 'medium', 'hard']
ratios = ['one', 'two']


base = 5*(10**-6)
#mask_bin_values = {4:5*(10**-6), 8:5*(10**-6)*4, 16:5*(10**-6)*16, 32:5*(10**-6)*64, 64:5*(10**-6)*256}
mask_bin_values = {'4':base, '8':base*4, '16':base*16, '32':base*64, '64':base*256}


#global lists, each element is an image
g_soft_flux_img_list = []
g_medium_flux_img_list = []
g_hard_flux_img_list = []
g_broad_flux_img_list = []

g_soft_img_list = []
g_medium_img_list = []
g_hard_img_list = []
g_broad_img_list = []

g_ratio_img_list = []




#CONSTRUCT MASK
masked_rows = []
masked_cols = []

read_in_images_to_list(ratio_global, bin_global)
construct_mask(rows_global, columns_global, ratio_global, bin_global)


#initialize to ensure array is float64 array	
temp_mask_array = zeros([rows_global, columns_global], float64)
temp_mask_array = transfer_vals_array(masked_rows, masked_cols)
mask_array = temp_mask_array + mask_array
draw_mask(mask_array, ratio_global, bin_global)



#OBSERVATION-WISE
obsid_index = 0
for observation in obsid:
	print 'Calculating OBSERVATION-WISE analysis for observation ' + str(observation) + ' ...'

	#MASK	
	masked_vals = []
	for k in range (0, len(masked_rows)):
		masked_vals.append(g_ratio_img_list[obsid_index][masked_rows[k], masked_cols[k]])


	print 'masked_vals'
	print len( masked_vals)
	

	#LOG & RENORMALIZE
	log_list = construct_log_list(masked_vals) #returns list of logged values of masked values 
	log_mean = calc_mean(log_list) #using logged values, finds mean
	log_stddev = calc_stddev(log_list) #using logged values, finds stddev
	
	renorm_list = construct_renorm_list(log_list,log_mean) #using mean of logged values, renormalizes
	renorm_mean = calc_mean(renorm_list)


	#ANALYSIS
	flagged_vals = []
	flagged_rows = []
	flagged_cols = []
	
	flag_vals(renorm_list, masked_rows, masked_cols, log_stddev, renorm_mean)
	
	flagged_vals = []
	for k in range (0, len(flagged_rows)):
		flagged_vals.append(g_ratio_img_list[obsid_index][flagged_rows[k], flagged_cols[k]])

		#10^[log10(mean) +log10(x/mean)]
		#value = 10**(renorm_list[k] + log_mean)
		#flagged_vals.append(value)



	print 'flagged_vals'
	print flagged_vals
	print 'flagged_rows'
	print flagged_rows
	print 'flagged_cols'
	print flagged_cols


	print '-----------------------------------------------'
	

	temp_flagged_vals_array = transfer_vals_array(flagged_rows, flagged_cols)
	flagged_vals_array = flagged_vals_array + temp_flagged_vals_array

	obsid_index = obsid_index + 1	

draw_mask(mask_array, ratio_global, bin_global)
draw_flagged_regions(flagged_vals_array, ratio_global, bin_global)

mask_array = zeros([rows_global, columns_global], float64)
flagged_vals_array = zeros([rows_global, columns_global], float64)



