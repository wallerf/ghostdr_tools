from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
import argparse

#Create Argparse instance
parser = argparse.ArgumentParser()
#Require a filename
parser.add_argument('filename', help='2D GHOST image')
#Add keyword option for summing individual exposures
parser.add_argument('-s', '--sum',action='store_true', help='Sum all exposure into one image')
#Add keyword option for showing cuts along x/y
parser.add_argument('-c', '--cuts',action='store_true', help='Plot cuts along x/y axis')
#Add keyword option for showing red chip
parser.add_argument('-r', '--red',action='store_true', help='Plot red chip')
#Add keyword option for showing red chip
parser.add_argument('-b', '--blue',action='store_true', help='Plot blue chip')

#Parse the args
args = parser.parse_args()

#Load the file
hdul = fits.open(args.filename)

#Get total number of extensions
nExten = len(hdul) - 1

#Get the optional keywords
combine_exposures=args.sum
show_cuts=args.cuts

#Create which chips to use
chips=[]
if args.red: chips.append('RED')
if args.blue: chips.append('BLUE')

if len(chips)==0:
	raise UserWarning("Please run code with at least -r or -b")


#Extract the x/y start/end indicies
#Are flipped
def getReg(headstr):
	tmp= headstr.split(',')
	y1 = tmp[0].split(':')[0].split('[')[1]
	y2 = tmp[0].split(':')[1]
	x1 = tmp[1].split(':')[0]
	x2 = tmp[1].split(':')[1].split(']')[0]
	return int(x1), int(x2), int(y1), int(y2)

#Function to get all the relevant section indicies
def getSecs(header):
	Secs = {}
	for key in ['CCDSEC', 'BIASSEC', 'AMPSEC', 'TRIMSEC', 'DATASEC', 'DETSEC']:
		if key in header.keys():
			Secs[key.split('SEC')[0]] = getReg(header[key])
	return Secs


#Dictionary to store images
imgs = {}
#loop through all the extensions
for nn in range(nExten):
	exten = nn+1
	head = hdul[exten].header
	#Determine camera type
	type = head['CAMERA']
	if type in chips:

		#Create entry into imgs dic if necessary
		if type not in imgs.keys():
			imgs[type] = {}

		#Get the section indicies
		Secs = getSecs(head)

		#Process the raw image
		raw_img = hdul[exten].data
		if raw_img is not None:
			raw_img = hdul[exten].data.astype(float)
			mean_overscan = 0.0
			#Get the Bias, if able
			if 'BIAS' in Secs.keys():
				x1,x2,y1,y2 = Secs['BIAS']
				overscan_data = raw_img[x1:x2, y1:y2]
				mean_overscan = np.mean(np.mean(overscan_data,axis=0))
			#Subtract the bias
			raw_img -= mean_overscan
			#Trim the data
			if 'TRIM' in Secs.keys():
				x1,x2,y1,y2 = Secs['TRIM']
				tmp = raw_img[x1:x2,y1:y2]
				raw_img = tmp[:,:]
			#Store the data based on the coordinates of the AMP on the detector
			if Secs['DET'] not in imgs[type]:
				imgs[type][Secs['DET']] = []
			imgs[type][Secs['DET']].append(raw_img)

#Close the file
hdul.close()

#Loop through the camera types and plot the images
for type in chips:

	#Create list of images
	limg = []
	#Loop through the images, and combine the exposures (if necessary)
	#Sort the keys to have them in a specific order to combine amps
	for sec in sorted(imgs[type].keys()):
		sumimg = np.zeros(imgs[type][sec][0].shape)
		if combine_exposures:
			for rimg in imgs[type][sec]:
				sumimg += rimg
		else:
			sumimg = imgs[type][sec][0]
		limg.append(sumimg)

	#Combine the amps
	Aimg=np.concatenate((limg[0], limg[1]), axis=1)
	Bimg=np.concatenate((limg[2], limg[3]), axis=1)
	img = np.concatenate((Aimg, Bimg), axis=0)

	#Plot the 2D data
	fig, ax = plt.subplots(1,1)
	ax.imshow(img, vmin=0, vmax=np.nanpercentile(img,80), origin='lower')
	ax.set_title(type)

	#Plot the cuts, if required
	if show_cuts:
		#Take centre of image for xcut
		xcut = int(img.shape[1]/2)
		#Get an initial guess of the y cut
		ycut = int(img.shape[0]/2)

		#Loop through increasing indicies to find maximum value close to centre of image
		yval=0
		keepGoing=True
		ind=0
		ii=0
		while keepGoing:
			if img[ycut+ii, xcut]>yval:
				yval = img[ycut+ii, xcut]
				ind+=ii
			elif ind>0 or ii>ycut/8:
				keepGoing=False
			ii+=1
		ycut+=ind

		#Plot the cuts on the image
		xmin,xmax=ax.get_xlim()
		ymin,ymax=ax.get_ylim()
		ax.plot([xcut,xcut],[ymin,ymax],':k')
		ax.plot([xmin,xmax],[ycut,ycut],':k')
		ax.set_xlim(xmin,xmax)
		ax.set_ylim(ymin,ymax)

		#Plot cut along y-axis
		plt.figure()
		xs = np.arange(0,img.shape[0])
		plt.plot(xs, img[:,xcut])
		plt.title(type)
		plt.ylabel('Counts')
		plt.xlabel('y-axis pixel')

		#Plot cut along x-axis
		plt.figure()
		ys = np.arange(0,img.shape[1])
		plt.plot(img[ycut,:],ys)
		plt.title(type)
		plt.ylabel('Counts')
		plt.xlabel('x-axis pixel')


plt.show()


