import astropy.io.fits as F
import numpy as np
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
from astropy.time import Time
import astropy.units as u

#Summary of code:
#This code finds pairs of stars that fall into given airmass/separation bins for testing dual acquistion, etc
#################


#Parameters of code:

#Bins of airmass; uples of [min,max)
AM_bins = [(1.0,1.25),(1.25,1.75),(1.75,2.25)]
#AM_bins = [(1.0,1.25),(1.25,1.75)]
#Bins of separation in arcmin [min,max)
#sep_bins = [(3.5,4.5), (6.0,7.0)]
sep_bins = [(1.5,2.5), (3.5,4.5), (6.0,7.0)]

#max magnitude of at least one star in the pair
maxMag = 10.0

applyBoth=True

#Minimum number of times during night the pair is visible
minNtime = 1


#Date/time information for getting AM

#Currently set up for getting in quarters of the night-ish
day = '2022-06-24'
#sunset = day+' 17:00:00'
#t1 = day+' 21:00:00'#First quarter
#t2 = day+' 00:00:00'#Second quarter
t3a = day+' 03:00:00'#Third quarter
t3 = day+' 04:00:00'#Third quarter
t3b = day+' 05:00:00'#Third quarter
t3c = day+' 06:00:00'#Third quarter
t3d = day+' 07:00:00'#Third quarter
sunrise = day+' 08:00:00'
#The following TIMES list is what is needed to compute airmass
times = [t3a, t3, t3b,t3c, t3d]
######################


#Fixed parameters

#Set location to Gemini South
GeminiSouth = EarthLocation.of_site('Gemini South')

#Set UTC offset from local time
utcoffset = 4.0*u.hour
###################

#Inner workings of code


#Function to get airmass at a given time
def getAM(coord, time, UT=False):
	time = Time(time)
	if not UT: time += utcoffset
	return coord.transform_to(AltAz(obstime=time,location=GeminiSouth)).secz


#Load pairs of data
data = F.open('Pairs.fits')[1].data


#Create an array to deterime the grid statistics
grid_stats = np.zeros((len(sep_bins),len(AM_bins)))


#Loop through pairs, and find which fit in the various bins


print("Possible pairs")
for ii, sep in enumerate(data['Separation']):
	#Check magnitude of one star in pair is sufficienlty bright
	magCheck=  data['Gaia_g_mag_A'][ii]<maxMag or data['Gaia_g_mag_B'][ii]<maxMag
	if applyBoth:
		magCheck=  data['Gaia_g_mag_A'][ii]<maxMag and data['Gaia_g_mag_B'][ii]<maxMag

	if magCheck:
		#make sure the separation fits in one of the bins
		for ss, (sepmin,sepmax) in enumerate(sep_bins):
			if sep>=sepmin and sep<sepmax:
				#Load the coordinates of the pair
				cooA = SkyCoord(ra=data['RA_A (J2000)'][ii], dec=data['Dec_A (J2000)'][ii], unit='degree')
				cooB = SkyCoord(ra=data['RA_B (J2000)'][ii], dec=data['Dec_B (J2000)'][ii], unit='degree')

				#loop through the times, and see if the pair falls in an AM bin
				Ntime =0
				atinds = []
				for tt,time in enumerate(times):
					amA = getAM(cooA, time, UT=False)
					amB = getAM(cooA, time, UT=False)
					AMmed = 0.5*(amA+amB)
					for aa, (AMmin, AMmax) in enumerate(AM_bins):
						if amA>=AMmin and amA<AMmax and amB>=AMmin and amB<AMmax:
							Ntime+=1
							atinds.append((aa,tt,AMmed))
				if Ntime> minNtime:
					for (aa,tt,AMmed) in atinds:
						AMmin,AMmax = AM_bins[aa]
						hr = times[tt].split(day)[1][:3]
						#This pair is in an AM bin and separation bin
						#print('\t',data['Name_A'][ii], data['Name_B'][ii], times[tt], 'AM:%.1f-%.1f'%(AMmin,AMmax), "Separation:%.2f'"%sep)
						print(data['Name_A'][ii], ';', data['RA_A (J2000)'][ii], ';', data['Dec_A (J2000)'][ii], ';', data['Gaia_g_mag_A'][ii], ';', data['BP-RP_A'][ii], ';', data['Type_A'][ii],'; Time-%shr, Sep-%.1f, AM-%.1f'%(hr, data['Separation'][ii], AMmed))
						print(data['Name_B'][ii], ';', data['RA_B (J2000)'][ii], ';', data['Dec_B (J2000)'][ii], ';', data['Gaia_g_mag_B'][ii], ';', data['BP-RP_B'][ii], ';', data['Type_B'][ii],'; Time-%shr, Sep-%.1f, AM-%.1f'%(hr, data['Separation'][ii], AMmed))
						grid_stats[ss,aa] +=1

#Examine statistics
print("Grid statistics")
for ss, (sepmin,sepmax) in enumerate(sep_bins):
	for aa, (AMmin, AMmax) in enumerate(AM_bins):
		n_match = grid_stats[ss,aa]
		print('\t','AM:%.1f-%.1f'%(AMmin,AMmax), "Sep:%.1f'-%.1f'"%(sepmin,sepmax), '%.0f matches'%n_match)



