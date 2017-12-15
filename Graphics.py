import numpy as np 
import matplotlib.pyplot as plt 
import matplotlib.colors as mcolors
import matplotlib.cm as cm
import re
import os

from integral_functions import *


#def GetDataSpectrum(fileSuperName,tempRange):
#    
#	filelist = []
#	for file in os.listdir(fileSuperName+'/Emission/collection'):
#		if file.endswith('spec'):
#			filelist.append(fileSuperName+'/Emission/collection/' + file)
#	data_list = []
#	temp_list = np.zeros(len(filelist))
##creating the colorbar
#	for i in range(len(filelist)):
#		file = filelist[i]
#		with open(file) as f:
#			line1 = f.readline()
#		temp_list[i] = float(re.search(r'(-?\d+\.\d+)',line1).group())/1000.0
#	nValues = temp_list
#	normalize = mcolors.Normalize(vmin=nValues.min(), vmax=nValues.max())
#	colormap = cm.jet
#
#
#
#
#
#	minTemps = abs(temp_list-tempRange)
#	print minTemps
#    
#	index1 = np.where(minTemps == min(minTemps))
#	print index1
#	#if len(index1)>1:
#	index1 = index1[1][0]
#	#else:
#	#	index1 = index1[0][0]
#	print('Closest temperature of {}'.format(temp_list[index1]))
#
#
#	for i in range(len(filelist)):
#		file = filelist[i]
#		with open(file) as f:
#			line1 = f.readline()
#		idx = 0
#	 	temp = float(re.search(r'(-?\d+\.\d+)',line1).group())/1000.0
#         
#	 	for j in range(len(temp_list)):
#	 		if temp_list[j] == temp:
#	 			idx = j
#	 #			print idx
#	 			break
#	 	data = np.loadtxt(file)
#	 	energy = data[:,0]
#	 	wavelength = (4.14E-15*3E8)/energy  *1E9  # in nm
#         
#	 	emissivity  = data[:,1]*4*np.pi*1.5112*10**26
#	#plt.semilogy(energy/1000,emissivity,label='T=%d[keV]'%(temp_list[i]/1000))
#
#	 	if temp == temp_list[index1]:
#			print len(emissivity)
#			break
#    
#        
#	return wavelength,emissivity   




def GetDataSpectrum(impurity,tempRange):
    
    rho,energy2D,emissivity2D,Te2D = loadDataNPZ(impurity,True)


    # Get the temperature I want
    minTemps = abs(Te2D-tempRange)
    index1 = np.where(minTemps == min(minTemps))
    index1 = index1[0][0]
    print('Closest temperature of {}'.format(Te2D[index1]))


    energy = energy2D[index1]
    wavelength = (4.14E-15*3E8)/energy  *1E9  # in nm       
    emissivity  = emissivity2D[index1]
    
    emissivity = [x for _,x in sorted(zip(wavelength,emissivity))]
    wavelength = [x for _,x in sorted(zip(wavelength,wavelength))]
    
    return wavelength,emissivity   




def Spectrum(AllTemps):
    
	f, ax0 = plt.subplots(1)
    
	NewWavelength = np.linspace(400,600,1000)


	for tempRange in AllTemps:

		wavelength,emissivity = GetDataSpectrum('D',tempRange)
		NewEmissivity = np.interp(NewWavelength,wavelength,emissivity)
        
		Totalemissivity = NewEmissivity
		for i in ['C','N','Mo']:
        
			wavelength,emissivity = GetDataSpectrum(i,tempRange)
			NewEmissivity = np.interp(NewWavelength,wavelength,emissivity)
    		    #plt.plot(wavelength,emissivity)
			Totalemissivity =Totalemissivity + NewEmissivity


		Totalemissivity = (Totalemissivity - Totalemissivity[0])/Totalemissivity[0]

		plt.plot(NewWavelength,Totalemissivity,linewidth = 2,label='Te = {}keV'.format(tempRange))
        
	ax0.set_xlabel('Wavelength [nm]')
	ax0.set_ylabel('Relative $photons/cm^3/s$')
	plt.show()

	#plt.ylim([0,5.0E14])
	ax0.set_ylim([-0.01,0.23])
    
	ax0.legend(loc='best').draggable()
    
    
    
	axi = ax0.twinx()
    
	x1 = [500,521.499,521.5,524.5,524.5001,600]
	y1 = [0,0,0.8,0.8,0,0]
    
	axi.plot(x1,y1,linestyle='--',color='k')
	axi.set_ylim([0,1.0])
	axi.set_ylabel(r'Filter function')
    
    
    
	ax0.set_xlim([515,535])
    
    

def ChargeState():
	filelist = []
	for file in os.listdir('.'):
		if file.endswith('spec'):
			filelist.append(file)
	data_list = []
	temp_list = np.zeros(len(filelist))

#creating the colorbar
	for i in range(len(filelist)):
		file = filelist[i]
		with open(file) as f:
			line1 = f.readline()
		temp_list[i] = float(re.search(r'(-?\d+\.\d+)',line1).group())/1000.0
	nValues = temp_list
	normalize = mcolors.Normalize(vmin=nValues.min(), vmax=nValues.max())
	colormap = cm.jet

	plt.figure()
	for i in range(len(filelist)):
		file = filelist[i]
		with open(file) as f:
			line1 = f.readline()
		idx = 0
	 	temp = float(re.search(r'(-?\d+\.\d+)',line1).group())/1000.0
	 	for j in range(len(temp_list)):
	 		if temp_list[j] == temp:
	 			idx = j
	 #			print idx
	 			break
	 	data = np.loadtxt(file)
	 	charge = data[:,0]
	 	population = data[:,1]
	# 	plt.plot(charge,population,label='T=%d[keV]'%(temp_list[i]/1000))
	 	plt.plot(charge,population,linewidth = 0.5,color = colormap(int(normalize(nValues[idx])*colormap.N)))
	scalarmappaple = cm.ScalarMappable(norm=normalize, cmap=colormap)
	scalarmappaple.set_array(nValues)
	clb = plt.colorbar(scalarmappaple)
	clb.ax.set_title(r'$T_e[keV]$')
	plt.xlabel('Charge State')
	plt.ylabel('Relative ion population')
	plt.show()


def Spectrum_star(element,nz_ne):
#Loading files
	filelist_Emissivity = []
	for file in os.listdir('.'):
		if file.startswith('Emissivity'):
			filelist_Emissivity.append(file)
	filelist_ChargeState = []
	for file in os.listdir('.'):
		if file.startswith('ChargeState'):
			filelist_ChargeState.append(file)


	temp_list = np.zeros(len(filelist_Emissivity))
#creating the colorbar
	for i in range(len(filelist_Emissivity)):
		file = filelist_Emissivity[i]
		with open(file) as f:
			line1 = f.readline()
		temp_list[i] = float(re.search(r'(-?\d+\.\d+)',line1).group())/1000.0
	nValues = temp_list
	normalize = mcolors.Normalize(vmin=nValues.min(), vmax=nValues.max())
	colormap = cm.jet

#Calculatint the average Z for each temperature
	Zavg_list = np.zeros(len(temp_list))
	for i in range(len(temp_list)):
		Te = temp_list[i]
		diff = 10000
		idx = 0
		for j in range(len(filelist_ChargeState)):
			file = filelist_ChargeState[j]
			with open(file) as f:
				line1 = f.readline()
			temp = float(re.search(r'(-?\d+\.\d+)',line1).group())/1000.0
			tmp = abs(temp/1000.0-temp_list[i])
			if tmp<diff:
				idx = j
		data = np.loadtxt(filelist_ChargeState[idx])
	 	charge = data[:,0]
	 	population = data[:,1]
	 	Zavg_list[i] = np.average(charge,weights = population)
	 
#Loading data given by Alex
	profiledata = np.loadtxt('/Users/munizhou/Dropbox (MIT)/SXR Info/FLYCHK/profile1dAll.txt',skiprows = 1)
	rhodata = profiledata[:,0]
	Tedata = profiledata[:,1]
	nedata = profiledata[:,3]
#Seraching for the corresponding ne and charge state to get the average Z
	ne_list = np.zeros(len(temp_list))
	rho_list = np.zeros(len(temp_list))
	for i in range(len(temp_list)):
		Te = temp_list[i]
		diff = 100000
		for j in range(len(Tedata)):
			tmp = abs(Tedata[j]/1000.0-temp_list[i])
			if tmp <diff:
				diff = tmp
				ne_list[i] = nedata[j]
				rho_list[i] = rhodata[j]
	ne_list = ne_list/10**20

	PhotonEnergy_list = []
	Emissivity_list = []
	plt.figure()
	for i in range(len(filelist_Emissivity)):
		file = filelist_Emissivity[i]
		with open(file) as f:
			line1 = f.readline()
		idx = 0
	 	temp = float(re.search(r'(-?\d+\.\d+)',line1).group())/1000.0
	 	for j in range(len(temp_list)):
	 		if temp_list[j] == temp:
	 			idx = j
	 #			print idx
	 			break
	 	data = np.loadtxt(file)
	 	energy = data[:,0]
	 	emissivity  = data[:,1]*4*np.pi*1.5112*10**26*ne_list[idx]*Zavg_list[idx]*nz_ne
	 	PhotonEnergy_list.append(energy)
	 	Emissivity_list.append(emissivity)
		plt.semilogy(energy/1000,emissivity,linewidth = 0.5,color = colormap(int(normalize(nValues[idx])*colormap.N)))
	scalarmappaple = cm.ScalarMappable(norm=normalize, cmap=colormap)
	scalarmappaple.set_array(nValues)
	clb = plt.colorbar(scalarmappaple)
	clb.ax.set_title(r'$T_e[keV]$')
	plt.xlabel('Photon Energy [keV]')
	plt.ylabel(r'Emissivity $photons/cm^3/s$')
	plt.ylim(ymin=10**8)
	plt.title(element)
	plt.show()
	PhotonEnergy_list = np.array(PhotonEnergy_list)
	Emissivity_list = np.array(Emissivity_list)
	outputfile = element+'_Emissivity'+'.npz'
	np.savez(outputfile,Te_list = temp_list,rho_list=rho_list,ne_list=ne_list,Zavg_list = Zavg_list,PhotonEnergy_list = PhotonEnergy_list,Emissivity_list=Emissivity_list)

