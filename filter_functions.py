import numpy as np
import scipy
import scipy.special
import re
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm


def S_pilatus(E, Ec): 
	''' Pilatus S threshold function.
	Formula: S = 1 - erfc(- \frac{E-E_c}{\sqrt(2) * E_w})
	'''
	Ew = 500 # Pilatus width
	S = 1 - 0.5*scipy.special.erfc(- (E -Ec)/(np.sqrt(2)* Ew))

	S = 1-S

	return S

def get_Si_transmission():
	''' Function to get the silicon transmission as a function of 
	energy. This data was obtained from 
	http://henke.lbl.gov/optical_constants/filter2.html
	using density 2.328 g/cm^3 and thickness=1000um
	'''
	Si_trans = open("Si_transmission_1000um.txt", "r")
	Si_trans_lines = Si_trans.readlines()
	data_Si = Si_trans_lines[2:]
	elements= [x.split() for x in data_Si]
	photon_energy = np.zeros(len(Si_trans_lines[2:]))
	Si_transmission = np.zeros(len(Si_trans_lines[2:]))
	for i in range(len(Si_trans_lines[2:])):
		photon_energy[i] = float(elements[i][0])
		Si_transmission[i] = float(elements[i][1])
	return (photon_energy, Si_transmission)

def get_Si_abs():
	''' Silicon absorption as a function of energy. 
	Simple operation on the transmission function above'''
	photon_energy, Si_transmission = get_Si_transmission()
	Si_absorption = 1 - Si_transmission
	return (photon_energy, Si_absorption)

def get_Be_transmission():
	''' Beryllium transmission as a function of energy. This was obtained 
	from http://henke.lbl.gov/optical_constants/filter2.html
	using a density 1.85 g/cm3 and thickness = 50 um. 
	'''
	Be_trans = open("Be_transmission_50um.txt", "r")
	Be_trans_lines = Be_trans.readlines()
	data_Be = Be_trans_lines[2:]
	elements= [x.split() for x in data_Be]
	photon_energy = np.zeros(len(Be_trans_lines[2:]))
	Be_transmission = np.zeros(len(Be_trans_lines[2:]))
	for i in range(len(Be_trans_lines[2:])):
		photon_energy[i] = float(elements[i][0])
		Be_transmission[i] = float(elements[i][1])  
    
	return (photon_energy, Be_transmission)


def get_Zmeter_abs():
	''' From Chen PPCF 2014 '''
    
	wavelength = [500,521.499,521.5,524.5,524.5001,600]
	absorption = [0,0,0.8,0.8,0,0]
    
	energy = (4.14E-15*3E8)/np.asarray(wavelength)  *1E9  # in eV
    
	absorption = [x for _,x in sorted(zip(energy,absorption))]
	energy = [x for _,x in sorted(zip(energy,energy))]
    

	return (energy, absorption)



def plot_Be_transmissivity():
	'''
	Convenience function for plotting. 
	'''
	photon_energy, Be_transmission = get_Be_transmission()
	fig=plt.figure()
	ax1=fig.add_subplot(111)
	ax1.plot(photon_energy,Be_transmission, 'r-')
	ax1.set_xlabel('Photon Energy [eV]')
	ax1.set_ylabel('Transmissivity')
	ax1.grid()
	plt.show()


def plot_Si_absorb():
	'''
	Convenience function for plotting.
	'''
	photon_energy, Si_absorption = get_Si_abs()
	fig=plt.figure()
	ax1=fig.add_subplot(111)
	ax1.plot(photon_energy,Si_absorption, 'r-')
	ax1.set_xlabel('Photon Energy [eV]')
	ax1.set_ylabel('Absoptivity')
	ax1.grid()
	plt.show()

def plot_pilatus_S(Ec=10000):
	'''
	Convenience function for plotting.
	Inputs: 
		Ec: threshold energy for the Pilatus 'filter' function.
		This can be a scalar or a list, e.g. Ec= [100,200,3000]
	'''
	Energy_range = np.linspace(100,25000,500)
	fig=plt.figure()
	ax1=fig.add_subplot(111)
	color=iter(cm.rainbow(np.linspace(0,1,len(Ec))))
	for E_in in Ec:
		c=next(color)
		S_pil = S_pilatus(Energy_range, E_in)
		ax1.plot(Energy_range, S_pil, c=c, label='E_c = %d eV' %(E_in))
	ax1.set_xlabel('Photon Energy [eV]')
	ax1.set_ylabel('S-curves')
	leg = ax1.legend()
	leg.draggable(True)
	ax1.grid()
	plt.show()