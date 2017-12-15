import numpy as np 
import matplotlib.pyplot as plt 
import matplotlib.colors as mcolors
import matplotlib.cm as cm
import re
import os
from matplotlib.pyplot import ion
ion()
import sys
sys.path.insert(0,'C:\Users\sciortino\Dropbox (MIT)\SXR Info\Integration')
sys.path.insert(0,'C:\Users\sciortino\Dropbox (MIT)\SXR Info\FLYCHK')
sys.path.insert(0,'C:\Users\sciortino\Dropbox (MIT)\SXR Info\FLYCHK\Mo_200_25000\Emission\collection')
from filter_functions import *

def Spectrum():
	filelist = []
	for file in os.listdir('.'):
		if file.endswith('.spec'):
			filelist.append(file)

	# If reading doesn't occur correctly, create file list alternatively
	if len(filelist) == 0:
		filelist = [os.path.join(root, name) 
			for root, dirs, files in os.walk('./')
			for name in files
			if name.endswith('.spec')]

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
	 	energy = data[:,0]
	 	emissivity  = data[:,1]*4*np.pi*1.5112*10**26
	#plt.semilogy(energy/1000,emissivity,label='T=%d[keV]'%(temp_list[i]/1000))
		plt.semilogy(energy/1000,emissivity,linewidth = 0.5,color = colormap(int(normalize(nValues[idx])*colormap.N)))
	scalarmappaple = cm.ScalarMappable(norm=normalize, cmap=colormap)
	scalarmappaple.set_array(nValues)
	clb = plt.colorbar(scalarmappaple)
	clb.ax.set_title(r'$T_e[keV]$')
	plt.xlabel('Photon Energy [keV]')
	plt.ylabel(r'Emissivity $photons/cm^3/s$')
	plt.ylim(ymin=10**8)
	plt.show()



def ChargeState():
	filelist = []
	for file in os.listdir('.'):
		if file.endswith('spec'):
			filelist.append(file)

	# If reading doesn't occur correctly, create file list alternatively
	if len(filelist) == 0:
		filelist = [os.path.join(root, name) 
			for root, dirs, files in os.walk('./')
			for name in files
			if name.endswith('spec')]

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

	# If reading doesn't occur correctly, create file list alternatively
	if len(filelist_Emissivity) <= 1:
		filelist_Emissivity = [os.path.join(root, name) 
			for root, dirs, files in os.walk('./')
			for name in files
			if name.startswith('Emissivity')]
	if len(filelist_ChargeState) <= 1:
		filelist_ChargeState = [os.path.join(root, name) 
			for root, dirs, files in os.walk('./')
			for name in files
			if name.startswith('ChargeState')]

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

	#Calculating the average Z for each temperature
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
	#profiledata = np.loadtxt('/Users/munizhou/Dropbox (MIT)/SXR Info/FLYCHK/profile1dAll.txt',skiprows = 1)
	profiledata = np.loadtxt('C:/Users/sciortino/Dropbox (MIT)/SXR Info/FLYCHK/profile1dAll.txt',skiprows = 1)

	rhodata = profiledata[:,0]
	Tedata = profiledata[:,1]
	nedata = profiledata[:,3]

	#Searching for the corresponding ne and charge state to get the average Z
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
	 			break
	 	data = np.loadtxt(file)
	 	energy = data[:,0]
	 	emissivity  = data[:,1]*4*np.pi*1.5112*10**26*ne_list[idx]*Zavg_list[idx]*nz_ne
	 	PhotonEnergy_list.append(energy)
	 	Emissivity_list.append(emissivity)
		plt.semilogy(energy/1000,emissivity,linewidth = 0.5,color = colormap(int(normalize(nValues[idx])*colormap.N)))

	scalarmappaple = cm.ScalarMappable(norm=normalize, cmap=colormap)
	scalarmappaple.set_array(nValues)
	clb = plt.colorbar(scalarmappaple,pad=0.15)

	clb.ax.set_title(r'$T_e[keV]$',fontsize=16)
	plt.xticks(fontsize=16)
	plt.yticks(fontsize=16)
	plt.xlabel('Photon Energy [keV]',fontsize=16)
	plt.ylabel(r'Emissivity $photons/cm^3/s$',fontsize=16)
	plt.ylim(ymin=10**8)
	plt.title(element)

	# Add Pilatus S-curves to the spectrum:
	# ax0 = plt.gca()
	# ax1 = ax0.twinx()
	# #plt.figure()

	# #Ec_range = np.linspace(6,24,10) *1000
	# Ec_range = [20000,22000]
	# for Ec in Ec_range:
	# 	ax1.plot(energy/1000, S_pilatus(energy, Ec), 'r--')
	# ax1.set_ylabel(r'Pixel Transmission',fontsize=16)

	plt.show()
	PhotonEnergy_list = np.array(PhotonEnergy_list)
	Emissivity_list = np.array(Emissivity_list)
	outputfile = element+'_Emissivity'+'.npz'
	np.savez(outputfile,Te_list = temp_list,rho_list=rho_list,ne_list=ne_list,Zavg_list = Zavg_list,PhotonEnergy_list = PhotonEnergy_list,Emissivity_list=Emissivity_list)


class DraggableColorbar(object):
    def __init__(self, cbar, mappable):
        self.cbar = cbar
        self.mappable = mappable
        self.press = None
        self.cycle = sorted([i for i in dir(plt.cm) if hasattr(getattr(plt.cm,i),'N')])
        self.index = self.cycle.index(cbar.get_cmap().name)

    def connect(self):
        """connect to all the events we need"""
        self.cidpress = self.cbar.patch.figure.canvas.mpl_connect(
            'button_press_event', self.on_press)
        self.cidrelease = self.cbar.patch.figure.canvas.mpl_connect(
            'button_release_event', self.on_release)
        self.cidmotion = self.cbar.patch.figure.canvas.mpl_connect(
            'motion_notify_event', self.on_motion)
        self.keypress = self.cbar.patch.figure.canvas.mpl_connect(
            'key_press_event', self.key_press)

    def on_press(self, event):
        """on button press we will see if the mouse is over us and store some data"""
        if event.inaxes != self.cbar.ax: return
        self.press = event.x, event.y

    def key_press(self, event):
        if event.key=='down':
            self.index += 1
        elif event.key=='up':
            self.index -= 1
        if self.index<0:
            self.index = len(self.cycle)
        elif self.index>=len(self.cycle):
            self.index = 0
        cmap = self.cycle[self.index]
        self.cbar.set_cmap(cmap)
        self.cbar.draw_all()
        self.mappable.set_cmap(cmap)
        self.mappable.get_axes().set_title(cmap)
        self.cbar.patch.figure.canvas.draw()

    def on_motion(self, event):
        'on motion we will move the rect if the mouse is over us'
        if self.press is None: return
        if event.inaxes != self.cbar.ax: return
        xprev, yprev = self.press
        dx = event.x - xprev
        dy = event.y - yprev
        self.press = event.x,event.y
        #print 'x0=%f, xpress=%f, event.xdata=%f, dx=%f, x0+dx=%f'%(x0, xpress, event.xdata, dx, x0+dx)
        scale = self.cbar.norm.vmax - self.cbar.norm.vmin
        perc = 0.03
        if event.button==1:
            self.cbar.norm.vmin -= (perc*scale)*np.sign(dy)
            self.cbar.norm.vmax -= (perc*scale)*np.sign(dy)
        elif event.button==3:
            self.cbar.norm.vmin -= (perc*scale)*np.sign(dy)
            self.cbar.norm.vmax += (perc*scale)*np.sign(dy)
        self.cbar.draw_all()
        self.mappable.set_norm(self.cbar.norm)
        self.cbar.patch.figure.canvas.draw()


    def on_release(self, event):
        """on release we reset the press data"""
        self.press = None
        self.mappable.set_norm(self.cbar.norm)
        self.cbar.patch.figure.canvas.draw()

    def disconnect(self):
        """disconnect all the stored connection ids"""
        self.cbar.patch.figure.canvas.mpl_disconnect(self.cidpress)
        self.cbar.patch.figure.canvas.mpl_disconnect(self.cidrelease)
        self.cbar.patch.figure.canvas.mpl_disconnect(self.cidmotion)