import numpy as np 
import matplotlib.pyplot as plt 
import matplotlib.colors as mcolors
import matplotlib.cm as cm
import scipy
import scipy.special
import re
import os
from filter_functions import *
import copy

def BellInversion(rho, n):
	'''
	Abel-like inversion for tokamak plasma geometry, following 
	R.E. Bell, "Inversion technique to obtain an emissivity profile 
	from tangential line-integrated hard x-ray measurements", RSI, 66, 558, 1995
	
	Parameters:
		rho: numpy.ndarray, e.g.(100x1)
			Array containing the vector representing rho_psi

		n: int
			Parameter referring to the cosine power to be assumed 
			for poloidally-asymmetric inversions. Set to 0 to assume 
			poloidal symmetry.

	Output: 
		M: numpy.ndarray
			Bell M matrix, see p.4 of Bell's paper.

	'''
	# First define a radius for each emission zone:
	RmajRho = getRealSpacePerpendicularCoord(rho)
	r_HFS = np.zeros(len(RmajRho)-1); 
	r_HFS[0] = RmajRho[0] - np.diff(RmajRho)[0]
	for i in np.arange(1,len(RmajRho)-1):
		r_HFS[i] = r_HFS[i-1] - np.diff(RmajRho)[i]

	r_T = np.concatenate((r_HFS[::-1],RmajRho))

	# Define Bell's "s" variables:
	s = np.zeros(len(r_T)-1)
	for j in range(len(r_T)-1):
		s[j] = (r_T[j]+r_T[j+1]) / 2.0

	# Define the length matrix L:
	L_Bell = np.zeros((len(r_T),len(s)))
	cos_theta = np.zeros((len(r_T),len(r_T)))
	for i in range(len(r_T)):
		for j in range(len(s)):
			if s[j]> r_T[i]:
				if i == j:
					L_Bell[i,j] = 2 * np.sqrt(s[j]**2 - r_T[i]**2)
				else:
					L_Bell[i,j] = 2 * np.sqrt(s[j]**2 - r_T[i]**2) - 2* np.sqrt(s[j-1]**2 - r_T[i]**2)

	cos_theta = [r_T[i]/r_T[j] for i in range(len(r_T)) for j in range(len(r_T)) ]
	cos_theta = np.reshape(cos_theta,(len(r_T),len(r_T)))

	# Define Bell M matrix:
	M_Bell = np.zeros((len(r_T),len(s)))
	for i in range(L_Bell.shape[0]):
		for j in range(L_Bell.shape[1]): 
			M_Bell[i,j] = L_Bell[i,j] * cos_theta[i,j]**n

	return M_Bell
	
def swapRhoAndY(rho,DetectedNew,Ec):

	rho_sorted = [x for _,x in sorted(zip(rho,rho))]

	if not Ec:
		DetectedNew_sorted = [x for _,x in sorted(zip(rho,DetectedNew))]
	else:
		DetectedNew_sorted = []
		for i_Ec in range(len(Ec)):
			DetectedNew_sorted.append([x for _,x in sorted(zip(rho,DetectedNew[i_Ec]))])

	return rho_sorted,DetectedNew_sorted

def arrangeRho(rho,energy2D,emissivity2D):

	# # Sort rho
	# rho_sorted = [x for _,x in sorted(zip(rho,rho))]

	# # Switch

	# energy2D = np.swapaxes(energy2D,0,1)
	# energy2D_sorted = []
	# for i in range(len(energy2D)):
	# 	energy2D_sorted.append([x for _,x in sorted(zip(rho,energy2D[i]))])
	# energy2D_sorted = np.swapaxes(energy2D_sorted,0,1)

	# rho = rho_sorted
	# energy2D = energy2D_sorted
	# emissivity2D = emissivity2D_sorted

	return rho,energy2D,emissivity2D

def loadDataNPZ(Impurity,LowYN):

	#Impurity e.g. 'Mo'

	if LowYN:
		fileName = '../FLYCHK/PlotsData/{}_Emissivity_0.1_20.npz'.format(Impurity)
	else:
		fileName = '../FLYCHK/PlotsData/{}_Emissivity.npz'.format(Impurity)

	print 'Reading data from {}'.format(fileName)

	aux = np.load(fileName)

	rho = aux['rho_list']
	energy2D = aux['PhotonEnergy_list']   #rhoxE
	emissivity2D = aux['Emissivity_list'] #rhoxE
	Te2D = aux['Te_list']

	rho,energy2D,emissivity2D = arrangeRho(rho,energy2D,emissivity2D)

	return rho,energy2D,emissivity2D,Te2D


def computeIntegral(EmissivityRadius,EnergyGrid,Ec,Bet_new,Sia_new,Zmetera_new,lowYN,PixelEnhancement,MSEnhancement):

	# Compute the photons absobed from that radius and divide by E for the dE integral

	if lowYN:
		integ = EmissivityRadius / EnergyGrid * Zmetera_new
	else:
		integ = EmissivityRadius / EnergyGrid * Bet_new * Sia_new

	# Multiply by pixel enhancement (more than 1 pixel)
	integ  = integ *PixelEnhancement 

	# Multiply by ms enhancement (more than 1 ms)
	integ  = integ *MSEnhancement


#	import pdb
#	pdb.set_trace()
     

	# Proceed with the cut-off S functions -> one integral per Ec

	Detected_Ec = []
	for i_Ec in range(len(Ec)):

		# Multiply by threshold and integrate

		integS = integ*S_pilatus(EnergyGrid, Ec[i_Ec]*1000)


		# Integrate

		integrResult = np.trapz(integS,x=EnergyGrid)
		Detected_Ec.append(integrResult)

	return Detected_Ec


def computeAllRadius(rho,newRho,energy2D,Ec,emissivity2D,Bet,Sia,en_Bet,en_Sia,en_Zmeter, Zmetera,lowYN,PixelEnhancement,MSEnhancement):

	Detected = []

	for i_rho in range(len(rho)):

		# Everything in the same Energy Grid (interpolate the energy grid for Berilium and Silicon to the new grid)

		EnergyGrid = energy2D[i_rho]
		EmissivityRadius = emissivity2D[i_rho]
		Bet_new = np.interp(EnergyGrid,en_Bet,Bet)
		Sia_new = np.interp(EnergyGrid,en_Sia,Sia)
		Zmetera_new = np.interp(EnergyGrid,en_Zmeter,Zmetera)
                
		Detected_Ec = computeIntegral(EmissivityRadius,EnergyGrid,Ec,Bet_new,Sia_new,Zmetera_new,lowYN,PixelEnhancement,MSEnhancement)

		Detected.append(Detected_Ec)

	# Convert Detected (Detected is now RhoxEc)
	DetectedNew = []

	for i_Ec in range(len(Ec)):

		Det_aux = []
		for i_rho in range(len(rho)):
			Det_aux.append(Detected[i_rho][i_Ec])

		DetectedNew.append(Det_aux)

	# Fix problem with radius SORT

	rho,DetectedNew = swapRhoAndY(rho,DetectedNew,Ec)

	# Interpolate to new rho

	DetectedNewRho = []
	for i_Ec in range(len(Ec)):
		DetectedNew_interp = np.interp(newRho,rho,DetectedNew[i_Ec])
		DetectedNewRho.append(DetectedNew_interp)

	# import pdb
	# pdb.set_trace()


	return DetectedNewRho


def getValuesMidplane(R,Z,Y):
    
    
    indexZ = [i for i,j in enumerate(Z) if j>=0][0]
    
    Ymidplane = []
    for i in range(len(Y)):
        Ymidplane.append(Y[i][indexZ])
    
    return Ymidplane


def getRealSpacePerpendicularCoord(rho,typeMapping):
	
    # typeMapping= 'LF','HF','BOTH'
    
    # --------------------------------------------------------
    # ------- READ ALL VALUES FROM EQUILIBRIUM RECONSTRUCTION
    # --------------------------------------------------------
    
    # MAJOR RADIUS: R

	f=open('Rvalues.txt',"r")
	lines=f.readlines()
	R=[]
	for x in lines:
	    R.append(float(x))
	f.close()

    # VERTICAL COORDINATE: Z
	f=open('Zvalues.txt',"r")
	lines=f.readlines()
	Z=[]
	for x in lines:
	    Z.append(float(x))
	f.close()

	# rho???????????????????????????????????????
	f=open('psiPolNormSqrt.txt',"r")
	lines=f.readlines()
	Psi2D=[]
	for x in lines:
		aux = x.split()
		vector = []
		for i in range(len(aux)):
			vector.append(float(aux[i]))

		Psi2D.append(vector)
	f.close()


    # --------------------------------------------------------
    # ------- CALCULATIONS
    # --------------------------------------------------------
    
     # ONLY CARE ABOUT MIDPLANE (Z=0)
	indexZ = [i for i,j in enumerate(Z) if j>=0][0]

	Psi1D = Psi2D[indexZ]


     # ONLY CARE ABOUT INSIDE OF LCFS
	Psi1D_new = []
	R_new = []
	for i in range(len(Psi1D)):
		if Psi1D[i]<=1.0:
			Psi1D_new.append(Psi1D[i])
			R_new.append(R[i])

	Psi1D = Psi1D_new
	R = R_new
    

     # ONLY CARE ABOUT LF or HF
	Psi1D_new = []
	R_new = []
	if typeMapping == 'LF':

		for i in range(len(Psi1D)):
			if R[i]>=1.63:
				Psi1D_new.append(Psi1D[i])
				R_new.append(R[i])

		Psi1D = Psi1D_new
		R = R_new
        
	elif typeMapping == 'HF':
		Psi1D_new = []
		R_new = []
		for i in range(len(Psi1D)):
			if R[i]<=1.63:
				Psi1D_new.append(Psi1D[i])
				R_new.append(R[i])

		Psi1D = Psi1D_new
		R = R_new
        
        
        

    # We have Psi(R) or R(Psi) -> I want R(Psi_user), which is R(rho)
 
	if typeMapping == 'HF':
		rho = rho[::-1]
    
	RMAJ = np.interp(rho,Psi1D,R)

	return RMAJ


def computeChordIntegral(rho,Emission):
    #EmissionNew = Emission

    R = 165.0  #cm
    a = 50.0   #cm
    
    r = getRealSpacePerpendicularCoord(rho,'LF') * 100.0 - R #cm  
    
    for i in range(len(r)):
        if r[i] < 0:
            r[i] = 0
    
    #xChord = 2.0*R*a*rho + a**2*rho**2
    xChord = np.sqrt(2.0*R*r + r**2)     #cm
    
    Emission_Diffchord = 2*np.trapz(Emission,x=xChord)
    
    lengthChord = 2*np.trapz(np.ones(len(xChord)),x=xChord)    
    
    return lengthChord,Emission_Diffchord
    


def computeGeometricIntegral(rho,Emission,EtendueCM,nameCounts,Ec):
    # In seconds!!!!
    
    lengthChord,Emission_Diffchord = computeChordIntegral(rho,Emission)
    
    Emission_chord = Emission_Diffchord*EtendueCM	#Detectable photons/s 
    
    
    Emission_chordMS = copy.deepcopy(Emission_chord)/1000.0
    
    print 'Ec = {} -> {} {}/pixel/ms'.format(str(Ec),str(Emission_chordMS),nameCounts)
    

    return Emission_chord  # In seconds!!!!



def computeUncertainty(photons):
    
    fractionalTe = 1/sqrt(photons)
    
    return fractionalTe


def plotProfiles(fileName):
      
    f, ax = plt.subplots(1)
    axi = ax.twinx()
    
    #Loading data given by Alex
    profiledata = np.loadtxt(fileName,skiprows = 1)
    rhodata = profiledata[:,0]
    Tedata = profiledata[:,1]/1000.0
    nedata = profiledata[:,3]*1.0E-20

    ax.plot(rhodata,Tedata,linewidth=2,color='b')
    axi.plot(rhodata,nedata,linewidth=2,color='r')
    ax.set_xlabel('rho')
    ax.set_ylabel('Electron Temperature (keV)')
    axi.set_ylabel('Electron density ($10^{20} m^{-3}$)')
    axi.set_ylim([0,5])
    
    
    plt.show()
