from __future__ import division
import numpy as np 
from filter_functions import *
from integral_functions import *

# MQ1 parameters:
Rmaj = 165.0  #cm
a = 50.0   #cm

def BellInversion(newRho, n):
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
	RmajRho = getRealSpacePerpendicularCoord(newRho)
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


#############################
def computeChordIntegral(newRho,R_T,Emission):
    ''' Compute emission integral for a chord having tangency
    radius R_T. 
    Inputs:
    	rho: np.ndarray
    		Radial coordinate given as rho_psi 
    	R_T: float [cm]
    		Value of R (major radius at emission surface) at the
    		point of tangency of the chord with the emission 
    		surface.
    	Emission: np.ndarray
    		Vector containing the emissivity as a function of rho. 
    		This may be the emissivity impinging on a pixel, or the 
    		one deriving from a difference of thresholds, etc.
	
	Outputs: 
		chord_length, chord_emission

    '''
    # R_T = getRealSpacePerpendicularCoord(newRho) * 100.0
    R = getRealSpacePerpendicularCoord(newRho) * 100.0 #cm  
    
    epsilon = []
 	for i in range(len(R)):
 		if R[i] > R_T:
 			epsilon.append(R[i]-R_T)

 			# TO BE COMPLETED


    # xChord = np.sqrt(2.0*Rmaj*r + r**2)     #cm
    
    # Emission_Diffchord = 2*np.trapz(Emission,x=xChord)
    
    # lengthChord = 2*np.trapz(np.ones(len(xChord)),x=xChord)    
    
    # return lengthChord,Emission_Diffchord
