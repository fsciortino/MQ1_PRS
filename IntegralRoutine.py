# import sys
# sys.path.append('/Users/pablorf/OMFIT-source/src/')
# from omfit_plot import *


import numpy as np 
import matplotlib.pyplot as plt 
from filter_functions import *
from integral_functions import *
import copy
import scipy.io
import math
from matplotlib.pyplot import ion
ion()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                               User INPUT
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Main parameters
Impurities       = ['Mo', 'C'] #['D','N','C','Mo']
LowYN            = False
Ec               = [2,7,14]#[7,8,9,10,11,12,13,14,15,16,17,18,19,20,21]#[2,5,10,15,20]    [-3,3] #     # in keV   [-3,3]
PixelEnhancement = 1.0#1.0    # More than 1 pixel?
MSEnhancement    = 1.0    # More than 1 ms?

# Etendue calculation
distancePinPixCM = 20              # in cm
AreaPixelCM2     = (8.3/487)**2    # in cm^2
AreaPinholeCM2   = 0.08*0.02       # in cm^2
cosines          = 1.0

# Neutrons calculation

#SXRInfoFolder                  = '/Users/pablorf/Dropbox (MIT)/'
SXRInfoFolder                  = 'C:\Users\sciortino\Dropbox (MIT)/'
ProbabilityNeutronsInteraction = 0.01
MaximumFluence                 = 1.0E15  # neutrons/cm^2
LifeMQ1                        = 3000.0  # 10s shots
NoCollimatedNeutrons           = 1.0E14  # neutrons/cm^2/s. In surface, Tinguely-> 10^18 n/m^2/s
FactorShielding                = 1E-5

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                           COMPUTATION
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Some parameters
fileNeutrons = SXRInfoFolder + 'SXR Info/Integration/DTNeutronYield2D_MQ1.mat'
newRho       = np.linspace(0,1,100)

EtendueCM    = AreaPinholeCM2/(4.0*np.pi*distancePinPixCM**2)*AreaPixelCM2*cosines**4.0
print 'Etendue Factor = {}cm'.format(str(EtendueCM))

col  = ['r','b','g','m','o','y','k','r','b','g','m','o','y','k','r','b','g','m','o','y','k']
lins = ['--','--','--','--','--','--','--','--','--','--','--','--','--','--','--','--','--','--','--','--','--','--','--','--']#'--','-.',':','-','--','-.',':','-','--','-.',':']

# Get materials data
en_Bet, Bet = get_Be_transmission()
en_Sia, Sia = get_Si_abs()
en_Zmeter, Zmetera = get_Zmeter_abs()

# Calculate Emissivity
rhoImpurity = []
energy2DImpurity = []
emissivity2DImpurity = []

DetectedImpurity = []
for Impurity in Impurities:

    # Get emissivity profiles (photons/cm^3/s, per impurity, radial position
    # (Te,ne mapping) and photon energy eV)
    rho,energy2D,emissivity2D,Te2D = loadDataNPZ(Impurity,LowYN)

    # Compute integral at each radial location (DETECTED photons/cm^3/s per
    # impurity per radial position)
    DetectedNewRho = computeAllRadius(rho,newRho,energy2D,Ec,emissivity2D,Bet,Sia,en_Bet,en_Sia,en_Zmeter, Zmetera,LowYN,PixelEnhancement,MSEnhancement)    

    print(Te2D[100])

    # SAVE
    outputfile = '{}_IntegralPhotons.npz'.format(Impurity)
    np.savez(outputfile,rho = newRho,Ec=Ec,Photons=DetectedNewRho)

    rhoImpurity.append(rho)
    energy2DImpurity.append(energy2D)
    emissivity2DImpurity.append(emissivity2D)
    DetectedImpurity.append(DetectedNewRho)

# Total and Difference
DetectedTOTAL = copy.deepcopy(DetectedImpurity[0])
for i_Imp in range(len(Impurities)-1):
    for i_Ec in range(len(Ec)):
        DetectedTOTAL[i_Ec] = DetectedTOTAL[i_Ec]+DetectedImpurity[i_Imp+1][i_Ec]

#DetectedTOTALDif = DetectedTOTAL[0] - DetectedTOTAL[1]
DetectedTOTALDiffs = [DetectedTOTAL[i+1] - DetectedTOTAL[i] for i in range(len(DetectedTOTAL)-1)]

# Geometric Integral
EmissionTotal = []
for i_Ec in range(len(Ec)):
    aux = computeGeometricIntegral(newRho,DetectedTOTAL[i_Ec],EtendueCM,'photons',Ec[i_Ec])
    EmissionTotal.append(aux/1000.0)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                            Neutrons
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Read Neutrons

mat = scipy.io.loadmat(fileNeutrons)
NeutronsR = []
for i in range(len(mat["R"])):
    NeutronsR.append(np.array(mat["R"].tolist())[i][0])
NeutronsZ = []
for i in range(len(mat["Z"])):
    NeutronsZ.append(np.array(mat["Z"].tolist())[i][0])
NeutronsY = []
for i in range(len(mat["Y"])):
    NeutronsYaux = []
    for j in range(len(mat["Y"])):
        NeutronsYaux.append(np.array(mat["Y"].tolist())[i][j] * 1.0E-6)     

    NeutronsY.append(NeutronsYaux)

# Map to midplane and to rho (and therefore corresponding RMAJ)
EmissionNeutronsR = getValuesMidplane(NeutronsR,NeutronsZ,NeutronsY)
RmajRho = getRealSpacePerpendicularCoord(newRho,'LF')
EmissionNeutrons = np.interp(RmajRho,NeutronsR,EmissionNeutronsR)

# Collimated neutrons per pixel (neutrons/pixel/ms)
print(' ')
print('Collimated (pinhole) neutrons:')
Neutrons_Collimated_Pixel = computeGeometricIntegral(newRho,EmissionNeutrons,EtendueCM,'neutrons',0) / 1000.0

# Neutrons from the surrondings, after shielding is applied (neutrons/cm^2/ms)
Flux_Uncollimated_CM2 = NoCollimatedNeutrons*FactorShielding / 1000.0

# Total flux of neutrons, accounting for collimated and surrondings (neutrons/cm^2/ms)
Flux_Total_CM2 = (Neutrons_Collimated_Pixel / AreaPixelCM2)  + Flux_Uncollimated_CM2
print('Flux Neutrons = {} neutrons/cm^2/ms'.format(str(Flux_Total_CM2)))

# Total neutrons on pixel, accounting for collimated and surrondings (neutrons/pixel/ms)
Neutrons_Uncollimated_Pixel = Flux_Uncollimated_CM2*AreaPixelCM2
Neutrons_Total_Pixel = Flux_Total_CM2*AreaPixelCM2

print('UnCollimated neutrons: {} neutrons/pixel/ms'.format(Neutrons_Uncollimated_Pixel))
print('Total neutrons per pixel: {} neutrons/cm^2/ms'.format(Neutrons_Total_Pixel))

# ~~~~~~~~~~~~~ Survivality
SurvivalSeconds = MaximumFluence / ( Flux_Total_CM2 )
ShotsNumber = SurvivalSeconds/10.0
PercentageLife = ShotsNumber / LifeMQ1 *100.0
print(' ')
print('Seconds survival = {}s'.format(str(SurvivalSeconds)))
print('Shots Survival (10s/shot) = {}'.format(str(int(ShotsNumber))))
print('Percentage of MQ1 life = {}%'.format(str(int(PercentageLife))))

# ~~~~~~~~~~~~~ Counts
FluxNeutronsInteract = ProbabilityNeutronsInteraction * copy.deepcopy(Flux_Total_CM2)
NeutronsInteractMS = ProbabilityNeutronsInteraction * copy.deepcopy(Neutrons_Total_Pixel)

# Compute neutrons that interact
print(' ')
print('Neutrons Interacting: {} neutrons/pixel/ms'.format(str(NeutronsInteractMS)))
#print('Flux Neutrons Interacting = {} neutrons/cm^2/ms'.format(str(FluxNeutronsInteract)))
print(' ')


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                            ELECTRON TEMPERATURE
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Ratios = []
TePredicted = []
for i_Ec in range(len(Ec)-1):
    
    RatioEc = EmissionTotal[i_Ec]/EmissionTotal[i_Ec+1]
    Ratios.append(RatioEc)
    TeEc    = -(Ec[i_Ec]-Ec[i_Ec+1])/math.log(RatioEc)
    TePredicted.append(TeEc)

print('Predicted Temperatures = {}'.format(TePredicted))

#Loading data given by Alex Creely
profiledata = np.loadtxt(SXRInfoFolder+'SXR Info/FLYCHK/profile1dAll.txt',skiprows = 1)
rhodata = profiledata[:,0]
Tedata = profiledata[:,1] / 1000.0

lengthChord,TeIntegrated = computeChordIntegral(rhodata,Tedata)

TeReal = TeIntegrated / lengthChord

print('Real Temperature = {}'.format(TeReal))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                   Plotting
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

f, ax0 = plt.subplots(1)
f, ax1 = plt.subplots(1)
f, ax2 = plt.subplots(1)

#fn = FigureNotebook(0,'RPLOT')
#fig, ax = fn.subplots(nrows=3)#,ncols=2)

axi = ax0.twinx()

for i_Imp in range(len(Impurities)):

    indexRho = 100
    rhoNow = rhoImpurity[i_Imp]
    energy2DNow = energy2DImpurity[i_Imp]
    emissivity2DImpurityNow = emissivity2DImpurity[i_Imp]
    
    ax0.semilogy(energy2DNow[indexRho]/1000.0,emissivity2DImpurityNow[indexRho],color=col[i_Imp],label='{} at $\\rho=${}'.format(Impurities[i_Imp],rhoNow[indexRho]))

    DetectedEc = DetectedImpurity[i_Imp]

    for i_Ec in range(len(Ec)):

        ax1.plot(newRho,EtendueCM*DetectedEc[i_Ec],color=col[i_Imp],linestyle=lins[i_Ec],label='{}, Ec = {}keV'.format(Impurities[i_Imp],str(Ec[i_Ec])))

DetectedEc = DetectedTOTAL
for i_Ec in range(len(Ec)):
    ax1.plot(newRho,EtendueCM*DetectedEc[i_Ec],color='k',linestyle=lins[i_Ec],linewidth=2,label='Total, Ec = {}keV'.format(str(Ec[i_Ec])))

axa = ax2
for i_Ec in range(len(Ec)):
    axi.plot(energy2D[0]/1000.0,S_pilatus(energy2D[0], Ec[i_Ec]*1000),'k',linestyle=lins[i_Ec])

axa.scatter(Ec,EmissionTotal)
axa.semilogy(Ec,EmissionTotal)
axa.set_ylabel('Detected $photons/pixel/ms$')
axa.set_xlabel('Ec (keV)')
axa.set_ylim([1,max(EmissionTotal)*1.2])

ax0.set_xlabel('Photon Energy [keV]')
ax0.set_ylabel('$photons/cm^3/s$')
ax0.set_xlim([0,25])
axi.set_ylabel('S-curves')

#ax0.legend(loc='best').draggable()
ax1.set_xlabel('$\\rho$')
ax1.set_ylabel('Detectable $photons/cm^3/s$')
ax1.legend(loc='best').draggable()
ax1.set_xlabel('$\\rho$')
ax1.set_ylabel('Detectable $photons/cm^3/s$')
ax1.legend(loc='best').draggable()

plt.show(block=False)
