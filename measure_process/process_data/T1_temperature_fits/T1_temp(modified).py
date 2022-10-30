import numpy as np
import matplotlib.pyplot as plt

import scipy
import scipy.constants 
from tqdm import tqdm
import sys

from measure_process.Init_file import data_dump
 #some of these values can be found in Relaxation of single-electron spin qubits in silicon, by Guido buckard

pi = np.pi
e = scipy.constants.value('atomic unit of charge')
ub = scipy.constants.value('Bohr magneton in eV/T')
g = 2 
m0 = scipy.constants.value('electron mass') 
kb = scipy.constants.value('Boltzmann constant in eV/K') 
ht = scipy.constants.value('reduced Planck constant in eV s')


rho = 2330 #silicon mass density 
vl = 9330 #speed sound lateral 
vt = 5420 #speed sound transverse 

PhiU = 8.77 #eV  , #deformation potential constant
PhiD = 5 #eV , deformation potential constant 

It = 1/vl**7 * (2/3*PhiD**2 + 4/15 *PhiD*PhiU + 2/35*PhiU**2) + 8/105 * PhiU**2 /vt**7   # slightly modified from original
Il = 1/vl**7 * (2/3*PhiD**2 + 4/5 *PhiD*PhiU + 2/7*PhiU**2) + 4/35 * PhiU**2 /vt**7 #modified from original


#Spin-Valley mixing with phonons
Evs = 300e-6  #valley splitting
r = 69e-9 #dipole moment, used the same for lateral and transversal phonons (how can i get a better value??)
def GammaValley_ph(omega):
	return  (ht**2/(2*pi))*(ht*omega)**5 * r**2 / (4*pi*rho*ht**6)* (It + Il)*e



#Spin-Valley mixing with Johnson
#Johnson noise parameters
l = 40e-9  #distance between gate and QD, this assumes the johnson noise is from the gates, not from 2DEG
R = 6e3 #?? #ohmic-to-ohmic resisntace? as measured it is ~152kOhm for SD1 and 98KoHm for SD2


Delta = 0.4e-9 #splitting at the anticrossing point (can i measure this better?)
Rk = 25.813e3  #quantum hall resistance
def GammaValley_jh(omega,T):
	return  2* R / Rk * omega * ht**2 * 3*r**2 / l**2 *1/np.tanh(omega*ht/(2*kb*T/1000)) #+ #2* R*3 / Rk * omega * ht**2 * 3*r**2 / l**2

S0=(0.1e-6)**2
l=50e-9 #roughly size of QD
def GammaValley_1_f(omega,T):
	return  2 * S0 * (T/1000)**2 /omega*3*r**2 / l**2#check if this has the correct units or maybe needs to be weighted by 
#Fast spin-valley-based quantum gates in Si with micromagnets -  HUANG


#bose distribution for number of phonons
def bose (omega, T):
    return  1/(np.exp(ht*omega/(kb*T/1000))-1)

#magnetic field applied externally, for qubit splitting
#B-fields Q=[0.5712, 0.581, 0.589, 0.59, 0.585, 0.578]
B = 0.578 #computed from the qubit splitting it should be 0.575 (16.3GHz)

omega = g*ub*B/ht


#%%plot energies vs B-field
plot = False
if plot:
    B = np.linspace(0,4,100)
    omega = g*ub*B/ht
    E0 = -(Evs + ht*omega)/2 
    E1 = -np.sqrt((Evs - ht*omega)**2+Delta**2) / 2 
    E2 = +np.sqrt((Evs - ht*omega)**2+Delta**2) / 2 
    E3 = (Evs + ht*omega)/2 
    E = np.zeros(4) 
    E = np.array([E0,E1,E2,E3])
    
    plt.figure()
    plt.title(f'E_vs = {Evs*1e6}ueV' )
    
    plt.plot(B,E0, label='E0')
    plt.plot(B,E1, label='E1')
    plt.plot(B,E2, label='E2')
    plt.plot(B,E3, label='E3')
    
    plt.xlabel('B-field')
    plt.ylabel('Energy')
    plt.legend()


#%% Simulation for fixed magnetic field



#Matrix of coefficients between states
a1 = -(Evs-ht*omega) / np.sqrt((Evs-ht*omega)**2+Delta**2)
a2 = -(Evs+ht*omega) / np.sqrt((Evs+ht*omega)**2+Delta**2) 

def SinGamma(x):
	return np.sqrt((1-x)/2)
def CosGamma(x):
	return np.sqrt((1+x)/2)

#how is matrix obtained???
valley_coupling = np.zeros((4,4)) 
valley_coupling[0,1] = -CosGamma(a1)*SinGamma(a2) - SinGamma(a1)*CosGamma(a2) 
valley_coupling[0,2] = SinGamma(a1)*SinGamma(a2) - CosGamma(a1)*CosGamma(a2) 
valley_coupling[0,3] = SinGamma(a2)*CosGamma(a2) 
valley_coupling[1,2] = SinGamma(a1)*CosGamma(a1) 
valley_coupling[1,3]= -CosGamma(a1)*CosGamma(a2) + SinGamma(a1)*SinGamma(a2)
valley_coupling[2,3] = SinGamma(a1)*CosGamma(a2) + CosGamma(a1)*SinGamma(a2) 


valley_coupling  =valley_coupling+ np.transpose(np.conjugate(valley_coupling))

# energy of states
E0 = -(Evs + ht*omega)/2 
E1 = -np.sqrt((Evs - ht*omega)**2+Delta**2) / 2 
E2 = +np.sqrt((Evs - ht*omega)**2+Delta**2) / 2 
E3 = (Evs + ht*omega)/2 
E = np.array([E0,E1,E2,E3])


# Getting Rates
index=np.array([0,1,2,3])
temp_res = 200
Gamma = np.zeros((4,4)) 
# original: zeros(4,4, temp_res) ; % ij means transition from i to j
rate = np.zeros((4,4,temp_res)) #ij means transition from i to j
rate_first_order = np.zeros((4,4,temp_res))
rate_second_order = np.zeros((4,4,temp_res))
rate_second_order_jh = np.zeros((4,4,temp_res))
rate_second_order_twophonon = np.zeros((4,4,temp_res))
#rate_orbit= zeros(1, temp_res)
rate_orbit = np.zeros((1,temp_res))# is that supposed to be a column??
Gammac = np.zeros((4, temp_res))
temperatures= np.linspace(115, 1200, temp_res)


#%% Calculate first order relaxation rates.



for kk, T in tqdm(enumerate(temperatures), total=len(temperatures), desc='1st order'):
    #these are the first order rates
    for ii in range(4):
        for jj in range(4):
            if jj != ii:
                if jj>ii:
                    rate_phonon = 2*pi/ht**2 * valley_coupling[ii,jj]**2 * GammaValley_ph(abs(E[ii]-E[jj])/ht) * bose(abs(E[ii]-E[jj])/ht, T)
                    rate_johnson = 2*pi/ht**2 * valley_coupling[ii,jj]**2 * GammaValley_jh(abs(E[ii]-E[jj])/ht, T) * bose(abs(E[ii]-E[jj])/ht, T) 
                    rate_1_f = 2*pi/ht**2 * valley_coupling[ii,jj]**2 * GammaValley_1_f(abs(E[ii]-E[jj])/ht, T) * bose(abs(E[ii]-E[jj])/ht, T) 
                else:
                    rate_phonon = 2*pi/ht**2 * valley_coupling[ii,jj]**2 * GammaValley_ph(abs(E[ii]-E[jj])/ht) * (1+bose(abs(E[ii]-E[jj])/ht, T)) 
                    rate_johnson = 2*pi/ht**2 * valley_coupling[ii,jj]**2 * GammaValley_jh(abs(E[ii]-E[jj])/ht, T) * (1+bose(abs(E[ii]-E[jj])/ht, T))                 
                    rate_1_f = 2*pi/ht**2 * valley_coupling[ii,jj]**2 * GammaValley_1_f(abs(E[ii]-E[jj])/ht, T) * (1+bose(abs(E[ii]-E[jj])/ht, T))                 
                
                

                rate_first_order[ii, jj, kk] = rate_phonon  + rate_johnson + rate_1_f
rate = rate_first_order   
#so far we've only computed the first order rates

#%% Calculate first order relaxation rates.

# def coeff1(omega1):
#     return valley_coupling[jj,m[0]] * valley_coupling[m[0],ii]/(E[m[0]]-E[ii]-np.sign(E[m[0]]-E[ii])*ht*omega1 + 0.5*ht*(ii+1)*Gammac[m[0], kk])
# def coeff2(omega1):
#     return valley_coupling[jj,m[1]] * valley_coupling[m[1],ii]/(E[m[1]]-E[ii]-np.sign(E[m[1]]-E[ii])*ht*omega1 + 0.5*ht*(ii+1)*Gammac[m[1], kk])
# def coeff(omega1):
#     return abs(coeff1(omega1)+coeff2(omega1))**2



# for kk, T in tqdm(enumerate(temperatures), total=len(temperatures),  desc='2nd order'):
#     tmp = sum(rate_first_order, 0)
#     Gammac[:, kk] = tmp[:,kk] + 1000 #;  %% 1 ms is the maximum lifetime limited by readout

#     omega_min = 100




#     # Calculate Raman, Orbach process
#     for ii in range(4):
#         for jj in range(4):
#             if jj != ii:
                
#                 m = index[np.where(index != ii)]
#                 m = m[np.where(m != jj)] 
#                 map = map_Raman([int2str(ii), int2str(jj)] ) #dafuq is this
                








# rate = rate_first_order    
#%% Solving rate equations

Rates_matrix = np.zeros((4,4,temp_res))
eigenvalues_all = np.zeros((4,4))
eigenvectors_all = np.zeros((4,4))
eigenvalues = np.zeros((temp_res,4))
eigenvectors = np.zeros((temp_res,4))

for ii in range(4):
    for jj in range(4):
        if jj != ii:
            Rates_matrix[ii,jj,:] =  rate[jj,ii,:]
        else:
            
            Rates_matrix[ii,ii,:] = -np.sum(rate[ii], axis=0)


for ii in range(temp_res):
   eigenvalues_all, eigenvectors_all = np.linalg.eig(np.reshape(Rates_matrix[:,:,ii],(4,4)))
   #there is something weird with the values signs and the first element, maybe global phase
   
   
   
   tmp = np.sort(abs(eigenvalues_all)) 
   tmp2 = eigenvectors_all[:,np.where(abs(eigenvalues_all) == tmp[1])]
   tmp2  = np.reshape(tmp2, (4))
   eigenvectors[ii,:] = abs(tmp2[:])
   
   eigenvalues[ii,:] = np.sort(abs(eigenvalues_all[:]))



#%% Plotting

y = 1/abs(eigenvalues[:,1]) #these are T1 times
x = temperatures


save = True
name = 'update'
if save:
    data_dic={'T':x, 'Q1': y}
    data_dump(data_dic, 'T1/fits',f'Q1_1p_{name}' )

plot=False
if plot:
    plt.figure()
    plt.semilogy(x, y, label = f'{Evs*1e6} ueV')
    plt.xlabel('Temperature [mK]')
    plt.xlabel('T1 [s]')
    plt.legend()
    
