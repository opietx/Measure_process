import numpy as np
import matplotlib.pyplot as plt

import scipy
import scipy.constants 
from tqdm import tqdm
import sys
if 'data_dump' not in dir():
    from measure_process.Init_file import data_dump


pi = np.pi
e = scipy.constants.value('atomic unit of charge')

ub = scipy.constants.value('Bohr magneton in eV/T')
g = 2 
m0 = scipy.constants.value('electron mass') 
kb = scipy.constants.value('Boltzmann constant in eV/K') 
ht = scipy.constants.value('reduced Planck constant in eV s')

omegaDebye = 13.81e13 
rho = 2200 
vl = 3750 
vt = 5900 

PhiU = 8.77 
PhiD = 5 

It = 16/105 * PhiU**2 
Il = 4*(PhiU**2/35 + 2*PhiD*PhiU/15 + PhiD**2/3) #is this 35 or 45???
mt = 0.19*m0


#%% Simulation for fixed magnetic field
B = 0.69 #compute from the qubit splitting
omega = g*ub*B/ht

#Matrix of coefficients between states
a1 = -(Evs-ht*omega) / np.sqrt((Evs-ht*omega)**2+Delta**2)
a2 = -(Evs+ht*omega) / np.sqrt((Evs+ht*omega)**2+Delta**2) 

def SinGamma(x):
	return np.sqrt((1-x)/2)
def CosGamma(x):
	return np.sqrt((1+x)/2)


valley_coupling = np.zeros((4,4)) 
valley_coupling[0,1] = -CosGamma(a1)*SinGamma(a2) - SinGamma(a1)*CosGamma(a2) 
valley_coupling[0,2] = SinGamma(a1)*SinGamma(a2) - CosGamma(a1)*CosGamma(a2) 
valley_coupling[0,3] = SinGamma(a2)*CosGamma(a2) 
valley_coupling[1,2] = SinGamma(a1)*CosGamma(a1) 
valley_coupling[1,3]= -CosGamma(a1)*CosGamma(a2) + SinGamma(a1)*SinGamma(a2)
valley_coupling[2,3] = SinGamma(a1)*CosGamma(a2) + CosGamma(a1)*SinGamma(a2) 


valley_coupling  =valley_coupling+ np.transpose(np.conjugate(valley_coupling))#' 

# energy of states
E = np.zeros(4) 
E[0] = -(Evs + ht*omega)/2 
E[1] = -np.sqrt((Evs - ht*omega)**2+Delta**2) / 2 
E[2] = +np.sqrt((Evs - ht*omega)**2+Delta**2) / 2 
E[3] = (Evs + ht*omega)/2 

# Getting Rates
index=np.array([0,1,2,3])
temp_res = 2
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

#Johnson noise parameters
l = 100e-9 
R = 10e3

Evs = 220e-6  #valley splitting
r = 69e-9 #dipole moment

Delta = 0.40e-9 #splitting at the anticrossing point
Rk = 25.813e3  #quantum resistance

def GammaValley_ph(omega):
	return  (ht**2/(2*pi))*(ht*omega)**5 * r**2 / (4*pi*rho*ht**6)* (It/vt**7 + 2*Il/vl**7)* e
def GammaValley_jh(omega):
	return  2* R / Rk * omega * ht**2 * 3*r**2 / l**2

def bose (omega, T):
    return  1/(np.exp(ht*omega/(kb*T/1000))-1)








for kk, T in tqdm(enumerate(temperatures), total=len(temperatures), desc='1st order'):
    #these are the first order rates
    for ii in range(4):
        for jj in range(4):
            if jj != ii:
                if jj>ii:
                    rate_phonon = 2*pi/ht**2 * valley_coupling[ii,jj]**2 * GammaValley_ph(abs(E[ii]-E[jj])/ht) * bose(abs(E[ii]-E[jj])/ht, T)
                    rate_johnson = 2*pi/ht**2 * valley_coupling[ii,jj]**2 * GammaValley_jh(abs(E[ii]-E[jj])/ht) * bose(abs(E[ii]-E[jj])/ht, T) 
                else:
                    rate_phonon = 2*pi/ht**2 * valley_coupling[ii,jj]**2 * GammaValley_ph(abs(E[ii]-E[jj])/ht) * (1+bose(abs(E[ii]-E[jj])/ht, T)) 
                    rate_johnson = 2*pi/ht**2 * valley_coupling[ii,jj]**2 * GammaValley_jh(abs(E[ii]-E[jj])/ht) * (1+bose(abs(E[ii]-E[jj])/ht, T))                 
                

                rate_first_order[ii, jj, kk] = rate_phonon  + rate_johnson
#so far we've only computed the first order rates

#%% Calculate first order relaxation rates.

def coeff1(omega1):
    return valley_coupling[jj,m[0]] * valley_coupling[m[0],ii]/(E[m[0]]-E[ii]-np.sign(E[m[0]]-E[ii])*ht*omega1 + 0.5*ht*(ii+1)*Gammac[m[0], kk])
def coeff2(omega1):
    return valley_coupling[jj,m[1]] * valley_coupling[m[1],ii]/(E[m[1]]-E[ii]-np.sign(E[m[1]]-E[ii])*ht*omega1 + 0.5*ht*(ii+1)*Gammac[m[1], kk])
def coeff(omega1):
    return abs(coeff1(omega1)+coeff2(omega1))**2



for kk, T in tqdm(enumerate(temperatures), total=len(temperatures),  desc='2nd order'):
    tmp = sum(rate_first_order, 0)
    Gammac[:, kk] = tmp[:,kk] + 1000 #;  %% 1 ms is the maximum lifetime limited by readout

    omega_min = 100




    # Calculate Raman, Orbach process
    for ii in range(4):
        for jj in range(4):
            if jj != ii:
                
                m = index[np.where(index != ii)]
                m = m[np.where(m != jj)] 
                map = map_Raman([int2str(ii), int2str(jj)] ) #dafuq is this
                








rate = rate_first_order    
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
if save:
    data_dic={'T':x, 'Q1': y}
    data_dump(data_dic, 'T1/fits','Q1_1p' )

plot=False
if plot:
    plt.figure()
    plt.semilogy(x, y, label = f'{Evs*1e6} ueV')
    plt.xlabel('Temperature [mK]')
    plt.xlabel('T1 [s]')
    plt.legend()
    
