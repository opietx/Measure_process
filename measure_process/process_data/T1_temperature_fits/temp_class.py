import numpy as np
import matplotlib.pyplot as plt

import scipy
import scipy.constants 
from tqdm import tqdm
import sys

from measure_process.Init_file import data_dump

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

Rk = 25.813e3  #quantum hall resistance


class T1_fits:
    def __init__(self, Evs, B, Delta, r, R_j, l_j, l_f, S0):
        self.Evs = Evs
        self.Delta = Delta
        self.r = r
        self.R_j = R_j
        self.l_j = l_j
        self.l_f = l_f
        self.S0 = S0
        
        self.temp_init = 15
        self.temp_final = 1500
        self.temp_res = 200
    
        self.B_init = 0.575
        self.B_final = 0.575
        self.B_res = 1
    
    
    @property
    def B(self):
        _B = np.linspace(self.B_init, self.B_final, self.B_res)
        # print(_temps)
        if len(_B) == 1:
            return _B[0] 
        else:
            return _B
    
    @property
    def omega(self):
        return g*ub*self.B/ht
    
    #Spin-Valley mixing with phonons
    def GammaValley_ph(self,omega ):
    	return  (ht**2/(2*pi))*(ht*omega)**5 * self.r**2 / (4*pi*rho*ht**6)* (It + Il)*e
    
    def GammaValley_jh(self,omega,T):
        return  2* self.R_j / Rk * omega * ht**2 * 3*self.r**2 / self.l_j**2 *1/np.tanh(omega*ht/(2*kb*T/1000))
    
    def GammaValley_1_f(self,omega,T):
        return  2 * (self.S0 * (T/1000)**2)/omega*3*self.r**2 / self.l_f**2#check if this has the correct units or maybe needs to be weighted by 
    #Fast spin-valley-based quantum gates in Si with micromagnets -  HUANG
    def GammaValley_power_law(self,omega,T):
        return 4*2 * self.S0 * (T/1000)**11 /omega*3*self.r**2 / self.l_f**2#check if this has the correct units or maybe needs to be weighted by 
    
    #bose distribution for number of phonons
    def bose (self,omega, T):
        return  1/(np.exp(ht*omega/(kb*T/1000))-1)   
    
    def plot_energy_levels(self):
        B = np.linspace(0,4,100)
        omega = g*ub*B/ht
        E0 = -(self.Evs + ht*omega)/2 
        E1 = -np.sqrt((self.Evs - ht*omega)**2+self.Delta**2) / 2 
        E2 = +np.sqrt((self.Evs - ht*omega)**2+self.Delta**2) / 2 
        E3 = (self.Evs + ht*omega)/2 
        
        plt.figure()
        plt.plot(B,E0, label='E0')
        plt.plot(B,E1, label='E1')
        plt.plot(B,E2, label='E2')
        plt.plot(B,E3, label='E3')
        plt.title(f'E_vs = {self.Evs*1e6}ueV' )
        plt.xlabel('B-field')
        plt.ylabel('Energy')
        plt.legend()
    
    @property
    def a1(self):
        return -(self.Evs-ht*self.omega) / np.sqrt((self.Evs-ht*self.omega)**2+self.Delta**2)
    
    @property
    def a2(self):
        return -(self.Evs+ht*self.omega) / np.sqrt((self.Evs+ht*self.omega)**2+self.Delta**2) 
    
    def SinGamma(self,x):
    	return np.sqrt((1-x)/2)
    def CosGamma(self,x):
    	return np.sqrt((1+x)/2)
    
    #how is matrix obtained???
    @property
    def valley_coupling(self):
        _vc = np.zeros((4,4)) 
        _vc[0,1] = -self.CosGamma(self.a1)*self.SinGamma(self.a2) - self.SinGamma(self.a1)*self.CosGamma(self.a2) 
        _vc[0,2] = self.SinGamma(self.a1)*self.SinGamma(self.a2) - self.CosGamma(self.a1)*self.CosGamma(self.a2) 
        _vc[0,3] = self.SinGamma(self.a2)*self.CosGamma(self.a2) 
        _vc[1,2] = self.SinGamma(self.a1)*self.CosGamma(self.a1) 
        _vc[1,3]= -self.CosGamma(self.a1)*self.CosGamma(self.a2) + self.SinGamma(self.a1)*self.SinGamma(self.a2)
        _vc[2,3] = self.SinGamma(self.a1)*self.CosGamma(self.a2) + self.CosGamma(self.a1)*self.SinGamma(self.a2) 
        return _vc+np.transpose(np.conjugate(_vc))
    
    @property
    def E(self):
        E0 = -(self.Evs + ht*self.omega)/2 
        E1 = -np.sqrt((self.Evs - ht*self.omega)**2+self.Delta**2) / 2 
        E2 = +np.sqrt((self.Evs - ht*self.omega)**2+self.Delta**2) / 2 
        E3 = (self.Evs + ht*self.omega)/2 
        return np.array([E0,E1,E2,E3])
    
    index=np.array([0,1,2,3])

    @property
    def temperatures(self):
        _temps = np.linspace(self.temp_init, self.temp_final, self.temp_res)
        # print(_temps)
        if len(_temps) == 1:
            return _temps[0] 
        else:
            return _temps
    
    
    
    
    temp_res = 200
    first_rate_johnson = np.zeros((4,4,temp_res))
    first_rate_phonon = np.zeros((4,4,temp_res))
    first_rate_1_f = np.zeros((4,4,temp_res))
    first_rate_power = np.zeros((4,4,temp_res))
    
    @property
    def first_rate(self):
        return self.first_rate_johnson + self.first_rate_1_f + self.first_rate_phonon
    
    
        
    def calculate_1st_rates(self):
        for kk, T in tqdm(enumerate(self.temperatures), total=self.temp_res, desc='1st order'):
        #these are the first order rates
            for ii in range(4):
                for jj in range(4):
                    arg=abs(self.E[ii]-self.E[jj])/ht
                    if jj != ii:
                        if jj>ii:
                            rate_phonon = 2*pi/ht**2 * self.valley_coupling[ii,jj]**2 * self.GammaValley_ph(arg) * self.bose(arg, T)
                            rate_johnson = 2*pi/ht**2 * self.valley_coupling[ii,jj]**2 * self.GammaValley_jh(arg, T) * self.bose(arg, T) 
                            rate_1_f = 2*pi/ht**2 * self.valley_coupling[ii,jj]**2 * self.GammaValley_1_f(arg, T) * self.bose(arg, T) 
                            rate_power = 2*pi/ht**2 * self.valley_coupling[ii,jj]**2 * self.GammaValley_power_law(arg, T) * self.bose(arg, T) 
                        else:
                            rate_phonon = 2*pi/ht**2 * self.valley_coupling[ii,jj]**2 * self.GammaValley_ph(arg) * (1+self.bose(arg, T)) 
                            rate_johnson = 2*pi/ht**2 * self.valley_coupling[ii,jj]**2 * self.GammaValley_jh(arg, T) * (1+self.bose(arg, T))                 
                            rate_1_f = 2*pi/ht**2 * self.valley_coupling[ii,jj]**2 * self.GammaValley_1_f(arg, T) * (1+self.bose(arg, T))                 
                            rate_power = 2*pi/ht**2 * self.valley_coupling[ii,jj]**2 * self.GammaValley_power_law(arg, T) * (1+self.bose(arg, T))                 
                        
                        
                        self.first_rate_johnson[ii, jj, kk] = rate_johnson
                        self.first_rate_phonon[ii, jj, kk] = rate_phonon
                        self.first_rate_1_f[ii, jj, kk] = rate_1_f
                        self.first_rate_power[ii, jj, kk] = rate_power
                    
        # return rate_phonon + rate_johnson + rate_1_f + rate_power


    @property
    def T1_phonon_1rate(self):
        return self.solve_rates(self.first_rate_phonon)

    @property
    def T1_johnson_1rate(self):
        return self.solve_rates(self.first_rate_johnson)
    @property
    def T1_1_f_1rate(self):
        return self.solve_rates(self.first_rate_1_f)
    @property
    def T1_power_1rate(self):
        return self.solve_rates(self.first_rate_power)
    @property
    def T1_total_1rate(self):
        return self.T1_phonon_1rate+self.T1_johnson_1rate+self.T1_1_f_1rate+self.T1_power_1rate

    
    
    def solve_rates(self, rate):
        Rates_matrix = np.zeros((4,4,self.temp_res))
        eigenvalues_all = np.zeros((4,4))
        eigenvectors_all = np.zeros((4,4))
        eigenvalues = np.zeros((self.temp_res,4))
        eigenvectors = np.zeros((self.temp_res,4))
            
        for ii in range(4):
            for jj in range(4):
                if jj != ii:
                    Rates_matrix[ii,jj,:] =  rate[jj,ii,:]
                else:
                    
                    Rates_matrix[ii,ii,:] = -np.sum(rate[ii], axis=0)    
                    
        for ii in range(self.temp_res):
           eigenvalues_all, eigenvectors_all = np.linalg.eig(np.reshape(Rates_matrix[:,:,ii],(4,4)))
           #there is something weird with the values signs and the first element, maybe global phase
                      
           tmp = np.sort(abs(eigenvalues_all)) 
           tmp2 = eigenvectors_all[:,np.where(abs(eigenvalues_all) == tmp[1])]
           tmp2  = np.reshape(tmp2, (4))
           eigenvectors[ii,:] = abs(tmp2[:])
           
           eigenvalues[ii,:] = np.sort(abs(eigenvalues_all[:]))
    
        return 1/abs(eigenvalues[:,1]) #these are T1 times
    
    def plot_rates(self):
        items_plot = [self.T1_phonon_1rate,self.T1_johnson_1rate,self.T1_1_f_1rate,self.T1_power_1rate,self.T1_total_1rate]
        for item in items_plot:
            plt.figure(55)
            plt.semilogy(self.temperatures, item,'.-')
            plt.xlabel('temps [mK]')
            plt.ylabel('T1 [s]')
            
        



    


a = T1_fits(300e-6,0.575,10e-9,0.75e-9,150e3,75e-9,55e-9,0.04e-6)
a.calculate_1st_rates()
a.plot_rates()



        
        
        
    
    
    
        