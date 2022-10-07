import numpy as np
import matplotlib.pyplot as plt


def FFT(x,y,title='',plot=True):
    yf = np.fft.fft(y)
    
    T =(x[1]-x[0])*1e-9
    N = len(yf)
    xf = np.fft.fftfreq(N, T)[:N//2]
    
    time= x[-1]*1e-9
    print(f'Frequency resolution = {round(1/time*1e-6,2)} MHz')
    print(f'BW = {round(N/(2*time)*1e-6,1)} MHz')
    plt.figure()
    plt.plot(xf/1e6, 2.0/N *abs(yf[0:N//2]), 'o-') 
    plt.xlabel('Frequency [MHz]')
    plt.ylabel('Amplitude')
    # plt.yscale('log')
    plt.title(title)
    plt.show()