import matplotlib.pyplot as plt
from core_tools.data.ds.data_set import load_by_uuid

#%%  #analyze resonance dips for boths sensors 
uuids = [1654610400421974076,1654614457460974076,1654615405396974076,1654617162067974076,1654619501694974076]#SD1
# uuids = [1654612345879974076,1654612793643974076,1654615872498974076,1654616613714974076,1654620017731974076]#SD2
t = [800,600,400,200,15]
ch='ch1'
plt.figure()
for ii, uuid in enumerate(uuids): 
    dic_data = {'x':None,'y':None}

    ds = load_by_uuid(uuid)
    x = ds[ch].average('y').x()
    
    
    I = ds['ch1'].average('y').y()
    Q = ds['ch2'].average('y').y()
    
    mod = np.sqrt(I**2 + Q**2)
    ang = np.rad2deg(np.arctan(I/Q))
    dic_data['x'] = x
    dic_data['y'] = mod
    data_dump(dic_data,'Resonator',str(t[ii]))
    plt.plot(x,mod,label=f'{t[ii]}')    

plt.xlabel('RF freq [Hz]')
plt.ylabel('Dig [mV]')
plt.legend()
plt.show()

#%% Plot IQ plane

ds = load_by_uuid(1654682272630974076)

# x = ds[ch].average('y').x()
# I = ds['ch3'].average('y').y()
# Q = ds['ch4'].average('y').y()

x = ds['ch3'].x()
I = ds['ch3'].y()
Q = ds['ch4'].y()


mod = np.sqrt(I**2 + Q**2)
ang = np.arctan(Q/I )+np.pi/2
f, ax1 = plt.subplots()
ax2=ax1.twinx() 
ax1.tick_params(axis='y', colors='blue')  
ax2.tick_params(axis='y', colors='red')  
plt1 = ax1.plot(x,np.rad2deg(ang), '-b',label='phase')   
ax1.set_ylabel('deg')
plt2 = ax2.plot(x,mod,'-r', label='magnitude')    
ax2.set_ylabel('dig[mV]')

ang = ang.max()
f, ax1 = plt.subplots()
ax1.grid()
ax1.plot(I, Q)

IQ_data  = I + 1j*Q
IQ_data *= np.exp(1j*ang)
ax1.plot(IQ_data.real,IQ_data.imag , 'r')

plt.show()

plt.figure()
plt.plot(x,I)
plt.plot(x,Q)
plt.plot(x,IQ_data.real)
plt.show()


#%% Analyze MW power

uuids = [1654627231928974076,1654627412181974076, 1654627601684974076, 1654628077774974076]
t = ['Off_11mK', '6dBm', '15dBm','18dBm']
ch='ch3'
plt.figure()
for ii, uuid in enumerate(uuids):    
    ds = load_by_uuid(uuid)
    x = ds[ch].average('y').x()
    
    
    #not sure which one of the two is better 
    y = np.sqrt(ds['ch3'].average('y').y()**2 + ds['ch4'].average('y').y()**2)
    # y = ds[ch].average('y').y()
    plt.plot(x,y,label=f'{t[ii]}')    

plt.xlabel('RF freq [Hz]')
plt.ylabel('Dig [mV]')
plt.legend()
plt.show()

#%% overlap peaks
#virtual plunger
uuids = [1654859913948974076,1654860953407974076,1654861890704974076,1654862323658974076,
         1654862679817974076,1654862989297974076,1654863425163974076,1654863913023974076,
         1654864124990974076,1654865136724974076,1654865364200974076,]
#physical plunger
t = [13,50,70, 100,125,150,175,200,225,250,275,300]
ch='ch1'
maxs_pos=[]

fig, axs = plt.subplots(1, 3)
for ii, uuid in enumerate(uuids):    
    ds = load_by_uuid(uuid)
    x = ds[ch].x()
    y = ds[ch].y()
    peak_pos = np.argpartition(y, 1)[:1]#this averages over the 'n'-minimum values of the plot
    peak_val = np.average(y[peak_pos])
    peak_pos = np.average(x[peak_pos])
    maxs_pos.append(peak_pos)
    axs[0].plot(x,y,'-',label=f'{t[ii]}mK')    
    axs[0].plot(peak_pos,peak_val,'o')
    axs[1].plot(x,y,'-',label=f'{t[ii]}mK')    
    axs[1].plot(peak_pos,peak_val,'o')


axs[0].set_xlabel('SD1_P [mV]')
axs[0].set_ylabel('Dig [mV]')
axs[0].set_title('Coulomb peaks')
axs[0].legend()
axs[1].set_xlabel('SD1_P [mV]')
axs[1].set_ylabel('Dig [mV]')
axs[1].set_title('Coulomb peaks')
axs[1].legend()


maxs_pos = np.array(maxs_pos)
axs[2].plot(t[:len(uuids)],maxs_pos,'-o')
axs[2].set_ylabel('Peak max position [mV]')
axs[2].set_xlabel('MXC temp [mK]')
plt.tight_layout()
plt.show()
#%%  analyze coulomb peaks
uuids = [1654629438122974076,1654629479143974076,1654629506014974076,1654629561860974076,1654629597536974076,
         1654629637295974076,1654629664782974076]
t = ['Off_11mK','0dBm','3dBm', '6dBm','9dBm','12dBm', '15dBm']
t = t[:len(uuids)]
ch='ch3'

maxs = []
Z = []
for ii, uuid in enumerate(uuids):    
    ds = load_by_uuid(uuid)
    Z.append(ds[ch].y())
    maxs.append(np.min(ds[ch].y())[0:200])

maxs = np.array(maxs)
Z = np.array(Z)    
x = ds[ch].x()
y = np.arange(len(t)+1)
X,Y = np.meshgrid(x,y)

if True:
    plt.figure()
    plt.xlabel('SDP [mV]')
    plt.ylabel('MW power [dBm]')
    plt.pcolormesh(X,Y,Z)
    plt.yticks(np.arange(len(t))+0.5,t)
    plt.colorbar()
    plt.show()





#make 1D cuts on peak and offpeak
if True:
    import matplotlib.pyplot as plt
    x,y = np.arange(0,18,3),Z[1:,200]
    plt.figure()
    ax1 = plt.subplot()
    ax1.plot(x,y,'og')
    m,b = np.polyfit(x, y, 1)
    ax1.plot(x,m*x+b,'-g', linewidth = 5, alpha=0.5, label='base')
    
    x,y = np.arange(0,18,3),maxs[1:]
    ax1.plot(x,y,'ob')
    m,b = np.polyfit(x, y, 1)
    ax1.plot(x,m*x+b,'-b', linewidth = 5, alpha=0.5, label='peak' )
    ax1.legend()
    ax1.set_xticks(np.arange(len(t)-1)*3)
    
    ax1.set_xlabel('MW_power[dBm]')
    ax1.set_ylabel('dig [mV]')

plt.show()
    


#%%
uuids = [1654629438122974076,1654633555502974076,1654633927948974076,1654634523928974076,1654634271191974076,1654635204091974076,
         1654635678754974076,1654636791307974076,1654637102995974076,1654637524350974076,1654637771354974076,1654637915752974076,
         1654638240740974076,1654638546572974076]

t = np.array([11,25,35,44,47,55,65,87,95,120,131,141.5,151,180])
ch='ch3'

maxs = []
Z = []
for ii, uuid in enumerate(uuids):    
    ds = load_by_uuid(uuid)
    Z.append(ds[ch].y())
    maxs.append(np.min(ds[ch].y()[0:200]))

Z = np.array(Z)    
x = ds[ch].x()

if True:
    import matplotlib.pyplot as plt
    x,y = t,Z[:,200]
    plt.figure()
    plt.plot(x,y,'om')
    m,b = np.polyfit(x, y, 1)
    plt.plot(x,m*x+b,'-m', linewidth = 5, alpha=0.5, label='base')
    
    x,y = t,maxs
    plt.plot(x,y,'or')
    m,b = np.polyfit(x, y, 1)
    plt.plot(x,m*x+b,'-r', linewidth = 5, alpha=0.5, label='peak' )
    plt.legend()
    # plt.xticks(,t[1:])
    plt.show()
    plt.xlabel('MX temperature [mK] ')
    plt.ylabel('dig [mV]')

#%% channel resistance vs temp
res = [152,152,152,151,151,151,152,153,153,154,157,159,164,165,165,166]
temps  = [800,600,450,300,250,225,200,175,150,125,100,75,50,37,30,20]
plt.figure()
plt.plot(temps,res, '-o')
plt.show()
plt.xlabel('MX temperature [mK] ')
plt.ylabel('Channel esistance [kOhm]]')

