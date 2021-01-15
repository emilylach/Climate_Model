# -*- coding: utf-8 -*-
"""
Created on Mon Feb 24 15:47:46 2020

@author: emlach
"""


import numpy as np
import matplotlib.pyplot as plt


surfaceP = 100000  # Pa
delta_p = 5000  # Pa

p_layer = np.zeros(21) #22 pressure interfaces range(0, 100000, 5000)

for P in range(0, 21):
    p_layer[P] = (100000-delta_p*P)

p_layer

S_0 = 1360
alpha_p = 0.3
sigma = 5.67e-8
g = 9.8 # m/s^2
c_p = 1004 # J/kg*K
emissivity_lw = [1, 0.3, 0.3, 0.3, 0.3, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1]
emissivity_sw = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]


Temp_time = np.zeros((17473, 21))
Temp_time[0][:] = 50
longwave_up = np.zeros((17473, 21))
longwave_down = np.zeros((17473, 21))
shortwave = np.zeros((17473, 21))
delta_q = np.zeros((17473, 21))


for time in range(0, 17472, 1): #this is hours
#for time in range(0,3,1):
    #longwave up 
    for layer_up in range(0,21,1):
        if layer_up == 0:
            through = 0
            emitted = sigma*Temp_time[time][layer_up]**4 # Wm^-2
            longwave_up[time][layer_up] = through + emitted 
        elif 1 <= layer_up <= 20:
            through = (1-emissivity_lw[layer_up])*longwave_up[time][layer_up-1]
            emitted = emissivity_lw[layer_up]*(sigma*Temp_time[time][layer_up]**4)
            longwave_up[time][layer_up] = through + emitted
        else:
            print('layering up is not working')
    # longwave down
    for layer_down in range(20, -1, -1):
        if layer_down == 20:
            through = 0
            emitted = emissivity_lw[layer_down]*(sigma*Temp_time[time][layer_down]**4)
            longwave_down[time][layer_down] = through + emitted
            shortwave[time][layer_down] = emissivity_sw[layer_down]*S_0/4*(1-alpha_p)
        elif 1 <= layer_down <= 19:
            through = (1-emissivity_lw[layer_down])*longwave_down[time][layer_down+1] 
            emitted = emissivity_lw[layer_down]*sigma*Temp_time[time][layer_down]**4 
            longwave_down[time][layer_down] = through + emitted
            shortwave[time][layer_down] = emissivity_sw[layer_down]*S_0/4*(1-alpha_p)
        elif layer_down == 0:
            through = (1-emissivity_lw[layer_down])*longwave_down[time][layer_down+1] 
            emitted = 0
            longwave_down[time][layer_down] = through + emitted
            shortwave[time][layer_down] = emissivity_sw[layer_down]*S_0/4*(1-alpha_p) # is this getting emitted down??
        else:
            print('layering down is not working')
    for layer in range(0, 21, 1): 
        if layer == 20:
            difference_down = 0 - longwave_down[time][layer]
            difference_up = longwave_up[time][layer-1] - longwave_up[time][layer]
            captured_longwave = 0
            delta_q[time][layer] = difference_down + difference_up + captured_longwave # Wm^-2
        if 1 <= layer <= 19:
            difference_down = longwave_down[time][layer+1] - longwave_down[time][layer]
            difference_up = longwave_up[time][layer-1] - longwave_up[time][layer]
            captured_longwave = 0
            delta_q[time][layer] = difference_down + difference_up + captured_longwave
        if layer == 0:
            difference_down = longwave_down[time][layer+1] - 0
            difference_up = 0 - longwave_up[time][0]
            captured_longwave = S_0/4*(1-alpha_p)
            delta_q[time][layer] = difference_down + difference_up + captured_longwave

    for temp_change in range(0,21,1):
        Temp_time[time+1][temp_change] = Temp_time[time][temp_change] + (g/(delta_p*c_p))* delta_q[time][temp_change]*3600 


x = []
for i in range(0,17000,1):
    x.append(i)





Temp_time[17472]

p_layer
fig, ax = plt.subplots()

ax.plot(Temp_time[17472], p_layer)
ax.set_yscale('log')
#ax.set_yticks([20, 300, 500])
plt.show()

plt.plot(x, Temp_time[0:17000,20], label = 'layer 20')
plt.plot(x, Temp_time[0:17000,19])
plt.plot(x, Temp_time[0:17000,18])
plt.plot(x, Temp_time[0:17000,17])
plt.plot(x, Temp_time[0:17000,0])
plt.legend()
plt.show()


# Temperature profile


plt.plot(Temp_time[17472], p_layer)
plt.yscale('log')
plt.title('Temperature decrease with Pressure')
plt.ylabel('Pressure, log scale (Pa)')
plt.xlabel('Temp, K')
plt.gca().invert_yaxis()
plt.show()



# Lapse rate of atmos w/o stratosphere

km_per_layer = 12/20 # assuming height of 12km 

lapse_rate = []
for i in range(0, 20, 1):
    lapse_rate.append((Temp_time[17472][i] - Temp_time[17472][i+1])/km_per_layer)
    
lapse_rate #K per km


plt.plot(lapse_rate, p_layer[1:21])
plt.yscale('log')
plt.title('Lapse rate')
plt.ylabel('Pressure, log scale(Pa)')
plt.xlabel('Lapse rate (K/km)')
plt.gca().invert_yaxis()
plt.plot()


# Potential temp of atmos w/o stratosphere


def theta(layer):
    return (Temp_time[17472][layer]*(p_layer[0]/p_layer[layer])**(0.286))

potential_temp = []
for i in range(0,20,1):
    potential_temp.append(theta(i))

plt.plot(potential_temp, p_layer[0:20])
plt.title('Potential Temp and Pressure')
plt.ylabel('Pressure, log scale(Pa)')
plt.xlabel('Potential Temp (K)')
plt.gca().invert_yaxis()
#plt.yscale('log')
plt.plot()

shortwave[17000]











