# -*- coding: utf-8 -*-
"""
Created on Mon Feb 24 21:57:33 2020

@author: emlach
"""

import numpy as np
import matplotlib.pyplot as plt

surfaceP = 100000  # Pa
delta_p = 5000  # Pa

p_layer = np.zeros(21)  # 22 pressure interfaces range(0, 100000, 5000)

for P in range(0, 21):
    p_layer[P] = (100000 - delta_p * P)

p_layer

S_0 = 1360
alpha_p = 0.3
sigma = 5.67e-8  # W m-2 K-4
g = 9.8  # m/s^2
c_p = 1004  # J/kg*K
R_d = 287
emissivity_lw = [1, 0.3, 0.3, 0.3, 0.3, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1]

emissivity_sw = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0015625, 0.003125, 0.00625, 0.0125, 0.025, 0.05]

Temp_time = np.zeros((17473, 21))
Temp_time[0][:] = 50
longwave_up = np.zeros((17473, 21))
longwave_down = np.zeros((17473, 21))
shortwave = np.zeros((17473, 21))
delta_q = np.zeros((17473, 21))
lapse_rate = np.zeros((17473, 21))

for time in range(0, 17472, 1):
    # for time in range(0,3,1):
    # longwave up
    for layer_up in range(0, 21, 1):
        if layer_up == 0:
            through = 0
            emitted = sigma * Temp_time[time][layer_up] ** 4  # Wm^-2
            longwave_up[time][layer_up] = through + emitted
        elif 1 <= layer_up <= 20:
            through = (1 - emissivity_lw[layer_up]) * longwave_up[time][layer_up - 1]
            emitted = emissivity_lw[layer_up] * (sigma * Temp_time[time][layer_up] ** 4)
            longwave_up[time][layer_up] = through + emitted
        else:
            print('layering up is not working')
    # longwave down + shortwave down
    for layer_down in range(20, -1, -1):
        if layer_down == 20:
            through = 0
            emitted = emissivity_lw[layer_down] * (sigma * Temp_time[time][layer_down] ** 4)
            longwave_down[time][layer_down] = through + emitted
            shortwave[time][layer_down] = emissivity_sw[layer_down] * S_0 / 4 * (1 - alpha_p)  # absorbed at L20
        elif 1 <= layer_down <= 19:
            through = (1 - emissivity_lw[layer_down]) * longwave_down[time][layer_down + 1]
            emitted = emissivity_lw[layer_down] * sigma * Temp_time[time][layer_down] ** 4
            longwave_down[time][layer_down] = through + emitted
            shortwave[time][layer_down] = (S_0 / 4 * (1 - alpha_p) - np.sum(shortwave[time][layer_down + 1:21])) * emissivity_sw[layer_down]
        elif layer_down == 0:
            through = (1 - emissivity_lw[layer_down]) * longwave_down[time][layer_down + 1]
            emitted = 0
            longwave_down[time][layer_down] = through + emitted
            shortwave[time][layer_down] = (S_0 / 4 * (1 - alpha_p) - np.sum(shortwave[time][layer_down + 1:21])) * emissivity_sw[layer_down]
        else:
            print('layering down is not working')
    for layer in range(0, 21, 1):
        if layer == 20:
            difference_down = 0 - longwave_down[time][layer]
            difference_up = longwave_up[time][layer - 1] - longwave_up[time][layer]
            captured_shortwave = shortwave[time][layer]
            delta_q[time][layer] = difference_down + difference_up + captured_shortwave  # Wm^-2
        if 1 <= layer <= 19:
            difference_down = longwave_down[time][layer + 1] - longwave_down[time][layer]
            difference_up = longwave_up[time][layer - 1] - longwave_up[time][layer]
            captured_shortwave = shortwave[time][layer]
            delta_q[time][layer] = difference_down + difference_up + captured_shortwave
        if layer == 0:
            difference_down = longwave_down[time][layer + 1] - 0
            difference_up = 0 - longwave_up[time][0]
            captured_shortwave = S_0 / 4 * (1 - alpha_p)
            delta_q[time][layer] = difference_down + difference_up + captured_shortwave

    for temp_change in range(0, 21, 1):
        Temp_time[time + 1][temp_change] = Temp_time[time][temp_change] + (g / (delta_p * c_p)) * delta_q[time][
            temp_change] * 3600
    for change in range(0, 20, 1):
        lapse_rate[time][change] = (Temp_time[time][change] - Temp_time[time][change+1])


rho = []
for density in range(21):
    rho.append(p_layer[density]/(R_d*Temp_time[17471][density]))

dz = []
for change in range(20):
    dz.append(((p_layer[change] - p_layer[change+1])/(rho[change]*g))/1000)


Temp_time[17000]

np.sum(shortwave[0][19:21])
shortwave[17000]

strat_temp = Temp_time[17471]

strat_temp

plt.plot(Temp_time[17472], p_layer)
# plt.yscale('log')
plt.title('Temperature decrease with Pressure')
plt.ylabel('Pressure, log scale (Pa)')
plt.xlabel('Temp, K')
plt.gca().invert_yaxis()
plt.show()

strat_lapse = lapse_rate[17471, :20]/dz

plt.plot(lapse_rate[17471, :20]/dz, p_layer[:20])
plt.yscale('log')
plt.title('Temperature decrease with Pressure')
plt.ylabel('Pressure, log scale (Pa)')
plt.xlabel('Temp, K')
plt.gca().invert_yaxis()
plt.show()



def theta(level):
    return Temp_time[17472][level] * (p_layer[0] / p_layer[level]) ** 0.286


potential_temp = []
for i in range(0, 20, 1):
    potential_temp.append(theta(i))

plt.plot(potential_temp, p_layer[0:20])
plt.title('Potential Temp and Pressure')
plt.ylabel('Pressure, log scale(Pa)')
plt.xlabel('Potential Temp (K)')
plt.gca().invert_yaxis()
# plt.yscale('log')
plt.plot()
