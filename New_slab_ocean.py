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
sigma = 5.67e-8
g = 9.8  # m/s^2
c_p = 1004  # J/kg*K
R_d = 287  # J/K*kg
#lapse_crit = 6.5  # delta degrees per km
lapse_dry = 9.8  # delta degrees per km
L = 2.5e6  # J/kg
c_w = 4218
rho_w = 1000
h_5 = 5
rho_sh = 1.25
C_e = 8e-3
V = 10
finalT_prof = [295.15249953, 293.63274284, 291.34483437, 289.01759418, 286.64231945, 284.15704909, 281.52610164, 278.70559272, 275.63756602, 272.24238259, 268.40648433, 263.96304841, 258.66225378, 252.07169314, 243.63753895, 238.85281145, 235.63567103, 232.84225561, 231.1298435,  231.88331856,237.90207178]
finalOLR = 237.98930041569122
h_changing = 5
forcing = 0


emissivity_lw = [1, 0.3, 0.3, 0.3, 0.3, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1]
emissivity_sw = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0015625, 0.003125, 0.00625, 0.0125, 0.025, 0.05]

Temp_time = np.zeros((52561, 21))
Temp_time[0][:] = 50
longwave_up = np.zeros((52561, 21))
longwave_down = np.zeros((52561, 21))
shortwave = np.zeros((52561, 21))
delta_q = np.zeros((52561, 21))
adjust_layer = np.zeros((52561, 21))
adjust_upperlayer = np.zeros((52561, 20))
layer_height = np.zeros((52561, 20))
midpoint = np.zeros((52561, 20))
delta_z = np.zeros((52561, 19))
lapse_crit = np.zeros((52561, 19))
dq_star = np.zeros((52561, 21))
temp_changes = np.zeros((52561, 21))
temp_at_layers = np.zeros((52561, 20))
SH = np.zeros(52561)

for time in range(0, 52560, 1):
    for temp_calc in range(19):
        temp_changes[time][temp_calc] = (temp_at_layers[time][temp_calc] - temp_at_layers[time][temp_calc + 1])
#for time in range(0, 10, 1):
    # longwave up
    for layer_up in range(0, 21, 1):
        if layer_up == 0:
            through = 0
            emitted = sigma * Temp_time[time][layer_up] ** 4  # Wm^-2
            longwave_up[time][layer_up] = through + emitted
        elif 1 <= layer_up <= 20:
            through = (1 - (emissivity_lw[layer_up] + forcing)) * longwave_up[time][layer_up - 1]
            emitted = (forcing + emissivity_lw[layer_up]) * (sigma * Temp_time[time][layer_up] ** 4)
            longwave_up[time][layer_up] = through + emitted
        else:
            print('layering up is not working')

    # longwave down
    for layer_down in range(20, -1, -1):
        if layer_down == 20:
            through = 0
            emitted = (forcing + emissivity_lw[layer_down]) * (sigma * Temp_time[time][layer_down] ** 4)
            longwave_down[time][layer_down] = through + emitted
            shortwave[time][layer_down] = emissivity_sw[layer_down] * S_0 / 4 * (1 - alpha_p)  # absorbed at L20
        elif 1 <= layer_down <= 19:
            through = (1 - (forcing + emissivity_lw[layer_down])) * longwave_down[time][layer_down + 1]
            emitted = (forcing + emissivity_lw[layer_down]) * sigma * Temp_time[time][layer_down] ** 4
            longwave_down[time][layer_down] = through + emitted
            shortwave[time][layer_down] = (S_0 / 4 * (1 - alpha_p) - np.sum(shortwave[time][layer_down + 1:21])) * emissivity_sw[layer_down]
        elif layer_down == 0:
            through = (1 - (forcing+emissivity_lw[layer_down])) * longwave_down[time][layer_down + 1]
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
            captured_shortwave = shortwave[time][layer]  # S_0 / 4 * (1 - alpha_p)
            delta_q[time][layer] = difference_down + difference_up + captured_shortwave

    for temp_change in range(0, 21, 1):  # Gives a temp value to next time step
        if temp_change == 0:
            Temp_time[time + 1][temp_change] = Temp_time[time][temp_change] + ((delta_q[time][temp_change] * 3600) / (rho_w * h_changing * c_w))
        else:
            Temp_time[time + 1][temp_change] = Temp_time[time][temp_change] + ((g * delta_q[time][temp_change] * 3600) / (delta_p * c_p))

#    SH[time] = rho_sh*c_p*C_e*V*(Temp_time[time+1][0]-Temp_time[time+1][1])
#    Temp_time[time+1][0] = Temp_time[time+1][0] - (SH[time] * 3600)/(c_w*rho_w*h_changing)
#    Temp_time[time+1][1] = Temp_time[time+1][1] + (SH[time]*g * 3600)/(delta_p * c_p)

    for actual in range(0, 20, 1):
        temp_at_layers[time + 1][actual] = (Temp_time[time+1][actual] + Temp_time[time+1][actual+1])/2

    for height in range(1, 20, 1):
        temporary = p_layer.copy()
        temporary[20] = 1
        layer_height[time][height] = (((R_d * Temp_time[time + 1][height]) / g) * np.log(temporary[height] / temporary[height + 1])) / 1000  # km
        midpoint[time][height] = (layer_height[time][height] / 2) / 1000 + np.sum(layer_height[time][:height])  # km

    for change in range(1, 19, 1):
        delta_z[time][change] = midpoint[time][change + 1] - midpoint[time][change]  # km

    for convection in range(1, 19, 1):
        temporary = p_layer.copy()
        temporary[20] = 1
        e_s_up = 0.6112*np.exp((17.67 * (Temp_time[time+1][convection+1] - 273.15))/((Temp_time[time+1][convection+1] - 273.15) + 243.5)) * 1000  # Pa
        e_s_lay = 0.6112*np.exp((17.67 * (Temp_time[time+1][convection] - 273.15))/((Temp_time[time+1][convection] - 273.15) + 243.5)) * 1000
        dq_star[time][convection] = ((0.622 * e_s_lay)/(temporary[convection] - e_s_lay)) - ((0.622 * e_s_up)/(temporary[convection] - e_s_up))
        alpha = (L * dq_star[time][convection]) / (c_p * (Temp_time[time + 1][convection] - Temp_time[time + 1][convection + 1]))
        lapse_crit[time][convection] = lapse_dry/(1+alpha)

        if ((Temp_time[time + 1][convection] - Temp_time[time + 1][convection + 1]) / delta_z[time][convection]) < lapse_crit[time][convection]:  # delta_z is the difference in height between midpoints
            adjust_layer[time][convection] = 0
            adjust_upperlayer[time][convection + 1] = 0

        elif ((Temp_time[time + 1][convection] - Temp_time[time + 1][convection + 1]) / delta_z[time][convection]) >= lapse_crit[time][convection]:
            adjust_upperlayer[time][convection + 1] = (Temp_time[time + 1][convection] + Temp_time[time + 1][convection + 1] - lapse_crit[time][convection] * delta_z[time][convection]) / 2
            adjust_layer[time][convection] = adjust_upperlayer[time][convection + 1] + lapse_crit[time][convection] * delta_z[time][convection]
            Temp_time[time+1][convection] = adjust_layer[time][convection]
            Temp_time[time+1][convection + 1] = adjust_upperlayer[time][convection + 1]
            temp_at_layers[time+1][convection] = adjust_layer[time][convection]
            temp_at_layers[time+1][convection+1] = adjust_upperlayer[time][convection + 1]


# at toa focing is difference in equilibrium drop in long wave fluxtwo equilibrium temperature 0.28

OriginalT = Temp_time.copy()
ocean5_noF = Temp_time.copy()  # F = 0, h = 5
ocean5_F = Temp_time.copy()  # F = 0.02, h = 5
ocean1_F = Temp_time.copy()
ocean2_F = Temp_time.copy()
ocean20_F = Temp_time.copy()
ocean50_F = Temp_time.copy()
ocean100_F = Temp_time.copy()
ocean150_F = Temp_time.copy()


plt.plot(range(52561), ocean1_F[:,0]-ocean5_noF[:,0], label='1 m')
plt.plot(range(52561), ocean2_F[:,0]-ocean5_noF[:,0], label='2 m')
plt.plot(range(52561), ocean5_F[:,0]-ocean5_noF[:,0], label='5 m')
plt.plot(range(52561), ocean20_F[:,0]-ocean5_noF[:,0], label='20 m')
plt.plot(range(52561), ocean50_F[:,0]-ocean5_noF[:,0], label='50 m')
plt.plot(range(52561), ocean100_F[:,0]-ocean5_noF[:,0], label='100 m')
plt.plot(range(52561), ocean150_F[:,0]-ocean5_noF[:,0], label='150 m')
plt.legend()
plt.xticks(range(0, 52561, 8760), range(7))
plt.xlabel('Number of Years')
plt.ylabel('Temperature Anamoly \u0394 T, K')
plt.title('Temperature Response at Ocean Depth')
plt.show()



plt.plot(Temp_time[26279].T, p_layer[:21], label='Temp with Feedback')
plt.plot(finalT_prof, p_layer[:21], label='Temp no Feedback')
plt.plot(Temp_time[26279]-finalT_prof, p_layer, label='\u0394 T')
#plt.yscale('log')
plt.title('Temp Difference with and without Feedback')
plt.ylabel('Pressure (Pa)')
plt.xlabel('Temp, K')
plt.legend()
plt.gca().invert_yaxis()
plt.show()

plt.plot(Temp_time[:,0],SH)
plt.show()


plt.plot(Temp_time[:,20],longwave_up[:,20])
plt.show()

delta_z[52559,0] = 5

env_lapse_rate = []
for i in range(0, 19, 1):
    env_lapse_rate.append((Temp_time[52559][i] - Temp_time[52559][i + 1])/delta_z[52559][i])


plt.plot(env_lapse_rate, p_layer[:19], label='Lapse Rate')
plt.yscale('log')
plt.title('Lapse Rate')
plt.ylabel('Pressure, log scale (Pa)')
plt.xlabel('Temp per km, K/km')
plt.legend()
plt.gca().invert_yaxis()
plt.show()
