
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

emissivity_lw = [1, 0.3, 0.3, 0.3, 0.3, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1]
emissivity_sw = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0015625, 0.003125, 0.00625, 0.0125, 0.025, 0.05]

Temp_time = np.zeros((17473, 21))
Temp_time[0][:] = 50
longwave_up = np.zeros((17473, 21))
longwave_down = np.zeros((17473, 21))
shortwave = np.zeros((17473, 21))
delta_q = np.zeros((17473, 21))
adjust_layer = np.zeros((17473, 21))
adjust_upperlayer = np.zeros((17473, 20))
layer_height = np.zeros((17473, 20))
midpoint = np.zeros((17473, 20))
delta_z = np.zeros((17473, 19))
lapse_crit = np.zeros((17473, 19))
dq_star = np.zeros((17473, 21))
temp_changes = np.zeros((17473, 21))
temp_at_layers = np.zeros((17473, 20))

for time in range(0, 17472, 1):
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
            through = (1 - emissivity_lw[layer_up]) * longwave_up[time][layer_up - 1]
            emitted = emissivity_lw[layer_up] * (sigma * Temp_time[time][layer_up] ** 4)
            longwave_up[time][layer_up] = through + emitted
        else:
            print('layering up is not working')

    # longwave down
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
            captured_shortwave = shortwave[time][layer]  # S_0 / 4 * (1 - alpha_p)
            delta_q[time][layer] = difference_down + difference_up + captured_shortwave

    for temp_change in range(0, 21, 1):  # Gives a temp value to next time step
        Temp_time[time + 1][temp_change] = Temp_time[time][temp_change] + (g / (delta_p * c_p)) * delta_q[time][temp_change] * 3600

    for actual in range(0, 20, 1):
        temp_at_layers[time + 1][actual] = (Temp_time[time+1][actual] + Temp_time[time+1][actual+1])/2

    for height in range(0, 20, 1):
        temporary = p_layer.copy()
        temporary[20] = 1
        layer_height[time][height] = (((R_d * Temp_time[time + 1][height]) / g) * np.log(temporary[height] / temporary[height + 1])) / 1000  # km
        midpoint[time][height] = (layer_height[time][height] / 2) / 1000 + np.sum(layer_height[time][:height])  # km

    for change in range(0, 19, 1):
        delta_z[time][change] = midpoint[time][change + 1] - midpoint[time][change]  # km

    for convection in range(0, 19, 1):
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




new_layer_height = np.zeros(20)
new_midpoint = np.zeros(20)
for height in range(0, 20, 1):
    temporary = p_layer.copy()
    temporary[20] = 1
    new_layer_height[height] = (((R_d * temp_at_layers[17471][height]) / g) * np.log(temporary[height] / temporary[height + 1])) / 1000  # km
    new_midpoint[height] = (new_layer_height[height] / 2) / 1000 + np.sum(new_layer_height[:height])  # km

new_delta_z = np.zeros(19)
for diff in range(19):
    new_delta_z[diff] = new_midpoint[diff + 1] - new_midpoint[diff]




# delta_z too large?

plt.plot(p_layer[:20], midpoint[17471])
plt.show()

strat_temp = [336.93376706, 319.48672169, 313.65003006, 307.46800462, 300.88884802, 295.18976712, 290.54842602, 285.67335107, 280.53508471, 275.09787867, 269.31770639, 263.13940603, 256.49244651, 251.24195937, 247.68451125, 244.45003439, 241.07938386, 238.08584361, 236.08643481, 236.39607759, 241.72042746]

plt.plot(Temp_time[17471], p_layer[:21], label='Temp with Convection')
#plt.plot(strat_temp, p_layer, label='Temp w/o Convection')
#plt.yscale('log')
plt.title('Temperature Profiles')
plt.ylabel('Pressure (Pa)')
plt.xlabel('Temp, K')
plt.legend()
plt.gca().invert_yaxis()
plt.show()

env_lapse_rate = []
for i in range(0, 19, 1):
    env_lapse_rate.append((Temp_time[17471][i] - Temp_time[17471][i + 1])/delta_z[17471][i])


plt.plot(env_lapse_rate, p_layer[:19], label='Lapse Rate w/Convection')
plt.plot(strat_lapse, p_layer[:20], label='Lapse Rate w/o Convection')
plt.yscale('log')
plt.title('Lapse Rate')
plt.ylabel('Pressure, log scale (Pa)')
plt.xlabel('Temp, K')
plt.legend()
plt.gca().invert_yaxis()
plt.show()

adj_env_lapse = []
blehh_env_lapse = []
for i in range(0, 19, 1):
    adj_env_lapse.append((temp_at_layers[17471][i] - temp_at_layers[17471][i + 1])/delta_z[17471][i])
    blehh_env_lapse.append((Temp_time[17471][i] - Temp_time[17471][i + 1]) / delta_z[17471][i])

np.subtract(blehh_env_lapse, adj_env_lapse)

temp_changes[17471, :19]/delta_z[17471, :19]
plt.plot(temp_changes[17471, :19]/delta_z[17471], p_layer[:19])
plt.yscale('log')
plt.gca().invert_yaxis()
plt.show()


plt.plot(env_lapse_rate, p_layer[:19], label='6.5')
plt.plot(adj_env_lapse, p_layer[:19], label='adjusted')
plt.yscale('log')
plt.title('Fixed vs Adjusted')
plt.ylabel('Pressure, log scale (Pa)')
plt.xlabel('Temp, K')
#plt.legend()
plt.gca().invert_yaxis()
plt.show()


def theta(level):
    return Temp_time[17471][level] * (p_layer[0] / p_layer[level]) ** 0.286


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

array = np.array([2, 2, 2, 2, 2, 2])

for i in range(5):
    temp = array.copy()
    temp[3] = 3
    print(temp[i])





'''        e_s_up = 0.6112*np.exp((17.67 * (Temp_time[time+1][convection+1] - 273.15))/((Temp_time[time+1][convection+1] - 273.15) + 243.5)) * 1000  # Pa
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
'''


















