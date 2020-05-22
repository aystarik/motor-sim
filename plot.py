#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
data = np.genfromtxt('PMSM.csv', delimiter=',', names=True)
d2 = np.genfromtxt('control.csv', delimiter=',', names=True)
fig = plt.figure()

ax1 = fig.add_subplot(111)
#ax1.plot(data['t'], data['Id'], color='purple', label='I_d')
#ax1.plot(data['t'], data['Iq'], color='brown', label='I_q')
#ax1.plot(data['t'], data['z'], color='r', label='I_q')
#ax1.plot(data['t'], data['Va'], color='r', label='I_q')
#ax1.plot(data['t'], data['Vb'], color='g', label='I_q')
#ax1.plot(data['t'], data['Vc'], color='b', label='I_q')
ax1.plot(data['t'], data['omega'], color='grey', label='omega')
ax1.plot(data['t'], data['theta'], color='black', label='theta')

ax1.plot(d2['t'], d2['va'], color='purple', label='bemf_a')
ax1.plot(d2['t'], d2['vb'], color='brown', label='bemf_b')
ax1.plot(d2['t'], d2['est_omega'], color='green', label='ia')
ax1.plot(d2['t'], d2['est_angle'], color='red', label='ib')
#ax1.plot(d2['t'], d2['freq'], color='red', label='omega_est')
#ax1.plot(d2['t'], d2['est_angle'], color='blue', label='theta_est')
#ax1.plot(d2['t'], d2['cntrl_angle'], color='grey', label='angle_err')
plt.show()
