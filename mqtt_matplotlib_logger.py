import matplotlib.pyplot as plt
plt.install_repl_displayhook()
import numpy as np
import yaqc
import time

daemon_port = 39016
xmin = -60 #seconds
poll_period = 1 #seconds
cached_vals = 200

data_cache = np.zeros((1, 2))
mqtt = yaqc.Client(daemon_port)


plt.ion()
t0  = time.time()
data0 = mqtt.get_measured()['bme680_humidity']
data_cache = np.append(data_cache, np.array([[0, data0]]), axis=0)

fig, ax = plt.subplots(1,1)
p_line, = ax.plot(*data_cache.T, marker='o', c='k', markersize=3)
plt.xlim(xmin, 0)
plt.ylim(0, 40)

t_prior = time.time()
while True:
    t_poll = time.time()
    if t_poll-t_prior>=poll_period:
        data_cache = np.append(data_cache, np.array([[0, mqtt.get_measured()['bme680_humidity']]]), axis=0)
        #print((time.asctime(), t_poll-t0, mqtt.get_measured()['bme680_humidity']))
        data_cache[:, 0] -= t_poll-t_prior
        p_line.set_data(*data_cache.T[int(xmin/poll_period):, :])
        plt.title(str(data_cache[-1][1]), fontsize=40)
        ax.relim()
        ax.autoscale_view()
        plt.draw_if_interactive()
        if len(data_cache)>cached_vals:
            np.delete(data_cache, 0, axis=0)
        t_prior = t_poll
    else:
        plt.pause(poll_period/5)
        time.sleep(poll_period/5)



