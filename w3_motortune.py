# attune store had bad w3 curve that only extended to ~1300 nm: try to expand it

import WrightTools as wt
import matplotlib.pyplot as plt
import pathlib
import attune
import numpy as np

apply = False
port = 39302

instrument_name = "w3"
arrangement = "sig"


bbo_folder = pathlib.Path(r"C:\\Users\\john\\bluesky-cmds-data\\2023-08-04 55418 motortune opa3 bbo motortune 3db48959")
grating_folder = pathlib.Path(r"C:\\Users\\john\\bluesky-cmds-data\\2023-08-04 55809 motortune opa3 grating motortune c3ab4a86")

d_bbo = wt.open(bbo_folder / "primary.wt5")
d_grating = wt.open(grating_folder/"primary.wt5")

d_bbo.print_tree()
d_grating.print_tree()

d_bbo.bring_to_front("daq_pyro_3")
#d_bbo.level(0, 0, 10)
#d_bbo.smooth(2)
d_grating.bring_to_front("daq_pyro_3")
#d_grating.level(0, 0, 10)
#d_grating.smooth(2)


instrument = attune.load(instrument_name)
fig, gs = wt.artists.create_figure()
ax = plt.subplot(gs[0])
ax.plot(instrument[arrangement]["bbo"].dependent, instrument[arrangement]["grating"].dependent)
ax.scatter(instrument[arrangement]["bbo"](1300), instrument[arrangement]["grating"](1300), color="red")
# ax.plot(new_instrument[arrangement]["grating"].dependent, new_instrument[arrangement]["bbo"].dependent+1)
ax.set_xlabel("grating")
ax.set_ylabel("bbo")
plt.grid()
plt.show()


'''



# now make a crude curve (to improve upon)
independent = np.arange(1100, 1620, 20)
grating = attune.Tune(independent, np.linspace(34.8, 42, independent.size))
bbo_dependent = []
for gi in grating.dependent:
    di = d.at(w3_grating_points=[gi, "mm"])
    bbo_i = di.w3_bbo_points[di.daq_pyro_3[:].argmax()]
    bbo_dependent.append(bbo_i)

bbo = attune.Tune(independent, bbo_dependent)

sig = attune.Arrangement("sig", dict(grating=grating, bbo=bbo))

new_instrument = instrument.as_dict()
new_instrument["arrangements"]["sig"] = sig.as_dict()
new_instrument = attune.Instrument(**new_instrument)

ax.scatter(new_instrument["sig"]["grating"].dependent, new_instrument["sig"]["bbo"].dependent, color="k")

if apply:
    wt.artists.savefig(folder / "custom_tune_extraction.png")
    attune.store(new_instrument)
else:
    plt.show()

'''
