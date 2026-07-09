# spectroscopy

Laser-lab code built on WrightTools / attune / yaqc.

Scripts:
- `beam_overlap.py` — interaction volume of crossed Gaussian beams, with a pyvista point-cloud view.
- `FWHM_1D.py` — Gaussian / skewed-Gaussian FWHM fits from `.wt5` scans behind an easygui front end.
- `sig_to_idl.py`, `sig_idl_arr_gen.py` — build an attune idler arrangement from a signal tune. `sig_idl_arr_gen.py` is WIP and won't run as-is (see repo notes).
- `w3_motortune.py` — extend a w3 OPA tuning curve.
- `mqtt_matplotlib_logger.py` — live humidity plot off a yaqc daemon.
- `raspicomm.py` — set SPI / I2C / GPIO device permissions on a Raspberry Pi.
- `liouville_pathways.py` — enumerate time-ordering sign combinations for wave-mixing pathways.

`notebooks/` — FWHM, matrix optics, laser layer propagation, columnated-beam overlap, and image-rotation variance experiments.
