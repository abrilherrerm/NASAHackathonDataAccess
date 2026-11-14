"""Genera datos SINTÉTICOS con el mismo formato que process.py, para probar la API
sin cuenta ni descarga. Se borran solos al correr process.py (mismos nombres)."""
import numpy as np
import pandas as pd
import xarray as xr

import config
from process import mexico_mask, save_yearly

b = config.MEXICO_BBOX
lat = np.arange(b["south"] + 0.05, b["north"], 0.1)
lon = np.arange(b["west"] + 0.05, b["east"], 0.1)
time = pd.date_range(config.START_DATE, config.END_DATE, freq="D")
rng = np.random.default_rng(0)
season = 1 + 2.5 * np.exp(-((time.day_of_year.to_numpy() - 210) / 50) ** 2)  # temporada de lluvias
wet = (1 - (lat - lat.min()) / np.ptp(lat))[:, None] * 1.5 + 0.2   # más lluvia al sur
data = rng.gamma(0.4, 1, (time.size, lat.size, lon.size)).astype("float32")
data *= season[:, None, None].astype("float32") * wet[None].astype("float32") * 6
da = xr.DataArray(data, coords={"time": time, "lat": lat, "lon": lon}, dims=("time", "lat", "lon"), name="precipitation")
save_yearly(mexico_mask(da))
