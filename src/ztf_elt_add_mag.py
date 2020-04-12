import numpy as np
import pandas as pd
import glob
from tqdm.auto import tqdm

# MSE imports
from ztf_data import load_ztf_det_all

ztf, mjd_unq = load_ztf_det_all()

file_paths = glob.glob('../data/ztf_elt/ztf_elt_*.h5')

for file_path in tqdm(file_paths):
    ztf_elt = pd.read_hdf(file_path)
    if 'mag_app' in ztf_elt.columns:
        continue
    loc = ztf_elt.columns.get_loc('uz')+1
    mag_app = ztf.mag_app.loc[ztf_elt.ztf_id].values
    ztf_elt.insert(loc=loc, column='mag_app', value=mag_app)
    ztf_elt.to_hdf(file_path, key='ztf_elt', mode='w')