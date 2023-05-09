import h5py
import numpy as np
from tqdm import tqdm

f = h5py.File("data/19F73_BiTiOx_all_oned.h5", 'r')


conds = list(f['exp'])

for cond in tqdm(conds):
    ion = f['exp'][cond].attrs['cations'][0]
    frac = f['exp'][cond].attrs['fracs'][0][0]

    q = f['exp'][cond]['0']['integrated_1d'][0]

    filename = f"{cond}_{ion}_{frac:.3f}_q"
    np.save(filename, q)


