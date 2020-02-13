import sys
import numpy as np
import h5py
from scipy.signal import savgol_filter, find_peaks
import laspy

from matplotlib.cm import magma_r
cmap = magma_r


# crude LAS file generation for inspection purposes

def exportlas(fn, var, pts, cmap):
    v = var - np.min(var)
    v /= v.max()
    rgb = cmap(v)
    rgb = rgb[:, :3]
    rgb *= 65535
    rgb = rgb.astype('uint')
    header = laspy.header.Header()
    header.data_format_id = 2
    f = laspy.file.File(fn, mode = 'w', header = header)
    f.header.scale = [0.001, 0.001, 0.001]
    f.header.offset = [pts[:,0].min(), pts[:,1].min(), pts[:,2].min()]
    f.x = pts[:, 0]
    f.y = pts[:, 1]
    f.z = pts[:, 2]
    if pts.shape[1] == 4:
        f.intensity = pts[:, 3]
    f.set_red(rgb[:, 0])
    f.set_green(rgb[:, 1])
    f.set_blue(rgb[:, 2])
    f.close()

fname = sys.argv[1]
f = h5py.File(fname, 'r')
m = f['Amplitude'].shape[0]
chunk = 100000000
k = 0
for i in range(0, m, chunk):
    k += 1
    sys.stdout.write('\b'*10+' %0.1f%%' % (100 * i/m))
    sys.stdout.flush()

    amp = f['Amplitude'][i:i+chunk]
    xyz = f['XYZ'][i:i+chunk]
    xyz = xyz[amp[:, 1] > 0, :]
    amp = amp[amp[:, 1] > 0, 0]
    mag = np.log10(amp - amp.min() + 1)
    sel, _ = find_peaks(savgol_filter(mag, 13, 2), height = np.percentile(mag, 50))
    amp = amp[sel]
    xyz = xyz[sel]
    c = np.log10(amp)
    c[c < np.percentile(c, 1)] = np.percentile(c, 1)
    c[c > np.percentile(c,99)] = np.percentile(c,99)
    exportlas('%s-%04d.las' % (fname[:-4], k), c, np.c_[xyz, amp], cmap)

sys.stdout.write('\b'*10+' %0.1f%%\n' % 100)
