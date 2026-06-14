"""
Example of usage:
python3 xdr_plotting.py lso_comp-s_lev4.0_20221102_090718_wave0656_pa023_obse_maparr.xdr
"""

import argparse
import datetime
import numpy as np
import astropy.units as u
from astropy.io import fits
from scipy.io import readsav
import sunpy.map 
import matplotlib.pyplot as plt

def rot(deg):
    # instrument frame -> helioprojective 
    a = np.deg2rad(deg)
    return np.array([[np.cos(a), -np.sin(a)],
                     [np.sin(a),  np.cos(a)]])


def build_header(rec):
    # build a 2-axis helioprojective WCS header
    d = rec['data']
    ny, nx = d.shape
    xc, yc = float(rec['xc']), float(rec['yc']) # arcsec, instrument frame
    dx, dy = float(rec['dx']), float(rec['dy']) # arcsec/pixel
    roll = float(rec['roll_angle']) # deg, FOV tilt to solar north
    R = rot(roll)
    t = datetime.datetime.strptime(rec['time'].decode().strip(), '%d-%b-%Y %H:%M:%S.%f')

    hdr = fits.Header()
    hdr['WCSAXES'] = 2
    hdr['CTYPE1'], hdr['CTYPE2'] = 'HPLN-TAN', 'HPLT-TAN'
    hdr['CUNIT1'], hdr['CUNIT2'] = 'deg', 'deg'
    hdr['CRPIX1'], hdr['CRPIX2'] = (nx + 1) / 2, (ny + 1) / 2
    hdr['CRVAL1'], hdr['CRVAL2'] = (R @ np.array([xc, yc])) / 3600.0
    hdr['CDELT1'], hdr['CDELT2'] = dx / 3600.0, dy / 3600.0
    hdr['PC1_1'], hdr['PC1_2'] = R[0]
    hdr['PC2_1'], hdr['PC2_2'] = R[1]
    hdr['RSUN_OBS'] = float(rec['rsun'])
    hdr['HGLT_OBS'] = float(rec['b0'])
    hdr['DATE-OBS'] = t.isoformat()
    return hdr


def make_map(xdr_path, frame=2):
    # load an XDR map array and return a north-up sunpy Map for one frame
    sav = readsav(xdr_path)
    maparr = sav['maparr']
    rec = maparr[frame]
    hdr = build_header(rec)
    return sunpy.map.Map(rec['data'], hdr).rotate(order=3, missing=0)


def main():
    p = argparse.ArgumentParser(description="Plot a CoMP-S XDR prominence frame.")
    p.add_argument('xdr', help="path to the ..._maparr.xdr file")
    p.add_argument('--frame', type=int, default=2, help="frame index")
    p.add_argument('--title', default='HI 656 nm Stokes I')
    p.add_argument('--save', default=None, help="output image path; shows if omitted")
    args = p.parse_args()

    m = make_map(args.xdr, frame=args.frame)

    fig = plt.figure(figsize=(8, 7))
    ax = fig.add_subplot(projection=m)
    m.plot(axes=ax, cmap='afmhot', clip_interval=(1, 99.5) * u.percent)
    m.draw_limb(axes=ax, color='white')
    m.draw_grid(axes=ax)
    plt.title(args.title)

    if args.save:
        fig.savefig(args.save, dpi=300, bbox_inches='tight')
        print(f"saved {args.save}")
    else:
        plt.show()


if __name__ == '__main__':
    main()