#!/usr/bin/env python

import argparse
import glob
import meshio
import numpy as np
import os
import pyvista as pv

parser = argparse.ArgumentParser()
parser.add_argument("--ndims", type=int, default=3, help="number of spatial dimensions", required=True)
parser.add_argument("--xvar_in", type=str, help="name of x-variable (as in mef-file)", required=True)
parser.add_argument("--yvar_in", type=str, help="name of x-variable (as in mef-file)", required=True)
parser.add_argument("--xvar_out", type=str, help="name of x-variable in out-jPDF", required=True)
parser.add_argument("--yvar_out", type=str, help="name of x-variable in out-jPDF", required=True)
parser.add_argument("-n", "--nbins", type=int, help="number of bins", required=True)
parser.add_argument("-i", "--infiles", metavar="PATTERN", help="mef infile names pattern (remember quotes for patterns, e.g. --infiles \"plt*.mef)", required=True)
args = parser.parse_args()
NDIM = args.ndims
XVAR_IN = args.xvar_in
YVAR_IN = args.yvar_in
XVAR_OUT = args.xvar_out
YVAR_OUT = args.yvar_out
NBINS = args.nbins
INFILES_PATTERN = args.infiles

def shape_from_header(h):
    start, stop, _, nfields = h.split()[-4:]
    nfields = int(nfields)
    start = np.array(start.split('(')[-1].replace(')','').split(','), dtype=int)
    stop = np.array(stop.replace('(', '').replace(')', '').split(','), dtype=int)
    shape = stop - start + 1
    return [shape[0], shape[1], shape[2], nfields]

def read_mef(filename):
    with open(filename, 'rb') as mef:
        time = mef.readline().decode('ascii')
        fields = mef.readline().decode('ascii').split()
        nfaces, verts_per_face = mef.readline().decode('ascii').split()
        nfaces = int(nfaces)
        verts_per_face = int(verts_per_face)

        header = mef.readline().decode('ascii')
        bin_data_shape = shape_from_header(header)
        bin_data_shape = (bin_data_shape[0], bin_data_shape[3])

        data = np.fromfile(mef, 'float64', np.prod(bin_data_shape))
        data = data.reshape(bin_data_shape, order='C')

        faces = np.fromfile(mef, 'int32', nfaces * verts_per_face)
        faces = faces.reshape((nfaces, verts_per_face), order='C') - 1

    ndims = NDIMS
    vertices = data[:, :ndims]
    point_data = {f: data[:, ndims + i] for i, f in enumerate(fields[ndims:])}

    m = meshio.Mesh(
        points=vertices,
        cells={'triangle': faces},
        point_data=point_data
    )
    mesh = pv.PolyData(m.points, np.hstack([np.full((len(faces), 1), 3), faces]))
    for key, val in m.point_data.items():
        mesh.point_data[key] = val

    return mesh, faces

mef_files = glob.glob(INFILES_PATTERN)

# Pass 1: find global min/max
xmin, xmax = np.inf, -np.inf   # global x min/max
ymin, ymax = np.inf, -np.inf   # global y min/max
for f in mef_files:
    mesh, faces = read_mef(f)
    mesh = mesh.triangulate()
    for key in list(mesh.point_data.keys()):
        mesh.cell_data[key + "_center"] = mesh.point_data[key][faces].mean(axis=1)
    xmin = min(xmin, mesh.cell_data[f"{XVAR_IN}_center"].min())
    xmax = max(xmax, mesh.cell_data[f"{XVAR_IN}_center"].max())
    ymin = min(ymin, mesh.cell_data[f"{YVAR_IN}_center"].min())
    ymax = max(ymax, mesh.cell_data[f"{YVAR_IN}_center"].max())

# Pass 2: accumulate histogram with fixed range
hist2d_total = None
for f in mef_files:
    print(f"[{f}]")
    mesh, faces = read_mef(f)
    mesh = mesh.triangulate()
    mesh = mesh.compute_cell_sizes(length=False, area=True, volume=False)
    areas = mesh.cell_data["Area"]
    for key in list(mesh.point_data.keys()):
        v = mesh.point_data[key]
        cell_values = v[faces].mean(axis=1)
        mesh.cell_data[key + "_center"] = cell_values
    x = mesh.cell_data[f"{XVAR_IN}_center"]
    y = mesh.cell_data[f"{YVAR_IN}_center"]
    hist2d, x_edges, y_edges = np.histogram2d(
        x, y,
        bins=NBINS,
        range=[[xmin, xmax], [ymin, ymax]],
        weights=areas
    )
    if hist2d_total is None:
        hist2d_total = hist2d
    else:
        hist2d_total += hist2d

hist2d_total /= hist2d_total.sum()
x_centers = 0.5 * (x_edges[:-1] + x_edges[1:])
y_centers = 0.5 * (y_edges[:-1] + y_edges[1:])

dirname = f"MEF_JPDFAverage_{XVAR_OUT}_{YVAR_OUT}"
os.makedirs(dirname, exist_ok=True)
x_filename = f"{dirname}/Pdf_{XVAR_OUT}_x.dat"
y_filename = f"{dirname}/Pdf_{YVAR_OUT}_x.dat"
jpdf_filename = f"{dirname}/Pdf_{XVAR_OUT}_{YVAR_OUT}.dat"
print(f"Saving {XVAR_OUT} to {x_filename}.")
np.savetxt(x_filename, x_centers)
print(f"Saving {YVAR_OUT} to {y_filename}.")
np.savetxt(y_filename, y_centers)
print(f"Saving jPDF to {jpdf_filename}.")
np.savetxt(jpdf_filename, hist2d_total)
print("Done.")
