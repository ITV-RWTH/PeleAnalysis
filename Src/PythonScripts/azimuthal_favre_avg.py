import yt
import numpy as np
from scipy.interpolate import RegularGridInterpolator
import os

# NOTE: *_cart -> cartesian, *_cyl -> cylindrical

N_THETA = 256
N_R = 512
DIR = "circ_average"
os.makedirs(DIR, exist_ok=True)
VAR_NAMES = [
    "rho_mean", 
    "rho_r_velocity_mean",
    "rho_t_velocity_mean",
    "rho_z_velocity_mean",
    "rho_z_velocity_times_Y(H2)_mean"
]
        

ds = yt.load("favreAvg_plt_all")
ad = ds.covering_grid(level=0, left_edge=ds.domain_left_edge, dims=ds.domain_dimensions)

# cart = [x, y, z]
dims_cart = ds.domain_dimensions
lo_cart = ds.domain_left_edge.d
hi_cart = ds.domain_right_edge.d
length_cart = hi_cart - lo_cart
dx = length_cart / dims_cart
x = lo_cart[0] + (np.arange(dims_cart[0]) + 0.5) * dx[0]   # 1D
y = lo_cart[1] + (np.arange(dims_cart[1]) + 0.5) * dx[1]   # 1D
z = lo_cart[2] + (np.arange(dims_cart[2]) + 0.5) * dx[2]   # 1D

# cyl = [r, theta, z]
r_max = min(abs(lo_cart[0]), abs(hi_cart[0]), abs(lo_cart[1]), abs(hi_cart[1]))   # max radius of cylinder still inside box
r = np.linspace(0, r_max, N_R)               # 1D
theta = np.linspace(0, 2 * np.pi, N_THETA)   # 1D
R, Theta = np.meshgrid(r, theta, indexing="ij")
X_cyl = R * np.cos(Theta)            # 2D
Y_cyl = R * np.sin(Theta)            # 2D

# INTERPOLATION:

rho_u_cart = ad[("boxlib", "rho_x_velocity_mean")]
u_interp = RegularGridInterpolator(
    (x, y, z),
    rho_u_cart,
    bounds_error=False,
    fill_value=None
)
rho_v_cart = ad[("boxlib", "rho_y_velocity_mean")]
v_interp = RegularGridInterpolator(
    (x, y, z),
    rho_v_cart,
    bounds_error=False,
    fill_value=None
)
for var_name in VAR_NAMES:
    is_cyl = var_name in ["rho_r_velocity_mean", "rho_t_velocity_mean"]
    if not is_cyl:
        vals_cart= ad[("boxlib", var_name)].d
        interp = RegularGridInterpolator(
            (x, y, z),
            vals_cart,
            bounds_error=False,
            fill_value=None
        )
    vals_mean = np.zeros((len(r), len(z)))
    
    # Do interpolation for one z-slice at a time:
    for i, z_val in enumerate(z):
        print(f"Var: {var_name} - Interp+circAvg at z = {z_val:.5f} ({i+1} | {len(z)})", end="\r")
        Z_cyl = np.full_like(X_cyl, z_val)   # 2D
        pts_cyl = np.column_stack((X_cyl.ravel(), Y_cyl.ravel(), Z_cyl.ravel()))
        if (var_name == "rho_r_velocity_mean"):
            u_cyl = u_interp(pts_cyl).reshape(R.shape)
            v_cyl = v_interp(pts_cyl).reshape(R.shape)
            vals_cyl = u_cyl * np.cos(Theta) + v_cyl * np.sin(Theta)   # radial vel
        elif (var_name == "rho_t_velocity_mean"):
            u_cyl = u_interp(pts_cyl).reshape(R.shape)
            v_cyl = v_interp(pts_cyl).reshape(R.shape)
            vals_cyl = -u_cyl * np.sin(Theta) + v_cyl * np.cos(Theta)   # tangential vel
        else:
            vals_cyl = interp(pts_cyl)
        vals_cyl = vals_cyl.reshape(R.shape)
        vals_mean[:, i] = vals_cyl.mean(axis=1)
    
    # SAVE:
    np.save(f"{DIR}/{var_name}.npy", vals_mean)
    print("\n")
np.save(f"{DIR}/r.npy", r)
np.save(f"{DIR}/z.npy", z)
