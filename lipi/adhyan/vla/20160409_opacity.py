import pickle
import numpy as np
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm


tau_file = "/data/Dropbox/20160409/opacity/tau_2d.pkl"
with open(tau_file, "rb") as f:
    data = pickle.load(f)

x_grid = data["x"]
y_grid = data["y"]
tau_grid = data["tau"]
line_LoS = data["line_LoS"]
line_dotted = data["line_dotted"]
line_solid = data["line_solid"]
line_dashed = data["line_dashed"]


fig, (ax, ax0) = plt.subplots(1, 2)
norm_tau = mcolors.LogNorm(vmin=1e-3, vmax=1e3)
cmap_tau = plt.get_cmap("bwr_r").reversed()
im = ax.pcolormesh(
    x_grid,
    y_grid,
    tau_grid,
    norm=norm_tau,
    cmap=cmap_tau,
    shading="nearest",
    rasterized=True,
)
fig.colorbar(im)
ax.axis("equal")
ax.set_xlabel("X [Mm]")
ax.set_ylabel("Z [Mm]")
ax.set_xlim(60, 100)
ax.set_ylim(10, 70)
ax.plot(line_LoS["x"], line_LoS["y"], "k-")
ax.plot(line_dotted["x"], line_dotted["y"], "k:")
ax.plot(line_solid["x"], line_solid["y"], "k-")
ax.plot(line_dashed["x"], line_dashed["y"], "k--")
ax.plot([75, 95], [35, 48], "g-")
ax0.imshow(
    tau_grid,
    origin="lower",
    aspect="auto",
    cmap="jet",
    norm=LogNorm(vmin=0.01, vmax=10),
)

plt.show()

x0 = np.arange(35, 75, 1)
y0 = np.arange(48, 95, 1)
xid0, xid1 = np.where(np.abs(x_grid - 35) == np.min(np.abs(x_grid - 35)))
yid0, yid1 = np.where(np.abs(y_grid - 48) == np.min(np.abs(y_grid - 48)))
print(xid0[0], xid1[0], yid0[0], yid1[0], tau_grid[31, 246])


plt.imshow(
    tau_grid,
    origin="lower",
    aspect="auto",
    cmap="jet",
    norm=LogNorm(vmin=0.01, vmax=10),
)
plt.show()
