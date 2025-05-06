# polar_plot_utils.py

import numpy as np
from matplotlib.path import Path
from matplotlib.patches import PathPatch
import matplotlib.pyplot as plt
from matplotlib import colors
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import cartopy.crs as ccrs

def plot_curvilinear_data(ax, data, lon, lat, title, vmin=None, vmax=None):
    norm = colors.TwoSlopeNorm(vmin=vmin, vcenter=0, vmax=vmax)

    pc = ax.pcolormesh(lon, lat, data,
                       transform=ccrs.PlateCarree(),
                       cmap='coolwarm',
                       shading='auto',
                       norm=norm)

    ax.set_extent([-180, 180, -90, -55], crs=ccrs.PlateCarree())
    ax.coastlines()
    ax.set_title(title)

    theta = np.linspace(0, 2 * np.pi, 100)
    circle = np.vstack([0.5 + 0.5 * np.cos(theta), 0.5 + 0.5 * np.sin(theta)]).T
    path = Path(circle)
    patch = PathPatch(path, transform=ax.transAxes, facecolor='none')
    ax.set_boundary(path, transform=ax.transAxes)

    gl = ax.gridlines(draw_labels=True,
                      x_inline=False, y_inline=False,
                      linestyle='--', linewidth=0.5, color='gray')
    gl.xlabels_top = gl.xlabels_right = False
    gl.xlocator = plt.FixedLocator(np.arange(-180, 181, 45))
    gl.ylocator = plt.FixedLocator(np.arange(-90, -54, 5))
    gl.xlabel_style = {'size': 8}
    gl.ylabel_style = {'size': 8}
    gl.longitude_formatter = LONGITUDE_FORMATTER
    gl.latitude_formatter = LATITUDE_FORMATTER

    cbar = plt.colorbar(pc, ax=ax, orientation='vertical', shrink=0.7, pad=0.05)

    return pc
