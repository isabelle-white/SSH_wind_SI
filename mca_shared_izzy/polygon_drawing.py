import matplotlib.pyplot as plt
from matplotlib.widgets import PolygonSelector
import numpy as np

class PolygonExtentDrawer:
    def __init__(self, ax):
        self.ax = ax
        self.selector = PolygonSelector(ax, self.onselect)
        self.verts = []

    def onselect(self, verts):
        self.verts = verts
        xs, ys = zip(*verts)
        extent = {
            "lon_min": min(xs),
            "lon_max": max(xs),
            "lat_min": min(ys),
            "lat_max": max(ys)
        }
        print("\n📦 Polygon Extents (Bounding Box):")
        print(f"Longitude: {extent['lon_min']:.2f} to {extent['lon_max']:.2f}")
        print(f"Latitude:  {extent['lat_min']:.2f} to {extent['lat_max']:.2f}")
        self.draw_polygon()

    def draw_polygon(self):
        x, y = zip(*self.verts)
        self.ax.plot(x + (x[0],), y + (y[0],), 'r--', lw=2)
        plt.draw()

# ---------------------------------------------
# Example setup: dummy lat/lon grid
# ---------------------------------------------
lon, lat = np.meshgrid(np.linspace(-180, 180, 360), np.linspace(-90, -60, 180))
data = np.random.rand(*lon.shape)

fig, ax = plt.subplots(figsize=(10, 5))
cs = ax.pcolormesh(lon, lat, data, cmap="coolwarm")
plt.colorbar(cs, ax=ax, orientation="vertical", label="Dummy Value")

ax.set_title("Draw a polygon (double-click to close)\nBounding box will be printed")
drawer = PolygonExtentDrawer(ax)

plt.show()

