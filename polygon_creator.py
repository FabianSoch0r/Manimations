import matplotlib.pyplot as plt
import numpy as np
import datetime

# --- Settings ---
GRID_STEP = 1.0      # grid spacing in data units (1 = one Manim unit)
SNAP = True           # snap to grid vertices
SHOW_GRID = True      # draw the grid
coords = []

fig, ax = plt.subplots(figsize=(10, 10))
ax.set_title("Left-click to add points, right-click to finish")

# Add point in middle for reference
ax.plot(0,0,'.',markersize=8,c='blue')

# Define visible range (centered around zero)
x_range = 20
y_range = 20
ax.set_xlim(-x_range / 2, x_range / 2)
ax.set_ylim(-y_range / 2, y_range / 2)
ax.set_aspect('equal')

# Draw integer grid and hide labels
if SHOW_GRID:
    ax.set_xticks(np.arange(-x_range / 2, x_range / 2 + GRID_STEP, GRID_STEP))
    ax.set_yticks(np.arange(-y_range / 2, y_range / 2 + GRID_STEP, GRID_STEP))
    ax.grid(True, linestyle='--', color='0.8')
    #ax.set_xticklabels([])  # remove digits
    #ax.set_yticklabels([])

line, = ax.plot([], [], 'r-', lw=2)
points, = ax.plot([], [], 'ro')

def snap_to_grid(x, y):
    gx = round(x / GRID_STEP) * GRID_STEP
    gy = round(y / GRID_STEP) * GRID_STEP
    return gx, gy

def onclick(event):
    if event.inaxes != ax:
        return

    if event.button == 1:  # left click
        x, y = event.xdata, event.ydata
        if SNAP:
            x, y = snap_to_grid(x, y)
        coords.append((x, y))
        update_plot()
    elif event.button == 3:  # right click
        save_coords()
        plt.close()

def save_coords():
    if coords:
        timestamp = datetime.datetime.now().strftime("%d%m%Y_%H%M%S")
        filename = f"{timestamp}.txt"
        np.savetxt(filename, np.array(coords), fmt="%.6f")
        print(f"Saved {len(coords)} coordinates to {filename}")
        print("Coordinates:", coords)

def update_plot():
    if len(coords) > 0:
        xs, ys = zip(*coords)
        line.set_data(xs, ys)  # close polygon visually
        points.set_data(xs, ys)
        fig.canvas.draw_idle()

fig.canvas.mpl_connect('button_press_event', onclick)
plt.show()
