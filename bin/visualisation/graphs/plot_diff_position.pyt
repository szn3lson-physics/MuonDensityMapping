import numpy as np
import matplotlib.pyplot as plt
from PIL import Image
import matplotlib.patheffects as patheffects
import matplotlib.patches as patches
import sys
import os

# ================= VISUALIZATION PARAMETERS =================
IMAGE_PATH = 'converted_black_white.png'
PIXEL_SIZE = 0.13333  

DETECTOR_X = 915      
DETECTOR_Y = 1116     

MAX_RADIUS_M = 40.0   
VIEW_RADIUS_M = 100.0 

# ---- HISTOGRAM CONTROL ----
PLOT_BIN_SIZE = 15    # Bin width for Matplotlib plot

# ---- C++ DATA FILES ----
FILE1_DETAILS = "D:\\MuonDensityMapping\\bin\\simulation\\simulation_data\\sim_physics_C.txt"
FILE2_MC_SUMMARY = "D:\\MuonDensityMapping\\bin\\simulation\\simulation_data\\sim_summary_C.txt"

# ================= STEP 1: LOAD MAP =================
print("1. Loading cave map...")
try:
    img = Image.open(IMAGE_PATH).convert('L')
    cave_map = np.array(img)
except Exception as e:
    print(f"[WARNING] Cannot load {IMAGE_PATH}: {e}")
    cave_map = None

# ================= STEP 2: LOAD C++ DATA =================
print("2. Loading data from text files...")

if not os.path.exists(FILE1_DETAILS) or not os.path.exists(FILE2_MC_SUMMARY):
    print("[ERROR] Missing text files. Make sure the C++ program has been executed.")
    sys.exit(1)

data_physics = np.loadtxt(FILE1_DETAILS, skiprows=1)
fwd_rock_3d = data_physics[:, 3]
bwd_rock_3d = data_physics[:, 6]

data_mc = np.loadtxt(FILE2_MC_SUMMARY, skiprows=1)
mc_counts = data_mc[:, 1]

# ================= STEP 3: DATA PROCESSING =================
print("3. Processing data for plots...")

front_empty = MAX_RADIUS_M - fwd_rock_3d
back_empty = MAX_RADIUS_M - bwd_rock_3d
total_empty = front_empty + back_empty

num_plot_bins = 180 // PLOT_BIN_SIZE
binned_flux_plot = np.zeros(num_plot_bins)
bin_centers = np.arange(0, 180, PLOT_BIN_SIZE) + (PLOT_BIN_SIZE / 2.0) - 0.5

for i in range(180):
    binned_flux_plot[i // PLOT_BIN_SIZE] += mc_counts[i]

# ================= STEP 4: DRAWING MATPLOTLIB =================
print("4. Generating plots...")

center = (DETECTOR_Y, DETECTOR_X)
max_radius_px = int(MAX_RADIUS_M / PIXEL_SIZE)
view_radius_px = int(VIEW_RADIUS_M / PIXEL_SIZE)
angles_180 = np.arange(180)

# Helpers for Plot 1
def draw_cave(ax):
    if cave_map is not None:
        im = ax.imshow(cave_map, cmap='gray')
        circle_mask = patches.Circle((center[1], center[0]), view_radius_px, transform=ax.transData)
        im.set_clip_path(circle_mask)
    ax.add_patch(plt.Circle((center[1], center[0]), view_radius_px, color='black', fill=False, linewidth=2, zorder=12))

def draw_infrastructure(ax):
    compass = {0: 'N', 45: 'NE', 90: 'E', 135: 'SE', 180: 'S', 225: 'SW', 270: 'W', 315: 'NW'}
    angles_to_draw = sorted(list(set(range(0, 360, 10)) | set(compass.keys())))
    
    text_radius_px = view_radius_px * 1.15 

    for a in angles_to_draw:
        rad = np.radians(a)
        dx, dy = np.sin(rad), -np.cos(rad)
        is_main = a in compass
        
        ax.plot([center[1], center[1] + dx * view_radius_px],
                [center[0], center[0] + dy * view_radius_px],
                color='red', linestyle='-' if is_main else ':', alpha=0.35 if is_main else 0.12, zorder=2)
        
        label_text = f'{a}° ({compass[a]})' if is_main else f'{a}°'
        txt = ax.text(center[1] + dx * text_radius_px, center[0] + dy * text_radius_px, 
                      label_text, color='black', fontsize=10 if is_main else 8, 
                      ha='center', va='center', fontweight='bold' if is_main else 'normal', zorder=10)
        txt.set_path_effects([patheffects.withStroke(linewidth=1.5, foreground='white')])

    ax.add_patch(plt.Circle((center[1], center[0]), max_radius_px, color='yellow', fill=False, linewidth=2, zorder=10))
    ax.plot(center[1], center[0], 'ro', markersize=8, zorder=20)


# Create figure
fig, axs = plt.subplots(1, 3, figsize=(22, 7))

# ---- Plot 1: Map (Circular Frame with Compass) ----
draw_cave(axs[0])
draw_infrastructure(axs[0])

margin = view_radius_px * 0.25 
lim_x = (center[1] - view_radius_px - margin, center[1] + view_radius_px + margin)
lim_y = (center[0] + view_radius_px + margin, center[0] - view_radius_px - margin)

axs[0].set_xlim(lim_x)
axs[0].set_ylim(lim_y)
axs[0].set_aspect('equal')
axs[0].axis('off')

# ZMIANA TUTAJ: pad=5 zamiast pad=20
axs[0].set_title(f'View around detector (Radius = {VIEW_RADIUS_M}m)', pad=5, fontsize=12)

# ---- Plot 2: Empty space ----
axs[1].plot(angles_180, front_empty, color='blue', linewidth=2, linestyle='--', alpha=0.7, label='Front direction')
axs[1].plot(angles_180, back_empty, color='red', linewidth=2, linestyle='--', alpha=0.7, label='Back direction (+180°)')
axs[1].plot(angles_180, total_empty, color='orange', linewidth=2.5, linestyle='-', label='TOTAL (Front + Back)')

axs[1].set_title(f'Empty space range (Max Radius = {MAX_RADIUS_M}m)', fontsize=12)
axs[1].set_xlabel('Indicated angle [degrees]')
axs[1].set_ylabel('Sum of empty space [meters]')
axs[1].set_xlim(0, 179)
axs[1].grid(True, alpha=0.3)
axs[1].legend(loc='upper right')

# ---- Plot 3: Monte Carlo Results ----
ax3 = axs[2]
ax3.bar(bin_centers, binned_flux_plot, width=PLOT_BIN_SIZE, color=(1, 0, 0, 0.3), edgecolor='red', linewidth=1.5, align='center')
ax3.set_title(f'Predicted flux histogram (bin_width = {PLOT_BIN_SIZE})', fontsize=12)
ax3.set_xlabel('Angle [degrees]')
ax3.set_ylabel('Particles (Normalized Flux)')
ax3.set_xlim(0, 180)

y_max_needed = np.max(binned_flux_plot) * 1.05
ax3.set_ylim(0, y_max_needed if y_max_needed > 0 else 1)
ax3.grid(True, linestyle='-', alpha=0.5)
ax3.tick_params(axis='both', direction='in', top=True, right=True)

plt.tight_layout()
plt.show()