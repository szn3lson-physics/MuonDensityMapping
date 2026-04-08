import numpy as np
import matplotlib.pyplot as plt
from PIL import Image
import matplotlib.patheffects as patheffects
import matplotlib.patches as patches
import sys

# ================= SIMULATION AND VISUALIZATION PARAMETERS =================
IMAGE_PATH = 'converted_black_white.png'
DATA_PATH = 'D:/MuonDensityMapping/output/all/zussamen_10/coin_34f.txt'
PIXEL_SIZE = 0.13333    

DETECTOR_X = 915      
DETECTOR_Y = 1116     

MAX_RADIUS_M = 40.0    
VIEW_RADIUS_M = 100.0  

# --- EXPOSURE PARAMETERS ---
BIN_SIZE_DEG = 6      
EXPONENT = 1    
BAR_COLOR = 'red'      
BAR_ALPHA = 0.6  # Set to 0.6 for better readability of the combined image      

# --- STATISTICAL PARAMETERS ---
SIGMA_MULT = 3  # Sigma multiplier (e.g., 1.0, 2.0, 3.0)

# ================= DATA LOADING =================
try:
    data_angles = np.loadtxt(DATA_PATH)
    print("Muon data loaded successfully.")
except:
    print("Data file not found. Generating random data for testing.")
    data_angles = np.random.normal(90, 20, 1500) % 180

full_data = np.concatenate([data_angles, (data_angles + 180) % 360])

# ================= CALCULATIONS AND SCALING =================
center = (DETECTOR_Y, DETECTOR_X)
max_radius_px = MAX_RADIUS_M / PIXEL_SIZE
view_radius_px = VIEW_RADIUS_M / PIXEL_SIZE

# Histogram
bins = np.arange(0, 361, BIN_SIZE_DEG)
counts, bin_edges = np.histogram(full_data, bins=bins)

# Statistics
mean_val = np.mean(counts)
std_val = np.sqrt(mean_val) # Poisson sigma = sqrt(N)
upper_limit = mean_val + (SIGMA_MULT * std_val)
lower_limit = max(0, mean_val - (SIGMA_MULT * std_val))

# Histogram scaling
available_space = view_radius_px - max_radius_px
max_c = np.max(counts) if np.max(counts) > 0 else 1
scale_ref = max(max_c, upper_limit)**EXPONENT
scale_factor = (available_space * 0.8) / scale_ref

# Cropping parameters (common for all plots)
margin = view_radius_px * 0.25
lim_x = (center[1] - view_radius_px - margin, center[1] + view_radius_px + margin)
lim_y = (center[0] + view_radius_px + margin, center[0] - view_radius_px - margin)

# ================= DRAWING FUNCTIONS =================
def apply_formatting(ax):
    """Ensures identical cropping for each layer."""
    ax.set_xlim(lim_x)
    ax.set_ylim(lim_y)
    ax.set_aspect('equal')
    ax.axis('off')

def draw_cave(ax):
    """Draws the cave map and borders."""
    try:
        img = Image.open(IMAGE_PATH).convert('L')
        cave_map = np.array(img)
        im = ax.imshow(cave_map, cmap='gray')
        circle_mask = patches.Circle((center[1], center[0]), view_radius_px, transform=ax.transData)
        im.set_clip_path(circle_mask)
    except Exception as e:
        print(f"[WARNING] Background file '{IMAGE_PATH}' not found. Drawing without map. Details: {e}")
    
    ax.add_patch(plt.Circle((center[1], center[0]), view_radius_px, color='black', fill=False, linewidth=2, zorder=12))

def draw_infrastructure(ax):
    """Draws the compass, detector axes, and its zone."""
    compass = {0: 'N', 45: 'NE', 90: 'E', 135: 'SE', 180: 'S', 225: 'SW', 270: 'W', 315: 'NW'}
    angles_to_draw = sorted(list(set(range(0, 360, 10)) | set(compass.keys())))
    text_radius_px = view_radius_px * 1.08 

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

def draw_statistics(ax):
    """Draws statistical circles (Mean +/- Sigma)."""
    if SIGMA_MULT == 0: return # Skips if sigma is set to 0
    
    stats_to_draw = [
        (lower_limit, f"-{SIGMA_MULT}σ", '--'),
        (mean_val, "Avg", '-'),
        (upper_limit, f"+{SIGMA_MULT}σ", '--')
    ]

    for val, label, lstyle in stats_to_draw:
        r_px = max_radius_px + (val**EXPONENT) * scale_factor
        circle = plt.Circle((center[1], center[0]), r_px, color='white', 
                            fill=False, linestyle=lstyle, linewidth=1.2, alpha=0.6, zorder=6)
        ax.add_patch(circle)
        
        txt_stat = ax.text(center[1], center[0] - r_px, f"{label} ({int(val)})", 
                           color='black', fontsize=9, ha='center', va='bottom', fontweight='bold', zorder=11)
        txt_stat.set_path_effects([patheffects.withStroke(linewidth=2, foreground='white')])

def draw_histogram(ax):
    """Draws bars representing muon counts."""
    for i in range(len(counts)):
        if counts[i] == 0: continue
        angle_deg = (bin_edges[i] + bin_edges[i+1]) / 2
        theta = np.radians(angle_deg)
        width = np.radians(BIN_SIZE_DEG) * 0.9 
        h_px = (counts[i]**EXPONENT) * scale_factor
        
        t_steps = np.linspace(theta - width/2, theta + width/2, 12)
        x_in = center[1] + max_radius_px * np.sin(t_steps)
        y_in = center[0] - max_radius_px * np.cos(t_steps)
        x_out = center[1] + (max_radius_px + h_px) * np.sin(t_steps)
        y_out = center[0] - (max_radius_px + h_px) * np.cos(t_steps)
        
        ax.fill(np.concatenate([x_in, x_out[::-1]]), 
                np.concatenate([y_in, y_out[::-1]]), 
                color=BAR_COLOR, alpha=BAR_ALPHA, edgecolor='darkred', linewidth=1, zorder=5)

# ================= MAIN PROGRAM (COMPOSITION) =================
print("Generating the complete image...")
fig, ax = plt.subplots(figsize=(12, 12))

# Combining all layers into one
draw_cave(ax)
draw_infrastructure(ax)
draw_statistics(ax)
draw_histogram(ax)

apply_formatting(ax)

plt.tight_layout()
plt.savefig('muon_map_complete.png', dpi=300, bbox_inches='tight', transparent=True)
print("File 'muon_map_complete.png' is ready.")
plt.show()


# ================= OPTIONAL: EXPORT TO SEPARATE LAYERS =================
# Uncomment the code block below if you want to generate 3 separate PNG files
# with perfectly matched cropping.

"""
def save_layer(filename, draw_functions):
    fig, ax = plt.subplots(figsize=(12, 12))
    for func in draw_functions:
        func(ax)
    apply_formatting(ax)
    plt.savefig(filename, dpi=300, bbox_inches='tight', transparent=True)
    plt.close(fig)
    print(f"Layer saved: {filename}")

print("\nGenerating separate layers (with transparency)...")

# If you want 100% histogram visibility on its layer, override alpha locally:
BAR_ALPHA = 1.0 

save_layer('layer_1_map.png', [draw_cave])
# In the second layer, we combine infrastructure and statistical circles
save_layer('layer_2_infrastructure.png', [draw_infrastructure, draw_statistics])
save_layer('layer_3_histogram.png', [draw_histogram])

print("Layers generated successfully!")
"""