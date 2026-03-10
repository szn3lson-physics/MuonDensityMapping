import numpy as np
import matplotlib.pyplot as plt
from PIL import Image
from scipy.ndimage import map_coordinates, convolve1d

# ================= PARAMETRY SYMULACJI =================
IMAGE_PATH = 'converted_black_white.png'
PIXEL_SIZE = 0.13333  # metra na piksel

# Pozycja detektora (w pikselach obrazu)
DETECTOR_X = 915      
DETECTOR_Y = 1116     

MAX_RADIUS_M = 40.0   
VIEW_RADIUS_M = 100.0 

# Parametry pionowe groty i detektora
CAVE_HEIGHT = 2.5
DETECTOR_Z = 1.0      # Wysokość detektora od podłogi w metrach

# Obliczenie odległości uderzenia w sufit i podłogę. 
# Tangens kąta pionowego (2.86 st) wynosi ok. 0.05
limit_down = DETECTOR_Z / 0.05                       # 20.0 metrów
limit_up = (CAVE_HEIGHT - DETECTOR_Z) / 0.05         # 30.0 metrów

ROCK_DENSITY = 2.65 

# ================= FUNKCJE POMOCNICZE =================

def load_and_process_image(path):
    img = Image.open(path).convert('L')
    arr = np.array(img)
    return arr > 128

def muon_attenuation(rock_length_m):
    mwe = rock_length_m * ROCK_DENSITY
    intensity = 1.0 * np.exp(-mwe / 15.0) + 0.01 * np.exp(-mwe / 150.0) 
    return intensity

# ================= GŁÓWNY SKRYPT =================

cave_map = load_and_process_image(IMAGE_PATH)
center = (DETECTOR_Y, DETECTOR_X)

max_radius_px = int(MAX_RADIUS_M / PIXEL_SIZE)
view_radius_px = int(VIEW_RADIUS_M / PIXEL_SIZE)

angles_deg = np.arange(360)
angles_rad = np.radians(angles_deg)

r_pixels = np.arange(max_radius_px)
r_grid, theta_grid = np.meshgrid(r_pixels, angles_rad)

# Zmiana kierunku rośnięcia kątów na odwrotny do wskazówek zegara (CCW)
dy = -np.cos(theta_grid) * r_grid
dx = -np.sin(theta_grid) * r_grid

y_coords = center[0] + dy
x_coords = center[1] + dx

ray_pixels = map_coordinates(cave_map, [y_coords, x_coords], order=0, cval=0)
dist_m_grid = r_grid * PIXEL_SIZE

# Zastosowanie osobnych limitów dla promieni biegnących w górę i w dół
empty_mask_up = (ray_pixels > 0) & (dist_m_grid <= limit_up)
empty_mask_down = (ray_pixels > 0) & (dist_m_grid <= limit_down)

empty_lengths_up = np.sum(empty_mask_up, axis=1) * PIXEL_SIZE
empty_lengths_down = np.sum(empty_mask_down, axis=1) * PIXEL_SIZE

# Detektor zbiera sygnał z całego przedziału pionowego naraz, więc uśredniamy pustą przestrzeń
empty_lengths_360 = (empty_lengths_up + empty_lengths_down) / 2.0
rock_lengths_360 = MAX_RADIUS_M - empty_lengths_360 

raw_muon_flux_360 = muon_attenuation(rock_lengths_360)

kernel = np.bartlett(23)
kernel /= kernel.sum() 
measured_flux_360 = convolve1d(raw_muon_flux_360, kernel, mode='wrap')

front_empty = empty_lengths_360[:180]
back_empty = empty_lengths_360[180:360]
total_empty = front_empty + back_empty

final_flux_180 = measured_flux_360[:180] + measured_flux_360[180:360]
angles_180 = np.arange(180)

# ================= RYSOWANIE WYKRESÓW =================
fig, axs = plt.subplots(1, 3, figsize=(22, 6))

# ---- Wykres 1: Mapa (Kadr kołowy CCW) ----
im = axs[0].imshow(cave_map, cmap='gray')
clip_patch = plt.Circle((center[1], center[0]), view_radius_px, transform=axs[0].transData)
im.set_clip_path(clip_patch)

axs[0].add_patch(plt.Circle((center[1], center[0]), view_radius_px, color='black', fill=False, linewidth=1.5))
axs[0].plot(center[1], center[0], 'ro', markersize=6, label='Detektor Mionowy')

circle_sim = plt.Circle((center[1], center[0]), max_radius_px, color='yellow', fill=False, linestyle='--', linewidth=1.5, alpha=0.7)
axs[0].add_patch(circle_sim)

compass = {0: 'N', 45: 'NW', 90: 'W', 135: 'SW', 180: 'S', 225: 'SE', 270: 'E', 315: 'NE'}
angles_to_draw = sorted(list(set(range(0, 360, 10)) | set(compass.keys())))

edge_radius = view_radius_px
text_radius = edge_radius * 0.90  

for a in angles_to_draw:
    rad = np.radians(a)
    dx_line = -np.sin(rad) * edge_radius
    dy_line = -np.cos(rad) * edge_radius
    is_main_compass = a in compass
    alpha_line = 0.5 if is_main_compass else 0.2
    line_style = '-' if is_main_compass else ':'
    
    axs[0].plot([center[1], center[1]+dx_line], [center[0], center[0]+dy_line], 
                color='red', linestyle=line_style, alpha=alpha_line)
    
    if is_main_compass:
        label_text = f'{a}° ({compass[a]})'
        font_weight = 'bold'; font_size = 9
    else:
        label_text = f'{a}°'
        font_weight = 'normal'; font_size = 7
        
    dx_text = -np.sin(rad) * text_radius
    dy_text = -np.cos(rad) * text_radius
    axs[0].text(center[1]+dx_text, center[0]+dy_text, label_text, 
                color='red', fontsize=font_size, ha='center', va='center', fontweight=font_weight)

axs[0].set_xlim(center[1] - view_radius_px, center[1] + view_radius_px)
axs[0].set_ylim(center[0] + view_radius_px, center[0] - view_radius_px)
axs[0].set_title(f'Widok wokół detektora (Promień kadru = {VIEW_RADIUS_M}m)')
axs[0].axis('off')

# ---- Wykres 2: Sumaryczna pusta przestrzeń ----
axs[1].plot(angles_180, front_empty, color='blue', linewidth=2, linestyle='--', alpha=0.7, label='Kierunek przedni')
axs[1].plot(angles_180, back_empty, color='red', linewidth=2, linestyle='--', alpha=0.7, label='Kierunek tylny (+180°)')
axs[1].plot(angles_180, total_empty, color='orange', linewidth=2.5, linestyle='-', label='SUMA (Przód + Tył)')

axs[1].set_title('Zasięg pustej przestrzeni (Asymetria: dół=20m, góra=30m)')
axs[1].set_xlabel('Kąt wskazywany [stopnie]')
axs[1].set_ylabel('Suma pustej przestrzeni [metry]')
axs[1].set_xlim(0, 179)
axs[1].grid(True, alpha=0.3)
axs[1].legend(loc='upper right')

# ---- Wykres 3: Przewidywany wynik pomiaru ----
ax3 = axs[2]
ax3.bar(angles_180, final_flux_180, width=1.0, color=(1, 0, 0, 0.3), edgecolor='red', linewidth=1.5, align='center')
ax3.set_title('Histogram przewidywanego strumienia (bin_width = 1)')
ax3.set_xlabel('Angle [stopnie]')
ax3.set_ylabel('Particles (Normalized Flux)')
ax3.set_xlim(0, 180)

y_max_needed = np.max(final_flux_180) * 1.05
ax3.set_ylim(0, y_max_needed)
ax3.grid(True, linestyle='-', alpha=0.5)
ax3.tick_params(axis='both', direction='in', top=True, right=True)

plt.tight_layout()
plt.show()