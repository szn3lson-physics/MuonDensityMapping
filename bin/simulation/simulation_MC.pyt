import numpy as np
from PIL import Image
from scipy.ndimage import map_coordinates, convolve1d
import os

# ================= PARAMETRY SYMULACJI =================
IMAGE_PATH = 'converted_black_white.png'
OUTPUT_FILE = 'simulated_muon_data.txt'
PIXEL_SIZE = 0.13333  # metra na piksel

# Pozycja detektora
DETECTOR_X = 915      
DETECTOR_Y = 1116     

MAX_RADIUS_M = 40.0   
ROCK_DENSITY = 2.65 

# Parametry Monte Carlo i tła elektroniki
EXPOSURE_EVENTS = 10000        # Oczekiwana całkowita liczba zarejestrowanych PRAWIDZIWYCH mionów
BACKGROUND_NOISE_PER_BIN = 0  # Oczekiwana liczba zdarzeń tła (szumu SiPM) na każdy 1 stopień

# ================= FUNKCJE POMOCNICZE =================

def load_and_process_image(path):
    img = Image.open(path).convert('L')
    arr = np.array(img)
    return arr > 128

def muon_attenuation(rock_length_m):
    """Zgrubne przybliżenie absorpcji mionów"""
    mwe = rock_length_m * ROCK_DENSITY
    intensity = 1.0 * np.exp(-mwe / 15.0) + 0.01 * np.exp(-mwe / 150.0) 
    return intensity

# ================= 1. MODELOWANIE FIZYKI ŚRODOWISKA =================
print("1. Skanowanie mapy jaskini i modelowanie absorpcji...")
cave_map = load_and_process_image(IMAGE_PATH)
center = (DETECTOR_Y, DETECTOR_X)
max_radius_px = int(MAX_RADIUS_M / PIXEL_SIZE)

angles_deg = np.arange(360)
angles_rad = np.radians(angles_deg)

r_pixels = np.arange(max_radius_px)
r_grid, theta_grid = np.meshgrid(r_pixels, angles_rad)

# Zgodnie z ustaleniami: obrót CCW, 0 st. = Północ
dy = -np.cos(theta_grid) * r_grid
dx = -np.sin(theta_grid) * r_grid

y_coords = center[0] + dy
x_coords = center[1] + dx

ray_pixels = map_coordinates(cave_map, [y_coords, x_coords], order=0, cval=0)
dist_m_grid = r_grid * PIXEL_SIZE

# Uwzględnienie asymetrii (detektor 1m od podłogi, 1.5m od sufitu)
tan_v = np.tan(np.radians(5.7))
limit_down = 1.0 / tan_v   # ~17.5 m (zaktualizowane do dokładnego kąta 5.7)
limit_up = 1.5 / tan_v     # ~26.3 m (zaktualizowane do dokładnego kąta 5.7)

empty_mask_up = (ray_pixels > 0) & (dist_m_grid <= limit_up)
empty_mask_down = (ray_pixels > 0) & (dist_m_grid <= limit_down)

empty_lengths_up = np.sum(empty_mask_up, axis=1) * PIXEL_SIZE
empty_lengths_down = np.sum(empty_mask_down, axis=1) * PIXEL_SIZE

empty_lengths_360 = (empty_lengths_up + empty_lengths_down) / 2.0
rock_lengths_360 = MAX_RADIUS_M - empty_lengths_360 

raw_muon_flux_360 = muon_attenuation(rock_lengths_360)

# Rozdzielczość pozioma (okno 22.6 stopnia)
kernel = np.bartlett(23)
kernel /= kernel.sum() # Tutaj normalizujemy kernel, bo ogólną skalę narzucimy zaraz w MC
measured_flux_360 = convolve1d(raw_muon_flux_360, kernel, mode='wrap')

# Złożenie wyników (nierozróżnialność przód/tył)
final_flux_180 = measured_flux_360[:180] + measured_flux_360[180:360]

# ================= 2. SYMULACJA MONTE CARLO =================
print("2. Generowanie zdarzeń Monte Carlo z uwzględnieniem szumu elektroniki...")

# Normalizujemy rozkład do 1, a następnie skalujemy przez pożądany czas ekspozycji
prob_distribution = final_flux_180 / np.sum(final_flux_180)
expected_muons_per_bin = prob_distribution * EXPOSURE_EVENTS

# Oczekiwany wynik w każdym stopniu to: Miony + Tło SiPM
expected_total_per_bin = expected_muons_per_bin + BACKGROUND_NOISE_PER_BIN

# Generowanie DYSKRETNEJ liczby zliczeń z rozkładu Poissona
mc_counts_per_bin = np.random.poisson(expected_total_per_bin)

# ================= 3. ZAPIS DO PLIKU (SUROWE DANE) =================
print("3. Tworzenie i tasowanie pliku z surowymi danymi (txt)...")
raw_data_stream = []

for angle in range(180):
    # Dodajemy dany kąt tyle razy, ile zostało fizycznie wylosowanych cząstek
    count = mc_counts_per_bin[angle]
    raw_data_stream.extend([angle] * count)

# Tasujemy listę, żeby symulowała ciągły zrzut danych w czasie (kolejność nie ma znaczenia dla histogramu)
np.random.shuffle(raw_data_stream)

# Zapis do pliku tekstowego (tylko 1 kolumna - kąt pomiaru)
with open(OUTPUT_FILE, 'w') as f:
    for angle in raw_data_stream:
        f.write(f"{angle}\n")

total_events = len(raw_data_stream)
print(f"\n[SUKCES] Wygenerowano plik: {OUTPUT_FILE}")
print(f"Całkowita liczba zapisanych zdarzeń (Miony + Tło): {total_events}")