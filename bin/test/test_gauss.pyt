import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import skew, kurtosis, norm

# ================= WCZYTYWANIE DANYCH =================
# Zakładamy, że masz plik z symulacji
try:
    data = np.loadtxt('D:/MuonDensityMapping/output/all/9_zussamen/coin__34f.txt')
except:
    # Generowanie danych dla testu, jeśli nie ma pliku
    data = np.random.normal(90, 30, 50000) % 180

# ================= 1. OBLICZANIE WYSOKOŚCI SŁUPKÓW =================
# Bierzemy 180 binów (po 1 stopniu każdy)
bins = np.arange(181)
counts, _ = np.histogram(data, bins=bins)

# ================= 2. STATYSTYKA =================
mu = np.mean(counts)
sigma = np.std(counts)
s = skew(counts)
k = kurtosis(counts) # "Fisher's definition" (kurtoza Gaussa = 0)

print(f"Statystyki wysokości słupków:")
print(f"Średnia: {mu:.2f}, Odchylenie std: {sigma:.2f}")
print(f"Skośność (Skewness): {s:.4f} (0 dla Gaussa)")
print(f"Kurtoza (Kurtosis): {k:.4f} (0 dla Gaussa)")

# ================= 3. RYSOWANIE HISTOGRAMU WYSOKOŚCI =================
fig, ax = plt.subplots(figsize=(10, 6))

# Histogram samych wysokości słupków
n, bins_h, patches = ax.hist(counts, bins=30, density=True, alpha=0.6, color='purple', edgecolor='black', label='Rozkład wysokości słupków')

# Nałożenie krzywej Gaussa dla porównania
xmin, xmax = ax.get_xlim()
x = np.linspace(xmin, xmax, 100)
p = norm.pdf(x, mu, sigma)
ax.plot(x, p, 'r--', linewidth=2, label=f'Gauss ($\mu={mu:.1f}, \sigma={sigma:.1f}$)')

# Opis statystyk na wykresie
stats_text = f"Skewness: {s:.2f}\nKurtosis: {k:.2f}"
ax.text(0.95, 0.95, stats_text, transform=ax.transAxes, fontsize=12,
        verticalalignment='top', horizontalalignment='right', 
        bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

ax.set_title('Histogram wysokości słupków (Test Gaussa)', fontsize=14, fontweight='bold')
ax.set_xlabel('Liczba zliczeń w jednym stopniu kąta', fontsize=12)
ax.set_ylabel('Częstość występowania (Gęstość)', fontsize=12)
ax.grid(True, linestyle='--', alpha=0.5)
ax.legend()

plt.tight_layout()
plt.show()