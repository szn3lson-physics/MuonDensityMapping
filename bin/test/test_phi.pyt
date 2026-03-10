import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import chisquare, ks_2samp
import os

# ================= NAZWY PLIKÓW =================
FILE_MC = 'D:/MuonDensityMapping/bin/simulation/simulated_muon_data.txt'
FILE_REAL = 'D:/MuonDensityMapping/output/all/9_zussamen/coin__34f.txt'

# ================= 1. WCZYTYWANIE DANYCH =================
data_mc = np.loadtxt(FILE_MC)
data_real = np.loadtxt(FILE_REAL)

# Zabezpieczenie przed kątami spoza przedziału
data_mc = data_mc[(data_mc >= 0) & (data_mc < 180)]
data_real = data_real[(data_real >= 0) & (data_real < 180)]

# ================= 2. BINOWANIE I NORMALIZACJA =================
bins = np.arange(181) # Krawędzie binów 0-180 (szerokość 1 stopień)
counts_mc, _ = np.histogram(data_mc, bins=bins)
counts_real, _ = np.histogram(data_real, bins=bins)

total_mc = np.sum(counts_mc)
total_real = np.sum(counts_real)

# Normalizacja - pole pod wykresem = 1
norm_mc = counts_mc / total_mc if total_mc > 0 else counts_mc
norm_real = counts_real / total_real if total_real > 0 else counts_real

# ================= 3. TESTY STATYSTYCZNE =================
# Test Kołmogorowa-Smirnowa - idealny do kształtów (działa na surowych danych)
ks_stat, ks_pval = ks_2samp(data_real, data_mc)

# Test Chi-Kwadrat - wymaga surowej statystyki (ilości zliczeń), żeby poprawnie
# wyliczyć błąd. Skalujemy MC do "wielkości" paczki Real Data, żeby porównać je bin po binie.
mask = counts_mc > 0
expected_counts_scaled = counts_mc * (total_real / total_mc)
chi2_stat, chi2_pval = chisquare(f_obs=counts_real[mask], f_exp=expected_counts_scaled[mask])

# Wypisywanie ("coutowanie") wyników w konsoli
print("\n========================================================")
print(" WYNIKI PORÓWNANIA HISTOGRAMÓW (Podobieństwo)")
print("========================================================")
print(f" -> Test Kołmogorowa-Smirnowa : {ks_pval * 100.0:.6f} %")
print(f" -> Test Chi-Kwadrat (Chi2)   : {chi2_pval * 100.0:.6f} %")
print("--------------------------------------------------------")
print(" Interpretacja:")
print(" Powyższe wartości (p-value) określają, z jakim")
print(" prawdopodobieństwem oba rozkłady pochodzą z tej")
print(" samej populacji (czyli z jakim prawdopodobieństwem")
print(" symulacja wiernie oddaje pomiar z detektora).")
print("========================================================\n")

# ================= 4. RYSOWANIE WYKRESU =================
fig, ax = plt.subplots(figsize=(12, 7))
bin_centers = bins[:-1] + 0.5

# Rysowanie histogramu Monte Carlo (Czerwony)
ax.hist(bin_centers, bins=bins, weights=norm_mc, 
        histtype='stepfilled', color=(1, 0, 0, 0.3), edgecolor='red', 
        linewidth=2, label='Monte Carlo')

# Rysowanie histogramu Prawdziwych Danych (Niebieski)
ax.hist(bin_centers, bins=bins, weights=norm_real, 
        histtype='stepfilled', color=(0, 0, 1, 0.3), edgecolor='blue', 
        linewidth=2, label='Prawdziwe Dane z Detektora')

# Ustawienia osi
ax.set_title('Porównanie Danych: Detektor vs Symulacja MC', fontsize=14, fontweight='bold')
ax.set_xlabel('Kąt pomiaru [stopnie]', fontsize=12)
ax.set_ylabel('Prawdopodobieństwo', fontsize=12)
ax.set_xlim(0, 180)

# Ustawienie limitu Y, by zmieścić legendę i najwyższe słupki
max_y = max(np.max(norm_mc), np.max(norm_real))
ax.set_ylim(0, max_y * 1.3)

ax.grid(True, linestyle='-', alpha=0.3)
ax.tick_params(axis='both', direction='in', top=True, right=True)

# Legenda w lewym górnym rogu
ax.legend(loc='upper left', fontsize=10, framealpha=1.0)

plt.tight_layout()
plt.show()