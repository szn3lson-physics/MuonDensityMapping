import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Wedge

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 7))

# ==========================================
# 1. WIDOK Z GÓRY (Przekrój poziomy X-Y) CCW - Kąt 45°
# ==========================================
ax1.set_title("Widok z góry: Płaszczyzna pozioma (Pomiar dla kąta 45°)", fontsize=13, fontweight='bold')
ax1.set_xlim(-30, 30)
ax1.set_ylim(-30, 30)
ax1.set_aspect('equal')
ax1.grid(True, linestyle='--', alpha=0.5)

# Oś N-S, E-W
ax1.axhline(0, color='black', linewidth=0.8)
ax1.axvline(0, color='black', linewidth=0.8)
ax1.text(0, 28, 'N (0°)', ha='center', va='bottom', fontweight='bold')
ax1.text(-28, 0, 'W (90°)', ha='right', va='center', fontweight='bold')
ax1.text(0, -28, 'S (180°)', ha='center', va='top', fontweight='bold')
ax1.text(28, 0, 'E (270°)', ha='left', va='center', fontweight='bold')
ax1.plot(0, 0, 'ko', markersize=8, label="Detektor")

# 45 stopni CCW to 90 + 45 = 135 stopni w układzie matplotlib.
angle_front = 90 + 45  
w_front = Wedge((0,0), 25, angle_front - 11.3, angle_front + 11.3, color='orange', alpha=0.5, label="Przedni stożek (+/- 11.3°)")
angle_back = angle_front + 180
w_back = Wedge((0,0), 25, angle_back - 11.3, angle_back + 11.3, color='red', alpha=0.4, label="Tylny stożek (nierozróżnialność)")

ax1.add_patch(w_front)
ax1.add_patch(w_back)

x_line = [np.cos(np.radians(angle_front))*30, np.cos(np.radians(angle_back))*30]
y_line = [np.sin(np.radians(angle_front))*30, np.sin(np.radians(angle_back))*30]
ax1.plot(x_line, y_line, 'k--', alpha=0.6, label="Oś pomiaru (45° i 225°)")
ax1.legend(loc='upper right')
ax1.set_xlabel("Szerokość [m]")
ax1.set_ylabel("Długość [m]")

# ==========================================
# 2. WIDOK Z BOKU (Przekrój pionowy w bezwzględnym Y)
# ==========================================
ax2.set_title("Widok z boku: Przekrój pionowy (Detektor na y = 1m)", fontsize=13, fontweight='bold')
ax2.set_xlim(0, 35)
ax2.set_ylim(-3.0, 5.5)
ax2.grid(True, linestyle='--', alpha=0.5)

# Rysowanie jaskini
ax2.axhline(2.5, color='dimgray', linewidth=4, label="Sufit groty (y = 2.5m)")
ax2.axhline(0.0, color='dimgray', linewidth=4, label="Podłoga groty (y = 0.0m)")

# Wypełnienie skały
ax2.fill_between([0, 35], 2.5, 5.5, color='gray', alpha=0.3, label="Skała nad/pod jaskinią")
ax2.fill_between([0, 35], -3.0, 0.0, color='gray', alpha=0.3)

# Detektor
ax2.plot(0, 1.0, 'ko', markersize=8, label="Detektor (y = 1.0m)")
ax2.axhline(1.0, color='black', linestyle=':', alpha=0.3, label="Oś pozioma detektora")

tan_v = np.tan(np.radians(5.7)) # ~0.1

# ================= GEOMETRIA PROMIENI =================

# 1. WNĘTRZE GROTY (Pusta przestrzeń - ucięta przez sufit i podłogę aż do 35m)
# Górna połowa (Pomarańczowa): od y=1.0 rośnie do y=2.5 przy x=15, potem zostaje na y=2.5
x_up_cave = [0, 15, 35]
y_up_cave = [1.0, 2.5, 2.5]
ax2.fill_between(x_up_cave, 1.0, y_up_cave, color='orange', alpha=0.6, label="Górna strefa w grocie")

# Dolna połowa (Żółta): od y=1.0 maleje do y=0.0 przy x=10, potem zostaje na y=0.0
x_down_cave = [0, 10, 35]
y_down_cave = [1.0, 0.0, 0.0]
ax2.fill_between(x_down_cave, y_down_cave, 1.0, color='gold', alpha=0.6, label="Dolna strefa w grocie")

# 2. WNĘTRZE SKAŁY (Pochłanianie - tylko fragmenty wystające poza obrys groty)
# Czerwony trójkąt nad sufitem (od x=15 do x=35)
x_up_rock = [15, 35]
y_up_ray = [2.5, 1.0 + 35 * tan_v]
ax2.fill_between(x_up_rock, 2.5, y_up_ray, color='red', alpha=0.4, label="Promienie w skale")

# Czerwony trójkąt pod podłogą (od x=10 do x=35)
x_down_rock = [10, 35]
y_down_ray = [0.0, 1.0 - 35 * tan_v]
ax2.fill_between(x_down_rock, y_down_ray, 0.0, color='red', alpha=0.4)

# Zewnętrzne krawędzie promieni 
ax2.plot([0, 35], [1.0, 1.0 + 35 * tan_v], color='darkorange', linestyle='-', linewidth=1.5)
ax2.plot([0, 35], [1.0, 1.0 - 35 * tan_v], color='darkgoldenrod', linestyle='-', linewidth=1.5)

# Linie oznaczające granice (uderzenie skrajnych promieni w podłogę i sufit)
ax2.axvline(10, color='blue', linestyle='--', linewidth=1.5)
ax2.text(9.5, -1.0, 'Początek uderzeń\nw podłogę (10m)', rotation=90, va='center', color='blue')

ax2.axvline(15, color='blue', linestyle='--', linewidth=1.5)
ax2.text(14.5, 3.5, 'Początek uderzeń\nw sufit (15m)', rotation=90, va='center', color='blue')

# Legenda w prawym górnym rogu
ax2.legend(loc='upper right', fontsize=9)
ax2.set_xlabel("Odległość od detektora [m]")
ax2.set_ylabel("Wysokość (y) [m]")

plt.tight_layout()
plt.show()