# MuonDensityMapping
MuonDensityMapping provides tools to create detailed density maps of rock masses using muon tomography. By analyzing cosmic muon trajectories, it reconstructs internal structures and topological variations, enabling visualization and analysis of geological formations.

# 📘 Muon Angular Data Analysis

This project processes muon detection data, calculates angles based on rotation and time, and produces histograms to visualize angular distributions.

## 📂 Files Overview

### `read.cpp`

Handles reading, processing, and filtering of `.log` files. Key functionalities include extracting fields, converting them to numbers, calculating angles, and filtering coincidence events.

#### Functions

`read_line(position, line_size, line)` – Extracts a semicolon-separated field from a line.
`read_int(position, line_size, line)` – Returns the integer value of a field.
`read_binary_string(position, line_size, line)` – Returns the field as a string.
`read_bit(position, line_size, line)` – Returns a single bit (0 or 1) from a field.
`output_data(file_data, file_angles)` – Reads a `.log` file, extracts fields 1, 2, and 13, and saves them to a new file.
`time_to_angle(file_angles)` – Computes the angle for each line based on rotation direction and time difference; appends as a new field.
`coindidence(file_angles, file_coincidence)` – Filters lines where field 1 equals `"000000"` and writes selected fields to output.

### Workflow

1. Read source `.log` file
2. Extract relevant fields using `read_line`
3. Compute angles with `time_to_angle`
4. Save processed data to `angles_*.txt` and `coin_*.txt`

## 📊 Histograms in ROOT

### `his1()` – 1D Histogram

Creates a 1D histogram for angles 0°–180°, styling: red line, semi-transparent fill, grid on X/Y axes.
Data read from `Output/test.txt`. Visualizes particle count vs. angle.

### `his2()` – 2D Angular Gaussian Map

Creates a 2D histogram representing angular Gaussian contributions.
Each angle generates contributions in two directions: θ and θ + 180°.
Intensity depends on distance from detector and angular difference. Detector is marked as a black dot at (0,0).
Data read from `Output/test.txt`.

## 🎯 Results

* `angles_*.txt` – Files with calculated angles for each event.
* `coin_*.txt` – Filtered files with events meeting `"000000"` criterion.
* `his1()` – 1D histogram showing counts vs. angle.
* `his2()` – 2D angular Gaussian map showing spatial distribution of particle directions.

Each generated muon is written to `test.txt` together with hour index, θ (deg), φ (deg), vector components (vx, vy, vz), and detector hit position (x, y). Example entry:
`hour: 12`
`theta_deg: 89.2`
`phi_deg: 41.5`
`vx: 0.74`
`vy: 0.67`
`vz: -0.03`
`x_cm: -1.5`
`y_cm: 4.2`
`---`



# 📘 Artificial simulation of Horizontal Muons
The simulation reproduces the angular distribution of muons expected in an underground environment with a nearby cavity that reduces the effective rock overburden in a specific azimuthal range.
## 🧭 Geometry and Physical Context
The detector is mounted horizontally and aligned along the North–South axis. It is only capable of detecting particles arriving horizontally. Because of the symmetry of the detector, it cannot distinguish whether a particle approached from direction φ or φ + 180°, so the azimuthal angle φ is restricted to the range: `0° ≤ φ ≤ 180°`. A significant feature of this simulated local topology is a nearby "empty" cave corridor covering the azimuthal range `30° – 50°`. This region has a much thinner rock overburden, so the muon flux coming from this direction is expected to be noticeably enhanced.
## 🎯 Monte Carlo Approach
The simulation generates one muon per hour over a configurable time period. Each muon is assigned a zenith angle θ (close to 90°, horizontal), an azimuth φ (with enhanced probability in the cave direction), a direction unit vector (vx, vy, vz), and an impact point on the detector plane (x, y).
### Zenith Angle θ
Since only horizontal muons can reach the detector, the zenith angle is drawn from the interval `θ ∈ [85°, 95°]`, allowing for small deviations from perfect horizontality.
### Azimuth φ: Two-Component Distribution
To represent the increased flux from the cave direction, the azimuth is drawn from a mixture model: with probability `p_cave = 0.5` the angle φ is sampled uniformly from the enhanced-flux region `30°–50°`, and with probability `1 - p_cave` it is sampled uniformly from the full range `0°–180°`. This produces a clear excess of muons in the cave direction. Mathematically: `φ = φ_cave with probability p_cave`, `φ = φ_uniform with probability 1 - p_cave`.
### Direction Vector
Each muon direction is computed using spherical coordinates: `vx = sinθ cosφ`, `vy = sinθ sinφ`, `vz = cosθ`. Because θ is close to 90°, the vz component stays near zero, ensuring horizontal trajectories.
### Detector Hit Position
The detector is a 10 cm × 20 cm rectangle. Impact coordinates are assigned uniformly across its surface: `x ∈ [−5 cm, +5 cm]`, `y ∈ [−10 cm, +10 cm]`.
## 📂 Output Format
Each generated muon is written to `test.txt` together with: the hour index, θ (deg), φ (deg), vector components (vx, vy, vz), and detector hit position (x, y). Example entry:  
`hour: 12`  
`theta_deg: 89.2`  
`phi_deg: 41.5`  
`vx: 0.74`  
`vy: 0.67`  
`vz: -0.03`  
`x_cm: -1.5`  
`y_cm: 4.2`  
`---`
## 🛠️ Extending the Simulation
Possible future improvements include: modeling energy loss in rock (dE/dx), depth-dependent transmission, realistic cosmic muon energy spectrum, 3D track visualization tools, and multi-detector coincidence simulation.
## 📄 Purpose
This simulation provides a reproducible, configurable tool for studying horizontal muon detection in cave environments and supports detector calibration, directionality studies, comparison with real data, and muon tomography of underground structures.
