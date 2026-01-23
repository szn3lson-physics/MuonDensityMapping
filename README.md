<h1 align="center">MuonDensityMapping</h1>

This project is part of thesis work "Investigation of muon tomography applications for imaging underground structures" done in University of Łódź by Kacper Dorszewski under supervision of prof. Tadeusz Wibig. The aim of the project is to adapt the CREDO-Maze software, which will process and visualise the obtained data, to the existing cosmic radiation measuring station system. The designed system enables imaging of internal rock mass structures. The idea is to analyse the directional absorption of secondary cosmic radiation muons. Field tests are currently being conducted in the Nagórzyckie Caves in Tomaszów Mazowiecki.

## Short characteristis of the detector
The muon detector is mounted horizontally and aligned along the North–South axis. It is only capable of detecting particles arriving horizontally. Because of the symmetry of the detector, it cannot distinguish whether a particle approached from direction φ or φ + 180°, so the azimuthal angle φ is restricted to the range: `0° ≤ φ ≤ 180°`.
<p align="center"><img width="789" height="576" alt="detektor_bkgless" src="https://github.com/user-attachments/assets/2766f0b1-8309-4071-813a-fbea22dc62d2" /></p>

The black boxes on the end are the scintillation detectors - created thanks to another project done in University of Łódź "Kosmos Widziany z Łodzi" as part of collaboration "CREDO-Maze" - more can be read in [CREDO-Maze Cosmic Ray Mini-Array for Educational Purposes](https://www.mdpi.com/2073-8994/14/3/500).

# ROOT practical remarks of usage on Windows OS

If you encounter the following error when trying to run `thisroot.ps1` in PowerShell:

```text
C:\root\bin\thisroot.ps1 : File C:\root\bin\thisroot.ps1 cannot be loaded because running scripts is disabled on this system. For more information, see about_Execution_Policies at 
https:/[go.microsoft.com/fwlink/?LinkID=135170](https://go.microsoft.com/fwlink/?LinkID=135170).
At line:1 char:1
+ C:\root\bin\thisroot.ps1
+ ~~~~~~~~~~~~~~~~~~~~~~~~
    + CategoryInfo          : SecurityError: (:) [], PSSecurityException
    + FullyQualifiedErrorId : UnauthorizedAccess
```
To fix this, run the following command in your PowerShell to allow executing local scripts:

`Set-ExecutionPolicy -ExecutionPolicy RemoteSigned -Scope CurrentUser`

And then proceed as usual - example:
`C:\root\bin\thisroot.ps1`
`root ROOT-Histograms\generate_all_histograms.C`

# 📊 2025 Summary & Results

<p align="center"><img width="945" height="603" alt="Screenshot 2026-01-18 at 16-58-18 D _MuonDensityMapping_results_all_3  December 17 12_coin_34f_hist_9 pdf - hist_9 pdf" src="https://github.com/user-attachments/assets/22286cca-0a2b-4177-810c-262c7dbc6e4c" /></p>

This section summarizes the muon flux measurements conducted throughout 2025. The following analysis is based on the cumulative dataset collected by pupilsof High School No. 1 in Tomaszów Mazowiecki.

The data visualized in the histogram represents the cumulative sum of all measurements collected throughout the research period. Each bin in the histogram represents an angular interval of 9 degrees, where the height of the bar indicates the number of particles detected from a specific direction. For instance, a prominent peak is observed in the 99–108° range, reaching approximately 365 counts.

The 9-degree interval was not chosen without reason; the detector can measure angles with an accuracy of 3 degrees. The 9-degree interval is simply the sum of 3 consecutive measurement intervals. Rebinning aims to increase accuracy, which depends on the number of counts, as the error $\sqrt{N}$ is reduced. 

From a statistical perspective, the mean value (Y) represents the average number of counts per bin. Given that particle detection follows a discrete distribution, the Poisson distribution was applied to analyze the data. The standard deviation (σ), derived as the square root of the mean, serves as a critical metric for distinguishing physical signals from background noise. In this analysis, 2σ threshold lines are utilized; these encompass approximately 95% of statistical fluctuations. Data points exceeding these boundaries are interpreted as significant signals, while values within the threshold are considered part of the statistical background noise.

The physical interpretation of the results are following:
* Increased particle counts in the 27–54°, 99–108°, and 126–135° intervals suggest areas of lower material density, such as geological voids or cavities.
* Significant drops in flux within the 72–81° and 117–126° ranges indicate higher matter density, likely due to the presence of metallic ores or dense rock formations that attenuate muon passage.

It should be noted that the first (0–9°) and last (171–180°) bins show approximately 50% lower counts due to specific detector geometry and edge-case triggering. However a technical solution for this characteristic is already being implemented.

Data collection will be still ongoing through 2026 till July. We belive that until this time we can achive at least 3σ in one direction, where is known that there is an empty space. This specific region corresponds to the entrance to the grottoes at 30-50 degrees from north (also visiable in the data set).

## Point Spread Function under development
What also was done in 2025 was a creation of PSF map. To generate it `his2()` is used. The programme is taking an angle $θ$ from which particle came (obtained from data set) and apply a uncertenty $θ \pm α$. The $α$ angle is multiple of three. The best results for histograms we obtain from 9 degrees and it have to be tested in 2026 if for PSF it's gonna be the right fit.

Firsty it was simulated for 5, 10 and 15 degrees to test how the program is working. The results are desplayed below.

<p align="center"><img width="624" height="612" alt="Tekst akapitu" src="https://github.com/user-attachments/assets/7f3d84d6-f780-46b3-a7a0-81f0c72b06e6" /></p>

In the end PSF map was applyed to the real data set from 26.11.2025 and generated this:

<p align="center"><img width="624" height="612" alt="Screenshot 2026-01-18 at 17-01-15 D _MuonDensityMapping_results_26_11_PSF_4_5_degree pdf - PSF_4_5_degree-1 pdf" src="https://github.com/user-attachments/assets/60822f03-ba32-488e-b289-8f5209671273" /></p>

For the folowing 2026 those actions are planned to do:
* apply to the map four cardinal directions North (N), East (E), South (S), and West (W) in corispondance how detector is placed in reality (0 angle should be on north)
* data normalisation to visually highlight directions where the sigma standard deviation is significant
* perform a simulation of the superposition of two maps created as a result of measurements taken at two different locations - possibly as `his3()`

Other tasks will be updated shortly after examination period.

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

Creates a 2D histogram representing angular distribution of muons.

Each angle generates contributions in two directions: θ and θ + 180°, because detector cannot tell the difference as is is measuring from both sides. The whole rotation of detector which is placed horizontally is 180 degree. Detector is marked as a black dot at (0,0). Data read from `Output/test.txt`.

## Results

* `angles_*.txt` – Files with calculated angles for each event.
* `coin_*.txt` – Filtered files with events meeting `"000000"` criterion.
* `his1()` – 1D histogram showing counts vs. angle.
* `his2()` – 2D angular PSF map showing spatial distribution of particle directions.


# 📘 Artificial simulation of Horizontal Muons with Monte Carlo Approach

This simulation reproduces the angular distribution of horizontal muons in an underground environment, specifically modeling an enhanced flux arriving from a nearby cavity located at an azimuth of 30∘–50∘. To achieve this, the system employs a Monte Carlo approach that generates one muon per hour, assigning each particle a zenith angle restricted to 85∘–95∘ to ensure a nearly horizontal trajectory. The azimuthal distribution uses a mixture model: while half the muons are distributed uniformly across the full 180∘ range, the other half are specifically sampled from the cave's coordinates to represent the reduced rock overburden.

The geometry assumes a 10 cm×20 cm horizontal detector aligned on a North–South axis, where each hit is assigned a precise impact coordinate and direction vector. By logging these parameters—including the hour index and xyz vector components—into a structured format, the simulation provides a controllable baseline for studying muon tomography. This setup allows researchers to calibrate physical detectors, distinguish local topological features from background noise, and develop future extensions such as depth-dependent energy loss and 3D track visualization.

`muon_generate.cpp` its a simulation program that generates a dataset of horizontal muon arrivals (technically program doesn't distinguish particles - just a coincidance between detectors). The simulation is based on Monte Carlo approach, using the Mersenne Twister engine to sample trajectories for a configurable number of hours. Each particle is assigned a zenith angle between 85∘ and 95∘ to maintain a horizontal orientation. The azimuthal distribution is determined by a probability parameter (pcave​), which directs a specific percentage of muons to originate from a narrow window between 30∘ and 50∘. This localized concentration is taken arbitrary just to represent a cave corridor or empty space from which we can detect the increased transmission of muons.

The code calculates these angles in radians using the <cmath> library and then converts the resulting azimuths back into degrees for storage. Each generated muon is recorded as a single angular entry in an output file named test.txt. This data reflects the expected physical symmetry of a horizontal detector that cannot distinguish between opposite directions, restricting the azimuthal range to 0∘–180∘. The program provides a direct numerical output of these events, allowing for the analysis of angular excesses caused by specific geological features.

Then this PSF map is created in order to test if visualisation is efficient for humans to read. The test dataset contain a sample of muons, that came more from one direction chosen arbitrary from 30 to 50 degrees. Results of this test can be seen below.

<p align="center"><img width="624" height="612" alt="Screenshot 2026-01-18 at 17-03-00 D _MuonDensityMapping_results_PSF_test_PSF_test pdf - PSF_test pdf" src="https://github.com/user-attachments/assets/0c7633a2-4475-4bc3-9862-b46f6e2efcb1" /></p>


Bibliography:
# Bibliography: Terrestrial and Underground Muon Measurements

This list contains publications and research papers focused on muon flux measurements at ground level, underground locations, and applications in muon tomography (muography).

### Muon Tomography & Structural Imaging (Pyramids, Mines, Volcanoes)
* **Borselli, D., et al. (2022).** *Three-dimensional muon imaging of cavities inside the Temperino mine (Italy)*. Scientific Reports, 12:22329. 
  - *Focus: 3D reconstruction of underground mine cavities using the MIMA hodoscope.*
* **Bross, A. D., et al. (2022).** *Tomographic Muon Imaging of the Great Pyramid of Giza*. Preprint arXiv:2202.08184.
  - *Focus: Proposed mission to use large-scale muon telescopes for full tomographic mapping of Khufu's Pyramid.*
* **ScanPyramids Collaboration (2017).** *Results of the analysis of the gas detectors*. Supplementary Materials for Nature findings.
  - *Focus: Discovery of internal voids in the Great Pyramid using gaseous tracking detectors.*

### Ground-Level Flux & Momentum Spectrum Measurements
* **Kremer, J., et al. (1999).** *Measurements of Ground-Level Muons at Two Geomagnetic Locations*. Physical Review Letters, Vol. 83, No. 21.
  - *Focus: Precise measurements of muon spectra and charge ratios at different geomagnetic latitudes.*
* **Pal, S., et al. (2012).** *Measurement of integrated flux of cosmic ray muons at sea level using the INO-ICAL prototype detector*. Journal of Cosmology and Astroparticle Physics, JCAP07(2012)033.
  - *Focus: Integrated flux data collected using Resistive Plate Chambers (RPC) at sea level.*
* **Mubashir, A., et al. (2023).** *Muon flux variations measured by low-cost portable cosmic ray detectors and their correlation with space weather activity*. Journal of Geophysical Research: Space Physics, 128.
  - *Focus: Real-time monitoring of muon flux using portable detectors across different global coordinates.*
* **Tanizaki, K. (2004).** *The Geomagnetic Latitude Effect on Atmospheric Muons at Ground Level Altitude*. Doctoral Thesis, Kobe University.
  - *Focus: Extensive study on how geomagnetic location affects muon intensity at the Earth's surface.*

### Atmospheric & Geomagnetic Studies
* **Bhattacharyya, D. P. (1970).** *The Dependence of Muon Intensity on the Geomagnetic Latitude*. Z. Physik 234, 17-22.
  - *Focus: Vertical momentum spectrum determined by counter-controlled neon hodoscope.*
* **Karmakar, N. L., Paul, A., & Chaudhuri, N. (1973).** *Measurements of Absolute Intensities of Cosmic-Ray Muons in the Vertical and Greatly Inclined Directions at Geomagnetic Latitudes 16 °N*. Il Nuovo Cimento Vol. 17 B, N. 1.
  - *Focus: Absolute integral intensities for vertical and highly inclined muon paths.*
* **Maghrabi, A. H., Alzahrani, S. A., & Alruhaili, A. S. (2023).** *The Role of Atmospheric Pressure, Temperature, and Humidity on Cosmic Ray Muons at a Low Latitude Station*. International Journal of Astronomy and Astrophysics, 13, 236-258.
  - *Focus: Environmental factors affecting cosmic ray detection at ground level.*

### Underground Measurements & Neutrino Backgrounds
* **Wolfendale, A. W., & Young, E. C. M. (1970).** *Cosmic Ray Muon Neutrino Intensities Below 1 GeV*. 
  - *Focus: Analysis of stopping muons deep underground related to neutrino interactions (Case-Wits Experiment).*


