# 🌊 HydroModel-X: Watershed Hydrological Modeling Suite (MATLAB)

*A comprehensive MATLAB toolkit for simulating watershed processes including snowmelt dynamics, evapotranspiration, and water balance calibration.*

---

## 📂 Repository Structure
HydroModel-X/
├── SnowRain_Partitioning.m # Temperature-based precipitation partitioning
├── PET_Soil_Reservoir.m # Potential Evapotranspiration + soil water balance
├── Model_Calibration.m # Parameter optimization (KGE/NSE metrics)
├── Validation_Plots.m # Observed vs. simulated discharge visualization
└── Data.mat # Sample dataset (P, Q, T timeseries)

---

## 🧰 Key Features
- **Process Modeling**
  - ✅ Snow/rain partitioning with adjustable temperature threshold
  - ✅ Degree-day snowmelt algorithm
  - ✅ Hamon's Potential Evapotranspiration (PET)
  - ✅ Dual-reservoir water balance (soil + groundwater)

- **Analysis Tools**
  - 🔍 Grid-search parameter calibration
  - 📊 Performance metrics (NSE, KGE)
  - ✂️ Split-sample validation

---

## 🚀 Quick Start
1. **Prerequisites**:
   - MATLAB R2018a or newer
   - Input data format in `Data.mat`:
     ```matlab
     D = [n×1 datetime]  % Dates
     P = [n×1 double]    % Precipitation (mm/day)
     Q = [n×1 double]    % Discharge (mm/day)
     T = [n×1 double]    % Temperature (°C)
     ```

2. **Run the model**:
   ```matlab
   load('Data.mat');
   run('SnowRain_Partitioning.m');  % Start with precipitation partitioning
---

🎓 Academic Context
Developed for:
🔬 Watershed Modeling Course @ University of Lausanne (UNIL), 2024
📚 Covers: Hydrological processes, model calibration, validation techniques
