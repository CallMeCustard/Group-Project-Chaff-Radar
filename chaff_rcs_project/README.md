    # Chaff RCS Simulation Project

This MATLAB-based project simulates the behaviour and radar scattering characteristics (RCS) of chaff clouds. It integrates physical modelling, 6DoF kinematics, Monte Carlo sampling, EM scattering approximations, and visualisation.

---

## 📁 Structure

- `sixdof/` – Contains 6DoF cloud generation logic based on aerodynamic equations.
- `cloud_utils/` – Shared functions, including plotting, colour mapping, and saving.
- `em_scattering/` – Contains RCS models (e.g. Rayleigh, Zhou-inspired, single/N dipole approximations).
- `class/ChaffCloud.m` – Main class to encapsulate chaff cloud data (positions, orientations, lengths, type).
- `main.m` – Driver script to run experiments and generate visualisations.

---

## 📌 Features

- **6DoF Cloud Generation**  
  Implements physics-based kinematics for chaff particles using realistic atmospheric parameters and aerodynamic interference coefficients. Includes orientation and RCS-effective visualisation.

- **Monte Carlo / Statistical Clouds**  
  Fast, non-kinematic cloud generation using Gaussian sampling, appropriate for scale testing.

- **EM Scattering Approximations**  
  - Zhou-inspired segmented dipole method  
  - Single/N dipole models  
  - Far-field occlusion & incoherent/coherent RCS calculations  
  - Basic Rayleigh model (planned)

- **Visualisation**  
  - 3D `scatter3` and `quiver3` plots with sin²(θ) or RCS heatmaps  
  - Temporal occlusion + RCS graphs  
  - Optional GIF export for time-resolved deployments

---

## ⚙️ How to Use

1. Run `main.m` to generate and visualise clouds.
2. Choose a cloud model (e.g. `'sixdof'`, `'mc-gaussian'`) in the `ChaffCloud` constructor.
3. Run scattering model functions on the generated `ChaffCloud` object (planned modular structure).
4. Export cloud data to `.mat` using the `.export()` method.

---

## 📈 Output Examples

- Time-resolved 3D chaff deployments
- Animated GIFs showing motion and coverage
- RCS plots across time or azimuth

---

## 🧠 Research References

- *Fast Algorithm for Full-Wave EM Scattering Analysis of Large-Scale Chaff Cloud With Arbitrary Orientation, Spatial Distribution, and Length* – Lee et al. (2023) [IEEE](https://ieeexplore.ieee.org/document/10778438)  
- *Experimental and Numerical Study of Chaff Cloud Kinetic Performance Under the Impact of High-Speed Airflow* – Zhang et al. (2018)

---

## 📌 Notes

- This is a research prototype. Accuracy is prioritised over real-time performance.
- More scattering methods and cloud models can be added easily using the `ChaffCloud` class.

