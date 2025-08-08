# Strain Engineering of Hexagonal Perovskites for Enhanced Oxide-Ion Conductivity

*An autonomous computational study demonstrating how mechanical strain can be engineered to enhance ion transport in next-generation energy materials.*

---

## 🧭 About This Project

### Real-World Problem

The development of energy technologies like solid oxide fuel cells (SOFCs) is limited by the need for materials with high oxygen-ion conductivity at lower temperatures (400–600°C). Materials like the hexagonal perovskite Ba₇Nb₄MoO₂₀ are promising candidates due to their unique layered structures, which contain intrinsic pathways for ion movement.

### Research Motivation

This research was driven by the need to systematically understand how **mechanical strain**, similar to what occurs when a thin film is grown on a substrate, affects the performance of these materials. By modeling how strain impacts ion migration, we can provide a roadmap for intentionally designing and fabricating materials with superior ionic conductivity.

### Significance of the Work

This study makes a key discovery: **compressive in-plane strain significantly *lowers* the energy barrier for oxygen-ion migration in Ba₇Nb₄MoO₂₀ by up to 15%**. This is a critical insight, as it provides a clear, actionable strategy for improving material performance. By carefully selecting substrates that induce compressive strain, it's possible to engineer more efficient solid electrolytes for future energy devices. The work provides atomic-scale insights into *why* this happens and shifts the material design paradigm from purely chemical adjustments to include mechanical and structural engineering.

---

## 🔎 Key Insights

- **Compressive strain boosts performance:** Applying compressive in-plane strain consistently lowers the oxygen-ion migration energy barrier by up to 0.14 eV, leading to a significant (~15%) enhancement in ionic conductivity. This is a counterintuitive finding compared to many other oxide materials.

- **Strain distorts the path for easier travel:** Strain induces specific structural distortions in the local atomic environment, creating more favorable pathways for oxygen ions to move through the crystal lattice.

- **A practical strategy for material design:** The results strongly suggest that epitaxial thin-film growth on substrates that induce compressive strain is an effective strategy for creating high-performance hexagonal perovskite conductors.

- **Mechanical tuning over chemical doping:** This work demonstrates that mechanically tuning a material's properties through strain can be a powerful alternative or supplement to traditional chemical doping.

---

## 📝 Methodology

This project used a fully automated workflow to investigate the effects of strain on oxygen-ion migration in Ba₇Nb₄MoO₂₀ using first-principles calculations.

- **Dataset and Systems:**  
  
  Atomic-scale models of Ba₇Nb₄MoO₂₀ were built based on known crystal structures. All structural data and calculation results were managed in a **SQLite3 database** using the **Atomic Simulation Environment (ASE)**. Various strain states (uniaxial and biaxial, up to ±3%) were systematically applied to simulate conditions relevant to experimental thin-film growth.

- **Workflow Management**  

  The entire computational pipeline—including structure relaxation, energy calculations, and migration barrier analysis—was orchestrated using **PerQueue**, a dynamic workflow manager. This enabled the automated, high-throughput execution of thousands of jobs on high-performance computing (HPC) resources. All calculation parameters are defined in the `project_config.cfg` file for full reproducibility.

- **Calculation Process:**  

  -   Oxygen vacancy formation energies were computed to assess thermodynamic stability under strain.
  -   The **Nudged Elastic Band (NEB)** method was used to calculate the energy barriers for oxygen-ion migration along key diffusion pathways.

- **Analytical Tools:**  

**Continuous Symmetry Measure (CSM)** analysis was used to quantify distortions in the local coordination of migrating ions, directly linking structural changes to variations in the migration barrier.

---

## 📊 Visualizations

**Formation Energies of Ba₇Nb₄MoO₂₀ Configurations**
![A plot showing the formation energies for three different configurations of the material.](https://github.com/armmorin/strained/blob/main/strained_perovskite_eof.pdf)

**Migration Barriers Under Different Strain Conditions**
![A plot showing the calculated migration barriers for the first configuration under various strain directions.](https://github.com/armmorin/strained/blob/main/neb_energies_fit.pdf)

---

## 🤔 Interpretation

The findings from this project strongly suggest that **strain engineering is a vital and practical tool for designing superior solid oxide ion conductors**. The atomic-scale insights gained here provide a clear design principle: by controlling mechanical strain during synthesis, we can purposefully create materials with enhanced ionic conductivity. This shifts the focus of materials design from relying solely on chemical composition to embracing the powerful interplay between a material's structure, mechanics, and function. Ultimately, this approach can accelerate the development of more efficient and durable energy technologies, such as next-generation solid oxide fuel cells.

---

## 🗂️ Reproducibility

- **Code Availability:**  

  The primary logic is contained in the [workflow.py](https://github.com/armmorin/strained/blob/main/workflow.py) script, which is managed by PerQueue. This script orchestrates the different calculation steps, including relaxation, strain application, and NEB calculations.

- **Data Accessibility:**  
  The complete dataset, including the initial structures, DFT calculation inputs/outputs, and the final SQLite3 database, is publicly available at the **DTU Data Repository**: [https://doi.org/10.11583/DTU.28581416.v1](https://doi.org/10.11583/DTU.28581416.v1).

- **Environment and Dependencies:**  
  
  Essential packages for this project were:
    -   Python (v3.11.3)
    -   VASP (v6.4.2)
    -   ASE (v3.22.1)
    -   PerQueue (v0.3.0a0)
    -   Pymatgen (v2024.5.1)
 
- **Instructions to Reproduce Key Results:**  

  After the workflow has completed, the `ce_nebs.py` script can be used to analyze the local coordination environments of the migrating ions and generate the final plots and data files.

---

## 🙏 References / Data Availability

-  **Main Publication:** For a detailed discussion of this work, please see the peer-reviewed article: [*Strain Engineering in Ba₇Nb₄MoO₂₀ Hexagonal Perovskites for Enhanced Oxygen-Ion Conductivity*](https://pubs.acs.org/doi/10.1021/acsmaterialslett.5c00520)

-  The workflow manager used is described in: [PerQueue: Dynamic Workflow Manager for Materials Discovery.](https://doi.org/10.1039/D4DD00134F)

-  The database and structures are also publicly available at: [*Dataset for" Strain Engineering in Ba₇Nb₄MoO₂₀ Hexagonal Perovskites for Enhanced Oxygen-Ion Conductivity"*](https://doi.org/10.11583/DTU.28581416.v1)
