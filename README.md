# Strain Engineering of Hexagonal Perovskites for Enhanced Oxide-Ion Conductivity

*Predicting oxygen-ion conductivities at different degrees of strain for a crystalline solid through autonomous workflows and Density Functional Theory (DFT).*

---

## 🧭 About This Project

### Real-World Problem

The challenge is to discover and optimize materials that demonstrate **high oxygen-ion conductivity at moderate to lower temperatures** (approximately 400–600 °C) to boost energy efficiency and device longevity. **Ba7Nb4MoO20** and similar hexagonal perovskite-derived materials are promising options owing to their **layered structures**. These structures include palmierite-like and perovskite blocks with inherent oxygen vacancies, which facilitate ion transport.

### Research Motivation

The research aims to **understand and improve oxygen-ion migration mechanisms in these complex materials**. It focuses on how **in-plane strain**, relevant to epitaxial thin-film growth, **impacts ion migration barriers and ionic conductivity**, an aspect that can be strategically managed by choosing different substrates during thin-film fabrication.

Hexagonal perovskites (HPs), such as **Ba7Nb4MoO20, exhibit anisotropic ionic conductivity**. The effects of strain on their crystal structure and oxygen-ion mobility remain underexplored. This study systematically investigates how uniaxial and biaxial strains, both compressive and tensile, influence oxygen migration energy barriers, vacancy formation energies, and local structural distortions, aiming to demonstrate strain engineering as a method for optimizing material performance.

### Significance of the Work

This study shows that in Ba7Nb4MoO20, unlike many other oxygen-ion conductors, **applying compressive in-plane strain *reduces* the migration energy barrier for oxygen ions by as much as 0.14 eV** (~15% decrease), leading to improved ionic conductivity. It identifies specific structural distortions and changes in the coordination environment at the transition states of migrating ions that depend on the strain direction, providing atomic-level insights into the ion migration pathways.

Additionally, the research offers practical advice for the epitaxial growth of HP thin films, indicating that **choosing substrates with certain lattice mismatches (especially inducing compressive strain) can enhance oxygen-ion conductivity**. This strain engineering approach could help design more efficient solid electrolytes for intermediate-temperature fuel cells and sensors.

By emphasizing nanoscale structural control through strain, **the study shifts the focus from solely chemical doping to include mechanical modulation of ion transport**, opening new possibilities for energy devices with better efficiency, stability, and scalability. It combines advanced Density Functional Theory calculations, Nudged Elastic Band techniques, and Continuous Symmetry Measure analyses to develop a comprehensive understanding of the mechano-chemical coupling that governs ionic conductivity in these materials.

In summary, this work addresses the critical energy materials challenge of low-temperature oxide ion conduction by demonstrating strain engineering as an effective method to improve oxygen transport in layered hexagonal perovskite oxides, with important implications for the advancement of next-generation solid-state energy conversion and storage technologies.

## 🔎 Key Insights

- **Compressive in-plane strain consistently lowers the oxygen-ion migration energy barrier** by up to 0.14 eV, resulting in nearly a **15% improvement in ionic conductivity**.

    - This behavior is distinct from many other oxide ion conductors where *tensile strain often enhances conductivity*, highlighting material-specific strain effects.

- The study reveals that **strain induces particular structural distortions and rearrangements in the local coordination environment around migrating oxygen ions**, which facilitate easier ion migration.

    - Continuous Symmetry Measure (CSM) analysis identifies specific polyhedral distortion patterns correlated with reduced migration barriers, providing atomic-scale mechanistic insights.

- The findings suggest that **epitaxial thin-film growth under carefully selected substrate-induced strain can be an effective strategy** to tailor and enhance oxygen-ion conductivity in hexagonal perovskite materials.

- This strain engineering approach offers a **promising avenue beyond chemical doping by mechanically tuning ionic transport properties**, potentially enabling more efficient and stable solid oxide fuel cells and related devices.


## 📝 Methodology

This project employed *state-of-the-art* workflow management tools to investigate oxygen-ion migration in Ba7Nb4MoO20 hexagonal perovskites under strain from first-principles.
The methodology can be summarized as follows:

- **Dataset and Systems:**  
  
  Atomic-scale models of Ba7Nb4MoO20 were constructed based on experimentally known crystal structures. The structural information was stored and managed in a SQLite3 database, which was accessed and manipulated using the Atomic Simulation Environment ([ASE](https://wiki.fysik.dtu.dk/ase/)) toolkit. 
  
  Various strain states were simulated by applying uniaxial (along a or b axes) and biaxial (a and b axes together) in-plane strains, representing both compressive and tensile regimes within ranges relevant to epitaxial thin-film growth (strain levels up to approximately ±3%).

- **Workflow Management**  

  The different stages of the computational workflow, including *structure relaxation*, *vacancy formation energy calculations*, and *Nudged Elastic Band* (NEB) migration barrier calculations, were efficiently orchestrated and managed using [PerQueue](https://asm-dtu.gitlab.io/perqueue/). This enabled automation, dynamic control, and high-throughput execution of the computational tasks across HPC resources.

  All of the calculation parameters can be found in the [project_config.cfg](https://github.com/armmorin/strained/blob/main/project_config.cfg), which PerQueue takes to add to all jobs in the workflow.

- **Calculation Process:**  

  - Oxygen vacancy formation energies were computed to assess the thermodynamic feasibility of vacancy creation under strain. 
  - Oxygen-ion migration barriers along key diffusion pathways were calculated using the Nudged Elastic Band (NEB) method to determine energy barriers for ion hopping between lattice sites.

- **Analytical Tools:**  

  Continuous Symmetry Measure (CSM) analysis quantified distortions in local polyhedral coordination environments around migrating ions, enabling correlation between structural changes and migration barrier variations.


## 📊 Visualizations

![Energies of formation for the three configurations of Ba7Nb4MoO20.](https://github.com/armmorin/strained/blob/main/strained_perovskite_eof.pdf)
![Migration barriers for the p1 configuration of Ba7Nb4MoO20 in the different strain directions.](https://github.com/armmorin/strained/blob/main/neb_energies_fit.pdf)

## 🤔 Interpretation

This work highlights **strain engineering as a key method for enhancing solid oxide ion conductors**, which is essential for low-temperature applications. Atomic-scale insights into local environment distortions help guide material design beyond Ba7Nb4MoO20. It stresses the importance of substrate choice and strain during thin-film growth, connecting fabrication processes to functional improvements. The focus shifts from just composition to nanoscale mechanical and structural adjustments of ionic pathways. These insights could accelerate the development of efficient, durable solid oxide fuel cells and energy technologies by managing ion conductivity through strain. Overall, incorporating strain engineering into material synthesis and device manufacturing offers the **potential for advanced oxide ion conductors with better conductivity and stability for clean energy solutions**.

## 🗂️ Reproducibility

- **Code Availability:**  

  The script that is managed by PerQueue is [Workflow](https://github.com/armmorin/strained/blob/main/workflow.py). It first finds the name of the structure one wishes to work with from the DB, and then depending on the range and the direction of the strain applied, generates an array of values that will be applied on the structures to be calculated. From here, the logic of the workflow kicks in, its steps are:
  
  1.  [Relaxation](https://github.com/armmorin/strained/blob/main/codes/relax.py)
  2.  [Apply Strain](https://github.com/armmorin/strained/blob/main/codes/apply_strain.py)
  3.  [preNEB](https://github.com/armmorin/strained/blob/main/codes/preneb.py)
  4.  [NEB](https://github.com/armmorin/strained/blob/main/codes/neb.py)
  5.  CINEB, it uses the same script as the previous step, but the 'climbing image' tag is enabled.

  At each step, the workflow manager retains some information to proceed, while other details are stored in the database for future use.

- **Data Accessibility:**  
  Indicate where the input data, structural models, and results can be accessed. Mention the SQLite3 database managed through Atomic Simulation Environment (ASE) if you share it, or provide instructions on how to regenerate data (e.g., from crystal structure files).

- **Environment and Dependencies:**  
  
  Essential packages for this project were:
  -  VASP v6.4.2 for the DFT calculations.
  -  A virtual environment running with Python v3.11.3 (but any version >3.10 should do).
  -  ASE v3.22.1
  -  MyQueue v24.1.0
  -  PerQueue v0.3.0a0 (Pre-release)
  -  Pymatgen v2024.5.1
  -  matplotlib v3.7.2
 
- **Instructions to Reproduce Key Results:**  

  After the workflow has finished with the calculation of the migration barriers, it is possible to jump in into every structure and analyze the evolution of their local environments, which happens through the [ce_nebs.py](https://github.com/armmorin/strained/blob/main/codes/ce_nebs.py) script. Here, we need to provide a few things like the name of the species at the A-site and the id of the migrating oxygen-ion, as well as a cutoff factor that defines the search limit for the metal-oxygen bonds.

- **Expected Outcomes:**  
  
  A CSV file will be generated once the `ce_nebs.py` script is run, together with plots of the evolution of local coordination environments for each of the structures.

## 🙏 References / Data Availability

-  For more information on this work, you can refer to the main article at: [Strain Engineering in Ba7Nb4MoO20 Hexagonal Perovskites for Enhanced Oxygen-Ion Conductivity](https://pubs.acs.org/doi/10.1021/acsmaterialslett.5c00520)

-  The database and structures are also publicly available at: [Dataset for "Strain Engineering in Ba7Nb4MoO20 Hexagonal Perovskites for Enhanced Oxygen-Ion Conductivity"](https://doi.org/10.11583/DTU.28581416.v1)
