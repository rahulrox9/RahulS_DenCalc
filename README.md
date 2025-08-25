# Density_Calculator
**Density Calculations of Icelandic Gabbroic Nodule Samples**
This project computes the density of carrier melt (scoria matrix glass), interstitial melt (nodule matrix), and mineral phases (plagioclase, olivine, clinopyroxene) in Icelandic gabbroic nodules.
Using oxide compositions, thermodynamic corrections, and modal proportions, it estimates both phase-specific and bulk nodule densities.

Author: Rahul Subbaraman <br>
Affiliation: University of Manchester <br>
Email: rahul.subbaraman@manchester.ac.uk <br>
Date: 25-08-2025

## How to Use
1. Install dependencies: pip install -r requirements.txt

2. Place your input CSV files inside the data/ folder:
   - Average_Phase_Comp.csv
   - OPAM_HS24_PT.csv

3. Run the calculation:
   python Density_Calc.py

Results will be written to:
   exports/Density.csv

## Data
- `Average_Phase_Comp.csv` : Major element oxide compositions of phases
- `OPAM_HS24_PT.csv` : Temperature (°C) and pressure (kbar) conditions
- `Modal_Propn.csv` : Modal proportions of the samples (excluding vesicles and normalised) from point counting

## Outputs 
- `exports/Density.csv`: Contains computed densities for liquid, matrix, minerals, and bulk nodules

## Workflow
1. Read Average_Phase_Comp.csv → create:
   - liq : carrier melt composition
   - mat : interstitial melt composition
   - xl : crystal compositions

2. Read OPAM_HS24_PT.csv → update liq, mat, xl with P-T columns

3. Calculate liquid densities:
   - liq and mat are passed to calculate_liquid_density()
   - H2O content is estimated as H2O = MgO * 0.052
   - `Density` dataframe is updated with mineral densities

4. Calculate mineral densities:
   - xl is passed to calculate_mineral_density()
   - `Density` dataframe is updated with mineral densities

5. Compute bulk density:
   - `Density` dataframe and `propn` dataframe (from Modal_Propn.csv) are passed to compute_bulk_density()

6. Save results:
   - `exports/Density.csv` contains liquid, matrix, mineral, and bulk nodule densities

## Notes
- H2O estimation method based on: Subbaraman et al., under revision, EPSL
- Tested with Python 3.12.7

License
MIT License

