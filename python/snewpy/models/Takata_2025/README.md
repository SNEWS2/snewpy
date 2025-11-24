# Takata_2025 1-D CCSN Simulation Models with Axion-like Particles

Neutrino data for 1-D core-collapse supernova models in the presence of axion-like particles.

The reference article is [**Progenitor dependence of neutrino-driven supernova explosions with the aid of heavy axionlike particles**](https://doi.org/10.1103/PhysRevD.111.103028) by Takata, T. _et al._, PhysRevD.111.103028, 2025.

## ALP parameters and progenitor masses simulated

| Progenitor mass [Msun] | Axion mass [MeV] | Axion-photon coupling [1e-10/GeV] |
| ---------------------- | ---------------- | --------------------------------- |
| 11.2 | 0, 40, 100, 150, 200, 300, 400, 600, 800 | 0, 2, 4, 6, 8, 10 |
| 20 | 0, 40, 100, 150, 200, 300, 400, 600, 800 | 0, 2, 4, 6, 8, 10 |
| 25 | 0, 40, 100, 150, 200, 300, 400, 600, 800 | 0, 2, 4, 6, 8, 10 |

## Data format

The data files are in the following format:

| Column(s) | Value(s) |
| --------- | -------- |
| 1         | Time [s] |
| 2         | Post-bounce time [s] |
| 3-5       | Number luminosity [1/s] |
| 6-8       | Luminosity [erg/s] |
| 9-11      | Mean energy [MeV] |
| 12-15     | RMS energy [MeV] |



