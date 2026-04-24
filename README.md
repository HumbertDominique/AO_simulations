# TM PSF-R simulations

The simulation is based on [OOMAO](https://github.com/cmcorreia/oomao/releases/tag/v0.1) for the "Astronomical Adaptive Optics Point Spread Function Reconstruction with Machine Learning Algorithms" master theis at HES-SO master.

# Quickstart

1. Clone the repository to your favorite work folder
2. Install MatLab
3. Configure the simulation in the _ao\_inputs.txt_.
    - No edgcases are handled
    - No input check is performed.
4. Run _PSF\_R\_AOSsim\_oversampled.m_
5. When the simulation is complete, some images will appear (Interaction matrix and lightfield)
6. A metadata file is saved. it contains
    1. The Simulation parameters
    2. The number of batch in which the data was generated. Note that the data is still temporaly coherent.
    3. The quantity of items in a batch
    4. The quantity of items in the last batch

## Misc.
- GPU support is not implemented
- Correctly setting the chunksize in important in order not to break MatLab when the RAM is full. The largest images are the wavefront. I do not know their size. The secon largest image that is reasonable to save is the lightfield. Is weighs $N_{lenslet}\cdot N_{lenslet}\cdot N_{px}\cdot N_{px}\cdot 8bits$. The pixel depth can be changed in the code (array prealocation and type conversion at the value assignment).
- _preprocessing\_template.ipynb_ is a basic example on how to reassemble the dataset.

## OOMAOSource code modifications
- Added HeNe source
- Commented most of the terminal output. Matlab will unfortunatley still print some. VSCode and the official matlab extension works fine. 