This Matlab code:
1) Reads the dielectric function of a material.
2) Lets you define the following parameters: particle radius, size of the simulation volume, mesh size, refractive index of the surrounding medium, energy range and steps.
3) Calculates the scattering, extinction and absorption cross-section of a spherical particle and lets you choose at which energy you want to calculate the scattered and internal fields.
4) Calculates the scattered and internal fields in a 3D box at the chosen energy using the formalism in Bohren & Huffman. In agreement with the book's system of reference, the incident plane wave propagates along the z direction, with the electric field polarized along the x-axis. The polar angle 'theta', defined as the angle between the scattered vector and the z-axis, spans from 0 degrees (forward scattering) to 180 degrees (back scattering), the azimuthal angle 'phi', defined as the angle between the x-axis and the projection of the scattered vector on the xy-plane, spans from 0 degrees to 360 degrees. The code needs the function "pin_andrea" (eq. 4.47).
5) Plots several 1D and 2D cuts across the simulation volume.
6) Does an error check by comparing the absorbed power computed from the internal fields with the one obtained from the absorption cross section. The error is minimized by decreasing the mesh size via the voxel side length "dx".
