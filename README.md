This MATLAB code:

1. Reads the dielectric function of a material (example files for Ag and Au are provided).

2. Lets the user define the particle radius, simulation-box size, mesh size, refractive index of the surrounding medium, and energy range.

3. Calculates the scattering, absorption, and extinction cross sections of a spherical particle as a function of photon energy, and lets the user select the energy at which the fields will be calculated.

4. Calculates the internal electric field inside the particle and the scattered electric field outside the particle in a 3D Cartesian box using the formalism of Bohren and Huffman.

   Geometry convention:
   - the incident plane wave propagates along +z;
   - the incident electric field is polarized along +x;
   - theta is the polar angle measured from +z;
   - phi is the azimuthal angle measured from +x in the xy-plane.

5. Plots several 1D and 2D cuts of the field-intensity distribution.

6. Performs a consistency check by comparing the absorbed power calculated from the internal fields with the absorbed power obtained from the absorption cross section. The agreement improves as the voxel side length `dx` is decreased.
