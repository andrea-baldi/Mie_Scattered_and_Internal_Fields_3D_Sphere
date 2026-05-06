% Andrea Baldi, 06/05/2026
%
% Mie-theory calculation for a spherical particle.
%
% The script first calculates the scattering, absorption, and extinction
% cross sections of the particle as a function of photon energy. The user
% is then prompted to select one energy from the plotted spectrum.
%
% At the selected energy, the script calculates the internal electric field
% inside the particle and the scattered electric field outside the particle
% in a 3D Cartesian box, following the formalism of Bohren and Huffman.
%
% Geometry convention:
% - the incident plane wave propagates along +z;
% - the incident electric field is polarized along +x;
% - theta is the polar angle measured from +z;
% - phi is the azimuthal angle measured from +x in the xy-plane.

clear all;
close all;

%% Input parameters

read = load('eVe1e2_Ag_JC.txt'); % Read the dielectric function of the particle [energy in eV, epsilon1, epsilon2]
E0 = 1; % Incident electric field amplitude in V/m
r = 40; % Radius of the particle in nm
Emin = 2; % Minimum energy in eV
Emax = 4; % Maximum energy in eV
Esteps = 200; % Number of energy steps
index = 1.333; % Refractive index of the medium
boxsize = 100; % Side length of the 3D calculation box in nm
dx   = 1;      % Voxel side length in nm

%% Calculate spectra and select energy

% Compute mesh so that mesh*dx = boxsize exactly, with mesh even
mesh = ceil(boxsize/dx);
if mod(mesh,2)==1
    mesh = mesh + 1;
end
dx = boxsize/mesh;    % adjust dx so mesh*dx == boxsize

% Fundamental constants
e = 1.60217646e-19; % Elementary charge in SI units
h = 6.626068e-34; % h in SI units
c = 2.99792458e8; % Speed of light in SI units
eps0 = 8.854187817e-12; % Vacuum permittivity in F/m

% Interpolate the experimental dielectric function
w = (Emin:(Emax-Emin)/Esteps:Emax)';
e1_read = interp1(read(:,1),read(:,2),w,'spline');
e2_read = interp1(read(:,1),read(:,3),w,'spline');

% Calculate the complex refractive index of the particle
n_read = (((e1_read.^2 + e2_read.^2).^(1/2) + e1_read)./2).^(1/2);
k_read = (((e1_read.^2 + e2_read.^2).^(1/2) - e1_read)./2).^(1/2);
etot_read = n_read + 1i*k_read;
m = etot_read ./ index; % Relative refractive index (page 100)
radius = r*1e-9; % Radius in meters
lambda = h*c./(e*w); % Converts energy in (eV) to wavelength in (m)
k = 2*pi*index./lambda; % Wavenumber 'k'
x = k .* radius;  % Size parameter outside the particle (page 86)
mx = m .* x; % Size parameter inside the particle

scaele = 0; % Initialise scattering element
extele = 0; % Initialise extinction element

for n=1:7
    jnx = sqrt(pi./(2.*x)).*besselj(n+0.5,x);
    jnminx = sqrt(pi./(2.*x)).*besselj(n-0.5,x);
    jnmx = sqrt(pi./(2.*mx)).*besselj(n+0.5,mx);
    jnminmx = sqrt(pi./(2.*mx)).*besselj(n-0.5,mx);
    hnx = sqrt(pi./(2.*x)).*besselh(n+0.5,x);
    hnminx = sqrt(pi./(2.*x)).*besselh(n-0.5,x);
    xjnxdiff = x.*jnminx - n.*jnx;
    mxjnmxdiff = mx.*jnminmx - n.*jnmx;
    xhnxdiff = x.*hnminx - n.*hnx;
    % Calculate the scattering coefficients from eq. 4.53
    an = (m.^2.*jnmx.*xjnxdiff-jnx.*mxjnmxdiff)./(m.^2.*jnmx.*xhnxdiff-hnx.*mxjnmxdiff);
    bn = (jnmx.*xjnxdiff-jnx.*mxjnmxdiff)./(jnmx.*xhnxdiff-hnx.*mxjnmxdiff);
    scaele = scaele + (2*n+1).*(an.*conj(an)+bn.*conj(bn));
    extele = extele + (2*n+1).*real(an+bn);
end

Csca = 2*pi./(k.^2).*scaele; % Scattering cross section, from eq. 4.61
Cext = 2*pi./(k.^2).*extele; % Extinction cross section, from eq. 4.62
Cabs = Cext-Csca;            % Absorption cross section
Qsca = Csca./(pi*radius^2);  % Scattering efficiency
Qext = Cext./(pi*radius^2);  % Extinction efficiency
Qabs = Cabs./(pi*radius^2);  % Absorption efficiency

% Plot the scattering, absorption, and extinction efficiencies and select
% the energy at which to calculate the fields
figure(1)
plot(w,Qsca,'b','linewidth',2)
hold on
plot(w,Qabs,'r','linewidth',2)
plot(w,Qext,'k','linewidth',2)
xlabel('Energy (eV)', 'FontSize', 10 );
ylabel(['Efficiencies of a ',num2str(radius*1E9),' nm radius sphere'],'FontSize',10);
legend('Scattering','Absorption','Extinction')
title('Select the energy');

% Select the energy at which to calculate the fields
[x0,y0]=ginput(1);
[~,indexplot]=min(abs(w-x0));
E = w(indexplot); % Energy at which we evaluate the fields in eV
lambdaplot = h*c./(e*E); % Wavelength at which we evaluate the fields
omega = 2 * pi * E * e / h; % Angular frequency in rad/s at which we evaluate the fields (E in eV, h in J*s)

%% Calculate internal and scattered fields

m = etot_read(indexplot) / index; % Relative refractive index m = N1/N
k = 2*pi*index / lambdaplot; % Wavenumber outside the particle
x = k * radius; % Size parameter outside the particle
k_int = m * k; % Wavenumber inside the particle
mx = m * x; % Size parameter inside the particle

nmax = round(x + 4 * x.^(1/3) + 2); % Maximum order of the vector spherical harmonics using the Wiscombe criterion

% Create an even, voxel-centered grid so that no grid point lies exactly at
% the origin
coordx = linspace(-boxsize/2 + dx/2, boxsize/2 - dx/2, mesh);
coordy = linspace(-boxsize/2 + dx/2, boxsize/2 - dx/2, mesh);
coordz = linspace(-boxsize/2 + dx/2, boxsize/2 - dx/2, mesh);

[X,Y,Z] = ndgrid(coordx, coordy, coordz);  % sizes: [nx,ny,nz]
X = X*1e-9;  Y = Y*1e-9;  Z = Z*1e-9; % convert to meters

R3D        = sqrt(X.^2 + Y.^2 + Z.^2); % distance from the center at any point
costheta3D = Z ./ R3D; % cosine of polar angle theta at any point
sinTheta3D = sqrt(1 - costheta3D.^2); % sine of polar angle theta at any point
phi3D      = atan2(Y, X); % azimuthal angle phi at any point
cosPhi3D   = cos(phi3D); % cosine of azimuthal angle phi at any point
sinPhi3D   = sin(phi3D); % sine of azimuthal angle phi at any point

% Initialize field-intensity array and Mie coefficient vectors
I_tot  = zeros(mesh,mesh,mesh);
En_all = zeros(nmax,1);
an_all = zeros(nmax,1);
bn_all = zeros(nmax,1);
cn_all = zeros(nmax,1);
dn_all = zeros(nmax,1);

for n=1:nmax
    
    % Prefactor, En (defined right after eq. 4.40)
    En_all(n) = 1i^n*E0*(2*n+1)/(n*(n+1));
    
    % Bessel functions and derivatives
    jnx = sqrt(pi./(2.*x)).*besselj(n+0.5,x);
    jnminx = sqrt(pi./(2.*x)).*besselj(n-0.5,x);
    jnmx = sqrt(pi./(2.*mx)).*besselj(n+0.5,mx);
    jnminmx = sqrt(pi./(2.*mx)).*besselj(n-0.5,mx);
    hnx = sqrt(pi./(2.*x)).*besselh(n+0.5,x);
    hnminx = sqrt(pi./(2.*x)).*besselh(n-0.5,x);
    xjnxdiff = x.*jnminx - n.*jnx;
    mxjnmxdiff = mx.*jnminmx - n.*jnmx;
    xhnxdiff = x.*hnminx - n.*hnx;
    
    % Calculate the internal-field coefficients at the selected wavelength
    % using equation (4.53) with \mu = \mu_l = 1
    cn_all(n) = (jnx*xhnxdiff-hnx*xjnxdiff)/(jnmx*xhnxdiff-hnx*mxjnmxdiff);
    dn_all(n) = (m*jnx*xhnxdiff-m*hnx*xjnxdiff)/(m^2*jnmx*xhnxdiff-hnx*mxjnmxdiff);
    
    % Calculate the scattering coefficients at the selected wavelength
    % using equation (4.53) with \mu = \mu_l = 1
    an_all(n) = (m^2*jnmx*xjnxdiff-jnx*mxjnmxdiff)/(m^2*jnmx*xhnxdiff-hnx*mxjnmxdiff);
    bn_all(n) = (jnmx*xjnxdiff-jnx*mxjnmxdiff)/(jnmx*xhnxdiff-hnx*mxjnmxdiff);
    
end

ix0 = mesh/2 + 1;   % first index at x ≥ 0
iy0 = mesh/2 + 1;   % first index at y ≥ 0

for px = ix0:mesh % only x ≥ 0
    
    if mod(px-ix0,10)==0 || px==ix0 || px==mesh
        fprintf('%.1f%% complete\n', 100*(px-ix0+1)/(mesh-ix0+1));
    end
    
    for py = iy0:mesh % only y ≥ 0
        
        for pz=1:mesh % all z values
            
            distance = R3D(px,py,pz);
            costheta = costheta3D(px,py,pz);
            sintheta = sinTheta3D(px,py,pz);
            phi      = phi3D(px,py,pz);
            cosphi   = cosPhi3D(px,py,pz);
            sinphi   = sinPhi3D(px,py,pz);
            
            % Initialize field components at this grid point
            E_sca_r = 0; % radial component of the scattered field
            E_sca_t = 0; % polar component of the scattered field
            E_sca_p = 0; % azimuthal component of the scattered field
            E_int_r = 0; % radial component of the internal field
            E_int_t = 0; % polar component of the internal field
            E_int_p = 0; % azimuthal component of the internal field
            
            if distance <= radius % Calculate the fields inside the particle
                
                for n=1:nmax
                    
                    % pi_n and tau_n, eq. 4.46 and 4.47
                    pin = pin_andrea(n,costheta);
                    pin_prev = pin_andrea(n-1,costheta);
                    taun = n*costheta*pin - (n+1)*pin_prev;
                    
                    % Bessel functions of the first kind from eq. 4.9, using the wavenumber inside the particle, k_int
                    jnkr     = sqrt(pi/(2*k_int*distance))*besselj(n+0.5,k_int*distance);
                    jnminkr  = sqrt(pi/(2*k_int*distance))*besselj(n-0.5,k_int*distance);
                    jnpluskr = sqrt(pi/(2*k_int*distance))*besselj(n+1.5,k_int*distance);
                    
                    % Bessel function's derivative on page 94 (eq. between 4.43 and 4.44)
                    djnkr_dkr = (n*jnminkr-(n+1)*jnpluskr)/(2*n+1);
                    
                    % Factor [\rho j_n(\rho)]' appearing in equations 4.50 for
                    % the Bessel function
                    dkrjnkr_dkr = jnkr+k_int*distance*djnkr_dkr;
                    
                    % Components of the vector spherical harmonics, equations 4.50
                    N1_eln_r = cosphi*n*(n+1)*sintheta*pin*jnkr/(k_int*distance);
                    N1_eln_t = cosphi*taun*dkrjnkr_dkr/(k_int*distance);
                    N1_eln_p = -sinphi*pin*dkrjnkr_dkr/(k_int*distance);
                    M1_oln_r = 0;
                    M1_oln_t = cosphi*pin*jnkr;
                    M1_oln_p = -sinphi*taun*jnkr;
                    
                    % Components of the internal fields, eq. 4.40
                    E_int_r = E_int_r + En_all(n)*(cn_all(n)*M1_oln_r - 1i*dn_all(n)*N1_eln_r);
                    E_int_t = E_int_t + En_all(n)*(cn_all(n)*M1_oln_t - 1i*dn_all(n)*N1_eln_t);
                    E_int_p = E_int_p + En_all(n)*(cn_all(n)*M1_oln_p - 1i*dn_all(n)*N1_eln_p);
                    
                end
                
                Iq = abs(E_int_r)^2+abs(E_int_t)^2+abs(E_int_p)^2;
                
            else % Calculate the fields outside the particle
                
                for n=1:nmax
                    
                    % pi_n and tau_n, eq. 4.46 and 4.47
                    pin = pin_andrea(n,costheta);
                    pin_prev = pin_andrea(n-1,costheta);
                    taun = n*costheta*pin - (n+1)*pin_prev;
                    
                    % Hankel functions, from combining equations 4.9 and 4.13,
                    % using the wavenumber outside the particle, k
                    hnkr     = sqrt(pi/(2*k*distance))*besselh(n+0.5,1,k*distance);
                    hnminkr  = sqrt(pi/(2*k*distance))*besselh(n-0.5,1,k*distance);
                    hnpluskr = sqrt(pi/(2*k*distance))*besselh(n+1.5,1,k*distance);
                    
                    % Hankel function's derivative on page 94 (eq. between 4.43 and 4.44)
                    dhnkr_dkr = (n*hnminkr-(n+1)*hnpluskr)/(2*n+1);
                    
                    % Factor [\rho h_n(\rho)]' appearing in equations 4.50 for
                    % the Hankel function
                    dkrhnkr_dkr = hnkr+k*distance*dhnkr_dkr;
                    
                    % Components of the vector spherical harmonics, equations 4.50
                    N3_eln_r = cosphi*n*(n+1)*sintheta*pin*hnkr/(k*distance);
                    N3_eln_t = cosphi*taun*dkrhnkr_dkr/(k*distance);
                    N3_eln_p = -sinphi*pin*dkrhnkr_dkr/(k*distance);
                    M3_oln_r = 0;
                    M3_oln_t = cosphi*pin*hnkr;
                    M3_oln_p = -sinphi*taun*hnkr;
                    
                    % Components of the scattered fields, eq. 4.45
                    E_sca_r = E_sca_r + En_all(n)*(1i*an_all(n)*N3_eln_r - bn_all(n)*M3_oln_r);
                    E_sca_t = E_sca_t + En_all(n)*(1i*an_all(n)*N3_eln_t - bn_all(n)*M3_oln_t);
                    E_sca_p = E_sca_p + En_all(n)*(1i*an_all(n)*N3_eln_p - bn_all(n)*M3_oln_p);
                    
                end
                
                Iq = abs(E_sca_r)^2+abs(E_sca_t)^2+abs(E_sca_p)^2;
                
            end
            
            % Use mirror symmetry in x and y to fill the remaining quadrants
            I_tot(px,        py,        pz) = Iq;
            I_tot(mesh+1-px, py,        pz) = Iq;
            I_tot(px,        mesh+1-py, pz) = Iq;
            I_tot(mesh+1-px, mesh+1-py, pz) = Iq;
            
        end
    end
end

%% Plot the field intensity in the xz plane

% Find the y = 0 index (closest to the center of the box)
[~, y_center] = min(abs(coordy));

% Extract the xz-slice from the 3D field matrix at y = 0
I_tot_xz = squeeze(I_tot(:, y_center, :)); % [x, z] slice

figure(2)
imagesc(coordx, coordz, I_tot_xz'/E0^2) % transpose to match axes
axis equal tight
title('Electric field intensity enhancement in the xz plane, |E/E_0|^2')
xlabel('x (nm)')
ylabel('z (nm)')
hold on
theta_circle = linspace(0,2*pi,300);
x_circle = r * cos(theta_circle);
z_circle = r * sin(theta_circle);
plot(x_circle, z_circle, 'k-', 'LineWidth', 1)
% Arrow settings
arrow_length = r * 0.4; % length of the arrow (adjust as needed)
arrow_linewidth = 2;
arrow_color = 'w'; % white arrows (visible on color map)
% Get top of the z-axis from coordz
z_top = min(coordz) + 0.1 * r; % slightly below top edge
x_origin = 0;
% Propagation direction arrow (downward along +z)
quiver(x_origin, z_top, 0, arrow_length, 0, ...
    'Color', arrow_color, 'LineWidth', arrow_linewidth, 'MaxHeadSize', 0.8)
% Polarization direction arrow (horizontal along +x)
quiver(x_origin, z_top, arrow_length, 0, 0, ...
    'Color', arrow_color, 'LineWidth', arrow_linewidth, 'MaxHeadSize', 0.8)
% Labels
text(x_origin, z_top + arrow_length + 2, 'k', ...
    'Color', arrow_color, 'FontSize', 12, 'HorizontalAlignment', 'center')
text(x_origin + arrow_length + 2, z_top, 'E', ...
    'Color', arrow_color, 'FontSize', 12, 'HorizontalAlignment', 'center')

%% Plot the field intensity along the x-axis and z-axis

[~, px_center] = min(abs(coordx)); % Find the index closest to x = 0
[~, py_center] = min(abs(coordy)); % Find the index closest to y = 0
[~, pz_center] = min(abs(coordz)); % Find the index closest to z = 0

I_zline = squeeze(I_tot(px_center, py_center, :)); % E along z-axis
I_xline = squeeze(I_tot(:, py_center, pz_center)); % E along x-axis

figure(3)
plot(coordz, I_zline / E0^2, '-ok')
hold on
xline(r, '--', 'Color', [0.5 0.5 0.5])
xline(-r, '--', 'Color', [0.5 0.5 0.5])
xlabel('z (nm)')
ylabel('|E/E_0|^2')
title('Field intensity enhancement along z-axis')

figure(4)
plot(coordx, I_xline / E0^2, '-ok')
hold on
xline(r, '--', 'Color', [0.5 0.5 0.5])
xline(-r, '--', 'Color', [0.5 0.5 0.5])
xlabel('x (nm)')
ylabel('|E/E_0|^2')
title('Field intensity enhancement along x-axis')

%% Check agreement between analytical and calculated power loss densities

inside_mask = (R3D <= radius); % inside the sphere

% Power loss density in the particle
P_loss = zeros(size(I_tot));
P_loss(inside_mask) = 0.5 * eps0 * omega * e2_read(indexplot) * I_tot(inside_mask);

% Voxel volume (m³)
dx = abs(coordx(2) - coordx(1)) * 1e-9;
dy = abs(coordy(2) - coordy(1)) * 1e-9;
dz = abs(coordz(2) - coordz(1)) * 1e-9;
dV = dx * dy * dz;

% Absorbed power inside the particle
P_total = sum(P_loss(:)) * dV;

% Absorbed power from absorption cross-section
P_total_analytical = 0.5 * eps0 * c * index * E0^2 * Cabs(indexplot);

% Comparison
fprintf('Simulated P_abs = %.4e W\n', P_total);
fprintf('Analytical P_abs from C_abs = %.4e W\n', P_total_analytical);
fprintf('Relative error = %.2f%%\n', 100 * abs(P_total - P_total_analytical) / P_total_analytical);

%% Local function

% Function pi_n as described in eq. 4.47 in Bohren and Huffman.
% Inputs are the order of the vector spherical harmonics, 'n', and the
% cosine of the polar angle, 'cos(theta)'
function p=pin_andrea(n,costheta)
if n==0
    p=0;
elseif n==1
    p=1;
else
    q=[0,1];
    for j=2:n
        q(j+1) = costheta*q(j)*(2*j-1)/(j-1)-j*q(j-1)/(j-1);
    end
    p=q(end);
end
end
