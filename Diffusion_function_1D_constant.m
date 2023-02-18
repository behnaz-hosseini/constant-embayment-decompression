%% This is the 1D decompression-diffusion function that is called upon by the main code (constant_embayment_decompression.m)
function [nodes, H2O_array, CO2_array] = Diffusion_function_1D_constant(P_i, P_f, lengthi, radius, H2O_i, rho, TC, dPdt)

TK = TC + 273.15; % Temperature [K], converted from input temperature [C]
R = 8.314; % Universal gas constant [J/mol K] 
H2O_molecular_mass = 18.015; % Molecular mass of H2O [g/mol]                               
CO2_molecular_mass = 44.01; % Molecular mass of CO2 [g/mol]
nnodes = 51; % Number of grid points
dx = lengthi / (nnodes - 1); % Set the node/grid spacing in the RE [um]    
time = 0; % Initialize time to use in progressing the main loop
dt = 0.05; % Change in time [s]

%% Initial CO2 calculation; construct the P_i [MPa] isobar using Liu et al. (2005) solubility model

Xw_isobar = 0:0.01:1;
    for i = 1:length(Xw_isobar)
    % At each value of P, H2O and CO2 concentrations at outlet are assumed
    % to be in equilibrium with gas in bubbles that are present in melt
    % outside the crystal, and thus can be calculated using the solubility
    % model of Liu et al. (2005) for rhyolitic melt
    Pw_isobar(i) = P_i * Xw_isobar(i); % Mole fraction of H2O in fluid multiplied by initial pressure
    Pc_isobar(i) = (1 - Xw_isobar(i)) * P_i; % Mole fraction of CO2 in fluid multiplied by initial pressure
    H2O_melt(i) = (354.94 * Pw_isobar(i)^0.5 + 9.623 * Pw_isobar(i) - 1.5223 * Pw_isobar(i)^1.5) / TK + 0.0012439 * Pw_isobar(i)^1.5 + Pc_isobar(i) * ((-1.084 * 10^-4) * Pw_isobar(i)^0.5 - (1.362 * 10^-5) * Pw_isobar(i)); % Total dissolved H2O content [wt %], eq. (2a) in Liu et al. (2005) or eq. (8) in Liu et al. (2007)
    CO2_melt(i) = Pc_isobar(i) * (5668 - 55.99 * Pw_isobar(i)) / TK + Pc_isobar(i) * (0.4133 * Pw_isobar(i)^0.5 + (2.041 * 10^-3) * Pw_isobar(i)^1.5); % CO2 content [ppm], eq. (2b) in Liu et al. (2005) or eq. (9) in Liu et al. (2007)
    end
    
% 1D linear data interpolation
H2Ov_i = interp1(H2O_melt', Xw_isobar', H2O_i, 'linear') * 100; % H2O in vapor [mol %]
CO2v_i = (1 - H2Ov_i / 100) * 100; % CO2 in vapor [mol %]
CO2_i = interp1(H2O_melt', CO2_melt', H2O_i, 'linear'); % Outputs initial CO2 [ppm]
Pw = H2Ov_i / 100 * P_i; % Outputs initial H2O in vapor phase [mol %]
Pc = CO2v_i / 100 * P_i; % Outputs initial CO2 in vapor phase [mol %]

%% Calculate initial moles of H2O and CO2 in bubble (Ideal Gas Law)

radius_meters = radius/(10^6); % Bubble radius (m) 
bubble_volume = 4/3 * pi * radius_meters.^3; % Volume of bubble (m^3)
P_i_Pa = P_i * 1000000; % initial pressure (Pa)

% Solve for moles H2O and CO2 using Ideal Gas Law
moles_total_i = (P_i_Pa * bubble_volume)/(R * TK); % Total moles of gas
moles_H2O_i = moles_total_i * (H2Ov_i/100);
moles_CO2_i = moles_total_i * CO2v_i/100;

%% Initialize all arrays
% CO2, H2O, X_mole_fraction, dD_CO2, dD_H2O, d_H2O, d_CO2, d2_H2O, d2_CO2, D_H2O, D_CO2 length arrays 

% CO2 and H2O concentration arrays (wt. % and ppm, respectively)
H2O_array = zeros(nnodes,1);
CO2_array = zeros(nnodes,1);
new_H2O_array = zeros(nnodes,1);
new_CO2_array = zeros(nnodes,1);

% For plotting the time evolution of the boundary nodes
% H2O_interface = zeros((P_i-P_f)/dPdt,1);
% CO2_interface = zeros((P_i-P_f)/dPdt,1);

% H2O concentration array (mol fraction of H2O, single oxygen basis -- for diffusivity of H2O) 
X_mole_fraction = zeros(nnodes,1);

% Initial concentrations in H2O and CO2 arrays
for k = 1:nnodes
    CO2_array(k) = CO2_i;
    H2O_array(k) = H2O_i;
end

% Derivative arrays
dD_CO2 = zeros(nnodes,1);
dD_H2O = zeros(nnodes,1);
d_H2O_array = zeros(nnodes,1);
d_CO2_array = zeros(nnodes,1);
d2_H2O_array = zeros(nnodes,1);
d2_CO2_array = zeros(nnodes,1);
D_H2O = zeros(nnodes,1);
D_CO2 = zeros(nnodes,1);

% Set up node positions array
nodes = zeros(nnodes,1);
nodes(1) = 0;
for i = 2:nnodes
    nodes(i) = nodes(i - 1) + dx;
end 

%% Diffusion of H2O and CO2

% Update H2O_moles_bubble to be value of moles_H2O_i
H2O_moles_bubble = moles_H2O_i;
CO2_moles_bubble = moles_CO2_i;

iter = 1;
P = P_i;

while P > P_f  
    P = P - dPdt * dt;
    m = -20.79 - 5030/TK - 1.4 * P/TK; % Calculate m for diffusivity of H2Ot equation

    % Calculate H2O and CO2 diffusivity constants  
    for i = 1:nnodes  
        % Convert H2O wt. % to H2O mole fraction (single oxygen basis)
        X_mole_fraction(i) = (H2O_array(i)/18.015)/(H2O_array(i)/18.015 + CO2_array(i)/440100 + (100 - H2O_array(i) - 0.0001*CO2_array(i))/32.49);
        D_H2O(i) = (X_mole_fraction(i) * exp(m) * (1 + exp(56 + m + X_mole_fraction(i) * (-34.1 + 44620/TK + (57.3 * P)/TK) - sqrt(X_mole_fraction(i)) * (0.091 + (4.77 * 10^6)/TK^2)))); % Eq 16 in Zhang et al. (2007)    
        D_CO2(i) = ((exp(-8.20 - (22963 + 2.005 * P)/TK + (-1.4262 + 2416.1/TK) * H2O_array(i)))*10^12); % Eq 29 in Zhang et al. (2007) 
    end
    % Update all the derivatives using finite differences for each term in
    % the continuity equation (Fick's second law with concentration
    % dependence of diffusivity)
    for i = 2:nnodes-1
        dD_H2O(i) = (D_H2O(i+1) - D_H2O(i-1))/(2*dx); % Finite difference of H2O diffusivity
        dD_CO2(i) = (D_CO2(i+1) - D_CO2(i-1))/(2*dx); % Finite difference of CO2 diffusivity
        d_H2O_array(i) = (H2O_array(i + 1) - H2O_array(i - 1))/(2*dx); % Finite difference of H2O concentration
        d_CO2_array(i) = (CO2_array(i + 1) - CO2_array(i - 1))/(2*dx); % Finite difference of CO2 concentration
        d2_H2O_array(i) = (H2O_array(i + 1) - 2 * H2O_array(i) + H2O_array(i - 1))/(dx^2);
        d2_CO2_array(i) = (CO2_array(i + 1) - 2 * CO2_array(i) + CO2_array(i - 1))/(dx^2);
    end
    % Diffusion happens
    for i = 2:nnodes-1
        % Diffuse H2O: combine finite difference terms from above into
        % Fick's second law with concentration-dependent diffusivity
        new_H2O_array(i) = (dD_H2O(i) * d_H2O_array(i) + D_H2O(i) * d2_H2O_array(i)) * dt + H2O_array(i);
        % Diffuse CO2: combine finite difference terms from above into
        % Fick's second law with concentration-dependent diffusivity
        new_CO2_array(i) = (dD_CO2(i) * d_CO2_array(i) + D_CO2(i) * d2_CO2_array(i)) * dt + CO2_array(i);   
    end
    for i = 2:nnodes-1
        H2O_array(i) = new_H2O_array(i);
        CO2_array(i) = new_CO2_array(i);
    end
% Increase time by time increment 
time = time + dt;

% Boundary conditions - vapor buffered 
% Pw_buffered = H2Ov_i/100 * P;                                          
% Pc_buffered = CO2v_i/100 * P; 
% H2O_array(nnodes) = (354.94 * Pw_buffered^0.5 + 9.623*Pw_buffered - 1.5223 * Pw_buffered^1.5)/TK + 0.0012439 * Pw_buffered^1.5 + Pc_buffered * ( (-1.084*10^-4) * Pw_buffered^0.5 - (1.362 * 10^-5) * Pw_buffered);
% CO2_array(nnodes) = Pc_buffered * (5668 - 55.99 * Pw_buffered)/TK + Pc_buffered * (0.4133 * Pw_buffered^0.5 + (2.041 * 10^-3) * Pw_buffered^1.5);

% Boundary condition - no flux
H2O_array(nnodes) = H2O_array(nnodes-1);
CO2_array(nnodes) = CO2_array(nnodes-1);

% Boundary condition at node 1 
% 1st and 2nd concentrations for H2O and CO2 arrays (mass fraction)
H2O_1_mass_fraction = H2O_array(1)/100;
H2O_2_mass_fraction = H2O_array(2)/100;
CO2_1_mass_fraction = CO2_array(1)/(10^6);
CO2_2_mass_fraction = CO2_array(2)/(10^6);

% 1st and 2nd concentrations for H2O and CO2 arrays (moles/m^3)
H2O_1_moles_per_meter_cubed = (H2O_1_mass_fraction * 1000 * rho)/H2O_molecular_mass;
H2O_2_moles_per_meter_cubed = (H2O_2_mass_fraction * 1000 * rho)/H2O_molecular_mass;
CO2_1_moles_per_meter_cubed = (CO2_1_mass_fraction * 1000 * rho)/CO2_molecular_mass;
CO2_2_moles_per_meter_cubed = (CO2_2_mass_fraction * 1000 * rho)/CO2_molecular_mass;

% H2O and CO2 flux (moles/(m^2 * s)) using Fick's first law
H2O_flux = D_H2O(1)/1e12 * (H2O_2_moles_per_meter_cubed - H2O_1_moles_per_meter_cubed)/(dx/10^6);
CO2_flux = D_CO2(1)/1e12 * (CO2_2_moles_per_meter_cubed - CO2_1_moles_per_meter_cubed)/(dx/10^6);

% Moles H2O and CO2 in bubble (moles)
H2O_moles_bubble = (H2O_flux * 4 * pi * radius_meters^2) * dt + H2O_moles_bubble; % Moles
CO2_moles_bubble = (CO2_flux * 4 * pi * radius_meters^2) * dt + CO2_moles_bubble; % Moles

% Calculate new H2O and CO2 vapor phases 
H2Ov = H2O_moles_bubble/(H2O_moles_bubble + CO2_moles_bubble)*100; % Mole percent of H2O in vapor
Xw = H2Ov/100; % Mole fraction of H2O in vapor
Pw = P * Xw; % Partial pressure of water in vapor
CO2v = CO2_moles_bubble/(H2O_moles_bubble + CO2_moles_bubble)*100;
Xc = CO2v/100; % Mole fraction of CO2 in vapor
Pc = P * Xc; % Partial pressure of CO2 in vapor

% Calculate new H2O and CO2 concentrations at node 1 (Liu et al., 2005)
H2O_node = (354.94 * Pw^0.5 + 9.623*Pw - 1.5223 * Pw^1.5)/TK + 0.0012439 * Pw^1.5 + Pc * ( (-1.084*10^-4) * Pw^0.5 - (1.362 * 10^-5) * Pw); 
CO2_node = Pc * (5668 - 55.99 * Pw)/TK + Pc * (0.4133 * Pw^0.5 + (2.041 * 10^-3) * Pw^1.5);

% Update H2O and CO2 concentrations at node 1 based on new vapor composition
H2O_array(1) = H2O_node; % Redundant
CO2_array(1) = CO2_node; % Redundant  

% %Storing the BCs from each timestep
% H2O_interface(iter) = H2O_node;
% CO2_interface(iter) = CO2_node;
% H2O_middle(iter) = H2O_array(nnodes);
% CO2_middle(iter) = CO2_array(nnodes);
% iter = iter + 1;
end
end
