% 1D diffusion of H2O and CO2 through a rhyolitic melt to one bubble during constant decompression 
% Optimized for decompression rate (dPdt), fragmentation pressure (P_f), and bubble radius (radius)

% Original MATLAB script created by Myers et al. (2018); this version modified by Hosseini et al. (in revision)

clear all % Clear any existing variables in the workspace
close all % Close all figures
clc % Clear command window

%% ASSUMPTIONS
% All diffusion occuring in one dimension from interior to outlet of embayment
% No volume change of bubble due to mass influx (i.e., melt in embayment is incompressible)
% Melt and bubble are initially in equilibrium 
% Isothermal
% Cartesian coordinates

%% ENTER VARIABLE INPUT PARAMETERS
% The following are all of the input parameters that can be adjusted for each sample.
              
P_i = 150; % Initial pressure [MPa]
lengthi = 90; % Distance from interior of embayment to glass-bubble interface [um]
H2O_i = 3; % Initial H2O concentration from melt inclusion or embayment interior plateau [wt. %]
rho = 2350; % Density of rhyolite [kg/m^3]-estimated density of obsidian from White and Harris (1997) 
TC = 780; % Temperature [deg C]
errorH2O = 0.32 ; % Analytical uncertainty in H2O measurement [wt. %]
errorCO2 = 22; % Analytical uncertainty in CO2 measurement [ppm]

%% ENTER MEASURED H2O AND CO2 CONCENTRATIONS
% The following arrays store the measured H2O and CO2 concentrations, as well as the distance along the 
% embayment at which measurements were made. Note: distances and measurements go from outlet/bubble to interior

% SAMPLE NAME: E3-1
MeasH2OCon = [2.01 2.13 2.27 2.29 2.32 2.33 2.28 2.06]; % Measured H2O concentration profile
MeasCO2Con = [27.28 39.00 50.74 58.56 74.19 81.99 81.97 89.67]; % Measured CO2 concentration profile
MeasDist = [5 15 25 35 45 55 65 75]; % Distance along profile, from outlet to interior

%% ENTER RANGE OF DPDT, FRAGMENTATION PRESSURES, AND BUBBLE RADII TO ITERATE THROUGH. MIDDLE NUMBER IS STEP SIZE.

% dPdt = (0.005:0.001:0.05); % dPdt to iterate through to find best fit to data [MPa/s]
% P_f = (0:5:50); % P_f to iterate through to find best fit to data [MPa]
% radius = (0.1:10:100.1); % Radii to iterate through to find best fit to data

dPdt = (0.008); % dPdt to iterate through to find best fit to data [MPa/s]
P_f = (30); % P_f to iterate through to find best fit to data [MPa]
radius = (0.1); % Radii to iterate through to find best fit to data

% Supplemental table from Myers et al. (2018): bubble radius 30-40 = 1 wt. %, 60-80 = 2 wt. %, 100-120 = 3 wt. %

%% CALCULATE DPDT, P_f, AND RADIUS THAT GENERATE THE LEAST MISFIT WITH THE MEASURED DATA.

misfit = zeros(length(dPdt),length(P_f),length(radius));

% Iterate through array of dPdt, P_f, and radius
for idpdt = 1:length(dPdt)
    disp(['Running Decompression ',num2str(idpdt),' of ',num2str(length(dPdt))])
    
    for ipress = 1:length(P_f)
        disp(['Running Fragmentation Pressure ',num2str(ipress),' of ',num2str(length(P_f))])
    
        for ir = 1:length(radius)
            disp(['Running Radius ',num2str(ir),' of ',num2str(length(radius))])
    
        % Call the 1D diffusion function (diffusion_function.m)
        [nodes, H2O_array, CO2_array] = diffusion_function(P_i,P_f(ipress),lengthi,radius(ir),H2O_i,rho,TC,dPdt(idpdt));

        % Analyze the output
        misfitModel = 0;
        for im = 1:length(MeasDist)
            
            % Iterate through measured H2O and CO2 concentrations and distances
            d = MeasDist(im); 
            H2O = MeasH2OCon(im);
            CO2 = MeasCO2Con(im);
            
            % Find the node that is closest to our measured distance
            % Tilde indicates that minimum value is discarded, and only the index of the minimum value is used
            [~,iN] = min(abs(nodes - d));
            H2O_model = H2O_array(iN);
            CO2_model = CO2_array(iN);
            
            misfitModel = misfitModel + ((((H2O_model - H2O) / errorH2O).^2) + ((CO2_model - CO2) / errorCO2).^2) / length(MeasCO2Con);
        end
        misfit(idpdt,ipress,ir) = misfitModel;
        end
    end
end

[~,loc] = min(misfit(:));
[IDPDT,IPF,IR] = ind2sub(size(misfit),loc);

[nodes, H2O_array, CO2_array] = diffusion_function(P_i,P_f(IPF),lengthi,radius(IR),H2O_i,rho,TC,dPdt(IDPDT));

disp(['Minimum misfit = ', num2str(min(misfit(:)))])
disp(['Decompression rate = ',num2str(dPdt(IDPDT))])
disp(['Fragmentation pressure = ',num2str(P_f(IPF))])
disp(['Radius = ',num2str(radius(IR))])

if min(misfit(:)) < 1
    disp(['Good fit achieved with parameters.'])
else
    disp(['Good fit not achieved with parameters.']) 
end

%% PLOT DIFFUSION PROFILES

figure(1)
clf
hold on

h = subplot(2,1,1);
meH2O = plot(MeasDist,MeasH2OCon,'s'); % Plot the measured H2O profile
error = ones(1,length(MeasH2OCon))*errorH2O; % Initiate an array containing the H2O analytical uncertainty
e = errorbar(MeasDist,MeasH2OCon,error,'s','MarkerEdgeColor','k','MarkerFaceColor','#EDB120');
e.Color = 'black';
e.MarkerSize = 10;
hold on
plot(nodes,H2O_array,'k-','LineWidth',3) % Plot the modeled H2O profile
hold on
xlim([0 lengthi])
ylim([0 max(MeasH2OCon)+1])
uistack(e,'top');
xlabel('Distance (\mum)','fontweight','bold','FontSize',12)
ylabel('H_{2}O (wt. %)','fontweight','bold','FontSize',12)

i = subplot(2,1,2);
meCO2 = plot(MeasDist,MeasCO2Con,'s'); % Plot the measured CO2 profile
error = ones(1,length(MeasCO2Con))*errorCO2; % Initiate an array containing the CO2 analytical uncertainty
e = errorbar(MeasDist,MeasCO2Con,error,'s','MarkerEdgeColor','k','MarkerFaceColor','#EDB120');
e.Color = 'black';
e.MarkerSize = 10;
hold on
plot(nodes,CO2_array,'k-','LineWidth',3) % Plot the modeled CO2 profile
hold on
xlim([0 lengthi])
ylim([0 max(MeasCO2Con)+200])
uistack(e,'top');
xlabel('Distance (\mum)','fontweight','bold','FontSize',12)
ylabel('CO_{2} (ppm)','fontweight','bold','FontSize',12)