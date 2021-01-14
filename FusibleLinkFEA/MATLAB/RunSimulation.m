clear all;
close all;

%variables for simulation accuracy
global params;


params.numElements = 30; %MUST BE AN EVEN NUMBER
params.timeStep = 1e-1; %seconds. Needs to be smaller 1e-2 to reduce radiation overshoot
params.simTime = 50; %seconds

params.numSteps = params.simTime / params.timeStep;

%Constants
params.stefBoltzConst = 5.67037e-8;%W⋅m−2⋅K−4, Stefan–Boltzmann Constant
params.gravity = 9.8; %m/s^2

%Electrical Condition
params.current = 45; %amps

%Fuse Properties
params.fuseBaselineDensity = 8900;%kg/m^3
params.fuseSpecificHeat = 385; %J/kg*C nickel = 456
params.fuseEmissivity = 0.1; %Nickel, not polished = 0.07. %Copper, ~0.1
params.fuseThermalExpansionCoeff = 1.4e-5; %1/C
params.meltingTemp = 1435 + 273.15; %Kelvin

%Fuse Dimensions
params.fuseWidth = 0.810 * 1e-3; %units of calculation are meters, but enter number in MM!
params.fuseThickness =  0.810 * 1e-3; %units of calculation are meters, but enter number in MM!
params.fuseLength = 400 * 1e-3; %m

%Air Characteristics
params.ambientTemp = 303.15; %K
params.airThermalConductivity = 0.026; %mw/mK
params.airSpecificHeat = 1; %kJ/kg
params.airDynamicViscosity = 1.81*10e-5; %kg/(m*s)
params.airDensity = 1.225 ;%kg/m^3
        
params.prandtlNumber = 0.7 ;

%Element Dimensions
params.elementLength = params.fuseLength / params.numElements;
params.firstElementLength = params.elementLength / 2; %with small first and last elements
params.lastElementLength = params.elementLength / 2; %With small first and last elements
params.midpointElement = params.numElements / 2; %for plotting midpoint temperature

%Progress Bar
progressBar = waitbar(0);

%Fuse Temps Initial Condition
fuseTemps = repmat ((params.ambientTemp), 1, params.numElements); %Initial Condition

%Von Neumann Stability
diffusivity = (GetThermalConductivity(fuseTemps(1)) / (GetDensity(fuseTemps(1)) * params.fuseSpecificHeat));
vonNeumannStability = (diffusivity * params.timeStep) / (params.firstElementLength ^ 2);
if vonNeumannStability > 0.5
    disp("Unstable");
    pause; 
end 

loggedTemps(1,:) = fuseTemps;
loggedRadiationCooling(1,:) = RadiationCooling(fuseTemps);
loggedResistiveHeating(1,:) = ResistiveHeating(fuseTemps);
loggedConductionCooling(1,:) = ConductionCooling(fuseTemps);
loggedConvectionCooling(1,:) = ConvectionCooling(fuseTemps);

for iter = 1 : params.numSteps

    dEnergy = ResistiveHeating(fuseTemps) + ConductionCooling(fuseTemps) - (ConvectionCooling(fuseTemps)) - ( 0.5 * RadiationCooling(fuseTemps));  %  %Conduction MUST be positive, sign is accounted for in the conduction cooling function! 
    dTemp = EnergyToTemp(dEnergy, fuseTemps);
    fuseTemps = fuseTemps + dTemp;
    loggedTemps(iter + 1, :) = fuseTemps;
    loggedConvectionCooling(iter + 1,:) = ConvectionCooling(fuseTemps);
    loggedRadiationCooling(iter + 1, :) = RadiationCooling(fuseTemps);
    loggedResistiveHeating(iter + 1, :) = ResistiveHeating(fuseTemps);
    loggedConductionCooling(iter + 1,:) = ConductionCooling(fuseTemps);
    
    DisplayPercentCompletion(iter, progressBar);
end


figure; %Temperatures of the middle element (50) during SimTime
plot((0:params.timeStep:params.simTime),loggedTemps(:,params.midpointElement),'o');
title('Midpoint Temp vs. Time');
xlabel('Time (s)');
ylabel('Temperature (K)');
%yline(params.meltingTemp, 'r-', 'Melting Point', 'LineWidth', 2);


figure; %Temperatures along fuse at end of SimTime
plot(loggedTemps(end,:),'o');
title(['Solution Profile at t = ' num2str(params.simTime) 's']);
xlabel('Element Number');
ylabel('Temperature (K)');


%Energy by heat type wrt time
figure; 
%ON SAME AXIS
% plot((0:params.timeStep:params.simTime),loggedResistiveHeating(:, params.midpointElement),'x');
% hold on;
% plot((0:params.timeStep:params.simTime),loggedRadiationCooling(:,params.midpointElement),'x');
% hold on;
% plot((0:params.timeStep:params.simTime),loggedConductionCooling(:, params.midpointElement),'o');
% hold on;
% plot((0:params.timeStep:params.simTime),loggedConvectionCooling(:,params.midpointElement),'o');
% title('Energy By Heat Transfer Type');
% xlabel('Time (s)');
% ylabel('Energy (J)');
% legend('Resitive Heating', 'Radiation','Conduction','Convection');
%ON DIFFERENT AXES
EnergyYs = cat(2,loggedResistiveHeating(:,params.midpointElement),loggedRadiationCooling(:,params.midpointElement),loggedConductionCooling(:,params.midpointElement),loggedConvectionCooling(:,params.midpointElement));
titles = {'Energy By Resistive Heating','Energy By Radiation Cooling','Energy By Conduction Cooling','Energy By Convection Cooling'};
 for i=1:4
     subplot(2,2,i)
     plot((0:params.timeStep:params.simTime),EnergyYs(:,i),'o')
     title(titles{i})
     ylabel('Energy (J)')
     xlabel('Time (s)')
 end

figure; %3D plot of Temp vs. Time at every element
surf(1:params.numElements,(0:params.timeStep:params.simTime),loggedTemps,'FaceLighting','flat','EdgeAlpha','0.3')
xlabel('Element');
ylabel('Time (s)');
zlabel('Temperature (K)');
title('Temperature as Function of Time and Element (Lateral Fuse Position)');

function dEnergy = ResistiveHeating(currentTemps)
    global params;
    
    for iter = 1 : params.numElements 
     if iter == 1 
          
         dEnergy(iter) = (((params.current ^ 2) * params.firstElementLength * GetElectricalResistivity(currentTemps(iter)))...
                     / ( GetCrossSectionalArea(currentTemps(iter)))) * params.timeStep;
               
     elseif iter == params.numElements
         
         dEnergy(iter) = (((params.current ^ 2) * params.lastElementLength * GetElectricalResistivity(currentTemps(iter)))...
                     / ( GetCrossSectionalArea(currentTemps(iter)))) * params.timeStep;
                  
     else
         
         dEnergy(iter) = (((params.current ^ 2) * params.elementLength * GetElectricalResistivity(currentTemps(iter)))...
                     / ( GetCrossSectionalArea(currentTemps(iter)))) * params.timeStep;
     end 
    end
end


function dEnergy = ConductionCooling(currentTemps)
   global params;
   
       for iter = 1 : params.numElements 
          if iter == 1 

             dEnergy(iter) = (((GetThermalConductivity(currentTemps(iter)) * GetCrossSectionalArea(currentTemps(iter))...
                           * (currentTemps(iter + 1) - (2 * currentTemps(iter)) + params.ambientTemp)) / params.firstElementLength) * params.timeStep);

          elseif iter == params.numElements


             dEnergy(iter) = (((GetThermalConductivity(currentTemps(iter)) * GetCrossSectionalArea(currentTemps(iter))...
                           * (currentTemps(iter - 1) - (2 * currentTemps(iter)) + params.ambientTemp)) / params.lastElementLength) * params.timeStep);

          else

             dEnergy(iter) = (((GetThermalConductivity(currentTemps(iter)) * GetCrossSectionalArea(currentTemps(iter))...
                           * (currentTemps(iter + 1) - (2 * currentTemps(iter)) + currentTemps(iter - 1))) / params.elementLength) * params.timeStep);
          end 
      end 
end 


function dEnergy = RadiationCooling(currentTemps)
   global params;
   
   for iter = 1 : params.numElements 
     if iter == 1 
          
         dEnergy(iter) = (params.fuseEmissivity * params.stefBoltzConst ...
                     * ((currentTemps(iter) ^ 4) - (params.ambientTemp ^ 4)) ...
                     * ((2 * GetFuseWidth(currentTemps(iter)) + (2 * GetFuseThickness(currentTemps(iter))) * params.firstElementLength))) ...
                     * params.timeStep; 
               
     elseif iter == params.numElements
         
         dEnergy(iter) = (params.fuseEmissivity * params.stefBoltzConst ...
                     * ((currentTemps(iter) ^ 4) - (params.ambientTemp ^ 4)) ...
                     * ((2 * GetFuseWidth(currentTemps(iter)) + (2 * GetFuseThickness(currentTemps(iter))) * params.lastElementLength))) ...
                     * params.timeStep; 
                  
     else
         
         dEnergy(iter) = (params.fuseEmissivity * params.stefBoltzConst ...
                     * ((currentTemps(iter) ^ 4) - (params.ambientTemp ^ 4)) ...
                     * ((2 * GetFuseWidth(currentTemps(iter)) + (2 * GetFuseThickness(currentTemps(iter))) * params.elementLength))) ...
                     * params.timeStep;
     end 
   end 
end 
 
function dEnergy = ConvectionCooling(currentTemps)
   global params;

   %OLD Rayleigh and Nusselt Number Derivations
%    RayleighNumbers = (GetGrasofNumbers(currentTemps) * params.prandtlNumber);
%   nusseltNumbers = (0.68 + ((0.67 * (RayleighNumbers .^(1/4))) / (1 + (0.492 / params.prandtlNumber) ^ (9/16)) ^ (4/9)));
   
 for iter = 1 : params.numElements 
     
     if iter == 1 
          
        nusseltNumber = ((4 / 3) * ((GetGrasofNumber(currentTemps(iter)) / 4) .^ (0.25)) * 0.499);
        convectionCoefficient = ((params.airThermalConductivity .* nusseltNumber) / GetFuseWidth(currentTemps(iter)));
        dEnergy(iter) = (convectionCoefficient * (2 * GetFuseWidth(currentTemps(iter)) * params.firstElementLength) * (currentTemps(iter) - params.ambientTemp) .* params.timeStep); %h*A(surface, both sides)*deltaT*time
               
     elseif iter == params.numElements
         
        nusseltNumber = ((4 / 3) * ((GetGrasofNumber(currentTemps(iter)) / 4) .^ (0.25)) * 0.499);
        convectionCoefficient = ((params.airThermalConductivity .* nusseltNumber) / GetFuseWidth(currentTemps(iter)));
        dEnergy(iter) = (convectionCoefficient * (2 * GetFuseWidth(currentTemps(iter)) * params.lastElementLength) * (currentTemps(iter) - params.ambientTemp) .* params.timeStep); %h*A(surface, both sides)*deltaT*time                 
     else
         
        nusseltNumber = ((4 / 3) * ((GetGrasofNumber(currentTemps(iter)) / 4) .^ (0.25)) * 0.499);
        convectionCoefficient = ((params.airThermalConductivity .* nusseltNumber) / GetFuseWidth(currentTemps(iter)));
        dEnergy(iter) = (convectionCoefficient * (2 * GetFuseWidth(currentTemps(iter)) * params.elementLength) * (currentTemps(iter) - params.ambientTemp) .* params.timeStep); %h*A(surface, both sides)*deltaT*time     end 
     end
 end
end
 
function GrasofNumber = GetGrasofNumber(temp)
    global params;
    
    GrasofNumber = (params.gravity * GetAirThermalExpansionCoeff(temp) * (temp - params.ambientTemp)...
                            * (GetFuseWidth(temp) ^ 3)) / (GetAirKinematicViscosity(temp) ^ 2);
end 

function airDensity = GetAirDensity(temp)
    global params;
    
    filmTemp = ((temp + params.ambientTemp) / 2); %calculated at the film temperature
    airDensity = (352.77 / filmTemp);
end

function airSpecificHeat = GetAirSpecificHeat(temp)
    global params;
    
    filmTemp = ((temp + params.ambientTemp) / 2); %calculated at the film temperature
    airSpecificHeat = ((-2e-8 * (filmTemp ^ 2)) + (2e-5 * filmTemp) + 0.9502);
end

function airKinematicViscosity = GetAirKinematicViscosity(temp)
    global params;

    filmTemp = ((temp + params.ambientTemp) / 2); %calculated at the film temperature
    airKinematicViscosity = (((6e-5 * (filmTemp ^ 2)) + (0.0627 * filmTemp) - 8.6444) * 10^-6);

end

function airThermalConductivity = GetAirThermalConductivity(temp)
    global params;
    
    filmTemp = ((temp + params.ambientTemp) / 2); %calculated at the film temperature
    airThermalConductivity = ((-2e-0));%%%%% * filmTemp ^ 2) + (0.0819 * filmTemp) + 2.9501);
end 

function airVolumetricExpansionCoeff = GetAirThermalExpansionCoeff(temp)
    global params;
    
    filmTemp = ((temp + params.ambientTemp) / 2); %calculated at the film temperature
    airVolumetricExpansionCoeff = 1 / filmTemp; %Ideal gas
end 

function fuseElectricalResistivity = GetElectricalResistivity(temp)
    fuseElectricalResistivity = 1.724e-8; % copper
   % nickel = fuseElectricalResistivity = ((0.0004 * temp) + 0.0085 ) * 10e-5; % Linear Approximation
%   fuseElectricalResistivity = (((-3e-7 * (temp ^ 2)) + (0.0008 * temp) - 0.1369) * 10e-6); %2nd order polynomial, fails after 1500K!!
end

function fuseThermalConductivity = GetThermalConductivity(temp)
    fuseThermalConductivity = 385; %copper
    % Nickel = fuseThermalConductivity = (0.00005 * (temp^2)) - (0.0782 * temp) + 87.622;
    %fuseThermalConductivity = (5 * 10e-4 * (temp ^ 2)) - (0.0782 * temp) + 87.622;
end

function crossSectionalArea = GetCrossSectionalArea(temp)
    global params;
    
    crossSectionalArea = params.fuseWidth * params.fuseThickness; %base cross sectional area
    crossSectionalArea = crossSectionalArea * (1 + (params.fuseThermalExpansionCoeff ^ 2 * (temp - params.ambientTemp))); %thermal expansion term. coefficient of thermal expansion varies with temperature (do this)
end

function fuseWidth = GetFuseWidth(temp)
    global params;
    
    fuseWidth = params.fuseWidth * (1 + (params.fuseThermalExpansionCoeff * (temp - params.ambientTemp)));
end

function fuseThickness = GetFuseThickness(temp)
    global params;
    
    fuseThickness = params.fuseThickness * (1 + (params.fuseThermalExpansionCoeff * (temp - params.ambientTemp))); 
end

function fuseDensity = GetDensity(temp)
    global params;
%     fuseDensity = repmat (8900, 1, params.numElements); %Initial Condition

    fuseDensity = params.fuseBaselineDensity / (1 + ((params.fuseThermalExpansionCoeff ^ 2) * (temp - params.ambientTemp))); 
end

function dTemp = EnergyToTemp(energy, fuseTemps)
    global params;

    for iter = 1 : params.numElements 
     if iter == 1 
          
             dTemp(iter) = energy(iter) / (GetDensity(fuseTemps(iter)) * GetCrossSectionalArea(fuseTemps(iter)) *  params.firstElementLength * params.fuseSpecificHeat); %Delta T = Energy/density*Volume*SpecificHeat

               
     elseif iter == params.numElements
         
             dTemp(iter) = energy(iter) / (GetDensity(fuseTemps(iter)) * GetCrossSectionalArea(fuseTemps(iter)) *  params.lastElementLength * params.fuseSpecificHeat); %Delta T = Energy/density*Volume*SpecificHeat
     else
         
             dTemp(iter) = energy(iter) / (GetDensity(fuseTemps(iter)) * GetCrossSectionalArea(fuseTemps(iter)) *  params.elementLength * params.fuseSpecificHeat); %Delta T = Energy/density*Volume*SpecificHeat

     end 
    end 
end


function fuseTemps = ApplyBoundaryConditions(fuseTemps)
   global params;
   
   fuseTemps(1) = params.ambientTemp;
   fuseTemps(params.numElements) = params.ambientTemp;
end

function DisplayPercentCompletion(iteration, progressBar)
    global params;
    
    completion = iteration / params.numSteps;
    percentComplete = round(completion * 100 - 0.5); %round down to not show 100% if not complete.
    waitbar(completion, progressBar, string(percentComplete) + "% Complete"); % show in progress bar
    disp(string(percentComplete) + "% Complete"); %print to log
end