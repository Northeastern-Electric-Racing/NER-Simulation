clear all;
close all;

%variables for simulation accuracy
global params;

params.numElements = 160; %MUST BE AN EVEN NUMBER
params.timeStep = 1e-4; %seconds. Needs to be smaller 1e-2 to reduce radiation overshoot
params.simTime = 20; %seconds

params.numSteps = params.simTime / params.timeStep;

%constants
params.stefBoltzConst = 5.67037e-8;%W⋅m−2⋅K−4, Stefan–Boltzmann Constant
params.gravity = 9.8; %m/s^2

%Electrical Condition
params.current = 20; %amps

%Fuse Properties
params.fuseBaselineDensity = 8900;%kg/m^3
params.fuseSpecificHeat = 385; %J/kg*C nickel = 456
params.fuseEmissivity = 0.7;%Nickel, not polished = 0.07
params.fuseThermalExpansionCoeff = 1.4e-5; %1/C
params.meltingTemp = 1435 + 273.15; %Kelvin

%Fuse Dimensions
params.fuseWidth = 0.437 * 1e-3; %units of calculation are meters, but enter number in MM!
params.fuseThickness =  0.437 * 1e-3; %units of calculation are meters, but enter number in MM!
params.fuseLength = 61 * 1e-3; %m

%Air Characteristics
params.ambientTemp = 303.15; %K
params.airThermalConductivity = 0.026; %mw/mK
params.airSpecificHeat = 1; %kJ/kg
params.airDynamicViscosity = 1.81*10e-5; %kg/(m*s)
params.airDensity = 1.225 ;%kg/m^3
        
params.prandtlNumber = 0.7 ;

params.elementLength = params.fuseLength / params.numElements;

params.midpointElement = params.numElements / 2; %for plotting midpoint temperature

%Progress Bar
progressBar = waitbar(0);

fuseTemps = repmat (params.ambientTemp, 1, params.numElements); %Initial Condition
fuseTemps = ApplyBoundaryConditions(fuseTemps);

loggedTemps(1,:) = fuseTemps;
loggedRadiationCooling(1,:) = RadiationCooling(fuseTemps);
loggedResistiveHeating(1,:) = ResistiveHeating(fuseTemps);
loggedConductionCooling(1,:) = ConductionCooling(fuseTemps);
loggedConvectionCooling(1,:) = ConvectionCooling(fuseTemps);

for iter = 1 : params.numSteps

    loggedConvectionCooling(iter + 1,:) = ConvectionCooling(fuseTemps);
    dEnergy = ResistiveHeating(fuseTemps) + ConductionCooling(fuseTemps); % - (ConvectionCooling(fuseTemps)) - RadiationCooling(fuseTemps) %  %Conduction MUST be positive, sign is accounted for in the conduction cooling function! 
    dTemp = EnergyToTemp(dEnergy, fuseTemps);
    fuseTemps = fuseTemps + dTemp;
    fuseTemps = ApplyBoundaryConditions(fuseTemps);
%     fuseTemps = ApplyContinuityConditions(fuseTemps);
    loggedTemps(iter + 1, :) = fuseTemps;
    loggedRadiationCooling(iter + 1, :) = RadiationCooling(fuseTemps);
    loggedResistiveHeating(iter + 1, :) = ResistiveHeating(fuseTemps);
    loggedConductionCooling(iter + 1,:) = ConductionCooling(fuseTemps);


    
    DisplayPercentCompletion(iter, progressBar);
end



figure; %Temperatures of the middle element (50) during SimTime
plot(loggedTemps(:,params.midpointElement),'o');
title('Midpoint Temp vs. Time');
xlabel('Time (0.01seconds)');
ylabel('Temperature, K');
%yline(params.meltingTemp, 'r-', 'Melting Point', 'LineWidth', 2);


figure; %Temperatures along fuse at end of SimTime
plot(loggedTemps(end,:),'o');
title('Solution Profile');
xlabel('Time (0.01seconds)');
ylabel('Temperature, K');

figure; %Temperatures along fuse at end of SimTime
plot(loggedResistiveHeating(:, params.midpointElement),'x');
hold on;
%plot(loggedRadiationCooling(:, params.midpointElement),'x');
hold on;
plot(loggedConductionCooling(:, params.midpointElement),'o');
%hold on;
%plot(loggedConvectionCooling(:, 50),'o');
title('Energy By Heat Transfer Type');
xlabel('Time (0.01seconds)');
ylabel('Energy, Joules');
legend( 'resistive heating', 'conduction'); %'Resistive Heating',



function dEnergy = ResistiveHeating(currentTemps)
    global params;
    
    for iter = 1 : params.numElements
       dEnergy(iter) = (((params.current ^ 2) * params.elementLength * GetElectricalResistivity(currentTemps(iter)))...
                     / ( GetCrossSectionalArea(currentTemps(iter)))) * params.timeStep;
    end
end

function dEnergy = ConductionCooling(currentTemps)
   global params;
   global leftConduction;
   global rightConduction;
   
   for iter = 1 : params.numElements 
       global leftConduction;
       
       if iter == 1 
          leftConduction(1) = 0; %There is no conduction from the left for the first element
       else
          leftConduction(iter) = (((GetThermalConductivity(currentTemps(iter)) * GetCrossSectionalArea(currentTemps(iter))...
                      * (currentTemps(iter - 1) - currentTemps(iter))) / params.elementLength) * params.timeStep); %
       end
   end
   
   for iter = 1 : (params.numElements) 
      global rightConduction; 
       
      if iter == params.numElements 
          rightConduction(params.numElements) = 0; %There is no conduction from the right for the last element
      else
          rightConduction(iter) = (((GetThermalConductivity(currentTemps(iter)) * GetCrossSectionalArea(currentTemps(iter))...
                      * (currentTemps(iter + 1) - currentTemps(iter))) / params.elementLength) * params.timeStep);
      end 
   end 
   
   dEnergy = leftConduction + rightConduction; 
   
   for iter = 1 : params.numElements
        if dEnergy(iter) > 0
            dEnergy(iter) = 0;
        end
   end    
end

function dEnergy = RadiationCooling(currentTemps)
   global params;
   
   for iter = 1 : params.numElements
       dEnergy(iter) = (params.fuseEmissivity * params.stefBoltzConst ...
                     * ((currentTemps(iter) ^ 4) - (params.ambientTemp ^ 4)) ...
                     * ((2 * GetFuseWidth(currentTemps(iter)) + (2 * GetFuseThickness(currentTemps(iter))) * params.elementLength))) ...
                     * params.timeStep; 

   end 
end 
 
function dEnergy = ConvectionCooling(currentTemps)
   global params;
   
%    RayleighNumbers = (GetGrasofNumbers(currentTemps) * params.prandtlNumber);
   nusseltNumbers = ((4 / 3) * ((GetGrasofNumbers(currentTemps) / 4) .^ (0.25)) * 0.499);
   %   nusseltNumbers = (0.68 + ((0.67 * (RayleighNumbers .^(1/4))) / (1 + (0.492 / params.prandtlNumber) ^ (9/16)) ^ (4/9)));
   convectionCoefficient = ((params.airThermalConductivity .* nusseltNumbers) ./ GetFuseWidth(currentTemps));
   dEnergy = (convectionCoefficient .* (2 .* GetFuseWidth(currentTemps) .* params.elementLength) .* (currentTemps - params.ambientTemp) .* params.timeStep); %h*A(surface, both sides)*deltaT*time
end
 
function GrasofNumbers = GetGrasofNumbers(currentTemps)
    global params;
    
    for iter = 1 : params.numElements
        GrasofNumbers(iter) = (params.gravity * GetAirThermalExpansionCoeff(currentTemps(iter)) * (currentTemps(iter) - params.ambientTemp)...
                            * (GetFuseWidth(currentTemps(iter)) ^ 3)) / (GetAirKinematicViscosity(currentTemps(iter)) ^ 2) ;
    end 
end 



function fuseTemps = ApplyContinuityConditions(fuseTemps)
    global params;
    
    for iter = 2 : params.midpointElement
        if fuseTemps(iter - 1) > fuseTemps(iter)
            fuseTemps(iter) = fuseTemps(iter - 1);
        end
    end 
    
    for iter = params.midpointElement : (params.numElements - 1)
        if fuseTemps(iter) < fuseTemps(iter + 1)
            fuseTemps(iter) = fuseTemps(iter + 1);
        end
    end 
    
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
    airThermalConductivity = ((-2e-0))%%%%% * filmTemp ^ 2) + (0.0819 * filmTemp) + 2.9501);
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

function fuseDensity = GetDensity(currentTemps)
    global params;
    fuseDensity = repmat (8900, 1, params.numElements); %Initial Condition

    %fuseDensity = params.fuseBaselineDensity / (1 + ((params.fuseThermalExpansionCoeff ^ 2) * (currentTemps - params.ambientTemp))); 
end

function dTemp = EnergyToTemp(energy, fuseTemps)
    global params;

    dTemp = energy ./ (GetDensity(fuseTemps) .* GetCrossSectionalArea(fuseTemps) .*  params.elementLength * params.fuseSpecificHeat); %Delta T = Energy/density*Volume*SpecificHeat
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