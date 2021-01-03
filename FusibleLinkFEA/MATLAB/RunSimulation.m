clear all;
close all;

%variables for simulation accuracy
global params;

params.numElements = 100;
params.timeStep = 1e-3; %seconds
params.simTime = 1; %seconds

params.numSteps = params.simTime / params.timeStep;

%constants
params.stefBoltzConst = 5.67037e-8;%W⋅m−2⋅K−4, Stefan–Boltzmann Constant
params.gravity = 9.8; %m/s^2

%Electrical Condition
params.current = 100; %amps

%Fuse Properties
params.fuseBaselineDensity = 8900;%kg/m^3
params.fuseSpecificHeat = 456; %J/kg*C
params.fuseEmissivity = 0.07;%Nickel, not polished
params.fuseElectricalConductivity = 1.43e7; %S/m
params.fuseThermalExpansionCoeff = 1.4e-5; %1/C

%Fuse Dimensions
params.fuseWidth = 0.01; %m
params.fuseThickness = 0.0075; %m
params.fuseLength = 10; %m

%Air Characteristics
params.ambientTemp = 30; %K
params.airThermalConductivity = 0.026; %mw/mK
params.airSpecificHeat = 1; %kJ/kg
params.airDynamicViscosity = 1.81*10e-5; %kg/(m*s)
params.airVolumetricExpansionCoeff = 1/300 ;%Approximated by 1/Tfilm (ABSOLUTE, in Kelvin)
params.airDensity = 1.225 ;%kg/m^3

%Calculated Values
% rayleightNumber = (gravity*airVolumetricExpansionCoeff*(TEMP - ambientTemp)...
%                  *(fuseWidth^3)*(airDensity^2)*airSpecificHeat)...
%                  /(airDynamicViscosity*airThermalConductivity);             
params.prandtlNumber = 0.7 ;

params.elementLength = params.fuseLength / params.numElements;

%Progress Bar
progressBar = waitbar(0);

fuseTemps = repmat (params.ambientTemp, 1, params.numElements); %Initial Condition
loggedTemps(1,:) = fuseTemps;
for iter = 1 : params.numSteps
    dEnergy = ResistiveHeating(fuseTemps) + ConvectionCooling(fuseTemps)+ ConductionCooling(fuseTemps);
    dEnergy = ApplyBoundaryConditions(dEnergy);
    dTemp = EnergyToTemp(dEnergy, fuseTemps);
    fuseTemps = fuseTemps + dTemp;
    loggedTemps(iter + 1, :) = fuseTemps;
    
    DisplayPercentCompletion(iter, progressBar);
end

function dEnergy = ResistiveHeating(currentTemps)
    global params;
    
    for iter = 1 : params.numElements
       dEnergy(iter) = ((params.current ^ 2) * params.elementLength)...
                     / (GetConductivity(currentTemps(iter)) * GetCrossSectionalArea(currentTemps(iter)))...
                     * params.timeStep;
    end
end

function dEnergy = ConvectionCooling(currentTemps)
    global params;
    dEnergy = zeros(1, params.numElements);
end

function dEnergy = ConductionCooling(currentTemps)
   global params;
   dEnergy = zeros(1, params.numElements);
end
 
function fuseElectricalConductivity = GetConductivity(temp)
    fuseElectricalConductivity = (5e-5 * temp ^ 2) - (0.0504 * temp) + 70.051; %Defined in Property Approximations
end

function crossSectionalArea = GetCrossSectionalArea(temp)
    global params;
    
    crossSectionalArea = params.fuseWidth * params.fuseThickness; %base cross sectional area
    crossSectionalArea = crossSectionalArea + crossSectionalArea * params.fuseThermalExpansionCoeff * (temp - params.ambientTemp); %thermal expansion term. coefficient of thermal expansion varies with temperature (do this)
end

function fuseDensity = GetDensity(currentTemps)
    global params;
  
    fuseDensity = params.fuseBaselineDensity ./ (1 + ((params.fuseThermalExpansionCoeff ^ 3) .* (currentTemps - params.ambientTemp))); 
end

function dTemp = EnergyToTemp(energy, fuseTemps)
    global params;

    dTemp = energy ./ (GetDensity(fuseTemps) .* GetCrossSectionalArea(fuseTemps) .*  params.elementLength .* params.fuseSpecificHeat); %Delta T = Energy/density*Volume*SpecificHeat
end

function dEnergy = ApplyBoundaryConditions(dEnergy)
   global params;
   
   dEnergy(1) = 0;
   dEnergy(params.numElements) = 0;
end

function DisplayPercentCompletion(iteration, progressBar)
    global params;
    
    completion = iteration / params.numSteps;
    percentComplete = round(completion * 100 - 0.5); %round down to not show 100% if not complete.
    waitbar(completion, progressBar, string(percentComplete) + "% Complete"); % show in progress bar
    disp(string(percentComplete) + "% Complete"); %print to log
end