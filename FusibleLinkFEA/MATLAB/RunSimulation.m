clear all;
close all;

%variables for simulation accuracy
global numElements;
global timeStep;
global simTime;
numElements = 100;
timeStep = 1e-3; %seconds
simTime = 1; %seconds

numSteps = simTime / timeStep;

%constants
stefBoltzConst = 5.67037e-8;%W⋅m−2⋅K−4, Stefan–Boltzmann Constant
gravity = 9.8; %m/s^2

%Electrical Condition
global current;
current = 100; %amps

%Fuse Properties
global fuseBaselineDensity;
global fuseThermalExpansionCoeff;
global fuseSpecificHeat;
fuseBaselineDensity = 8900;%kg/m^3
fuseSpecificHeat = 456; %J/kg*C
fuseEmissivity = 0.07;%Nickel, not polished
fuseElectricalConductivity = 1.43e7; %S/m
fuseThermalExpansionCoeff = 1.4e-5; %1/C

%Fuse Dimensions
global fuseWidth;
global fuseThickness;
fuseWidth = 0.01; %m
fuseThickness = 0.0075; %m
fuseLength = 10; %m

%Air Characteristics
global ambientTemp;
ambientTemp = 30; %K
airThermalConductivity = 0.026; %mw/mK
airSpecificHeat = 1; %kJ/kg
airDynamicViscosity = 1.81*10e-5; %kg/(m*s)
airVolumetricExpansionCoeff = 1/300 ;%Approximated by 1/Tfilm (ABSOLUTE, in Kelvin)
airDensity = 1.225 ;%kg/m^3

%Calculated Values
% rayleightNumber = (gravity*airVolumetricExpansionCoeff*(TEMP - ambientTemp)...
%                  *(fuseWidth^3)*(airDensity^2)*airSpecificHeat)...
%                  /(airDynamicViscosity*airThermalConductivity);             
prandtlNumber = 0.7 ;

global elementLength; 
elementLength = fuseLength / numElements;

global fuseTemps;
fuseTemps = repmat (ambientTemp, 1, numElements); %initialCondition
loggedTemps(1,:) = fuseTemps;
for iter = 1 : numSteps
    dEnergy = ResistiveHeating(fuseTemps) + ConvectionCooling(fuseTemps)+ ConductionCooling(fuseTemps);
    dEnergy = ApplyBoundaryConditions(dEnergy);
    dTemp = EnergyToTemp(dEnergy);
    fuseTemps = fuseTemps + dTemp;
    loggedTemps(iter + 1, :) = fuseTemps; 
end

function dEnergy = ResistiveHeating(currentTemps)
    global numElements;
    global current;
    global elementLength;
    global timeStep;
    for iter = 1:numElements
       dEnergy(iter) = ((current ^ 2) * elementLength) / GetConductivity(currentTemps(iter))...
                       * GetCrossSectionalArea(currentTemps(iter)) * timeStep;
    end
end

function dEnergy = ConvectionCooling(currentTemps)
    global numElements;
    dEnergy = zeros(1, numElements);
end

function dEnergy = ConductionCooling(currentTemps)
   global numElements;
   dEnergy = zeros(1, numElements);
end
 
function fuseElectricalConductivity = GetConductivity(temp)
    fuseElectricalConductivity = (5e-5 * temp ^ 2) - (0.0504 * temp) + 70.051; %Defined in Property Approximations
end

function crossSectionalArea = GetCrossSectionalArea(temp)
    global fuseWidth;
    global fuseThickness;
    global ambientTemp;
    global fuseThermalExpansionCoeff;
    
    crossSectionalArea = fuseWidth * fuseThickness; %base cross sectional area
    crossSectionalArea = crossSectionalArea + crossSectionalArea * fuseThermalExpansionCoeff * (temp - ambientTemp); %thermal expansion term. coefficient of thermal expansion varies with temperature (do this)
end

function fuseDensity = GetDensity(currentTemps)
    global fuseBaselineDensity;
    global fuseThermalExpansionCoeff;
    global ambientTemp;
   
    fuseDensity = fuseBaselineDensity ./ ((fuseThermalExpansionCoeff ^ 3) .* (currentTemps - ambientTemp)); 
end

function dTemp = EnergyToTemp(energy)
    global elementLength;
    global fuseSpecificHeat;
    global fuseTemps;
    
    dTemp = energy .* GetDensity(fuseTemps) .* GetCrossSectionalArea(fuseTemps) .*  (elementLength * fuseSpecificHeat) ; %Delta T = Energy/density*Volume*SpecificHeat
end

function dEnergy = ApplyBoundaryConditions(dEnergy)
   global numElements;
   
   dEnergy(1) = 0;
   dEnergy(numElements) = 0;
end



