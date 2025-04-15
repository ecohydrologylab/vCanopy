function [Constants,EnergyOptions,LeafBoundaryLayer] = callConstantsInput(Options)

%% Constants
Constants.R = 8.314/1000; % Universal gas constant [kJ mol-1 K-1]
Constants.convert = 1E6/(2.35E5); % [u mole J-1] Convert radiation from [W m-2] to [u mol m-2 s-1]
Constants.Boltzman = 5.67E-8; % Stefan Boltzman's constant [J K-1]
Constants.Cp = 29.3; % Specific heat capacity of air [J K-1 mol-1]
Constants.soilK = 3.9; % Soil thermal conductivity from Drewry
Constants.vonk = 0.41; % von Karman constant
Constants.rhoW = 1000; % density of water [kg m-3]
Constants.airdensity = 1.2923*1000/28.97; % [mol/m3]

LeafBoundaryLayer.cForced = 4.322*1E-3; % 1.2035*1E-3 for coniferous shoot
LeafBoundaryLayer.cFree = 1.6361*1E-3; % 0.86691*E-3 for coniferous shoot
LeafBoundaryLayer.leafDimension = 0.08; % Leaf dimension [m] 

EnergyOptions.isotropy = Options.isotropy; % Uniform overcast sky
EnergyOptions.switch = Options.LEBSwitch; % Energy balance switch in each layer On: 1; Off: 0
EnergyOptions.mixingSwitch = Options.mixingSwitch; % MicroEnvironment On: 1; Off: 0
EnergyOptions.fixedMicro = Options.fixedMicro; % This will fix the microenvironment/longwave but allow the leaf energu balance to work
EnergyOptions.longwaveSwitch = Options.longwaveSwitch; 
EnergyOptions.aNIRLeaf = Options.aNIR; % PAR absorptivity of leaf [-]
EnergyOptions.aPARLeaf = Options.aPAR; % NIR absorptivity of leaf [-]
EnergyOptions.epsilonLeaf = Options.epsilonLeaf; % Emissivity leaf [-]
EnergyOptions.aPARSoil = 0.8; % PAR absorptivity of soil [-]
EnergyOptions.aNIRSoil = 0.8; % NIR absorptivity of soil [-]
EnergyOptions.epsilonSoil = 0.9; % Emissivity soil [-]
EnergyOptions.LEfactor = 1.0; % Stomata on one side of leaf [-]
EnergyOptions.LWfactor = 2.0; % Emission from two sides of leaf [-]
EnergyOptions.Hfactor = 1.0; % Sensible heat loss from two sides of leaf [-]
end