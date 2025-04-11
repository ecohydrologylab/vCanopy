%% Function to get model input parameters
function [Photosynthesis,Stomata,Weather,Constants,LeafBoundaryLayer,...
    EnergyOptions] = callInputData(fileName)

%% Constants
Constants.R = 8.314/1000; % Universal gas constant [kJ mol-1 K-1]
Constants.convert = 4.57;%1E6/(2.35E5); % [u mole J-1] Convert radiation from [W m-2] to [u mol m-2 s-1]
Constants.Boltzman = 5.67E-8; % Stefan Boltzman's constant [J K-1]
Constants.Cp = 29.3; % Specific heat capacity of air [J K-1 mol-1]
Constants.soilK = 3.9; % Soil thermal conductivity from Drewry
Constants.vonk = 0.41; % von Karman constant
Constants.rhoW = 1000; % density of water [kg m-3]
Constants.airdensity = 1.2923*1000/28.97; % [mol/m3]

LeafBoundaryLayer.cForced = 4.322*1E-3; % 1.2035*1E-3 for coniferous shoot
LeafBoundaryLayer.cFree = 1.6361*1E-3; % 0.86691*E-3 for coniferous shoot
LeafBoundaryLayer.leafDimension = 0.08; % [m] Nikolov 1995

EnergyOptions.switch = 1; % Energy balance switch in each layer
EnergyOptions.aNIRLeaf = 0.23; % Drewry 2010
EnergyOptions.aPARLeaf = 0.80; % Drewry 2010
EnergyOptions.longwaveSwitch = 1 ;
EnergyOptions.fixedMicro = 1;
EnergyOptions.epsilonLeaf = 0.94; % Emissivity leaf [-]
EnergyOptions.LEfactor = 1.0; % Stomata on one side of leaf [-]
EnergyOptions.LWfactor = 2.0; % Emission from two sides of leaf [-]
EnergyOptions.Hfactor = 1.0; % Sensible heat loss from two sides of leaf [-]

%% Photosynthesis and stomatal parameters, and weather data
inputData = readtable(fileName,'ReadVariableNames',true); % Read input data in table format
Photosynthesis = inputData(1:end,{'plant','leafID','vpr','vcmax25','jmax25','vpmax25', ...
    'rd25','sco25','theta','gbs','alpha','x'});
Stomata = inputData(1:end,{'slope','intercept'});
Weather = inputData(1:end,{'PAR','NIR','long','ea','tAir','wind','ca','O2',...
    'controlTemp','pressure'});
clear inputData

Stomata.model(:) = "BB";
if Stomata.model(1) == "BB"
    %% Ball Berry Model
    Stomata.slope = zeros(length(Stomata.slope),1) + 3;% 7;%
    Stomata.intercept = zeros(length(Stomata.slope),1) + 0.08;% 0.008;%
elseif Stomata.model(1) == "BM"
    %% Belinda Model
    Stomata.slope = zeros(length(Stomata.slope),1) + 0.44;%
    Stomata.intercept = zeros(length(Stomata.slope),1) + 0.07;%
end
Stomata.fsv = zeros(length(Stomata.slope),1) + 1;
Stomata.s = zeros(length(Stomata.slope),1) + 0.71; % ratio of stomatal density on abaxial and adaxian side
Stomata.Mb = 0.5.*(1 + Stomata.s).^2./(1 + Stomata.s.^2);

%% Calculating the absoprbed radiation from incident radiation
% Weather.PAR = Weather.PAR*EnergyOptions.aPARLeaf; % Absorbed PAR
% Weather.NIR = Weather.NIR*EnergyOptions.aNIRLeaf; % Absorbed NIR
Weather.skyemissivity =  0.926 + 1.006*(0.001*Weather.ea./...
    (Weather.tAir+273.15)).^(1.067); % Sky emittivity modified from Brutsaert (1982)
Weather.long = EnergyOptions.LWfactor*EnergyOptions.epsilonLeaf*...
    Weather.skyemissivity.*Constants.Boltzman.*(273.15 + Weather.tAir).^4; % Absorbed longwave from sky

Weather = table2struct(Weather);
Stomata = table2struct(Stomata);
Photosynthesis = table2struct(Photosynthesis);
end