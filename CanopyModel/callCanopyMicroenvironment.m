function [Canopy] = callCanopyMicroenvironment(EnergyOptions,Constants,Weather,Soil...
    ,Canopy,Error,CanopyIteration)

% Compute canopy wind and scalar profile
% Canopy = callCanopyScalarProfile2(EnergyOptions,Weather,Canopy,Soil,Constants); % Compute temperature profile
Canopy = callCanopyScalarProfile5('temperature',Weather,Canopy,Soil,Constants); % Compute temperature profile
Canopy = callCanopyScalarProfile5('ea',Weather,Canopy,Soil,Constants); % Compute vapour profile
Canopy = callCanopyScalarProfile5('ca',Weather,Canopy,Soil,Constants); % Compute CO2 profile
Canopy = callCanopyRHProfile(Canopy); % Compute RH profile

if EnergyOptions.switch == 1
    Canopy.windProfile = Error.relax*Canopy.windProfile + ...
        (1-Error.relax)*CanopyIteration.windProfile;
    Canopy.tAirProfile = Error.relax*Canopy.tAirProfile +...
        (1-Error.relax)*CanopyIteration.tAirProfile;
    Canopy.eaProfile = Error.relax*Canopy.eaProfile + ...
        (1-Error.relax)*CanopyIteration.eaProfile;
    Canopy.caProfile = Error.relax*Canopy.caProfile + ...
        (1-Error.relax)*CanopyIteration.caProfile;

    Canopy.sunTleaf = (Error.relax*Canopy.sunTleaf + ...
        (1-Error.relax)*CanopyIteration.sunTleaf);
    Canopy.shadeTleaf = (Error.relax*Canopy.shadeTleaf + ...
        (1-Error.relax)*CanopyIteration.shadeTleaf);
else
    Canopy.sunTleaf = Canopy.sunTleaf;
    Canopy.shadeTleaf = Canopy.shadeTleaf;
end
end