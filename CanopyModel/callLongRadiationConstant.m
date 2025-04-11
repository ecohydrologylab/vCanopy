function [Soil,Canopy] = callLongRadiationConstant(Constants,EnergyOptions,Weather,Soil,Canopy)


% Based on the same principal as the leaf level paper
% Top leaf longwave: sky + leaf just below
% Middle leaf longwave: leaf just above + leaf just below
% Lowest leaf longwave: leaf just above  + soil

emissivitySoil = EnergyOptions.epsilonSoil;
leafAbsorptivity = EnergyOptions.epsilonLeaf; % Leaf absorptivity Changed from 0.96
soilAbsorptivity = EnergyOptions.epsilonSoil; % Changed from 0.8

coeff = [0.9258    1.0059    1.0677];
skyEmissivity = coeff(1) +  coeff(2)*(Canopy.eaProfile/1000./...
    (Canopy.tAirProfile+273.15)).^coeff(3); % ea is used in kPa
layeredLongWave = skyEmissivity.*Constants.Boltzman.*(Canopy.tAirProfile+273.15).^4; % Long wave radiation [W m-2]
Soil.longEmitted = emissivitySoil*Constants.Boltzman*(Soil.temperature+273.15).^4;

Canopy.longEmitted = leafAbsorptivity*Constants.Boltzman*(Canopy.sunTleaf+273.15).^4.*Canopy.sunFraction ...
    + leafAbsorptivity*Constants.Boltzman*(Canopy.shadeTleaf+273.15).^4.*Canopy.shadeFraction;
Canopy.longEmitted(1) = soilAbsorptivity*Constants.Boltzman*(Soil.temperature+273.15)^4;

incidentTop = [0,layeredLongWave(3:end-1),Weather.long]; % indicent on top of canopy to incident on the last canopy layer (ie., i = 2)
incidentBottom = [0,Soil.longEmitted,layeredLongWave(2:end-2)]; % incident on bottom of the layer so not accounting the top layer and the additional 0 at the start
Canopy.longAbsorbed = leafAbsorptivity.*Canopy.deltaLAI.*(incidentTop+incidentBottom);
Soil.longAbsorbed = soilAbsorptivity*layeredLongWave(2);

Canopy.sunLongAbsorbed = Canopy.longAbsorbed.*Canopy.sunFraction;
Canopy.shadeLongAbsorbed = Canopy.longAbsorbed.*Canopy.shadeFraction;
Canopy.longUp = Canopy.longAbsorbed*NaN;


end