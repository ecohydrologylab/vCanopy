function [Soil,Canopy] = callLongRadiationAttenuation(Constants,EnergyOptions,Weather,Soil,Canopy)

%% Input

emissivityLeaf = EnergyOptions.epsilonLeaf;
emissivitySoil = EnergyOptions.epsilonSoil;

%% Compute longwave radiation attenuation
% Canopy.tauD = exp(-Canopy.kD.*Canopy.omega.*Canopy.deltaLAI);
Canopy.longEmitted = (1-Canopy.tauD).*(emissivityLeaf*Constants.Boltzman*...
    (Canopy.sunTleaf+273.15).^4.*Canopy.sunFraction ...
    + emissivityLeaf*Constants.Boltzman*(Canopy.shadeTleaf+273.15).^4.*Canopy.shadeFraction);
Canopy.longEmitted(1) = emissivitySoil*Constants.Boltzman*(Soil.temperature+273.15)^4;

Canopy.scatteredDown = [Canopy.longEmitted(2:end) Weather.long]; % Shift as emitted is bottom of layer
Canopy.scatteredUp = Canopy.longEmitted;
Canopy.scatteredAbsorbed = zeros(1,size(Canopy.z,2));
Options.wavelength = 'long';
Options.orientation = 'diffuse';
Options.isotropy = EnergyOptions.isotropy;
Options.kDVariation = EnergyOptions.kDVariation;
Canopy.zenith = Weather.zenith;

[Canopy] = callScatteredRadiation(EnergyOptions,Canopy,Options);

Canopy.longAbsorbed = Canopy.scatteredAbsorbed;
Canopy.longUp = Canopy.scatteredUp;

%% Error

Error.longAbsolute = Weather.long+ sum(Canopy.longEmitted(2:end))+sum(Canopy.longEmitted)- ...
    sum(Canopy.longAbsorbed)-Canopy.longUp(end); % Longwave radiation

%% Compute air and soil radiation
Soil.longAbsorbed = Canopy.longAbsorbed(1);
Soil.longEmitted = Canopy.longEmitted(1);

%% Long absorbed sun-shade model

Canopy.sunLongAbsorbed = Canopy.sunFraction.*Canopy.longAbsorbed;
Canopy.shadeLongAbsorbed = Canopy.shadeFraction.*Canopy.longAbsorbed;

% Canopy.sunLongEmitted = Canopy.sunFraction.*Canopy.longEmitted;
% Canopy.shadeLongEmitted = Canopy.shadeFraction.*Canopy.longEmitted;
%% Remove duplicate fields in structure

Canopy = rmfield(Canopy,{'down','up','scatteredUp','scatteredDown', ...
    'scatteredAbsorbed','absorbed','transmitted','reflected', ...
    'zenith'});

end