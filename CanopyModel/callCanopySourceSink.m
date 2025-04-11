%% Coumpute canopy profile of souces and sinks
function [Canopy] = callCanopySourceSink(LeafMassFlux,LeafEnergyFlux,Canopy,Soil)

% Photosynthesis flux [J m-2 ground area s-1]
nLayers = length(Canopy.z);
Canopy.aNet = ...
    ([LeafMassFlux(1:nLayers).aNet].*Canopy.sunFraction.*Canopy.deltaLAI+ ...
    [LeafMassFlux(nLayers+1:end).aNet].*Canopy.shadeFraction.*Canopy.deltaLAI);

% Latent heat flux [J m-2 ground area s-1]
Canopy.LE = ...
    ([LeafEnergyFlux(1:nLayers).LE].*Canopy.sunFraction.*Canopy.deltaLAI+ ...
    [LeafEnergyFlux(nLayers+1:end).LE].*Canopy.shadeFraction.*Canopy.deltaLAI);

% Sensible heat flux [J m-2 ground area s-1]
Canopy.H = ...
    ([LeafEnergyFlux(1:nLayers).H].*Canopy.sunFraction.*Canopy.deltaLAI+ ...
    [LeafEnergyFlux(nLayers+1:end).H].*Canopy.shadeFraction.*Canopy.deltaLAI);

%% Adding Soil flux
Canopy.aNet(1) = -Soil.CO2Flux;
Canopy.LE(1) = Soil.LE;
Canopy.H(1) = Soil.H;

end

