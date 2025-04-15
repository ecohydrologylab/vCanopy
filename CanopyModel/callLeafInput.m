%% Function to assemble leaf input from canopy model, convert from per ground to leaf area
function [LeafWeather,CanopyLayer] = callLeafInput(Weather,Canopy,savedtLeaf)

LeafWeather = table();
LeafWeather.PAR = [Canopy.sunPARAbsorbed./Canopy.deltaLAI./Canopy.sunFraction,...
    Canopy.shadePARAbsorbed./Canopy.deltaLAI./Canopy.shadeFraction]';
LeafWeather.NIR = [Canopy.sunNIRAbsorbed./Canopy.deltaLAI./Canopy.sunFraction,...
    Canopy.shadeNIRAbsorbed./Canopy.deltaLAI./Canopy.shadeFraction]';
LeafWeather.long = [Canopy.longAbsorbed./Canopy.deltaLAI,...
    Canopy.longAbsorbed./Canopy.deltaLAI]';
LeafWeather.long(isnan(LeafWeather.long)) = 0;
LeafWeather.fraction = [Canopy.sunFraction.*logical(Canopy.deltaLAI),...
    Canopy.shadeFraction.*logical(Canopy.deltaLAI)]';


CanopyLayer = table();
CanopyLayer.tauD = [Canopy.tauD,Canopy.tauD]';
CanopyLayer.deltaLAI = [Canopy.deltaLAI,Canopy.deltaLAI]';
CanopyLayer.savedtLeaf = savedtLeaf';
CanopyLayer.tLeafLag = [Canopy.sunTleafLag.*(Canopy.sunFraction>0),...
    Canopy.shadeTleafLag.*(Canopy.shadeFraction>0)]';
CanopyLayer.MPLA(:) = 1./Canopy.SLA;
CanopyLayer.tLeafLag(isnan(CanopyLayer.tLeafLag) | CanopyLayer.tLeafLag <= 0) = 0;
CanopyLayer.MPLA(CanopyLayer.tLeafLag <= 0) = 0;

LeafWeather.wind = [Canopy.windProfile(1:end-1),Canopy.windProfile(1:end-1)]';
LeafWeather.tAir = [Canopy.tAirProfile(1:end-1),Canopy.tAirProfile(1:end-1)]';
LeafWeather.ca = [Canopy.caProfile(1:end-1),Canopy.caProfile(1:end-1)]';
LeafWeather.O2 = LeafWeather.ca*0 + Weather.O2;
LeafWeather.pressure = LeafWeather.ca*0 + Weather.pressure;
LeafWeather.ea = [Canopy.eaProfile(1:end-1),Canopy.eaProfile(1:end-1)]';
LeafWeather.deltaT(:) = Weather.deltaT; 

CanopyLayer = table2struct(CanopyLayer);
LeafWeather = table2struct(LeafWeather);

end