function [LeafState] = callBoundaryLayerConductance(LeafBoundaryLayer,Stomata,Weather,LeafState,LeafMassFlux)

leafDimensions = LeafBoundaryLayer.leafDimension; % Leaf width/needle diameter [m]
tAirKelvin = Weather.tAir + 273.15; % Air temperature [Kelvin]
convert = Weather.pressure/(8.309*tAirKelvin); % Convertion of boundary layer conductance from [m s-1] to [moles m-2 s-1] Nikolov et al 1995 (Page 212)
tLeafKelvin = LeafState.tLeaf + 273.15; % Leaf temperature [Kelvin]
callEs = @(temperature) (0.611*exp(17.502*temperature/(240.97+temperature)))*1000; %% Clausius–Clapeyron relation Saturation Vapour Pressure [Pa]
LeafState.ei = callEs(LeafState.tLeaf);
% callEs = @(temperature) (0.611*exp(17.502*temperature/(240.97+temperature)))*1000; % Saturation Vapour Pressure [Pa]
% callEa = @(temperature) (Weather.RH/100*callEs(temperature)); % Atmospheric Vapour Pressure [Pa]
% Weather.ea = callEa(Weather.tAir); % Atmospheric vapour pressure [Pa]

%% Forced convection conductance
LeafState.gbForced = LeafBoundaryLayer.cForced*tAirKelvin^0.56*((tAirKelvin + 120.0)...
    * (Weather.wind / leafDimensions / Weather.pressure))^0.5; % Nikolov et al 1995 Boundary layer forced conductance [m s-1] 

%% Free convection conductance
TDifference = (tLeafKelvin / (1.0 - 0.378 * LeafState.eb / Weather.pressure)) -...
    (tAirKelvin / (1.0 - 0.378 * Weather.ea / Weather.pressure)); % Virtual temperature difference [Kelvin] Nikolov et al 1995 (Page 211)
LeafState.gbFree = LeafBoundaryLayer.cFree*tLeafKelvin^0.56*((tLeafKelvin + 120.0)...
    / Weather.pressure)^0.5*(abs(TDifference) / leafDimensions)^0.25; % Nikolov et al 1995 Boundary layer free conductance [m s-1]

%% Maximum of two conductances
LeafState.gb = max(LeafState.gbFree, LeafState.gbForced); % Boundary layer conductance for vapour [m s-1]

%% Compute leaf boundary layer conductance for vapour
LeafState.gb = LeafState.gb * convert; % Boundary layer conductance for vapour [mol m-2 s-1]
LeafState.g = LeafState.gs*Stomata.Mb*LeafState.gb/(LeafState.gs + Stomata.Mb*LeafState.gb); % Total leaf vapour conductance [mol m-2 sec-1]
LeafState.eb = (LeafState.gs / convert * LeafState.ei + LeafState.gb  / convert * Weather.ea)...
    / (LeafState.gs / convert + LeafState.gb / convert); % Vapour Pressure at leaf for boundary layer conductance [Pa]

LeafState.cb = Weather.ca - 1.37 * LeafMassFlux.aNet / (Stomata.Mb*LeafState.gb); % Leaf boundary layer concentration of CO2 [ppm]

end