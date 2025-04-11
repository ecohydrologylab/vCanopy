%% Boundary layer Conductance Module
temperatureKelvin = Weather.tAir + 273.15;       % Air temperature [Kelvin]
convert = Weather.pressure./(8.309*temperatureKelvin); % Convertion of boundary layer conductance from [m s-1] to [moles m-2 s-1] Nikolov et al 1995 (Page 212)
leafTemperatureKelvin = LeafState.tLeaf + 273.15; % Leaf temperature [Kelvin]

% Forced convection conductance
leafDimensions = LeafBoundaryLayer.leafDimension; % Leaf width/needle diameter [m]
LeafState.gbForced = LeafBoundaryLayer.cForced.*temperatureKelvin.^0.56.*((temperatureKelvin + 120.0)...
    .* (Weather.wind ./ leafDimensions ./ Weather.pressure)).^0.5; % Boundary layer forced conductance [m s-1]

% Free convection conductance
TDifference = (leafTemperatureKelvin ./ (1.0 - 0.378 * LeafState.eb ./ Weather.pressure)) -...
    (temperatureKelvin ./ (1.0 - 0.378 * Weather.ea ./ Weather.pressure)); % Virtual temperature difference [Kelvin] Nikolov et al 1995 (Page 211)
LeafState.gbFree = LeafBoundaryLayer.cFree.*leafTemperatureKelvin.^0.56.*((leafTemperatureKelvin + 120.0)...
    ./ Weather.pressure).^0.5.*(abs(TDifference) ./ leafDimensions).^0.25; % Boundary layer free conductance [m s-1]

% Maximum of two conductances
LeafState.gb = max(LeafState.gbFree, LeafState.gbForced); % Boundary layer conductance for vapour [m s-1]

% Compute leaf boundary layer conductance for vapour
LeafState.gb = LeafState.gb .* convert; % Boundary layer conductance for vapour [mol m-2 s-1]
