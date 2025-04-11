function [Constants,EnergyOptions,LeafBoundaryLayer,Weather,Soil,Canopy,Photosynthesis,Stomata] = callInputCanopy(Options,WeatherData)

addpath('./CanopyModel')

% This code checks the model setup
callModelSetupChecks

[Weather] = callWeatherInput(Options,WeatherData);
[Soil,Canopy] = callSoilCanopyAirInput(Weather,Options); % Only model inputs no outputs
[Constants,EnergyOptions,LeafBoundaryLayer] = callConstantsInput(Options);
[Photosynthesis,Stomata] = callPhotosynthesisStomataInput(Options,Canopy);

% Correction to wind
Weather.wind(Weather.wind<0.1) = 0.1;
d0 = 2/3*Canopy.hN;
z0 = 0.13*Canopy.hN;
vonK = 0.41;
if ismember('uStar', Weather.Properties.VariableNames)
    uStar = Weather.uStar;
    uStar(isnan(uStar)) = 0.41.*Weather.wind(isnan(uStar))./(log((10-d0(isnan(uStar)))./z0(isnan(uStar))));
else
    uStar = 0.41.*Weather.wind./(log((10-d0)./z0));
end
uCorrection = (uStar./vonK).*log(10./Canopy.hN);
Weather.wind = Weather.wind - uCorrection;
% correctin for wind
Weather.wind(Weather.wind < 0.4) = NaN;
ind = find(isnan(Weather.wind));
hour = unique(Weather.hour);
U_mean = nanmean(reshape(Weather.wind,...
    [length(hour),length(Weather.wind)/length(hour)]),2);
for i = 1:1:length(ind)
    wHour = Weather.hour(ind(i));
    hIndex = find(wHour == hour);
    Weather.wind(ind(i)) = U_mean(hIndex);
end

% Applying photosynthesis decay
callPhotosynthesisDecay

% Assign specific leaf area
callBiomassCalculation

Weather = table2struct(Weather);
Soil = table2struct(Soil);
Canopy = table2struct(Canopy);
end