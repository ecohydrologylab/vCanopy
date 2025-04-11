function [Constants,Weather,Soil,Air,Canopy,Photosynthesis,Stomata] = callInputCanopyDrewrySeason(Options,WeatherData)

addpath('./CanopyModel')

% This code checks the model setup
callModelSetupChecks

[Weather] = callWeatherInput(Options,WeatherData);
[Soil,Canopy,Air] = callSoilCanopyAirInput(Weather,Options);
[Constants] = callConstantsInput;
[Photosynthesis,Stomata] = callPhotosynthesisStomataInput(Options,Canopy);

end