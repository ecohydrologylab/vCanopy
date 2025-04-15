%% Function to create canopy input

function [Weather] = callWeatherInput(Options,WeatherData)

Weather = WeatherData; % Loading the WeatherData
deltaT = diff(unique(Weather.hour)); deltaT = [deltaT(1);deltaT];
Weather.deltaT = repmat(deltaT,[length(Weather.hour)/length(deltaT),1]); % Time interval between two simulation points [h]
Weather.PAR = WeatherData.Rg*0.45; % 45% of incoming SW radiation is PAR [W m-2]
Weather.NIR = WeatherData.Rg*0.55; % 55% of incoming SW radiation is PAR [W m-2]

%% Changes based on future weather
if Options.deltaTair ~= 0 && Options.deltaRH == 0
% Changes the Tair so that RH doesn't change
% 1) Calculate the current RH
% 2) Increase the Tair
% 3) Calculate the es again to maintain same RH
    RH = WeatherData.ea./(611.0*exp(17.502*WeatherData.tAir./ ...
        (240.97+WeatherData.tAir)))*100 + Options.deltaRH; % Atmospheric RH [%]    
    Weather.tAir = WeatherData.tAir + Options.deltaTair; % Atmospheric Tair [Celsius]
    Weather.ea = RH./100.*(611.0*exp(17.502*Weather.tAir./ ...
    (240.97+Weather.tAir))); % Atmospheric water vapour pressure [Pa]
elseif Options.deltaTair == 0 && Options.deltaRH ~= 0
    % Changes the Ea based on desired RH change
    RH = WeatherData.ea./(611.0*exp(17.502*WeatherData.tAir./ ...
        (240.97+WeatherData.tAir)))*100 + Options.deltaRH; % Atmospheric RH [%]
    Weather.ea = RH./100.*(611.0*exp(17.502*WeatherData.tAir./ ...
        (240.97+WeatherData.tAir))); % Atmospheric water vapour pressure [Pa]
elseif Options.deltaTair ~= 0 && Options.deltaRH ~= 0
    % Changes the Ea and Tair corresponding to each other    
    RH = WeatherData.ea./(611.0*exp(17.502*WeatherData.tAir./ ...
    (240.97+WeatherData.tAir)))*100 + Options.deltaRH; % Atmospheric RH [%]
    Weather.tAir = WeatherData.tAir + Options.deltaTair; % Atmospheric Tair [Celsius]
    Weather.ea = RH./100.*(611.0*exp(17.502*Weather.tAir./ ...
        (240.97+Weather.tAir))); % Atmospheric water vapour pressure [Pa]
end

Weather.ca = WeatherData.ca + Options.deltaCO2; % Atmospheric [CO2]

% Import the GDD if avaiable in the input
if ismember('GDD', WeatherData.Properties.VariableNames)
    Weather.GDD = WeatherData.GDD;
end

% Calculate the Solar Zenith when not available in the input
index = isnan(Weather.zenith);
h = 2*pi*(Weather.hour-12)./24;
tilda = asin(-sin(23.45.*pi/180).*cos(2*pi*(Weather.julian + 10)/365));
zenith = acosd(sind(Weather.latitude).*sin(tilda) + cosd(Weather.latitude).*cos(tilda).*cos(h));
Weather.zenith(index) = zenith(index);

% Initialise the direct and diffused fraction for PAR and NIR
Weather.diffuseFraction = nan(height(Weather),1); % Diffuse fraction
Weather.PARDirect = nan(height(Weather),1); % Direct PAR
Weather.PARDiffuse = nan(height(Weather),1); % Diffuse PAR
Weather.NIRDirect = nan(height(Weather),1); % Direct NIR
Weather.NIRDiffuse = nan(height(Weather),1); % Diffuse NIR

Weather = callLongWaveIncoming(Weather); % Compute longwave fraction of radiation
Weather = callDiffuseFraction(Weather,Options); % Compute diffuse fraction of PAR

if sum(sum(isnan(table2array(Weather))))
    warning("Check the Weather data some data is missing or is 'NaN'");
end


end