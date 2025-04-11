%% Function to create canopy input

function [Weather] = callWeatherInputDrewrySeason(Options,WeatherData)

Weather = table();

%% Carefullt change these things based on future weather
Weather.temperature = WeatherData.Ta + Options.deltaTair;
Weather.RH = WeatherData.Ea*1000./(611.0*exp(17.502*WeatherData.Ta./ ...
    (240.97+WeatherData.Ta)))*100 + Options.deltaRH;
Weather.ea = Weather.RH./100.*(611.0*exp(17.502*Weather.temperature./ ...
    (240.97+Weather.temperature)));
Weather.ca = Weather.temperature*0+420 + Options.deltaCO2;
Weather.wind = (1+Options.deltaWind/100).*WeatherData.U; 
WeatherData.Rg = (1+Options.deltaRg/100).*WeatherData.Rg;
%%

Weather.julian = WeatherData.DOY;
Weather.year = WeatherData.Year;
Weather.hour = WeatherData.Hour;
Weather.latitude = Weather.hour*0+WeatherData.latitude;
Weather.longitude = Weather.hour*0+WeatherData.longitude;
Weather.longitudeStandard = Weather.hour*0-82.5;
deltaT = diff(Weather.hour); deltaT = [deltaT(1);deltaT];
Weather.deltaT = deltaT;


if ismember('GDD', WeatherData.Properties.VariableNames)
    Weather.GDD = WeatherData.GDD;
end

if ~ismember('PAR', WeatherData.Properties.VariableNames)
    Weather.PAR = WeatherData.Rg*0.5;
else
    WeatherData.PAR = (1+Options.deltaRg/100).*WeatherData.PAR;
    Weather.PAR = WeatherData.PAR;
end
if ~ismember('NIR', WeatherData.Properties.VariableNames)
    Weather.NIR = WeatherData.Rg*0.5;
else
    WeatherData.NIR= (1+Options.deltaRg/100).*WeatherData.NIR;
    Weather.NIR = WeatherData.NIR;
end

Weather.long = WeatherData.LW_in;
Weather.O2 = Weather.hour*0+210;
Weather.zenith = WeatherData.Zenith;

index = isnan(Weather.zenith);
h = 2*pi*(Weather.hour-12)./24;
tilda = asin(-sin(23.45.*pi/180).*cos(2*pi*(Weather.julian + 10)/365));
zenith = acosd(sind(Weather.latitude).*sin(tilda) + cosd(Weather.latitude).*cos(tilda).*cos(h));
Weather.zenith(index) = zenith(index);


Weather.diffuseFraction = nan(height(Weather),1); % Diffuse fraction
Weather.PARDirect = nan(height(Weather),1); % Direct PAR
Weather.PARDiffuse = nan(height(Weather),1); % Diffuse PAR
Weather.NIRDirect = nan(height(Weather),1); % Direct NIR
Weather.NIRDiffuse = nan(height(Weather),1); % Diffuse NIR

% Weather(timeLoop,:) = callZenith(Weather(timeLoop,:)); % Compute zenith angle
Weather = callLongWaveIncoming(Weather); % Compute longwave fraction of radiation


if ~ismember('PARdirect', WeatherData.Properties.VariableNames) && ... 
    ~ismember('PARdiffuse', WeatherData.Properties.VariableNames) && ...
    ~ismember('NIRdirect', WeatherData.Properties.VariableNames) && ...
    ~ismember('NIRdiffuse', WeatherData.Properties.VariableNames)
    Weather = callDiffuseFraction(Weather,Options); % Compute diffuse fraction of PAR
else
    WeatherData.PARdirect = (1+Options.deltaRg/100).*WeatherData.PARdirect;
    WeatherData.NIRdirect= (1+Options.deltaRg/100).*WeatherData.NIRdirect;
    WeatherData.PARdiffuse = (1+Options.deltaRg/100).*WeatherData.PARdiffuse;
    WeatherData.NIRdiffuse = (1+Options.deltaRg/100).*WeatherData.NIRdiffuse;

    Weather.PARDirect = WeatherData.PARdirect;
    Weather.PARDiffuse = WeatherData.PARdiffuse;
    Weather.NIRDirect = WeatherData.NIRdirect;
    Weather.NIRDiffuse = WeatherData.NIRdiffuse;
    Weather.diffuseFraction = Weather.PARDiffuse./(Weather.PAR);
    Weather.diffuseFraction(isnan(Weather.diffuseFraction)) = 1;
end

if sum(sum(isnan(table2array(Weather))))
    error("Check the Weather data some data is missing or is 'NaN'");
end

end