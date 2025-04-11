%% Function to create canopy input

function [Weather] = callWeatherInput(Options,WeatherData)

Weather = WeatherData; % Loading the WeatherData
deltaT = diff(unique(Weather.hour)); deltaT = [deltaT(1);deltaT];
Weather.deltaT = repmat(deltaT,[length(Weather.hour)/length(deltaT),1]);
Weather.PAR = WeatherData.Rg*0.45;
Weather.NIR = WeatherData.Rg*0.55;

%% Changes based on future weather
if Options.deltaTair ~= 0 && Options.deltaRH == 0
    % Changes the Tair so that RH doesn't change
    RH = WeatherData.ea./(611.0*exp(17.502*WeatherData.tAir./ ...
        (240.97+WeatherData.tAir)))*100 + Options.deltaRH;    
    Weather.tAir = WeatherData.tAir + Options.deltaTair; 
    Weather.ea = RH./100.*(611.0*exp(17.502*Weather.tAir./ ...
    (240.97+Weather.tAir)));
elseif Options.deltaTair == 0 && Options.deltaRH ~= 0
    % Changes the Ea based on desired RH change
    RH = WeatherData.ea./(611.0*exp(17.502*WeatherData.tAir./ ...
        (240.97+WeatherData.tAir)))*100 + Options.deltaRH;
    Weather.ea = RH./100.*(611.0*exp(17.502*WeatherData.tAir./ ...
        (240.97+WeatherData.tAir)));
elseif Options.deltaTair ~= 0 && Options.deltaRH ~= 0
    % Changes the Ea and Tair corresponding to each other    
    RH = WeatherData.ea./(611.0*exp(17.502*WeatherData.tAir./ ...
    (240.97+WeatherData.tAir)))*100 + Options.deltaRH;
    Weather.tAir = WeatherData.tAir + Options.deltaTair;
    Weather.ea = RH./100.*(611.0*exp(17.502*Weather.tAir./ ...
        (240.97+Weather.tAir)));
end

% RH = Weather.ea./(611.0*exp(17.502*Weather.tAir./ ...
%     (240.97+Weather.tAir)))*100; 

% if Options.deltaTair ~= 0 && Options.deltaRH == 0
%     % Changes the Tair only without affecting the Ea
%     Weather.tAir = WeatherData.tAir + Options.deltaTair; 
% elseif Options.deltaTair == 0 && Options.deltaRH ~= 0
%     % Changes the Ea based on desired RH change
%     RH = WeatherData.ea./(611.0*exp(17.502*WeatherData.tAir./ ...
%         (240.97+WeatherData.tAir)))*100 + Options.deltaRH;
%     Weather.ea = RH./100.*(611.0*exp(17.502*WeatherData.tAir./ ...
%         (240.97+WeatherData.tAir)));
% elseif Options.deltaTair ~= 0 && Options.deltaRH ~= 0
%     % Changes the Ea and Tair corresponding to each other
%     Weather.tAir = WeatherData.tAir + Options.deltaTair;
%     RH = WeatherData.ea./(611.0*exp(17.502*WeatherData.tAir./ ...
%         (240.97+WeatherData.tAir)))*100 + Options.deltaRH;
%     Weather.ea = RH./100.*(611.0*exp(17.502*WeatherData.tAir./ ...
%         (240.97+WeatherData.tAir)));
% end

Weather.ca = WeatherData.ca + Options.deltaCO2;

if ismember('GDD', WeatherData.Properties.VariableNames)
    Weather.GDD = WeatherData.GDD;
end

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


% if ~ismember('PARdirect', WeatherData.Properties.VariableNames) && ...
%     ~ismember('PARdiffuse', WeatherData.Properties.VariableNames) && ...
%     ~ismember('NIRdirect', WeatherData.Properties.VariableNames) && ...
%     ~ismember('NIRdiffuse', WeatherData.Properties.VariableNames)
Weather = callDiffuseFraction(Weather,Options); % Compute diffuse fraction of PAR
% else
%     WeatherData.PARdirect = (1+Options.deltaRg/100).*WeatherData.PARdirect;
%     WeatherData.NIRdirect= (1+Options.deltaRg/100).*WeatherData.NIRdirect;
%     WeatherData.PARdiffuse = (1+Options.deltaRg/100).*WeatherData.PARdiffuse;
%     WeatherData.NIRdiffuse = (1+Options.deltaRg/100).*WeatherData.NIRdiffuse;
%
%     Weather.PARDirect = WeatherData.PARdirect;
%     Weather.PARDiffuse = WeatherData.PARdiffuse;
%     Weather.NIRDirect = WeatherData.NIRdirect;
%     Weather.NIRDiffuse = WeatherData.NIRdiffuse;
%     Weather.diffuseFraction = Weather.PARDiffuse./(Weather.PAR);
%     Weather.diffuseFraction(isnan(Weather.diffuseFraction)) = 1;
% end

if sum(sum(isnan(table2array(Weather))))
    warning("Check the Weather data some data is missing or is 'NaN'");
end


end