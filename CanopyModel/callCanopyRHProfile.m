%% Function to compute canopy RH profile
function [Canopy] = callCanopyRHProfile(Canopy)
% Saturation Vapour Pressure [Pa]
% totalPressure = (Constants.Pressure + Canopy.eaProfile(1:end))/1000;
callEs = @(temperature) (0.611*exp(17.502*temperature./(240.97+temperature)))*1000;
Canopy.RHProfile = Canopy.eaProfile./callEs(Canopy.tAirProfile)*100;

% index = Canopy.RHProfile >= 100;
% Canopy.RHProfile(index) = 99.9;
% ESProfile = 611.*exp(17.502.*Canopy.tAirProfile./(240.97+Canopy.tAirProfile));
% Canopy.eaProfile = Canopy.RHProfile/100.*ESProfile;

end