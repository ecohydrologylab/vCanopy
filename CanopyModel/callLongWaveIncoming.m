%% Function to compute incoming longwave radiationin sunlight

function [Weather] = callLongWaveIncoming(Weather)

nanIndex = isnan(Weather.long);
% skyEmissivity = 1.0105*((Weather.ea*10/1000)./(Weather.tAir+273)).^(1/151.2986); % ea is used in kPa
% coeff = [0.9258    1.0059    1.0677];
% skyEmissivity = coeff(1) +  coeff(2)*(Weather.ea/1000./(Weather.tAir+273.15)).^coeff(3); % ea is used in kPa
skyEmissivity = 1.72.*(Weather.ea/1000./(Weather.tAir+273.15)).^(1/7);
boltzmann = 5.6697E-8; % Stefan-Boltzmann constant W m-2 K-4
longWave = skyEmissivity.*boltzmann.*(Weather.tAir+273.15).^4; % Long wave radiation [W m-2]
Weather.long(nanIndex) = longWave(nanIndex);
end