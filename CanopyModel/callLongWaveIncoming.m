function [Weather] = callLongWaveIncoming(Weather)
% Function to compute incoming longwave radiationin
nanIndex = isnan(Weather.long);
skyEmissivity = 1.72.*(Weather.ea/1000./(Weather.tAir+273.15)).^(1/7);
boltzmann = 5.6697E-8; % Stefan-Boltzmann constant W m-2 K-4
longWave = skyEmissivity.*boltzmann.*(Weather.tAir+273.15).^4; % Long wave radiation [W m-2]
Weather.long(nanIndex) = longWave(nanIndex);
end