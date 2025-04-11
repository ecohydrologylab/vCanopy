%% Function to compute zenith angle
function [Weather] = callZenith(Weather)

standardMeridians = -345:15:345;
[minValue,closestIndex] = min(abs(Weather.longitude-standardMeridians));
closestValue = standardMeridians(closestIndex);
Weather.longitudeStandard = closestValue;
Weather.longitudeStandard = -82.5;
longitudeCorrection = (Weather.longitude-Weather.longitudeStandard)/15; % Longitude correction [h]
phi = 279.575+0.9856*Weather.julian; % [degrees]
timeEquation = (-104.7*sind(phi)+596.2*sind(2*phi)+4.3*sind(3*phi)- ...
    -12.7*sind(4*phi)-429.3*cosd(phi)-2.0*cosd(2*phi)+19.3*cosd(3*phi))/3600; % Equation of time [h]

hNoon = 12-longitudeCorrection-timeEquation; % Solar noon hour [h]
% hNoon = 13;
eta = 15*(Weather.hour-hNoon); % Hour angle [degrees]
% delta = asind(0.39785*sind(278.97+0.9856*Weather.julian+1.9165* ...
%     sind(356.6+0.9856*Weather.julian))); % Solar declination angle
%     [degrees] Campbell and Norman
% delta = -asind(0.39779*cos((0.98565*((Weather.julian-1.0)+10.0) + ...
%     1.914* sin(0.98565*((Weather.julian-1.0)-2.0)*pi/180.0))*pi/180.0)); % Solar declination angle [degrees]

delta = asind(-sind(23.45).*cosd(360.*(Weather.julian+10)./365)); %Spitters

Weather.zenith = acosd(sind(Weather.latitude).*sind(delta)+cosd(Weather.latitude).* ...
    cosd(delta)*cosd(eta)); % Solar zenith angle [degrees]

% if Weather.zenith >= 89
%     Weather.zenith = 89;
% end

end