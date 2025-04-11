%% Function to compute the fraction of diffuse radiation in sunlight

function [Weather] = callDiffuseFraction(Weather,Options)
solarConstant = 1370.0; % [W m-2]
elevation = 90 - Weather.zenith; % Solar elevation [degrees]

Weather.diffuseFraction = ones(length(elevation),1);

zenithIndex = Weather.zenith > 85;%(elevation < 1 & elevation > 90);
So = solarConstant*(1.0+0.033*cosd(360*Weather.julian./365)).*sind(elevation); % [W m-2]
rho = 0.847-1.61*sind(elevation)+1.04*(sind(elevation)).^2.0;
K = (1.47-rho)./1.66;

SW_So = (Weather.PAR + Weather.NIR)./ So;

index1 = SW_So <= 0.22;
Weather.diffuseFraction(index1) = 1;

index2 = SW_So > 0.22 & SW_So <= 0.35;
Weather.diffuseFraction(index2) = 1 - 6.4*(SW_So(index2) - 0.22).^2;

index3 = SW_So > 0.35 & SW_So <= K;
Weather.diffuseFraction(index3) = 1.47 - 1.66*SW_So(index3);

index4 = SW_So > K;
Weather.diffuseFraction(index4) = rho(index4);

Weather.diffuseFraction(zenithIndex) = 1;
Weather.diffuseFraction(Weather.diffuseFraction<0) = 1;

% timePoints = length(unique(Weather.hour));
% zenith_Array = reshape(Weather.zenith,[timePoints,length(Weather.zenith)/timePoints]);
% diffFrac_Array = reshape(Weather.diffuseFraction,[timePoints,length(Weather.diffuseFraction)/timePoints]);

% for dLoop = 1:1:size(diffFrac_Array,2)
%     % Early noon
%     zenith80Index = find(zenith_Array(1:round(timePoints/2),dLoop) < 80);
%     zenith90Index = find(zenith_Array(1:round(timePoints/2),dLoop) >= 90);
%     smootheningIndex = zenith90Index(end):1:zenith80Index(1);
%     if length(smootheningIndex) > 2
%         diffFrac_Array(smootheningIndex(2:end-1),dLoop) = 0.8;
% %         interp1(zenith_Array(smootheningIndex([1,end]),dLoop),...
% %             diffFrac_Array(smootheningIndex([1,end]),dLoop),zenith_Array(smootheningIndex([2:end-1]),dLoop),...
% %             'spline');
%     end
%     % After noon
%     zenith80Index = find(zenith_Array(round(timePoints/2)+1:timePoints,dLoop) < 80);
%     zenith90Index = find(zenith_Array(round(timePoints/2)+1:timePoints,dLoop) >= 90);
%     smootheningIndex = round(timePoints/2)+[zenith80Index(end):1:zenith90Index(1)];
%     if length(smootheningIndex) > 2
%         diffFrac_Array(smootheningIndex(2:end-1),dLoop) = 0.8;
% %         interp1(zenith_Array(smootheningIndex([1,end]),dLoop),...
% %             diffFrac_Array(smootheningIndex([1,end]),dLoop),zenith_Array(smootheningIndex([2:end-1]),dLoop),...
% %             'spline');
%     end
% end

% x = reshape(Weather.diffuseFraction,[48,length(Weather.diffuseFraction)/48]);
% D90_Zenith = [0;diff(zenith90Index)];
% zz90 = reshape(D90_Zenith,[48,length(D90_Zenith)/48]);
% D_Zenith = [0;diff(zenith80Index)];
% zz = reshape(D_Zenith,[48,length(D_Zenith)/48]);
% earlyMorningIndex = find(D_Zenith == -1);
% earlyMorning_DFvalues = 
% earlyMorningArray = [earlyMorningIndex-3,earlyMorningIndex-2,earlyMorningIndex-1];


% Weather.diffuseFraction = diffFrac_Array(:);

Weather.PARDirect = (1-Options.diffuseFactor/100).*(1-Weather.diffuseFraction).*Weather.PAR;
Weather.PARDiffuse = Weather.PAR - Weather.PARDirect;
Weather.NIRDirect = (1-Options.diffuseFactor/100).*(1-Weather.diffuseFraction).*Weather.NIR;
Weather.NIRDiffuse = Weather.NIR - Weather.NIRDirect;

end