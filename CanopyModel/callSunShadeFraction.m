%% Compute sun and shade fractions
function [Canopy] = callSunShadeFraction(Weather,Canopy)

%% Compute sun fraction
Canopy.sunFraction = zeros(1,length(Canopy.z)); % Initialize sun fraction

if Weather.zenith <= 90
    kBeam = sqrt(Canopy.x^2 + tand(Weather.zenith).^2)/ ...
        (Canopy.x+1.774*(Canopy.x+1.182).^(-0.733)); % Fraction of leaf area projected on the horizontal
    
% if length(Canopy.deltaLAI) == 2     
%     % Top layer
%     LAIIndex = length(Canopy.deltaLAI);
% %     Canopy.sunFraction(1,LAIIndex) = (1-exp(-kBeam*Canopy.omega* ...
% %         sum(Canopy.deltaLAI(1,LAIIndex:end))))/ ...
% %         (kBeam*Canopy.omega*Canopy.deltaLAI(1,LAIIndex)); % Sun fraction
%     Canopy.sunFraction(1,LAIIndex) = exp(-kBeam*Canopy.omega...
%         *Canopy.deltaLAI(1,LAIIndex)); % Sun fraction Drewry  
% else
    % Top layer
    LAIIndex = length(Canopy.deltaLAI);
    Canopy.sunFraction(1,LAIIndex) = (1-exp(-kBeam*Canopy.omega* ...
        sum(Canopy.deltaLAI(1,LAIIndex:end))))/ ...
        (kBeam*Canopy.omega*Canopy.deltaLAI(1,LAIIndex)); % Sun fraction
%     Canopy.sunFraction(1,LAIIndex) = exp(-kBeam*Canopy.omega* ...
%         sum(Canopy.deltaLAI(1,LAIIndex:end))); % Sun fraction Drewry       
% end

    % All other layers
    for LAIIndex = length(Canopy.deltaLAI)-1:-1:1
        Canopy.sunFraction(1,LAIIndex) = (...
            exp(-kBeam*Canopy.omega*sum(Canopy.deltaLAI(1,LAIIndex+1:end)))- ...
            exp(-kBeam*Canopy.omega*sum(Canopy.deltaLAI(1,LAIIndex:end))))/ ...
            (kBeam*Canopy.omega*Canopy.deltaLAI(1,LAIIndex)); % Sun fraction
%         Canopy.sunFraction(1,LAIIndex) = exp(-kBeam*Canopy.omega* ...
%             sum(Canopy.deltaLAI(1,LAIIndex:end))); % Sun fraction Drewry
    end 
end

Canopy.sunFraction(Canopy.sunFraction<0.01) = 0.0;

%% Compute shade fraction
Canopy.shadeFraction = 1-Canopy.sunFraction; % Shade fraction

%% Remove values for layers without leaves
index = find(Canopy.deltaLAI(1,:) == 0);
Canopy.sunFraction(1,index) = 0;
Canopy.shadeFraction(1,index) = 1;

%% Overwrite for big Leaf model: Make all leaves sunlit
if length(Canopy.z) == 2
    Canopy.sunFraction = [0,1];
    Canopy.shadeFraction = [1,0];
end

end