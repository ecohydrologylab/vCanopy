%% Function to attenuate scattered radiation

function [Canopy] = callScatteredRadiation(EnergyOptions,Canopy,Options)

incident = sum(Canopy.scatteredDown+Canopy.scatteredUp);

Error.count = 0;
Error.countMax = 30;
Error.threshold = 0.1;
Error.absolute = 1;

while(Error.absolute >= Error.threshold && Error.count <= Error.countMax)
    
    % Down attenuation
    Canopy.down = Canopy.scatteredDown;
    [Canopy] = callDownRadiation(EnergyOptions,Canopy,Options,length(Canopy.deltaLAI));
    Canopy.scatteredAbsorbed = Canopy.scatteredAbsorbed+Canopy.absorbed;
    Canopy.scatteredDown = [Canopy.transmitted(2:end) 0]; % Shift as transmitted is bottom of layer
    Canopy.scatteredUp = Canopy.scatteredUp+Canopy.reflected;
    Error.absolute = abs(incident-sum(Canopy.scatteredAbsorbed)- ...
        sum(Canopy.scatteredDown+Canopy.scatteredUp));
    
    % Up attenuation
    Canopy.up = Canopy.scatteredUp;
    [Canopy] = callUpRadiation(EnergyOptions,Canopy,Options,length(Canopy.deltaLAI));
    Canopy.scatteredAbsorbed = Canopy.scatteredAbsorbed+Canopy.absorbed;
    Canopy.scatteredUp = Canopy.transmitted;
    Canopy.scatteredUp(end) = Canopy.scatteredUp(end)+Canopy.up(end); % Add unintercepted at top
    Canopy.scatteredDown = Canopy.scatteredDown+[Canopy.reflected(2:end),0]; % Shift as reflected is top of layer
    Error.absolute = abs(incident-sum(Canopy.scatteredAbsorbed)- ...
        sum(Canopy.scatteredDown+Canopy.scatteredUp));

    % Error analysis
    Error.absolute = abs(incident-sum(Canopy.scatteredAbsorbed)- ...
        Canopy.scatteredUp(end));
    Error.count = Error.count+1;
    
end

if (Error.count >= Error.countMax)
    disp(strcat('Scattered radiation not converged at canopy ID = ',num2str(Canopy.ID), ...
        'for',Options.wavelength,'and for',Options.orientation))
end

end