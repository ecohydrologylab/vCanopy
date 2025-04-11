function [Photosynthesis,Stomata] = callPhotosynthesisStomataInput(Options,Canopy)

%% Photosynthesis
% Photosynthesis = table;
leafID = 0:1:size(Canopy.deltaLAI,2)-1;
Photosynthesis.leafID = [1+leafID/100,2+leafID/100]';
Photosynthesis.plant = Options.plant*ones(2*size(Canopy.deltaLAI,2),1);

Photosynthesis.vcmax25 = 50*ones(2*size(Canopy.deltaLAI,2),1); %
Photosynthesis.rd25 = 0.05*ones(2*size(Canopy.deltaLAI,2),1);
Photosynthesis.vpr = 80*ones(2*size(Canopy.deltaLAI,2),1);
Photosynthesis.jmax25 = Photosynthesis.vcmax25*Options.Jmax_Vcmax; %
Photosynthesis.vpmax25 = 110*ones(2*size(Canopy.deltaLAI,2),1);
Photosynthesis.gbs = 0.003*ones(2*size(Canopy.deltaLAI,2),1);
Photosynthesis.leafDimensions = 0.08*ones(2*size(Canopy.deltaLAI,2),1); % Leaf width/needle diameter [m]
%% Stomata values
Stomata.model = "BB";
if Stomata.model == "BB"
    %% Ball Berry Model
    Stomata.slope = 3.23; %7; %4.53; % 7;% 3;% 3.13;% 3;% 
    Stomata.intercept = 0.06; % 0.008; %0.017; %0.008;% 0.08;%0.078;% 0.08;% 
elseif Stomata.model == "BM"
    %% Belinda Model
    Stomata.slope = 0.44;%
    Stomata.intercept = 0.07;%
elseif Stomata.model == "BBL2"
    %% Ball-Berry-Leunning Model
    Stomata.slope = 3.39;% 7;%
    Stomata.intercept = 0.088;% 0.008;%  
    Stomata.D0 = 1.5;% 0.008;%
end
Stomata.fsv = 1;
Stomata.s = 0.71; % ratio of stomatal density on abaxial and adaxian side
Stomata.Mb = 0.5*(1 + Stomata.s)^2/(1 + Stomata.s^2);

Photosynthesis.theta = 0.76*ones(2*size(Canopy.deltaLAI,2),1);
Photosynthesis.alpha = 0*ones(2*size(Canopy.deltaLAI,2),1);
Photosynthesis.waterStress = ones(2*size(Canopy.deltaLAI,2),1);
Photosynthesis.x = 0.4*ones(2*size(Canopy.deltaLAI,2),1);
Photosynthesis.sco25 = 2590*ones(2*size(Canopy.deltaLAI,2),1);
Photosynthesis.kn = Options.kn; % Exponential nitrogen decay
Photosynthesis.PcSeasonalVariation = Options.PcSeasonalVariation; % Decay in photosynthetic parameters
Photosynthesis.linearDTD = Options.linearDTD;
Photosynthesis.Jmax_Vcmax = Options.Jmax_Vcmax;
Photosynthesis.nitrogenDecay = Options.nitrogenDecay;
end
