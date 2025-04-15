function [Photosynthesis,Stomata] = callPhotosynthesisStomataInput(Options,Canopy)

%% Photosynthesis
leafID = 0:1:size(Canopy.deltaLAI,2)-1;
Photosynthesis.leafID = [1+leafID/100,2+leafID/100]';
Photosynthesis.plant = Options.plant*ones(2*size(Canopy.deltaLAI,2),1);

Photosynthesis.vcmax25 = 50*ones(2*size(Canopy.deltaLAI,2),1); % Maximum rubisco carboxylation rate @25 ◦C
Photosynthesis.rd25 = 0.05*ones(2*size(Canopy.deltaLAI,2),1); % Ratio of respiratio rate to maximum rubisco carboxylation rate @25 ◦C
Photosynthesis.vpr = 80*ones(2*size(Canopy.deltaLAI,2),1); % maximum PEP regeneration rate in the mesophyll
Photosynthesis.jmax25 = Photosynthesis.vcmax25*Options.Jmax_Vcmax; % Maximum whole chain electron transport rate @25 ◦C
Photosynthesis.vpmax25 = 110*ones(2*size(Canopy.deltaLAI,2),1); % Maximum PEP carboxylation rate @25 ◦
Photosynthesis.gbs = 0.003*ones(2*size(Canopy.deltaLAI,2),1); % CO2 conductance of bundle sheath
Photosynthesis.leafDimensions = 0.08*ones(2*size(Canopy.deltaLAI,2),1); % Leaf width/needle diameter [m]

%% Stomata values
Stomata.model = "BB";
if Stomata.model == "BB"
    %% Ball Berry Model
    Stomata.slope = 3.23; % Seasonal average BB slope 
    Stomata.intercept = 0.06; % Seasonal average BB intercept
elseif Stomata.model == "BM"
    %% Belinda Model
    Stomata.slope = 0.44; % Fitted in Ball 1987 data (Shared by author)
    Stomata.intercept = 0.07; % Fitted in Ball 1987 data (Shared by author)
end
Stomata.fsv = 1;
Stomata.s = 0.71; % Ratio of stomatal density on adaxial and abaxial side
Stomata.Mb = 0.5*(1 + Stomata.s)^2/(1 + Stomata.s^2); % Correction factor for the boundary layer conductance

Photosynthesis.theta = 0.76*ones(2*size(Canopy.deltaLAI,2),1); % Curvature parameter
Photosynthesis.alpha = 0*ones(2*size(Canopy.deltaLAI,2),1); % Fraction of O2 evolving in bundle sheath
Photosynthesis.waterStress = ones(2*size(Canopy.deltaLAI,2),1); % Water stress:  Not limiting
Photosynthesis.x = 0.4*ones(2*size(Canopy.deltaLAI,2),1); % Fraction of electron transport rate partitioned to PEP
Photosynthesis.sco25 = 2590*ones(2*size(Canopy.deltaLAI,2),1); % Rubisco CO2-Oxygen Specificity @ 25C
Photosynthesis.kn = Options.kn; % Exponential nitrogen decay
Photosynthesis.PcSeasonalVariation = Options.PcSeasonalVariation; % Decay in photosynthetic parameters
Photosynthesis.linearDTD = Options.linearDTD; % Rate of linear decline of photosynthetic capacity
Photosynthesis.Jmax_Vcmax = Options.Jmax_Vcmax; % Ratio of Jmax/Vcmax
Photosynthesis.nitrogenDecay = Options.nitrogenDecay; % Vertical exponential decline rate of nitrogen content
end
