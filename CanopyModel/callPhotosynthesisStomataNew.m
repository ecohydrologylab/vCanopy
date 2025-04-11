function [Photosynthesis,Stomata] = callPhotosynthesisStomata(Canopy)

%% Photosynthesis
Photosynthesis = table;
Photosynthesis.leafID = 1+(1:1:2*length(Canopy.deltaLAI(1,:)))'/100.0;

% Exponential change starts
Photosynthesis.vpr = 80*ones(2*length(Canopy.deltaLAI(1,:)),1); 
Photosynthesis.vcmax25 = 65*ones(2*length(Canopy.deltaLAI(1,:)),1);
Photosynthesis.jmax25 = 350*ones(2*length(Canopy.deltaLAI(1,:)),1);
Photosynthesis.vpmax25 = 110*ones(2*length(Canopy.deltaLAI(1,:)),1);
Photosynthesis.rd25 = 0.05*ones(2*length(Canopy.deltaLAI(1,:)),1);
% Exponential change stops

Photosynthesis.theta = 0.76*ones(2*length(Canopy.deltaLAI(1,:)),1);
Photosynthesis.O2 = 350*ones(2*length(Canopy.deltaLAI(1,:)),1);
Photosynthesis.gbs = 0.003*ones(2*length(Canopy.deltaLAI(1,:)),1);
Photosynthesis.alpha = 0*ones(2*length(Canopy.deltaLAI(1,:)),1);
Photosynthesis.waterStress = 0*ones(2*length(Canopy.deltaLAI(1,:)),1);
Photosynthesis.x = 0.4*ones(2*length(Canopy.deltaLAI(1,:)),1);

%% Stomata values
Stomata = table;
Stomata.slope = 3*ones(2*length(Canopy.deltaLAI),1);
Stomata.intercept = 0.08*ones(2*length(Canopy.deltaLAI),1);

end
