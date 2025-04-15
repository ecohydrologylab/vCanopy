function [Photosynthesis] = callModifyPhotosyntheticCapacity(orgPhotosynthesis,Weather,Canopy)

Photosynthesis = orgPhotosynthesis;
for yLoop = 1:1:length(Photosynthesis.years)
    if Photosynthesis.PcSeasonalVariation == 1
        %% Simple linear decay per day
        deltaPC = -max((Weather.julian - Photosynthesis.decayStart(yLoop)),0)*Photosynthesis.linearDTD;
    else
        deltaPC = 0;
    end
    
    Photosynthesis.vcmax25 = orgPhotosynthesis.vcmax25+deltaPC;
    Photosynthesis.jmax25 = Photosynthesis.vcmax25*orgPhotosynthesis.Jmax_Vcmax; % From Miner et. al. 2018
end

%% Nitrogen Decay ======================
CLAI = Canopy.LAI;
if length(CLAI) == 2
    CLAI(end) = 0;
end
CLAI = [CLAI,CLAI];

if Photosynthesis.nitrogenDecay == 1
    kn = Photosynthesis.kn;%
    Photosynthesis.vpr = Photosynthesis.vpr.*exp(-kn*CLAI)';
    Photosynthesis.vcmax25 = Photosynthesis.vcmax25.*exp(-kn*CLAI)';
    Photosynthesis.jmax25 = Photosynthesis.jmax25.*exp(-kn*CLAI)';
    Photosynthesis.vpmax25 = Photosynthesis.vpmax25.*exp(-kn*CLAI)';
end

end