Initialize = table();

% Use previous canopy results to initialize the canopy and leaf
if sum(isnan(PreviousCanopy.tAirProfile(2:end))) == 0 || ...
        sum(isnan(PreviousCanopy.eaProfile(2:end))) == 0
    Canopy.tAirProfile = PreviousCanopy.tAirProfile;
    Canopy.eaProfile = PreviousCanopy.eaProfile;
    Canopy.caProfile = PreviousCanopy.caProfile;
    
    %% Dummy Initializations
    gs = 0.1*ones(2*length(Canopy.z),1);
    ci = 0.7*[Canopy.caProfile(2:end),Canopy.caProfile(2:end)]';
    cb = 0.8*[Canopy.caProfile(2:end),Canopy.caProfile(2:end)]';
    cbs = 10*[Canopy.caProfile(2:end),Canopy.caProfile(2:end)]';
    aNet = 0.1*[Canopy.caProfile(2:end),Canopy.caProfile(2:end)]';
    GammaStar = 0.0143*ones(2*length(Canopy.z),1);
    tLeaf = [Canopy.tAirProfile(2:end),Canopy.tAirProfile(2:end)]';
    eb = 1.1*[Canopy.eaProfile(2:end),Canopy.eaProfile(2:end)]';
    
    sunIndex = logical(Canopy.sunFraction);
    shadeIndex = logical(Canopy.shadeFraction);
    previousSunIndex = logical(PreviousCanopy.sunAnet);
    previousShadeIndex = logical(PreviousCanopy.shadeAnet);
    
    Initialize.gs = [PreviousCanopy.sunGs.*sunIndex,PreviousCanopy.shadeGs.*shadeIndex]';
    Initialize.ci = [PreviousCanopy.sunCi.*sunIndex,PreviousCanopy.shadeCi.*shadeIndex]';
    Initialize.cb = [PreviousCanopy.sunCb.*sunIndex,PreviousCanopy.shadeCb.*shadeIndex]';
    Initialize.cbs = [PreviousCanopy.sunCbs.*sunIndex,PreviousCanopy.shadeCbs.*shadeIndex]';
    Initialize.aNet = [PreviousCanopy.sunAnet.*sunIndex,PreviousCanopy.shadeAnet.*shadeIndex]';
    Initialize.GammaStar = [PreviousCanopy.sunGammaStar.*sunIndex,PreviousCanopy.shadeGammaStar.*shadeIndex]';
    Initialize.eb = [PreviousCanopy.sunEb.*sunIndex,PreviousCanopy.shadeEb.*shadeIndex]';
    Initialize.tLeaf = [PreviousCanopy.sunTleaf.*sunIndex,PreviousCanopy.shadeTleaf.*shadeIndex]';
    
    % Update the variables
    leafIndex = [sunIndex,shadeIndex]';
    previousLeafIndex = [previousSunIndex,previousShadeIndex]';
    
    Initialize.gs = Initialize.gs.*(1-abs(leafIndex-previousLeafIndex)) + ...
        abs(leafIndex-previousLeafIndex).*gs;
    Initialize.ci = Initialize.ci.*(1-abs(leafIndex-previousLeafIndex)) + ...
        abs(leafIndex-previousLeafIndex).*ci;
    Initialize.cb = Initialize.cb.*(1-abs(leafIndex-previousLeafIndex)) + ...
        abs(leafIndex-previousLeafIndex).*cb;
    Initialize.cbs = Initialize.cbs.*(1-abs(leafIndex-previousLeafIndex)) + ...
        abs(leafIndex-previousLeafIndex).*cbs;
    Initialize.aNet = Initialize.aNet.*(1-abs(leafIndex-previousLeafIndex)) + ...
        abs(leafIndex-previousLeafIndex).*aNet;
    Initialize.GammaStar = Initialize.GammaStar.*(1-abs(leafIndex-previousLeafIndex)) + ...
        abs(leafIndex-previousLeafIndex).*GammaStar;
    Initialize.eb = Initialize.eb.*(1-abs(leafIndex-previousLeafIndex)) + ...
        abs(leafIndex-previousLeafIndex).*eb;
    Initialize.tLeaf = Initialize.tLeaf.*(1-abs(leafIndex-previousLeafIndex)) + ...
        abs(leafIndex-previousLeafIndex).*tLeaf;
    
    
else % When previous canopy was faulty use default values for canopy initialization
    Canopy.tAirProfile = linspace(Soil.temperature, ...
        Weather.tAir,length(Canopy.z)+1);
    Canopy.eaProfile = linspace(Weather.ea+100,Weather.ea,length(Canopy.z)+1);
    Canopy.caProfile = linspace(Weather.ca+5,Weather.ca,length(Canopy.z)+1);
    
    Initialize.gs = 0.1*ones(2*length(Canopy.z),1);
    Initialize.ci = 0.7*[Canopy.caProfile(2:end),Canopy.caProfile(2:end)]';
    Initialize.cb = 0.8*[Canopy.caProfile(2:end),Canopy.caProfile(2:end)]';
    Initialize.cbs = 10*[Canopy.caProfile(2:end),Canopy.caProfile(2:end)]';
    Initialize.aNet = 0.1*[Canopy.caProfile(2:end),Canopy.caProfile(2:end)]';
    Initialize.GammaStar = 0.0143*ones(2*length(Canopy.z),1);
    Initialize.tLeaf = [Canopy.tAirProfile(2:end),Canopy.tAirProfile(2:end)]';
    Initialize.eb = 1.1*[Canopy.eaProfile(2:end),Canopy.eaProfile(2:end)]';
    
end
