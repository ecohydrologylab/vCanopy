PreviousCanopy.sunGs = CanopyOut(hLoop).sunGs;
PreviousCanopy.sunCi = CanopyOut(hLoop).sunCi;
PreviousCanopy.sunCb = CanopyOut(hLoop).sunCb;
PreviousCanopy.sunCbs = CanopyOut(hLoop).sunCbs;
PreviousCanopy.sunAnet = CanopyOut(hLoop).sunAnet;
PreviousCanopy.sunTleaf = CanopyOut(hLoop).sunTleaf;
PreviousCanopy.sunEb = CanopyOut(hLoop).sunEb;
PreviousCanopy.sunGammaStar = CanopyOut(hLoop).sunGammaStar;

PreviousCanopy.shadeGs = CanopyOut(hLoop).shadeGs;
PreviousCanopy.shadeCi = CanopyOut(hLoop).shadeCi;
PreviousCanopy.shadeCb = CanopyOut(hLoop).shadeCb;
PreviousCanopy.shadeCbs = CanopyOut(hLoop).shadeCbs;
PreviousCanopy.shadeAnet = CanopyOut(hLoop).shadeAnet;
PreviousCanopy.shadeTleaf = CanopyOut(hLoop).shadeTleaf;
PreviousCanopy.shadeEb = CanopyOut(hLoop).shadeEb;
PreviousCanopy.shadeGammaStar = CanopyOut(hLoop).shadeGammaStar;
PreviousCanopy.tempLeafLag1 = NaN*x;

PreviousCanopy.caProfile = CanopyOut(hLoop).caProfile;
PreviousCanopy.tAirProfile = CanopyOut(hLoop).tAirProfile;
PreviousCanopy.eaProfile = CanopyOut(hLoop).eaProfile;
% PreviousCanopy.longAbsorbed = CanopyOut.longAbsorbed;

if hLoop < length(Weather)
    Canopy(hLoop+1).tAirProfileLag = CanopyOut(hLoop).tAirProfile;
    Canopy(hLoop+1).sunTleafLag = CanopyOut(hLoop).sunTleaf;
    Canopy(hLoop+1).tLeafProfileLag = CanopyOut(hLoop).tLeafProfile;
    Canopy(hLoop+1).shadeTleafLag = CanopyOut(hLoop).shadeTleaf;
    Canopy(hLoop+1).airQProfileLag = CanopyOut(hLoop).airQProfile;
end