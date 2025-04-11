hLoop = 1;
x = linspace(0.1,1,length(Canopy(hLoop).z))*NaN;

PreviousCanopy.sunGs = 0.01*ones(1,length(x));
PreviousCanopy.sunCi = 0.7*Weather(hLoop).ca.*x;
PreviousCanopy.sunCb = 0.8*Weather(hLoop).ca.*x;
PreviousCanopy.sunCbs = 10*Weather(hLoop).ca.*x;
PreviousCanopy.sunAnet = 0.1*Weather(hLoop).ca.*x;
PreviousCanopy.sunTleaf = Weather(hLoop).tAir*ones(1,length(x));
PreviousCanopy.sunEb = 1.1*Weather(hLoop).ea.*x;
PreviousCanopy.sunGammaStar = 0.0143*x;


PreviousCanopy.shadeGs = 0.01*ones(1,length(x));
PreviousCanopy.shadeCi = 0.7*Weather(hLoop).ca.*x;
PreviousCanopy.shadeCb = 0.8*Weather(hLoop).ca.*x;
PreviousCanopy.shadeCbs = 10*Weather(hLoop).ca.*x;
PreviousCanopy.shadeAnet = 0.1*Weather(hLoop).ca.*x;
PreviousCanopy.shadeTleaf = Weather(hLoop).tAir*ones(1,length(x));
PreviousCanopy.shadeEb = 1.1*Weather(hLoop).ea.*x;
PreviousCanopy.shadeGammaStar = 0.0143*x;

PreviousCanopy.caProfile = NaN*linspace(Weather(hLoop).ca,Weather(hLoop).ca,length(Canopy(hLoop).z)+1);
PreviousCanopy.tAirProfile = NaN*linspace(Soil(hLoop).temperature,...
    Weather(hLoop).tAir,length(Canopy(hLoop).z)+1);
PreviousCanopy.eaProfile = NaN*linspace(Weather(hLoop).ea+200,Weather(hLoop).ea,...
    length(Canopy(hLoop).z)+1);
