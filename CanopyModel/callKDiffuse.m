%% Function to compute kDiff from x and LAI

function [kD,tauD] = callKDiffuse(sunZenith,x,omega,deltaLAI,Options)

%% Input

if strcmp(Options.direction,'down')
    if strcmp(Options.isotropy,'SOC')
        callPDF = @(sunZenith,zenith) (6/7*(1+2*cosd(zenith)).*sind(zenith).*cosd(zenith));
    elseif strcmp(Options.isotropy,'UOC')
        callPDF = @(sunZenith,zenith) (2.*sind(zenith).*cosd(zenith));
    elseif strcmp(Options.isotropy,'Clear')
        callPDF = @(sunZenith,zenith) ((6/(3+4*(sind(sunZenith)+cosd(sunZenith))))...
            *(1+2*cosd(zenith-sunZenith)).*sind(zenith).*cosd(zenith));
    end
    callKbeam = @(x,zenith) (x^2 + tand(zenith).^2).^0.5 / ...
        (x+1.774*(x+1.182).^(-0.733));
    callTauD = @(zenith,LAI) exp(-callKbeam(x,zenith).*omega.*LAI).*callPDF(sunZenith,zenith);
    
elseif strcmp(Options.direction,'up')
    callPDF = @(zenith) (2.*sind(zenith).*cosd(zenith));
    callKbeam = @(x,zenith) (x^2 + tand(zenith).^2).^0.5 / ...
        (x+1.774*(x+1.182).^(-0.733));
    callTauD = @(zenith,LAI) exp(-callKbeam(x,zenith).*omega.*LAI).*callPDF(zenith);
end

%% Compute tauD and kD

for deltaLAILoop = 1:1:length(deltaLAI)
    if strcmp(Options.direction,'down')
        deltaLAIIndex = [deltaLAILoop:1:length(deltaLAI)];
    elseif strcmp(Options.direction,'up')
        deltaLAIIndex = [1:1:deltaLAILoop];
    end
    tauD(deltaLAILoop) = integral(@(zenith) ...
        callTauD(zenith,sum(deltaLAI(1,deltaLAIIndex))),0,90).*pi/180; % Tranmission coefficient (fraction of radiation transmitted)
    kD(deltaLAILoop) = -log(tauD(deltaLAILoop))./ ...
        (omega .*sum(deltaLAI(1,deltaLAIIndex))); % Extinction coefficient
end

kD(isinf(kD)) = 0;
kD(isnan(kD)) = 0;
% % kD = kD*0 + 0.6;
% kB = @(psi) sqrt(x.^2 + tand(psi).^2)./(x+1.774.*(x+1.182).^(-0.733));
% tauB = @(psi) exp(-kB(psi).*omega.*deltaLAI);
% calltauD = @(psi) 2.*tauB(psi).*sin(psi).*cos(psi);
% tauD = integral(calltauD,0,pi/2);
% kD = -log(tauD)/deltaLAI;

end