%% Function to compute wind profile within the canopy
function [Canopy] = callCanopyWindProfile(Canopy,Weather,Soil)

dz = diff(Canopy.z);

%% Initialize canopy wind profile [m s-1]
windProfile = linspace(Soil.wind,Weather.wind,length(Canopy.z)+1);
dz = [dz(1) dz];
windProfileOld = windProfile;

%% Computational variables
relax = 0.5; % Relaxation coefficient for successive over relaxation
sumSquareError = 1; % Initialize sum of squares error between previous and current wind profile
minError = 0.00001; % Sum of squares error tolerance between previous and current wind profile
errorLoop = 1; % Initialize iterative loop
maxLoop = 500; % Maximum iterations for convergence

%% Set up finite difference matrix
aMatrix = zeros(length(Canopy.z)+1); % Initialize A matrix
bVector = zeros(length(Canopy.z)+1,1); % Initiaize b vector

%% Compute iterative finite difference model
while (errorLoop < maxLoop) && (sumSquareError > minError)
    
    dU = diff(windProfileOld);
    Km = Canopy.lm^2*abs(dU./dz);
    Km = [Km(1),Km];    
    dKm_dz = diff(Km)./dz;
    dKm_dz = [dKm_dz(1),dKm_dz];    
    
    % Bottom boundary node
    aMatrix(1,1) = 1;
    bVector(1) = Soil.wind;   
    
    % Internal nodes
    for zLoop = 2:1:length(windProfileOld)-1
          aMatrix(zLoop,zLoop-1) = -Km(zLoop)/dz(zLoop)^2 ;
          aMatrix(zLoop,zLoop) = 2*Km(zLoop)/dz(zLoop)^2 + dKm_dz(zLoop)/dz(zLoop) ...
              + 0.5*Canopy.cd*Canopy.LADProfile(zLoop)*abs(windProfileOld(zLoop));
          aMatrix(zLoop,zLoop+1) = -Km(zLoop)/dz(zLoop)^2 - dKm_dz(zLoop)/dz(zLoop);          
    end

    % Top boundary node
    aMatrix(end,end) = 1;
    bVector(end) = Weather.wind;
    
    %% Compute solution
    windProfile = (aMatrix\bVector)';
    index = isnan(windProfile);
    windProfile(index) = 0;
    
    %% Picard loop convergence check
    errorLoop = errorLoop+1;
    sumSquareError = sum((windProfile-windProfileOld).^2)/length(windProfile);
    
    %% Successive over relaxation
    windProfileOld = (1-relax)*windProfileOld+relax*windProfile;
    
end
Canopy.windProfile = windProfile;
end