%% Input
modelSetupError = 0;

if (Options.model(1) == 1 || Options.model(1) == 2) && Options.canopyPoints(1) > 2
    str = 'User Input Error: Big leaf model and Sun-Shade leaf model are Single layer models.\n%s Please, change the number of cnaopy layers (canopyLayers) to 1; deactivate scalar transport, (scalarMixing) to 0';
    modelSetupError = 1;
end
if (Options.model(1) == 3 || Options.model(1) == 4) && Options.canopyPoints(1) < 3
    str = 'User Input Error: You are using Multi layer models.\n%s Please, change the keep the number of cnaopy layers (canopyLayers) > 1';
    modelSetupError = 1;
end
if Options.model(1) == 3 && Options.mixingSwitch(1) == 1
    str = 'User Input Error: You are using Multi layer models with well mixed assumption.\n%s Please, deactivate scalar transport, (scalarMixing) to 0';
    modelSetupError = 1;
end
if Options.model(1) == 4 && Options.mixingSwitch(1) == 0
   str = 'User Input Error: You are using Multi layer models with scalar transport.\n%s Please, activate scalar transport, (scalarMixing) to 1';
    modelSetupError = 1;
end

if Options.LEBSwitch(:) ~= 1 
    warning("Leaf energy balance is deactivated")
elseif Options.nitrogenDecay(:) ~= 1 
    warning("Vertical leaf nitrogen variation is deactivated")
elseif Options.PcSeasonalVariation(:) ~= 1
    warning("Seasonal variation in leaf photosynthetic capacity is deactivated")
end

if modelSetupError == 1
   error(str) 
end
