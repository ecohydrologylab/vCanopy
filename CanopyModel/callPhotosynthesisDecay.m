% Plant ageing effect
Photosynthesis.years = unique(Weather.year);
hour = unique(Weather.hour);
decayStartDOY = 185;
Photosynthesis.decayStart = [];

if Photosynthesis.PcSeasonalVariation == 1 

    Photosynthesis.linearDTD = Options.linearDTD; % Daily reduction in plant Vcmax
    
    if isnan(Weather.GDD(1))
        disp("Applying photosynthesis decay from predecided DOY")
        Photosynthesis.decayStart = decayStartDOY*ones(1,length(Photosynthesis.years));
        for yLoop = 1:1:length(Photosynthesis.years)
            yearIndex = find(Weather.year == Photosynthesis.years(yLoop))';
            Tair = reshape(Weather.tAir(yearIndex),[length(hour),length(yearIndex)/length(hour)]);
            Tmin = min(Tair); Tmax = max(Tair);
            gdd = max(0.5*(Tmax+Tmin)-10,0);
            offset = max(decayStartDOY - min(Weather.julian(yearIndex)),-1);
            cumgdd = [-fliplr(cumsum(fliplr(gdd(1:offset+1)))),cumsum(gdd(offset+2:end))];
            Weather.GDD(yearIndex) = reshape(repmat(Options.GDDpre + cumgdd,[length(hour),1]),[],1);
        end
    else
        for yLoop = 1:1:length(Photosynthesis.years)
            yearindex = find(Weather.year == Photosynthesis.years(yLoop));
            doy = Weather.julian(yearindex);
            gdd = Weather.GDD(yearindex);
            decayDelay = doy(gdd >= Options.GDDpre);
            if ~isempty(decayDelay)
                Photosynthesis.decayStart(yLoop) = decayDelay(1);
                disp(strcat("Applying photosynthesis decay from GDD of ",num2str(Options.GDDpre)," for ",num2str(Photosynthesis.years(yLoop))))
            else
                warning(strcat("Did not find a DOY for GDD of ",num2str(Options.GDDpre),", so not applying photosynthesis decay for ",...
                    num2str(Photosynthesis.years(yLoop))))
                Photosynthesis.decayStart(yLoop) = 9999;
            end
        end
    end
    
else
    
    Photosynthesis.linearDTD = 0;%
    Photosynthesis.decayStart = zeros(1,length(Photosynthesis.years));%
    
end