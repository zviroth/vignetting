close all
clear all
% contrasts = [1 20 100];
contrasts = [1 100];
% contrasts = [20 100];
interpMethod = 'linear';
normResp = 1;
numPhases=1;
numFreqs = 20;

freqs = linspace(1,150,numFreqs);
freqs = logspace(-1.5,1,numFreqs);
freqs = logspace(-1.5,-0.5,numFreqs);
freqs = logspace(-1.8,-0.8,numFreqs);
freqs = logspace(-2.3,-1.3,numFreqs);

phases = linspace(0,pi,numPhases+1);
phases = phases(1:end-1);

% create the stimulus
gratingSize = 256;
apertureSize = 64;

numOri = 4;
gratingOrientations = linspace(0, pi, numOri+1);
gratingOrientations = gratingOrientations(1:numOri);

rfSize = 15;

rfY = (gratingSize+1)/2; rfX= 58+(gratingSize+1)/2;%bottom
rfY = (gratingSize+1)/2; rfX= (gratingSize+1)/2;%bottom
% rfY = 58 + (gratingSize+1)/2; rfX= (gratingSize+1)/2;


% extract the response from RF
rf = mkDisc([gratingSize gratingSize], rfSize, [rfX rfY]);


% construct quad frequency filters
numOrientations = 4;
bandwidth = 1;
dims = [gratingSize gratingSize];
numLevels = maxLevel(dims,bandwidth);
[freqRespsImag, freqRespsReal, pind] = makeQuadFRs(dims, numLevels, numOrientations, bandwidth);



% initizlize grating
gratings = [];

for icontrast=1:length(contrasts)
    contrast = contrasts(icontrast);
    % loop over aperture sizes
    outerAperture = mkDisc([gratingSize gratingSize], apertureSize, [gratingSize gratingSize]/2, 0);
    % loop over stimulus orientations
    for iori = 1:length(gratingOrientations)
        for ifreq=1:length(freqs)
            for iphase = 1:numPhases
                %         gratings(icontrast,:,:,iOri,ifreq) = outerAperture .* mkSine(gratingSize, gratingFreq, gratingOrientations(iOri), contrast/100);
                %mkSine uses freq as cycle length in pixels!!
%                 gratings(:,:,icontrast,iori,ifreq,iphase) = outerAperture .* mkSine(gratingSize, freqs(ifreq), gratingOrientations(iori), contrast/100,phases(iphase));
                gratings(:,:,icontrast,iori,ifreq,iphase) = outerAperture .* mkSine(gratingSize, 1/freqs(ifreq), gratingOrientations(iori), contrast/100,phases(iphase));
            end
        end
    end
end


for icontrast=1:length(contrasts)
    % loop over stimulus orientations
    for iori = 1:length(gratingOrientations)
        for ifreq=1:length(freqs)
            for iphase = 1:numPhases
                % build pyramid for oriented grating
                [pyr, pind] = buildQuadBands(squeeze(gratings(:,:,icontrast,iori,ifreq,iphase)), freqRespsImag, freqRespsReal);
                for iLev = 1:numLevels
                    % loop over levels and orientations of the pyramid
                    % initialize output
                    temp = zeros(gratingSize, gratingSize);
                    for orientation = 1:numOrientations
                        if normResp
                            nEnergies = normEnergies(pyr,pind,numOrientations,0.1);
                            thisBand = abs(accessSteerBand(nEnergies,pind,numOrientations,iLev,orientation));
                        else
                            thisBand = abs(accessSteerBand(pyr, pind, numOrientations,iLev, orientation)).^2;
                        end
                        temp = temp + thisBand;
                    end
                    sumBandsContrastOriFreqPhase{iLev}(:,:,icontrast,iori,ifreq,iphase) = temp;
                end
                
            end
        end
    end
end

%sum over orientations, to get SF tuning

for iLev=1:numLevels
    sumBandsContrastOriFreq{iLev} = mean(sumBandsContrastOriFreqPhase{iLev},6);%mean over phases
    sumBandsContrastFreq{iLev} = mean(sumBandsContrastOriFreq{iLev},4);%mean over orientation
end

for iLev = 1:numLevels
    levMean(iLev) =  mean(sumBandsContrastOriFreq{iLev}(:));
    levMax(iLev) = max(sumBandsContrastOriFreq{iLev}(:));
end
[maxVal,whichLev] = max(levMean); whichLev
% [maxVal,whichLev] = max(levMax); whichLev

%%
%now focusing only on ONE LEVEL

%get RF response for each orientation, contrast, and frequency
for icontrast=1:length(contrasts)
    
    %frequency tuning per contrast
    temp = squeeze(sumBandsContrastOriFreq{whichLev}(:,:,icontrast,:,:,:)) .* repmat(rf, [1 1 numOri length(freqs)]);
    temp(temp == 0) = NaN;
    rfFreq(icontrast,:) = squeeze(nanmean(nanmean(nanmean(temp,1),2),3));
    temp = squeeze(sumBandsContrastOriFreq{whichLev}(:,:,icontrast,:,:,:));
    temp(temp == 0) = NaN;
    v1Freq(icontrast,:) = squeeze(nanmean(nanmean(nanmean(temp,1),2),3));
end


%%
interpFactor = 10;
interpFreqs = logspace(log10(freqs(1)), log10(freqs(end)), interpFactor*length(freqs));

interpRfFreq = interp1(freqs,rfFreq',interpFreqs,interpMethod,'extrap')';
interpV1Freq = interp1(freqs,v1Freq',interpFreqs,interpMethod,'extrap')';
if length(contrasts)==1
    interpRfFreq = interpRfFreq';
    interpV1Freq = interpV1Freq';
end

rows=length(contrasts);
cols=4;
for icontrast=1:length(contrasts)
    
    [fwhmx, halfMax, index1, index2] = findFWHM(interpFreqs,interpRfFreq(icontrast,:));
    disp(sprintf('RF: Contrast=%i: FWHM=%i', contrasts(icontrast), rad2deg(fwhmx)));
    
    
    
    
    figure(1)
    
    % stimulus image
    subplot(rows,cols,1 + (icontrast-1)*cols); cla
    imagesc(squeeze(gratings(:,:,icontrast,1,1,1)));
    temp = squeeze(gratings(:,:,icontrast,2,1,1));
    temp = temp + rf./2;
    imagesc(temp);
    axis image ; axis off;
    hold on
    colormap(gray);
    caxis([-1 1])
    
    %RF frequency tuning
    subplot(rows,cols,2 + (icontrast-1)*cols); cla
    %     plot(interpFreqs, interpRfFreq(icontrast,:), 'k.-', 'markerfacecolor', 'b', 'markersize',5);
    semilogx(interpFreqs, interpRfFreq(icontrast,:), 'k.-', 'markerfacecolor', 'b', 'markersize',5);
    line(interpFreqs([index1 index2]), [halfMax halfMax], 'color', 'red', 'linewidth', 2);
    xaxis([interpFreqs(1) interpFreqs(end)]);
    %yaxis([-100 100]);
    axis square
    xlabel('Stimulus frequency');
    ylabel('Simulated response amplitude');
    drawPublishAxis('xTick', [interpFreqs(1) interpFreqs(end)],'yLabelOffset', -6/64,'xLabelOffset', -6/64,'labelFontSize',10);
    
    if icontrast==1
        drawPublishAxis('xTick', [interpFreqs(1) interpFreqs(end)],'yLabelOffset', -6/64,'xLabelOffset', -6/64,'labelFontSize',10,'titleStr','white RF');
    else
        drawPublishAxis('xTick', [interpFreqs(1) interpFreqs(end)],'yLabelOffset', -6/64,'xLabelOffset', -6/64,'labelFontSize',10);
    end
    
    
    
    [fwhmx, halfMax, index1, index2] = findFWHM(interpFreqs,interpV1Freq(icontrast,:));
    disp(sprintf('V1: Contrast=%i: FWHM=%i', contrasts(icontrast), rad2deg(fwhmx)));
    
    
    
    %V1 frequency tuning
    subplot(rows,cols,3 + (icontrast-1)*cols); cla
    %     plot(interpFreqs, interpV1Freq(icontrast,:), 'k.-', 'markerfacecolor', 'b', 'markersize',5);
    semilogx(interpFreqs, interpV1Freq(icontrast,:), 'k.-', 'markerfacecolor', 'b', 'markersize',5);
    line(interpFreqs([index1 index2]), [halfMax halfMax], 'color', 'red', 'linewidth', 2);
    xaxis([interpFreqs(1) interpFreqs(end)]);
    %yaxis([-100 100]);
    axis square
    xlabel('Stimulus frequency');
    ylabel('Simulated response amplitude');
    
    drawPublishAxis('xTick', [freqs(1) freqs(end)],'yLabelOffset', -6/64,'xLabelOffset', -8/64,'labelFontSize',10);
    
    if icontrast==1
        drawPublishAxis('xTick', [freqs(1) freqs(end)],'yLabelOffset', -6/64,'xLabelOffset', -8/64,'labelFontSize',10,'titleStr','entire FOV');
    else
        drawPublishAxis('xTick', [freqs(1) freqs(end)],'yLabelOffset', -6/64,'xLabelOffset', -8/64,'labelFontSize',10);
    end
end
set(gcf,'position',[2 5 35 15]);
%
rows=1;
cols=3;
figure(2); clf

subplot(rows,cols,1)
% plot(interpFreqs, zscore(interpRfFreq,0,2), 'markersize',5);
semilogx(interpFreqs, zscore(interpRfFreq,0,2), 'markersize',5);
title('z-scored RF SF tuning');
xlabel('spatial frequency');

subplot(rows,cols,2)
% plot(interpFreqs, zscore(interpV1Freq,0,2), 'markersize',5);
semilogx(interpFreqs, zscore(interpV1Freq,0,2), 'markersize',5);
title('z-scored v1 SF tuning');
xlabel('spatial frequency');

%     plot(interpGratingOrientations, interpData, 'k.-', 'markerfacecolor', 'b', 'markersize',5);
% line(interpGratingOrientations([index1 index2]), [halfMax halfMax], 'color', 'red', 'linewidth', 2);
set(gcf,'position',[100 300 1100 250]);






% function [fwhmx, halfMax, index1, index2] = findFWHM(X,Y)
% % Find the half max value.
% halfMax = (min(Y) + max(Y)) / 2;
% index1 = find(Y >= halfMax, 1, 'first');
% index2 = find(Y >= halfMax, 1, 'last');
% fwhm = index2-index1 + 1; % FWHM in indexes.
% fwhmx = X(index2) - X(index1);
%
% end