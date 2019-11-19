close all
clear all
% contrasts = [1 20 100];
contrasts = [1 100];
% contrasts = [20 100];
interpMethod = 'linear';
normResp = 1;
numPhases=2;
phases = linspace(0,pi,numPhases+1);
phases = phases(1:end-1);
whichLev = 5;

% rfSize = 90;

% create the stimulus
gratingSize = 257;
apertureSize = 64;
% apertureSize = 4;
gratingFreq = 16;
gratingDirection = pi;
numOri = 16;
gratingOrientations = linspace(0, pi, numOri);

rfSize = 20;
rfY = (gratingSize+1)/2; rfX= 40+(gratingSize+1)/2;

% construct quad frequency filters
numOrientations = 6;
bandwidth = 1;
dims = [gratingSize gratingSize];
numLevels = maxLevel(dims,bandwidth);
[freqRespsImag, freqRespsReal, pind] = makeQuadFRs(dims, numLevels, numOrientations, bandwidth);


rows=length(contrasts);
cols=3;
    % initizlize grating
    gratings = [];
    
for icontrast=1:length(contrasts)
    contrast = contrasts(icontrast);
%     if ieNotDefined('contrast'), contrast = 100;end%[20 100];end

    

    % loop over aperture sizes
    outerAperture = mkDisc([gratingSize gratingSize], apertureSize, [gratingSize gratingSize]/2, 0);
    % loop over stimulus orientations
    for iOri = 1:length(gratingOrientations)
        for iphase = 1:numPhases
        gratings(icontrast,:,:,iOri) = outerAperture .* mkSine(gratingSize, gratingFreq, gratingOrientations(iOri), contrast/100);
        gratings(icontrast,:,:,iOri,iphase) = outerAperture .* mkSine(gratingSize, gratingFreq, gratingOrientations(iOri), contrast/100,phases(iphase));

        end
    end

    % loop over stimulus orientations
    for iOri = 1:length(gratingOrientations)
        for iphase = 1:numPhases
            % build pyramid for oriented grating
            [pyr, pind] = buildQuadBands(squeeze(gratings(icontrast,:,:,iOri,iphase)), freqRespsImag, freqRespsReal);
            for iLev = whichLev %1:numLevels
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
                sumBandsPhase{iLev}(:,:,iOri,iphase) = temp;
            end
        end
        sumBands{iLev} = squeeze(mean(sumBandsPhase{iLev},4));%mean over phases
    end
    
    % extract the response from RF
    %rf = mkDisc([gratingSize gratingSize], 16, [24 129], 0);
    %rf = mkDisc([gratingSize gratingSize], 32, [32 129], 0);
%     rf = mkDisc([gratingSize gratingSize], rfSize, [(gratingSize+1)/2 (gratingSize+1)/2]);
%     rf = mkDisc([gratingSize gratingSize], rfSize, [50+(gratingSize+1)/2 (gratingSize+1)/2]);
    rf = mkDisc([gratingSize gratingSize], rfSize, [rfX rfY]);
    
%     rf = mkDisc([gratingSize gratingSize], gratingSize, [128 128], 0);
    
    for iLev = whichLev %1:numLevels
        temp = sumBands{iLev} .* repmat(rf, [1 1 numOri]);
        temp(temp == 0) = NaN;
        rfResp(iLev,:) = squeeze(nanmean(nanmean(temp,1),2));
        temp = sumBands{iLev};
        temp(temp == 0) = NaN;
        v1Resp(iLev,:) = squeeze(nanmean(nanmean(temp,1),2));%mean across entire image
    end
    
    data(icontrast,:) = rfResp(whichLev,:);
    v1Data(icontrast,:) = v1Resp(whichLev,:);
    
end
%%
interpFactor = 100;
interpGratingOrientations = linspace(0, pi, interpFactor*numOri);
    
interpData = interp1(gratingOrientations,data',interpGratingOrientations,interpMethod)';

for icontrast=1:length(contrasts)
    contrast = contrasts(icontrast);


%     interpData = interp1(gratingOrientations,data(icontrast,:),interpGratingOrientations,interpMethod);
    
    x = 1:length(interpData(icontrast,:));
    % Find the half max value.
    halfMax = (min(interpData(icontrast,:)) + max(interpData(icontrast,:))) / 2;
    % Find where the data first drops below half the max.
    index1 = find(interpData(icontrast,:) >= halfMax, 1, 'first');
    % Find where the data last rises above half the max.
    index2 = find(interpData(icontrast,:) >= halfMax, 1, 'last');
    fwhm = index2-index1 + 1; % FWHM in indexes.
    % OR, if you have an x vector
    fwhmx = interpGratingOrientations(index2) - interpGratingOrientations(index1);
    disp(sprintf('Contrast=%i: FWHM=%i', contrast, rad2deg(fwhmx)));
    
    % open a new figure
    % smartfig('fig5', 'reuse');
    figure(1)
    
    % stimulus image
    subplot(rows,cols,1 + (icontrast-1)*cols); cla
    imagesc(squeeze(gratings(icontrast,:,:,1,1)));
    
    temp = squeeze(gratings(icontrast,:,:,2,1));
temp = temp + rf./2;
imagesc(temp);


    
    axis image ; axis off;
    hold on
    colormap(gray);
    caxis([-1 1])
    % circle([gratingSize/2 64], rfSize, 100);
%     rf = mkDisc([gratingSize gratingSize], rfSize, [(gratingSize+1)/2 (gratingSize+1)/2]);
%     mglGluDisk((gratingSize+1)/2, (gratingSize+1)/2, rfSize);
%  mglGluDisk(120,120,50);

    %
    subplot(rows,cols,2 + (icontrast-1)*cols); cla
%     thisResp = rfResp(whichLev,:);
%     plot(gratingOrientations, thisResp, '-.k', 'markeredgecolor', 'w', 'markerfacecolor', 'b', 'markersize',10);
    plot(interpGratingOrientations, interpData(icontrast,:), 'k.-', 'markerfacecolor', 'b', 'markersize',5);
    line(interpGratingOrientations([index1 index2]), [halfMax halfMax], 'color', 'red', 'linewidth', 2);
    xaxis([0 pi]);
    %yaxis([-100 100]);
    axis square
    xlabel('Stimulus orientation');
    ylabel('Simulated response amplitude');
    drawPublishAxis('xTick', [0 pi],'yLabelOffset', -18/64,'xLabelOffset', -8/64,'labelFontSize',10);
    

    subplot(rows,cols,3 + (icontrast-1)*cols); cla
    plot(gratingOrientations, v1Data(icontrast,:), 'k.-', 'markerfacecolor', 'b', 'markersize',5);
    
    xaxis([0 pi]);
    %yaxis([-100 100]);
    axis square
    xlabel('Stimulus orientation');
    ylabel('Simulated response amplitude');
    drawPublishAxis('xTick', [0 pi],'yLabelOffset', -18/64,'xLabelOffset', -8/64,'labelFontSize',10);
%     figure(2)
%     hold all
%     plot(interpGratingOrientations, interpData, 'k.-', 'markersize',5);
% %     plot(interpGratingOrientations, interpData, 'k.-', 'markerfacecolor', 'b', 'markersize',5);
%     line(interpGratingOrientations([index1 index2]), [halfMax halfMax], 'color', 'red', 'linewidth', 2);
end

%%
rows=1;
cols=2;
figure(2); clf
subplot(rows,cols,1)
plot(interpGratingOrientations, zscore(interpData,0,2), 'markersize',5);
subplot(rows,cols,2)
% plot(gratingOrientations, zscore(v1Data,0,2), 'markersize',5);
plot(gratingOrientations, v1Data, 'markersize',5);
%     plot(interpGratingOrientations, interpData, 'k.-', 'markerfacecolor', 'b', 'markersize',5);
% line(interpGratingOrientations([index1 index2]), [halfMax halfMax], 'color', 'red', 'linewidth', 2);

