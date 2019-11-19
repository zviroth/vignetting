close all
clear all
% contrasts = [1 20 100];
contrasts = [1 100];
% contrasts = [20 100];
interpMethod = 'linear';
normResp = 1;
numPhases=1;
freqs = linspace(1,100,3);
% freqs= [10 16 24];
phases = linspace(0,pi,numPhases+1);
phases = phases(1:end-1);
whichLev = 5;

% rfSize = 90;

% create the stimulus
gratingSize = 257;
apertureSize = 64;
% apertureSize = 4;
% gratingFreq = 16;
gratingDirection = pi;
numOri = 16;
gratingOrientations = linspace(0, pi, numOri);

rfSize = 20;
rfY = (gratingSize+1)/2; rfX= 40+(gratingSize+1)/2;
% extract the response from RF
rf = mkDisc([gratingSize gratingSize], rfSize, [rfX rfY]);
        
        
% construct quad frequency filters
numOrientations = 6;
bandwidth = 1;
dims = [gratingSize gratingSize];
numLevels = maxLevel(dims,bandwidth);
[freqRespsImag, freqRespsReal, pind] = makeQuadFRs(dims, numLevels, numOrientations, bandwidth);


rows=length(contrasts);
cols=4;
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
                gratings(:,:,icontrast,iori,ifreq,iphase) = outerAperture .* mkSine(gratingSize, freqs(ifreq), gratingOrientations(iori), contrast/100,phases(iphase));
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
                    sumBandsContrastOriFreqPhase{iLev}(:,:,icontrast,iori,ifreq,iphase) = temp;
                end
                sumBandsContrastOriFreq{iLev} = mean(sumBandsContrastOriFreqPhase{iLev},6);%mean over phases
            end
            sumBandsContrastOri{iLev} = mean(sumBandsContrastOriFreq{iLev},5);%mean over freq
        end
    end
end

%sum over orientations, to get SF tuning
for icontrast=1:length(contrasts)
    for ifreq=1:length(freqs)
        sumBandsContrastFreq{iLev} = mean(sumBandsContrastOriFreq{iLev},4);%mean over orientation
    end
end

%%
%now focusing only on ONE LEVEL

%get RF response for each orientation, contrast, and frequency
for icontrast=1:length(contrasts)
    %orientation tuning per contrast
    temp = squeeze(sumBandsContrastOri{whichLev}(:,:,icontrast,:,:,:)) .* repmat(rf, [1 1 numOri]);
    temp(temp == 0) = NaN;
    rfOri(icontrast,:) = squeeze(nanmean(nanmean(temp,1),2));
    temp = squeeze(sumBandsContrastOri{whichLev}(:,:,icontrast,:,:,:));
    temp(temp == 0) = NaN;
    v1Ori(icontrast,:) = squeeze(nanmean(nanmean(temp,1),2));
    
    %frequency tuning per contrast
    temp = squeeze(sumBandsContrastFreq{whichLev}(:,:,icontrast,:,:,:)) .* repmat(rf, [1 1 length(freqs)]);
    temp(temp == 0) = NaN;
    rfFreq(icontrast,:) = squeeze(nanmean(nanmean(temp,1),2));
%     temp = squeeze(sumBandsContrastFreq{whichLev}(:,:,icontrast,:,:,:));
%     temp(temp == 0) = NaN;
%     v1Freq(icontrast,:) = squeeze(nanmean(nanmean(temp,1),2));
end
    


%     % loop over stimulus orientations
%     for iori = 1:length(gratingOrientations)
%         for ifreq=1:length(freqs)
%             
%             for iLev = whichLev %1:numLevels
%                 temp = squeeze(sumBandsContrastOri{iLev}(:,:,icontrast,:,:,:)) .* repmat(rf, [1 1 numOri]);
%                 temp(temp == 0) = NaN;
%                 rfResp(iLev,:) = squeeze(nanmean(nanmean(temp,1),2));
%                 
%                 temp = sumBandsContrastOri{iLev};
%                 temp(temp == 0) = NaN;
%                 v1Resp(iLev,:) = squeeze(nanmean(nanmean(temp,1),2));%mean across entire image
%             end
%             
% %             data(icontrast,:) = rfResp(whichLev,:);
% %             v1Data(icontrast,:) = v1Resp(whichLev,:);
%         end
%     end
% end

%%
interpFactor = 10;
interpGratingOrientations = linspace(0, pi, interpFactor*numOri);
interpFreqs = linspace(freqs(1), freqs(end), interpFactor*length(freqs));
    
interpRfOri = interp1(gratingOrientations,rfOri',interpGratingOrientations,interpMethod)';
interpRfFreq = interp1(freqs,rfFreq',interpFreqs,interpMethod)';
for icontrast=1:length(contrasts)

    x = 1:length(interpRfOri(icontrast,:));
    % Find the half max value.
    halfMax = (min(interpRfOri(icontrast,:)) + max(interpRfOri(icontrast,:))) / 2;
    % Find where the data first drops below half the max.
    index1 = find(interpRfOri(icontrast,:) >= halfMax, 1, 'first');
    % Find where the data last rises above half the max.
    index2 = find(interpRfOri(icontrast,:) >= halfMax, 1, 'last');
    fwhm = index2-index1 + 1; % FWHM in indexes.
    % OR, if you have an x vector
    fwhmx = interpGratingOrientations(index2) - interpGratingOrientations(index1);
    disp(sprintf('Contrast=%i: FWHM=%i', contrasts(icontrast), rad2deg(fwhmx)));
    

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

    % RF orientation tuning
    subplot(rows,cols,2 + (icontrast-1)*cols); cla
    plot(interpGratingOrientations, interpRfOri(icontrast,:), 'k.-', 'markerfacecolor', 'b', 'markersize',5);
    line(interpGratingOrientations([index1 index2]), [halfMax halfMax], 'color', 'red', 'linewidth', 2);
    xaxis([0 pi]);
    %yaxis([-100 100]);
    axis square
    xlabel('Stimulus orientation');
    ylabel('Simulated response amplitude');
    drawPublishAxis('xTick', [0 pi],'yLabelOffset', -18/64,'xLabelOffset', -8/64,'labelFontSize',10);
    
    % V1 orientation tuning
    subplot(rows,cols,3 + (icontrast-1)*cols); cla
    plot(gratingOrientations, v1Ori(icontrast,:), 'k.-', 'markerfacecolor', 'b', 'markersize',5);
    xaxis([0 pi]);
    %yaxis([-100 100]);
    axis square
    xlabel('Stimulus orientation');
    ylabel('Simulated response amplitude');
    drawPublishAxis('xTick', [0 pi],'yLabelOffset', -18/64,'xLabelOffset', -8/64,'labelFontSize',10);

    
    %RF frequency tuning
    subplot(rows,cols,4 + (icontrast-1)*cols); cla
    plot(interpFreqs, interpRfFreq(icontrast,:), 'k.-', 'markerfacecolor', 'b', 'markersize',5);
%     line(interpFreqs([index1 index2]), [halfMax halfMax], 'color', 'red', 'linewidth', 2);
    xaxis([interpFreqs(1) interpFreqs(end)]);
    %yaxis([-100 100]);
    axis square
    xlabel('Stimulus frequency');
    ylabel('Simulated response amplitude');
    drawPublishAxis('xTick', [interpFreqs(1) interpFreqs(end)],'yLabelOffset', -18/64,'xLabelOffset', -8/64,'labelFontSize',10);
end

%%
rows=1;
cols=3;
figure(2); clf
subplot(rows,cols,1)
plot(interpGratingOrientations, zscore(interpRfOri,0,2), 'markersize',5);
subplot(rows,cols,2)
% plot(gratingOrientations, zscore(v1Data,0,2), 'markersize',5);
plot(gratingOrientations, zscore(v1Ori,0,2), 'markersize',5);

subplot(rows,cols,3)
plot(interpFreqs, zscore(interpRfFreq,0,2), 'markersize',5);
%     plot(interpGratingOrientations, interpData, 'k.-', 'markerfacecolor', 'b', 'markersize',5);
% line(interpGratingOrientations([index1 index2]), [halfMax halfMax], 'color', 'red', 'linewidth', 2);

