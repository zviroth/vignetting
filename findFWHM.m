function [fwhmx, halfMax, index1, index2] = findFWHM(X,Y)
% Find the half max value.
halfMax = (min(Y) + max(Y)) / 2;
index1 = find(Y >= halfMax, 1, 'first');
index2 = find(Y >= halfMax, 1, 'last');
fwhm = index2-index1 + 1; % FWHM in indexes.
fwhmx = X(index2) - X(index1);

end