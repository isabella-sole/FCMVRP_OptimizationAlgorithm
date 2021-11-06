function lh = updateSalesmanPlot(xopt,idxs,stopsLat,stopsLon)
% Plotting function for tsp_intlinprog example

%   Copyright 2014-2018 The MathWorks, Inc. 

%lh : zero matrix of length equal to x_optimization 
% if ( lh ~= zeros(size(lh)) ) % First time through lh is all zeros
%     delete(lh) % Get rid of unneeded lines
% end

% Find index of value equal to one in the decision variable x_tsp
segments = find(round(xopt)); % Indices to trips in solution

% Loop through the trips then draw them
Lat = zeros(3*length(segments),1);
Lon = zeros(3*length(segments),1);
for ii = 1:length(segments)
    % Nella matrice degli indici idxs trova il segmento ii definedno come
    % punto d'inizio del segmento il punto nella colonna 1 di indxs e come
    % fine del segmento il punto nella colonna 2 di indxs
    start = idxs(segments(ii),1);
    stop = idxs(segments(ii),2);
    % Separate data points with NaN's to plot separate line segments
    % Lat e Lon sono inizializzate tre volte piu grandi perche poi va a
    % salavare per ciascun segmento punto di start, stop e NaN
    Lat(3*ii-2:3*ii) = [stopsLat(start); stopsLat(stop); NaN];
    Lon(3*ii-2:3*ii) = [stopsLon(start); stopsLon(stop); NaN];  
    % (3*ii-2:3*ii) crea un intervallo di tre righe dopo NaN per salbvare i
    % sucessivi tre valori nel vettore [stopsLat(start); stopsLat(stop); NaN];
end
% lh = plot(Lat,Lon,'k:','LineWidth',2);
lh = [Lat, Lon];
% drawnow; % Add new lines to plot