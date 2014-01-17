%% GouyCost - example beam path optimization cost function
function cost = gouyCost(path)
    % this is an example cost function for use in optimizePath or
    % chooseComponents. It takes a beam path as an argument, and produces a
    % number which represents the cost. The optimization functions will then
    % try to minimize this cost function by moving components and/or switching
    % out components.

    %% Beam size cost
    % beam size minima only happen on lenses or at waists
    zdomain = [[path.components.z] path.getWaists().'];
    qs = path.qPropagate(zdomain);

    % this is the smallest radius in the beam path
    smallestRadius = min([qs.beamWidth]);

    % we don't want radii smaller than this
    badRadius = 0.0005;

    % strongly penalize small radii
    radiusCost = exp(-(smallestRadius/badRadius)^3);

    %% Gouy phase cost
    % also penalize gouy phase that is not close to 90 degress
    angleCost = cos(path.gouySeparation('WFS A','WFS B')*pi/180)^2;

    %% Total cost
    % total cost is the sum of the two components
    cost = radiusCost + angleCost;
end