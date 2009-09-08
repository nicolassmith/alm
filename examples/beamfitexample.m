% beam width fitting

% data from kate's elog:

%             diameter (um)
% hole       x             y
% ----     ------        ------
%   3     3550 ± 30     3480 ± 10
%  19     3490 ± 20     3360 ± 10
%  35     3460 ± 20     3310 ± 10

close all
clear classes

zScan = [3 19 35] * 0.0508; % hole spacing is 2 inches
xWidth = [3550 3490 3460] * 1e-6/2;
yWidth = [3480 3360 3310] * 1e-6/2;

goo = beamPath;

goo.seedWaist(1700e-6,2);

goox = goo.fitBeamWidth(zScan,xWidth);
gooy = goo.fitBeamWidth(zScan,yWidth);

zdomain = 0:.01:4;%linspace(min(zScan),max(zScan),100);
figure(1)
clf
hold on
plot(zScan,xWidth,'bo')
goox.plotBeamWidth(zdomain,'b')
goox.plotBeams(zdomain,.5e-3,'b')
plot(zScan,yWidth,'ro')
gooy.plotBeamWidth(zdomain,'r')
gooy.plotBeams(zdomain,-.5e-3,'r')
hold off