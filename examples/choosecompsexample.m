% Example script using chooseComponents method

disp(' chooseComponents example script.')
disp(' ')

focalLengthList = [-.75;.5;1.75;-2;-1;2;3;1;2.5;1.25];
lensList = component.lens(focalLengthList);

goo = beamPath;


% put down the lenses which we want to swap for lenses in our list.
% The positions are used as initial conditions when optimizing the mode 
% overlap.
goo.addComponent(component.lens(1.25,.75,'lens1'));
goo.addComponent(component.lens(1.75,3.25,'lens2'));

goo.seedWaist(.2e-3,0); % define input beam
goo.targetWaist(.4e-3,5); % define target beam

zdomain = -1:.01:6;

figure(266)
hold on
orighandle=goo.plotBeamWidth(zdomain);
goo.plotComponents(zdomain)

[pathList,overlapList] = goo.chooseComponents(...
                'lens1',lensList,[0.5 3],...  % choose lens1 from the list,
                'lens2',lensList.duplicate,[3.5 4],... %duplicate the list, this allows
                ...                                    %  the same component to be chosen more than once
                'target',[4.5,6]... % we can also allow the target waist position to vary while optimizing the overlap
                ,'-vt',.25); % set the minimum initial overlap to 0.25, if a combination of components
                             % has an overlap less than this, it will be skipped without trying to optimize the lens positions
           
% note about duplicating the list:
% If you have a box of lenses, such that you can't use the same lens twice,
% you can pass the same list to the function and it will make sure that
% each lens is used only once. If you can order as many lenses as you want,
% then duplicate the list, which makes an array of new component objects which
% are not linked to the originals.

% now let's cut out all solutions with modematching < 99%

pathList = pathList(overlapList >= 0.99);

% make an array with the combined position sensitivity of the components

sensitivityList = pathList.positionSensitivity;

% now sort in increasing sensitivity

[sensitivityList,sortIndex] = sort(sensitivityList);
pathList = pathList(sortIndex);

% plot the best solution
newhandle=pathList(1).plotBeamWidth(zdomain,'r');
pathList(1).plotComponents(zdomain,'r*')
pathList(1).plotBeams(zdomain,'k')

hold off
legend([orighandle newhandle],'Original Beam Path','Optimized Beam Path')
ylabel('Beam Width (m)')
xlabel('Propagation axis (m)')


% print the component list to the command window
disp(' ')
disp(' Optimized Path Component List:')
display(pathList(1).components)