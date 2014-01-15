% ---------- Example script for using a la mode mode matching utilities ------------

% create a new beam path object
PSLypath = beamPath; 

% add components to the beam path 
PSLypath.addComponent(component.lens(.3559,.37,'SL1'));
%                  lens syntax:     (focal length,z position,string label)
PSLypath.addComponent(component.lens(.175,.84,'SL2'));
PSLypath.addComponent(component.lens(.2542,.1277,'SLnew'));
PSLypath.addComponent(component.flatMirror((46-24+2)*0.0254,'M3'))
%                  flat mirror:           (z position,string label)
PSLypath.addComponent(component.flatMirror((56-24+2)*0.0254,'M4'))
% flat mirrors don't change modematching but they let you know if you're going to be putting
% stuff on top of eachother.

% The other useful component part is curved mirror, to make a curved mirror:
% component.curvedMirror(radius of curvature,z position,label)

% define "input beam" but it doesn't have to be at the input, it can be anywhere in the beam path
PSLypath.seedWaist(132e-6,-.137);
% seedWaist syntax:(waist width,z position)

% define the beam you are trying to match into, the target.
PSLypath.targetWaist(371e-6,(52+10-24+2)*.0254);
% targetWaist syntax:(waist width,z position)

% slide components to optimize mode overlap.
PSLypath = PSLypath.optimizePath('SL1',[5 39]*0.0254,'SL2',[.8 1]);
% optimizePath syntax           (component name,[(lower bound) (upper bound)],another component name,...)
% you can choose to optimize as many components as you would wish. Result is sensitive to initial conditions
% which are defined by the z position of the components before running optimize path.
% You can make it unbounded on either or both sides by using inf.
% if a component is not named, it will stay put.

%% after y path optimized

% duplicate the optimized beampath in order to work with the components exclusive to the x path
PSLxpath = PSLypath.duplicate;
% If you just did PSLxpath = PSLypath; you would just have two names for the same object, changing one would
% change the other. (think pointers)

%the x path has a different starting waist than the y path
PSLxpath.seedWaist(118e-6,-.066);

% add a cylindrical lens
PSLxpath.addComponent(component.lens(.752,.3,'CL1'));

% optimize the position of the cylindrical lens
PSLxpath = PSLxpath.optimizePath('CL1',[.25 .29]);

% the targetOverlap method calculates the mode overlap assuming it's doing an x and y integral,
% because these are actually the two dimensions of the same beam we square root and multiply them together.
modematch = sqrt(PSLxpath.targetOverlap*PSLypath.targetOverlap);
disp(['modematching = ',num2str(modematch)])

%% plot

% define plotting domain
zdomain = -.2:.01:42*0.0254;

figure(661)
subplot(2,1,1)
hold on % right now all the plot commands act like the matlab plot command and will overwrite the existing
        % figure unless you turn hold on

% The plot commands actually plots two traces, the top and bottom of the beam.
% The output of the plot commands returns the plot handle of the top so when we make the legend we don't 
% have to put a label on the top and bottom of the beam.
yplot = PSLypath.plotBeamWidth(zdomain,'b');
xplot = PSLxpath.plotBeamWidth(zdomain,'r');
PSLxpath.plotComponents(zdomain,0,'r*');
axis tight

legend([yplot xplot],'Y','X') % if we didn't use handles we would need to do 
                              %    legend('Y top','Y bottom','X top','X bottom') or something.
ylabel('Beam width (m)')
grid on
hold off

subplot(2,1,2)
hold on
PSLypath.plotGouyPhase(zdomain,'wrap','b');
PSLxpath.plotGouyPhase(zdomain,'wrap','r');
PSLxpath.plotComponents(zdomain,0,'r*');
axis tight
grid on
hold off

ylabel('Gouy Phase (degrees)')
xlabel('axial distance from MOPA aperture (m)')
