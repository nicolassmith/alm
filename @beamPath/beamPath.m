classdef beamPath < handle
    % -- beamPath --
    % A beam path object consists of a few things:
    %
    %   * A seed beam. This defines the "input beam" of your system.
    %           When calculating the beam properties somewhere in a beam
    %           path, this is the origin of the beam which is propagated 
    %           around the system. It is defined by a beamq object which is
    %           stored in the beamPath.seedq property, also the location
    %           of the seed beam in the beam path is defined by the 
    %           beamPath.seedz property. This beam can be defined anywhere
    %           in the beam path and need not be at the "beginning."
    %
    %   * A target beam. This defines the "output beam" of your system.
    %           This is the beam which you would like to match into.
    %           It is defined similarly to the seed beam by the properties
    %           beamPath.targetq and beamPath.targetz. It is also used by 
    %           the various optimization routines when trying to maximize 
    %           the mode overlap with the seed beam.
    %
    %   * Components. These are the various optical components in the beam
    %           path which affect the propagation of beams in the path.
    %           It is a vector array of objects in the component class.
    %           Each component has a certain transfer matrix and position,
    %           these properties are stored in the component class.
    %           The beamPath class stores the array of components in the 
    %           beamPath.components array. The ordering of the elements in 
    %           this array is always by increasing z position. This means
    %           that the array index of components in the component array
    %           may change if the z positions are changed. The component
    %           class also has a 'label' property, which allows a component
    %           to be indexed unambiguously. See beamPath.component below.
    %
    % The beamPath class also has methods which facilitate in calculating 
    %     beam propagation and laying out your beam path.
    %
    % Methods:
    %  - These methods are used for making beamPath objects:
    %     beamPath - the beam path constructor, to create a blank beam path
    %         object.
    %     <a href="matlab:help beamPath.duplicate">duplicate</a> - creates a copy of the calling beamPath object.
    %     <a href="matlab:help beamPath.branchPath">branchPath</a> - creates a new beamPath object identical to the calling
    %         object, but has a seed beam defined as a propagated beam from the 
    %         original beamPath object.
    %  - These methods allow manipulation of the components in a beam path:
    %     <a href="matlab:help beamPath.addComponent">addComponent</a> - adds a component to the beam path.
    %     <a href="matlab:help beamPath.deleteComponent">deleteComponent</a> - deletes a component from the beam path.
    %     <a href="matlab:help beamPath.moveComponent">moveComponent</a> - changes the position of a component in the beam path.
    %     <a href="matlab:help beamPath.replaceComponent">replaceComponent</a> - replaces a component in the beam path with another
    %         component.
    %     <a href="matlab:help beamPath.component">component</a> - returns the requested component object, usually referenced
    %         by label.
    %     <a href="matlab:help beamPath.findComponentIndex">findComponentIndex</a> - returns the list index of the requested component.
    %  - These methods are for calculating beam propagation in a beam path:
    %     <a href="matlab:help beamPath.seedWaist">seedWaist</a> - defines the seed beam as a waist of a given size at a given
    %         location.
    %     <a href="matlab:help beamPath.seedWaistR">seedWaistR</a> - defines the seed beam as a waist of a given size and radius at a given
    %         location.
    %     <a href="matlab:help beamPath.targetWaist">targetWaist</a> - defines the target beam as a waist of a given size at a given
    %         location.
    %     <a href="matlab:help beamPath.targetWaistR">targetWaistR</a> - defines the target beam as a waist of a given size and radius at a given
    %         location.
    %     <a href="matlab:help beamPath.qPropagate">qPropagate</a> - returns a beamq object at a desired position given an 
    %         input beam, usually the seed beam.
    %     <a href="matlab:help beamPath.getTransferMatrix">getTransferMatrix</a> - returns the ABCD transfer matrix between two points 
    %         in the beam path.
    %     <a href="matlab:help beamPath.gouyPhase">gouyPhase</a> - calculates the accumulated gouy phase from the input beam
    %         (usually the seed beam) to a given location.
    %     <a href="matlab:help beamPath.targetOverlap">targetOverlap</a> - calculates mode overlap of the seed and target beams.
    %  - These methods are for analyzing properties of a beam path:
    %     <a href="matlab:help beamPath.gouySeparation">gouySeparation</a> - calculates the accumulated gouy phase between two
    %         components in the beam path.
    %     <a href="matlab:help beamPath.eigenMode">eigenMode</a> - calculates the eigenmode of a part of the beampath.
    %     <a href="matlab:help beamPath.angleSensitivity">angleSensitivity</a> - calculates the sensitivity of components to misalignments.
    %     <a href="matlab:help beamPath.lateralSensitivity">lateralSensitivity</a> - calculates the sensitivity of components to lateral motion.
    %     <a href="matlab:help beamPath.positionSensitivity">positionSensitivity</a> - calculates the sensitivity of components to axial motion 
    %         (motion along the beam axis).
    %     <a href="matlab:help beamPath.sensitivitySummary">sensitivitySummary</a> - displays a summary of the sensitivity of components to 
    %         angle, lateral position, and axial position motion.
    %     <a href="matlab:help beamPath.getWaists">getWaists</a> - returns an array of the locations of beam waists.
    %  - These methods are for making plot representations of a beam path:
    %     <a href="matlab:help beamPath.plotBeamWidth">plotBeamWidth</a> - plots the beam width calculated from the seed beam.
    %     <a href="matlab:help beamPath.plotGouyPhase">plotGouyPhase</a> - plots the accumulated gouy phase from the seed beam.
    %     <a href="matlab:help beamPath.plotComponents">plotComponents</a> - plots the components in a beam path.
    %     <a href="matlab:help beamPath.plotBeams">plotBeams</a> - plots information about the seed and target beams.
    %     <a href="matlab:help beamPath.plotSummary">plotSummary</a> - makes a plot with a summary of the beam path.
    %  - These methods are for manipulating the components in a path to optimize
    %        the mode overlap to the target beam.
    %     <a href="matlab:help beamPath.optimizePath">optimizePath</a> - adjusts the positions of components in the path in order
    %         to maximize the mode overlap.
    %     <a href="matlab:help beamPath.chooseComponents">chooseComponents</a> - chooses components from lists supplied in the input 
    %         arguments in order to maximize mode overlap.
    %  - These methods are for fitting the seed beam to measurements of beam width.
    %     <a href="matlab:help beamPath.fitBeamWidth">fitBeamWidth</a> - returns a beamPath object identical to the calling object, 
    %         but with a seed beam which is chosen to best match the supplied beam
    %         width values.
    properties (Dependent)
        components;
    end
    properties
        seedq;
        seedz;
        targetq;
        targetz;
    end
    properties (Hidden, SetAccess = private)
        components_raw; %this is the component list that can be accessed by the sortComponents method
    end
    methods (Hidden)
        function losses = lossFunc(pathobj,zVec,varargin)
            path2 = pathobj.duplicate;
            path2.batchMove(zVec);

            losses = 1-path2.targetOverlap;
            path2.delete;
        end
        function cost = applyCostFunc(pathobj,costFunc,zVec,varargin)
            path2 = pathobj.duplicate;
            path2.batchMove(zVec);
            
            cost = costFunc(path2);
            path2.delete;
        end
        function batchMove(pathobj,moveVec)
            comps = pathobj.components; % this component array will not reorder itself and get indicies mixed up
            
            compCount = length(comps);
            moveCount = length(moveVec);
            
            if moveCount>compCount
                pathobj.targetz = moveVec(compCount+1);
            end
            for j = 1:compCount;
                if comps(j).z == moveVec(j)
                    continue
                end
                comps(j).z = moveVec(j);
            end
        end
        function sortComponents(pathobj)
            complist = [pathobj.components_raw];

            zlist = [complist.z];
            [zsorted,zindex] = sort(zlist);
%             if length(complist) < 1 || all(complist == complist(zindex))
%                 return
%             end
            pathobj.components = complist(zindex);
        end
        function error = beamFitError(pathobj,waistParams,zPred,widthPred,lambda)
            
            pathdup = pathobj.duplicate;
            pathdup.seedq=beamq.beamWaistAndZ(waistParams(1),waistParams(2),lambda);
            qout = pathdup.qPropagate(zPred);
            
            widthPred = reshape(widthPred,1,numel(widthPred));

            width = [qout.beamWidth];
            
            error = sum((widthPred-width).^2);
        end
    end
    methods
        % methods for creating beampaths
        function pathObjOut = beamPath
            pathObjOut.seedq = beamq;
            pathObjOut.seedz = [];
            pathObjOut.targetq = beamq;
            pathObjOut.targetz = [];
            
            pathObjOut.components_raw = component;            
            pathObjOut.components_raw(1)=[]; %empty the component list
        end
        function pathObjOut = duplicate(pathobj)
            % -- beamPath.duplicate --
            % Creates a new beampath with the same properties as the original.
            % Example:
            % path1copy = path1.duplicate;
            pathObjOut = beamPath;
            pathObjOut.components = pathobj.components.duplicate;
            pathObjOut.seedq = pathobj.seedq.duplicate;
            pathObjOut.seedz = pathobj.seedz;
            pathObjOut.targetq = pathobj.targetq.duplicate;
            pathObjOut.targetz = pathobj.targetz;
        end
        function pathObjOut = branchPath(pathobj,zlink)
            % -- beamPath.branchPath --
            % create a new beamPath object which is identical to the calling object,
            % but has a seed beam which is calculated from the seed beam in the first
            % object but located in another place.
            % 
            % this can be useful when you want to make a new path which has the same
            % beam as the original path at a given location, and then allows you to alter
            % the components in the new path without changing the beam at the location you chose.
            %
            % syntax: path2 = path1.branchPath(zlink)
            % zlink is the position of the beam you would like to be the seed beam of the new
            % path.
            qlink = pathobj.qPropagate(zlink);
            pathObjOut = pathobj.duplicate;

            pathObjOut.seedq = qlink;
            pathObjOut.seedz = zlink;
        end
        % property access methods
        function pathobj = set.components(pathobj,comps) %#ok<MCHV2>
            pathobj.components_raw = comps;
        end
        function compArray = get.components(pathobj)
            pathobj.sortComponents;
            compArray = pathobj.components_raw;
        end
        % overridden builtin methods
        function display(pathobj)
            if length(pathobj)>1
                sizeobj = size(pathobj);
                disp(' ')
                disp(['  ' num2str(sizeobj(1)) 'x' num2str(sizeobj(2)) ' <a href="matlab:help beamPath">beamPath</a> object: ' inputname(1)])
                disp(' ')
                return
            end
            
            disp(' ')
            disp(['  <a href="matlab:help beamPath">beamPath</a> object: ' inputname(1)])
            disp(' ')
            disp(['  Contains ' num2str(length(pathobj.components)) ' components.'])
            disp(['      <a href="matlab:' inputname(1) '.components">' inputname(1) '.components</a>'])
            disp(' ')
            disp('  Beams:')

            seedzString = ['          Location:    <a href="matlab:' inputname(1) '.seedz">' inputname(1) '.seedz</a>'];
            seedqString = ['          q parameter: <a href="matlab:' inputname(1) '.seedq">' inputname(1) '.seedq</a>'];
            
            if ~isempty(pathobj.seedq.q)
                disp('      Seed beam;')
                seedzString = [seedzString ' = ' num2str(pathobj.seedz) ' m.'];
                seedqString = [seedqString ' = ' num2str(pathobj.seedq.q) ' m.'];
            else
                disp('      Seed beam not defined;')
            end
            
            disp(seedzString)
            disp(seedqString)


            targetzString = ['          Location:    <a href="matlab:' inputname(1) '.targetz">' inputname(1) '.targetz</a>'];
            targetqString = ['          q parameter: <a href="matlab:' inputname(1) '.targetq">' inputname(1) '.targetq</a>'];
            
            if ~isempty(pathobj.targetq.q)
                disp('      Target beam;')
                targetzString = [targetzString ' = ' num2str(pathobj.targetz) ' m.'];
                targetqString = [targetqString ' = ' num2str(pathobj.targetq.q) ' m.'];
            else
                disp('      Target beam not defined;')
            end
            
            disp(targetzString)
            disp(targetqString)

            if ~isempty(pathobj.seedq.q) && ~isempty(pathobj.targetq.q)
                disp(' ')
                disp(['  Mode overlap with target beam: ' num2str(pathobj.targetOverlap) '.'])
            end
            disp(' ')
        end
        % these methods allow manipulation of the components in a beam path
        function addComponent(pathobj,newComponent) %idea - add component by name?
            % -- beamPath.addComponent --
            % add a component object to a beampath.
            % Example:
            % mylens = component.lens(2,0,'mylens');
            % path1.addComponent(mylens);
            
            listSize = size(newComponent);
            if all(listSize~=1)
                error('When adding components, make list a Nx1 array')
            end
            
            if listSize(2)>listSize(1)
                newComponent = newComponent.';
            end
            
            pathobj.components = [pathobj.components;newComponent];
        end
        function deleteComponent(pathobj,delLabel)
            % -- beamPath.deleteComponent --
            % remove a component object from a beampath
            % Example:
            % path1.deleteComponent('mylens');
            delIndex = pathobj.findComponentIndex(delLabel);
            
            pathobj.components = pathobj.components(1:end~=delIndex);
            
        end
        function moveComponent(pathobj,componentLabel,displacement,isabsolute) 
            % -- beamPath.moveComponent --
            % Change the z parameter of a component in a beam path by a given displacement.
            % Example:
            % path1.moveComponent('lens1',0.5)
            % This will move 'lens1' 0.5m in the positive z direction.
            % To move a component to an absolute position use the extra argument 'absolute.'
            % Example:
            % path1.moveComponent('lens1',2.5,'absolute')
            % This will move 'lens1' to the position z = 2.5m
            if nargin < 4
                isabsolute = '';
            end
            componentIndex=pathobj.findComponentIndex(componentLabel);
            
            componentToMove = pathobj.components(componentIndex);
            zstart = componentToMove.z;
            
            if strcmp(isabsolute,'absolute')
                displacement = displacement - zstart;
            end
            
            componentToMove.z = zstart + displacement;
        end
        function replaceComponent(pathobj,componentLabel,newComponent)
            % -- beamPath.replaceComponent --
            % Remove a component and replace it with another, the new component inherits the
            % z position and label of the previous component being removed.
            % Example:
            % newlens = component.lens(1)
            % path1.replaceComponent('lens1',newlens)
            % This removes the old 'lens1' component and adds the new component at the 
            % position where 'lens1' was. The name 'lens1' is inherited by the new component.
            componentIndex=pathobj.findComponentIndex(componentLabel);
            
            newComponent.z = pathobj.components(componentIndex).z;
            newComponent.label = pathobj.components(componentIndex).label;
            
            pathobj.deleteComponent(componentIndex);
            pathobj.addComponent(newComponent);
        end
        function returnComponent = component(pathobj,componentLabel)
            % -- beamPath.component --
            % Allows access to component in a beam path by use of the component label.
            % Example:
            % To find out the z position of 'mylens' which might be the 
            % third lens in the component list, one could either do:
            % 
            % path1.components(3).z
            % or
            % path1.component('mylens').z  %%% <-- notice the use of the 
            %                                      singular word 'component' in this case
            % which would yield the same result.
            % However, the component index may change if new components are 
            % added, the indexing is always in order of increasing z position,
            % so this method allows one to access the desired component
            % unambiguously.
            componentIndex=pathobj.findComponentIndex(componentLabel);
            
            returnComponent=pathobj.components(componentIndex);
        end
        function indexout = findComponentIndex(pathobj,componentLabel)
            % -- beamPath.findComponentIndex --
            % Returns the index of a component in the component list
            % Example:
            % Suppose 'goodlens' is the fourth component in the component list 
            % of beampath path1.
            % path1.findComponentIndex('goodlens')
            % This statement would return the number 4.
            if isnumeric(componentLabel)
                indexout = componentLabel;
                componentCount = length(pathobj.components);
                if indexout>componentCount
                    error(['Sorry, the component index (' num2str(indexout)...
                        ') is greater than the number of components. ('...
                        num2str(componentCount) ')']);
                end
                return
            end
            
            pathLabels = {pathobj.components.label};
            
            listCompare = find(strcmpi(componentLabel,pathLabels),2);
            
            if isempty(listCompare)
                error(['Couldn''t find component with label: ' componentLabel])
            elseif length(listCompare) > 1
                warning('beamPath:ComponentLabelNotUnique',...
                    ['Found multiple components with label: '...
                    componentLabel '. Using only the first instance.'])
            end

            indexout = listCompare(1);
        end
        % these methods are for calculating beam propagation given a beam path
        function seedWaist(pathobj,waistSize,waistPos,lambda)
            % -- beamPath.seedWaist --
            % Sets the seed beam of the beam path to a waist of given size and position.
            % Example:
            % path1.seedWaist(w0,position,lambda)
            % This will set the seed beam to be a waist at z = position
            % with beam width = w0. The wavelength is set by lambda, if omitted,
            % the default value is 1064nm.
            if nargin<4
                lambda = 1064e-9;
            end
            
            pathobj.seedq = beamq.beamWaistAndZ(waistSize,0,lambda);
            pathobj.seedz = waistPos;
            
        end
        function seedWaistR(pathobj,waistSize, waistR, waistPos, lambda)
            % -- beamPath.seedWaistR --
            % Sets the seed beam of the beam path with a waist and radius at the position.
            % Example:
            % path1.seedWaistR(w0,R0,position,lambda)
            % This will set the seed beam to be a beam waist w0  and radius of curvature R0
            % at z = position. The wavelength is set by lambda, if omitted,
            % the default value is 1064nm.
            if nargin<5
                lambda = 1064e-9;
            end
            
            pathobj.seedq = beamq.beamWaistAndR(waistSize,waistR,lambda);
            pathobj.seedz = waistPos;
            
        end
        function targetWaist(pathobj,waistSize,waistPos,lambda)
            % -- beamPath.targetWaist --
            % Sets the target beam of the beam path to a waist of given size and position.
            % Example:
            % path1.targetWaist(w0,position,lambda)
            % This will set the target beam to be a waist at z = position
            % with beam width = w0. The wavelength is set by lambda, if omitted,
            % the default value is 1064nm.
            if nargin<4
                lambda = 1064e-9;
            end
            
            pathobj.targetq = beamq.beamWaistAndZ(waistSize,0,lambda);
            pathobj.targetz = waistPos;
        end
        function targetWaistR(pathobj,waistSize,waistR, waistPos,lambda)
            % -- beamPath.targetWaistR --
            % Sets the target beam of the beam path to a waist of given size and position.
            % Example:
            % path1.targetWaist(w0,R0,position,lambda)
            % This will set the target beam to be a waist at z = position
            % with beam waist size = w0 and radius R0. The wavelength is set by lambda, if omitted,
            % the default value is 1064nm.
            if nargin<5
                lambda = 1064e-9;
            end
            
            pathobj.targetq = beamq.beamWaistAndR(waistSize,waistR,lambda);
            pathobj.targetz = waistPos;
        end
        function qout = qPropagate(pathobj,zdomain,qin,zqin)
            % -- beamPath.qPropagate --
            % Propagate a beam q parameter through your beam path.
            % Example:
            % qout = path1.qPropagate(z,qin,zqin)
            % Notes:
            % If z is a vector of postitions, qout will be a vector of beamq objects.
            % If qin and zqin are omitted, the seedq and seedz properties of 
            %    the path object will be used.
            if nargin<3
                qin = pathobj.seedq;
                zqin = pathobj.seedz;
            end
            
            if isempty(qin.q) || isempty(zqin)
                error('Initial beam is undefined. Either define a seed beam, or use more input arguments')
            end
            
            zlength = length(zdomain);
            qout(zlength,1) = beamq;
            
            for j = 1:zlength;
                transferM = pathobj.getTransferMatrix(zqin,zdomain(j));
                qout(j) = qin.transform(transferM);
                qin = qout(j);
                zqin = zdomain(j);
            end
        end
        function Mout = getTransferMatrix(pathobj,z1,z2) %needs better support of components with length
            % -- beamPath.getTransferMatrix --
            % Calculate and return the transfer (ABCD) matrix through the 
            % beam path between two points.
            % Example:
            % path1.getTransferMatrix(z1,z2)
            % This returns the ABCD matrix between z1 and z2.
                        
            % handle the case where z1>z2, flip them then invert the matrix at the end
            invertPower = 1;
            if z1 > z2
                z1new = z2;
                z2 = z1;
                z1 = z1new;
                invertPower = -1;
            end
            
            compCount = length(pathobj.components);
            if compCount > 0
                zlist = [pathobj.components.z];
                % the endpoints are just to the + side of a component that happens to be there
                zindex = z1 < zlist & zlist <= z2;
            end
            
            if compCount > 0 && any(zindex) % short circuiting is used here
                comps = pathobj.components;
                
                comps = comps(zindex);
                % make a new beampath which only has the components between z1 and z2
                compPath = beamPath;
                compPath.components = comps;
                
                % zstarts and zends define the start and enpoints of the propagators that will be added
                zstarts = [z1,[comps.z]];
                zends = [[comps.z],z2];
                
                %this block here is to deal with components with a 'length' parameter
                params = {comps.parameters};
                for j=1:length(params)
                    paramNames = fieldnames(params{j});
                    
                    listCompare = find(strcmpi('length',paramNames),1);
                    if ~isempty(listCompare)
                        % this line makes the next propagator start after the end of the component
                        zstarts(j+1) = zstarts(j+1) + params{j}.length;
                    end
                end
                dz = zends - zstarts;

                for j = 1:length(dz)
                    compPath.addComponent(component.propagator(dz(j),zstarts(j))) % this is a slow way to do it, needs optimization
                end

                fullPathComponent = compPath.components.combine();
                Mout = fullPathComponent.M;
                compPath.delete; %clear the temporary path object
            else
                % if there are no components in the path we are looking at,
                % we know what the matrix is. This is not necessary, but is done for speed.
                
                Mout = [1 z2-z1;0 1];
            end
            %invert if it needs it
            Mout = Mout^(invertPower);
        end
        function [gPhase,qout] = gouyPhase(pathobj,zdomain,qin,zqin)
            % -- beamPath.gouyPhase --
            % Returns the accumulated gouy phase (in degrees) from zqin to zdomain, given
            % initial beam qin. If no qin, zqin are given, then the seed beam is used
            % Example:
            % [gPhase,qout] = path1.gouyPhase(zdomain,qin,zqin)
            % 
            % it also returns a beamq array (qout) if you don't want to do the propogation
            % calculation twice.
            % (reference: Erden, Ozaktas [1997])
            % see also beamPath.qPropagate
            
            if nargin<3
                qin = pathobj.seedq;
                zqin = pathobj.seedz;
            end
            
            if isempty(qin.q) || isempty(zqin)
                error('Initial beam is undefined. Either define a seed beam, or use more input arguments.')
            end
            
            zlength = length(zdomain);
            
            qout(zlength,1) = beamq;
            gPhase = zeros(zlength,1);
            gAccum = 0;
            
            for j = 1:zlength;
                transferM = pathobj.getTransferMatrix(zqin,zdomain(j));
                A = transferM(1,1);
                B = transferM(1,2);
                win = qin.beamWidth;
                rin = qin.radiusOfCurvature;
                lambda = qin.lambda;
                
                tanPhiA = lambda*B; % see reference
                tanPhiB = (A + B/rin )*pi*win^2;
                
                gAccum = gAccum + 180/pi * atan2(tanPhiA,tanPhiB);
                gPhase(j) = gAccum;
                
                qout(j) = qin.transform(transferM);
                qin = qout(j);
                zqin = zdomain(j);
            end
        end
        function overlapFrac = targetOverlap(pathobj)
            % -- beamPath.targetOverlap --
            % This calculates the result of propogating the 'seed' beam through
            % the beampath to the 'target.' Here the mode overlap is calculated and
            % returned, where 1.0 is perfect overlap.
            % Example:
            % path1.targetOverlap()
            
            % handle array of beamPaths
            if numel(pathobj)>1
                arraySize = size(pathobj);
                overlapFrac = zeros(arraySize);
                for jj = 1:arraySize(1)
                    for kk = 1:arraySize(2)
                        overlapFrac(jj,kk) = pathobj(jj,kk).targetOverlap;
                    end
                end
                return
            end
            
            if isempty(pathobj.targetz) || isempty(pathobj.targetq.q)
                error('Can''t do mode overlap, no target beam is defined.')
            end
            
            ztarget = pathobj.targetz;
            qAtTarget = pathobj.qPropagate(ztarget);
            overlapFrac = overlap(pathobj.targetq,qAtTarget);
        end
        % these are for analyzing properties of the beam path.
        function gPhase = gouySeperation(varargin)
            %compatibility with previous versions
            warning(['The function gouySeperation is misspelled, ' ...
                     'correct spelling: gouySeparation.']);
            gPhase = gouySeparation(varargin{:});
        end
        function gPhase = gouySeparation(pathobj,compLabel1,compLabel2,dontWrap)
            % -- beamPath.gouySeparation --
            % Returns the accumulated gouy phase between two components in the beam path.
            % syntax: gouyPhase = path1.gouySeparation('label1','label2')
            % 
            % By default this function will return values in the range -90 to +90, 
            % to have the function return the gouy phase unwrapped, give the string
            % 'nowrap' as the third argument
            if nargin<4
                dontWrap = '';
            end
            pathArraySize = size(pathobj);
            gPhase = zeros(pathArraySize);
            for jj = 1:pathArraySize(1)
                for kk = 1:pathArraySize(2)
                    z1 = pathobj(jj,kk).component(compLabel1).z;
                    z2 = pathobj(jj,kk).component(compLabel2).z;

                    g1 = pathobj(jj,kk).gouyPhase(z1);
                    g2 = pathobj(jj,kk).gouyPhase(z2);
                    
                    gPhase(jj,kk) = g2-g1;
                    
                    if ~strcmpi(dontWrap,'nowrap')
                        gPhase(jj,kk) = mod(gPhase(jj,kk)+90,180)-90;
                    end    
                end
            end
        end
        function qObjOut = eigenMode(pathobj,z1,z2,lambda)
            %  -- beamPath.eigenMode --
            % Returns a beamq object that represents the eingenmode of propagation
            % through a part of a beam path. 
            %
            % syntax: eigen = path1.eigenMode(z1,z2,lambda)
            % eigen is the returned beamq object. z1 is the starting point, and
            % z2 is the end point. Thus if a beam path has a q value given by
            % eigen at z1, it will have the same q value at z2. If lambda is not
            % specified, it is assumed to be 1064nm.
            % 
            
            if nargin<4
                lambda=1064e-9;
            end
            
            % get ABCD matrix
            m = pathobj.getTransferMatrix(z1,z2);
            
            A=m(1,1);
            B=m(1,2);
            C=m(2,1);
            D=m(2,2);
            if (A-D)^2 + 4*B*C <0
                qout=(A-D)/2/C + sqrt((A-D)^2 +4*B*C)/2 /abs(C);
                % C could be negative and we'd like to choose a physical
                % solution where the imaginary part is positive.
            else
                error('No stable eigenmode inside this resonator.');
            end
            
            qObjOut=beamq(qout,lambda);
        end
        function sensitivity = angleSensitivity(pathobj,varargin)
            % -- beamPath.angleSensitivity -
            % syntax: sensitivity = path1.angleSensitivity('mirror 1','mirror 2')
            % (if no arguent is given, the sensitivity is calculated
            % for all componenets in the beam path)
            % Returns an array of values representing the
            % sensitivity of the beampath to angular motion of the
            % components. The value returned is the matrix element
            % which connects the TEM00 mode to the 01 mode, divided
            % by the misalignment. So for a component having a
            % misalignment matrix element of A, the power remaining
            % in the TEM00 mode is P=1-(A*Theta)^2 for misalignment
            % Theta (in radians).
            % 
            % options:
            % To calculate only the combined sensitivity of the entire
            % path, by adding all sensitivities in quadrature, use the
            % argument '-c'.
            % Note: the only types of objects for which this returns a
            % non-zero value are curved mirrors and flat mirrors.
            
            arraySize = size(pathobj);
            
            if numel(pathobj)>1  % deals with arrays of beampaths
                sensitivity = zeros(arraySize);
                for jj = 1:arraySize(1)
                    for kk = 1:arraySize(2)
                        sensitivity(jj,kk) = pathobj(jj,kk).angleSensitivity('-c');
                    end
                end
                return
            end
            
            if isempty(pathobj.targetz) || isempty(pathobj.targetq.q)
                error('Can''t do sensitivity, no target beam is defined.')
            end
            
            lvargin = length(varargin);
            
            flagIndex = find(strncmp(varargin,'-',1),2); % find the option flag
            if length(flagIndex)>1
                error('Please only use the option flag "-" once in arguments.')
            end
            
            combined = 0;
            if flagIndex
                if any(varargin{flagIndex}=='c') || any(varargin{flagIndex}=='C') % identify combined flag
                    combined = 1; 
                end
                varargin = {varargin{1:end~=flagIndex}}; % remove options from the arguments
                lvargin = length(varargin);
            end
            
            if lvargin == 0
                complist = pathobj.components;
            else
                complist(lvargin,1) = component;
                for jj = 1:lvargin
                    complist(jj) = pathobj.component(varargin{jj});
                end
            end
            
            numComps = length(complist);
            
            sensitivity = zeros(1,numComps);
            
            for kk = 1:numComps
                if strcmp(complist(kk).type,'curved mirror') || strcmp(complist(kk).type,'flat mirror')
                    %calculate
                    S = pathobj.getTransferMatrix(complist(kk).z,pathobj.targetz);
                    
                    w0 = pathobj.targetq.waistSize;
                    divAngle = pathobj.targetq.divergenceAngle;
                    
                    % the factor of 2 is because the beam angle is twice the mirror angle
                    sensitivity(kk) = 2* sqrt( (S(1,2)/w0)^2 + (S(2,2)/divAngle)^2 );
                else
                    sensitivity(kk) = 0;
                end
            end
            
            if combined
                sensitivity = sqrt(sum(sensitivity.^2));
            end
        end
        function sensitivity = lateralSensitivity(pathobj,varargin)
            % -- beamPath.lateralSensitivity --
            % syntax: sensitivity = path1.lateralSensitivity('mirror 1','mirror 2')
            % (if no arguent is given, the sensitivity is calculated
            % for all componenets in the beam path)
            % Returns an array of values representing the
            % sensitivity of the beampath to lateral motion of the
            % components. The value returned is the matrix element
            % which connects the TEM00 mode to the 01 mode, divided
            % by the lateral displacement. So for a component having a
            % lateral displacement matrix element of A, the power remaining
            % in the TEM00 mode is P=1-(A*DeltaY)^2 for lateral displacement
            % DeltaY (in meters).
            % 
            % options:
            % To calculate only the combined sensitivity of the entire
            % path, by adding all sensitivities in quadrature, use the
            % argument '-c'.
            %
            % Note: the only objects which return a non-zero value are
            % curved mirrors and lenses.
            
            arraySize = size(pathobj);
            
            if numel(pathobj)>1 % deals with arrays of beampaths
                sensitivity = zeros(arraySize);
                for jj = 1:arraySize(1)
                    for kk = 1:arraySize(2)
                        sensitivity(jj,kk) = pathobj(jj,kk).lateralSensitivity('-c');
                    end
                end
                return
            end
            
            if isempty(pathobj.targetz) || isempty(pathobj.targetq.q)
                error('Can''t do sensitivity, no target beam is defined.')
            end
            
            lvargin = length(varargin);
            
            flagIndex = find(strncmp(varargin,'-',1),2); % find the option flag
            if length(flagIndex)>1
                error('Please only use the option flag "-" once in arguments.')
            end
            
            combined = 0;
            if flagIndex
                if any(varargin{flagIndex}=='c') || any(varargin{flagIndex}=='C') % identify combined flag
                    combined = 1; 
                end
                varargin = {varargin{1:end~=flagIndex}}; % remove options from the arguments
                lvargin = length(varargin);
            end
            
            if lvargin == 0
                complist = pathobj.components;
            else
                lvargin = length(varargin);
                complist(lvargin,1) = component;
                for jj = 1:lvargin
                    complist(jj) = pathobj.component(varargin{jj});
                end
            end
            
            numComps = length(complist);
            
            sensitivity = zeros(1,numComps);
            
            for kk = 1:numComps
                if strcmp(complist(kk).type,'curved mirror')
                    %calculate
                    S = pathobj.getTransferMatrix(complist(kk).z,pathobj.targetz);
                    
                    w0 = pathobj.targetq.waistSize;
                    divAngle = pathobj.targetq.divergenceAngle;
                    f = complist(kk).parameters.ROC / 2;
                    
                    sensitivity(kk) = sqrt( (S(1,2)/w0)^2 + (S(2,2)/divAngle)^2 )/abs(f);
                elseif strcmp(complist(kk).type,'lens')
                    %calculate
                    S = pathobj.getTransferMatrix(complist(kk).z,pathobj.targetz);
                    
                    w0 = pathobj.targetq.waistSize;
                    divAngle = pathobj.targetq.divergenceAngle;
                    f = complist(kk).parameters.focalLength;
                    
                    sensitivity(kk) = sqrt( (S(1,2)/w0)^2 + (S(2,2)/divAngle)^2 )/abs(f);
                else
                    sensitivity(kk) = 0;
                end
            end
            
            if combined
                sensitivity = sqrt(sum(sensitivity.^2));
            end
        end
        function sensitivity = positionSensitivity(pathobj,varargin)
            % -- position sensitivity --
            % syntax: sensitivity = path1.positionSensitivity('mirror 1','mirror 2')
            % (if no arguent is given, the sensitivity is calculated
            % for all componenets in the beam path)
            % Returns an array of values representing the
            % sensitivity of the beampath to axial position motion of the
            % components (motion along the z axis). 
            % The value returned is the matrix element
            % which connects the TEM00 mode to the TEM20+TEM02 mode, divided
            % by the axial displacement. So for a component having a
            % axial displacement matrix element of A, the power remaining
            % in the TEM00 mode is P=1-(A*DeltaZ)^2 for axial displacement
            % DeltaZ (in meters).
            % 
            % options:
            % To calculate only the combined sensitivity of the entire
            % path, by adding all sensitivities in quadrature, use the
            % argument '-c'.
            %
            % Note: the only objects which return a non-zero value are
            % curved mirrors, flat mirrors and lenses.
            
            arraySize = size(pathobj);
            
            if numel(pathobj)>1 % deals with arrays of beampaths
                sensitivity = zeros(arraySize);
                for jj = 1:arraySize(1)
                    for kk = 1:arraySize(2)
                        sensitivity(jj,kk) = pathobj(jj,kk).positionSensitivity('-c');
                    end
                end
                return
            end
            
            if isempty(pathobj.targetz) || isempty(pathobj.targetq.q)
                error('Can''t do sensitivity, no target beam is defined.')
            end
            
            lvargin = length(varargin);
            
            flagIndex = find(strncmp(varargin,'-',1),2); % find the option flag
            if length(flagIndex)>1
                error('Please only use the option flag "-" once in arguments.')
            end
            
            combined = 0;
            if flagIndex
                if any(varargin{flagIndex}=='c') || any(varargin{flagIndex}=='C') % identify combined flag
                    combined = 1; 
                end
                varargin = {varargin{1:end~=flagIndex}}; % remove options from the arguments
                lvargin = length(varargin);
            end
            
            if lvargin == 0
                complist = pathobj.components;
            else
                lvargin = length(varargin);
                complist(lvargin,1) = component;
                for jj = 1:lvargin
                    complist(jj) = pathobj.component(varargin{jj});
                end
            end
            
            numComps = length(complist);
            
            sensitivity = zeros(1,numComps);
            
            for kk = 1:numComps
                if strcmp(complist(kk).type,'curved mirror') || strcmp(complist(kk).type,'flat mirror')
                    %calculate
                    
                    fudge = 1e-9;
                    
                    qOut = pathobj.qPropagate(complist(kk).z,pathobj.targetq,pathobj.targetz);
                    qIn = pathobj.qPropagate(complist(kk).z-fudge,pathobj.targetq,pathobj.targetz);
                    
                    qIn = qIn.transform([1 fudge ; 0 1]);
                    
                    A = (qOut.q/qIn.q).^2;
                    
                    zRO = qOut.rayleighRange;
                    
                    alpha= (1i*.5*(1+real(A))./zRO+.5./zRO .*imag(A));
                    
                    sensitivity(kk) = abs(alpha);
                elseif strcmp(complist(kk).type,'lens')
                    %calculate
                    
                    fudge = 1e-9;
                    
                    qOut = pathobj.qPropagate(complist(kk).z,pathobj.targetq,pathobj.targetz);
                    qIn = pathobj.qPropagate(complist(kk).z-fudge,pathobj.targetq,pathobj.targetz);
                    
                    qIn = qIn.transform([1 fudge ; 0 1]);
                    
                    zRO = qOut.rayleighRange;
                    
                    f = complist(kk).parameters.focalLength;
                    
                    B = qOut.q/f * ( 1 + qOut.q/qIn.q );
                    
                    alpha= (1i*.5*(real(B))./zRO+.5./zRO .*imag(B));
                    
                    sensitivity(kk) = abs(alpha);
                else
                    sensitivity(kk) = 0;
                end
            end
            
            if combined
                sensitivity = sqrt(sum(sensitivity.^2));
            end
        end
        function sensitivitySummary(pathobj,varargin)
            % -- beamPath.sensitivitySummary --
            % syntax: path1.sensitivitySummary('mirror 1','mirror 2')
            % (if no arguent is given, the sensitivity is calculated
            % for all componenets in the beam path)
            %
            % Displays a summary of the beam path components
            % sensitivity to misalignments, lateral displacements and
            % axial displacements.
            
            if isempty(pathobj.targetz) || isempty(pathobj.targetq.q)
                error('Can''t do sensitivity, no target beam is defined.')
            end
            
            if nargin<2
                complist = pathobj.components;
            else
                lvargin = length(varargin);
                complist(lvargin,1) = component;
                for jj = 1:lvargin
                    complist(jj) = pathobj.component(varargin{jj});
                end
            end
            
            labelColumn = {'';...
                           'Ang. Sensitivity';...
                           'Lat. Sensitivity';...
                           'Pos. Sensitivity'}; %#ok<NASGU>
                       
            unitsColumn = {'';...
                           '(1/rad)';...
                           '(1/m)';...
                           '(1/m)'}; %#ok<NASGU>
            
            numComps = length(complist);
            compColumns = cell(4,numComps);
            
            for kk = 1:numComps
                compColumns{1,kk} = complist(kk).label;
                compColumns{2,kk} = num2str(pathobj.angleSensitivity(complist(kk).label));
                compColumns{3,kk} = num2str(pathobj.lateralSensitivity(complist(kk).label));
                compColumns{4,kk} = num2str(pathobj.positionSensitivity(complist(kk).label));
            end
            
            dispstring = evalc('disp([labelColumn,compColumns,unitsColumn])');
            dispstring = dispstring(dispstring~=''''); % remove quotes from output
            
            disp(' ');
            disp([' Sensitivity summary for beam path: ' inputname(1)])
            disp(' ');
            disp(dispstring);
            disp(' ');
        end
        function [zList,qList] = getWaists(pathobj,z1,z2)
            % -- beamPath.getWaists --
            % Returns an array of the position of all waists in the beam path
            % between two points.
            % syntax:
            % waistPositions = path1.getWaists(z1,z2)
            % waistPostions is an array of the positions of all waists between z1
            % and z2. If no arguments are given, then all waists in the beam 
            % path are returned.
            %
            % A second optional output argument may be requested which contains 
            % an array of beamq objects for each of the waists found:
            % [waistPostions,waistBeamQs] = path1.getWaists(z1,z2)
            
            % if z limits aren't defined, make them -inf..inf
            if nargin<2
                z1=-inf;
                z2=inf;
            end
            
            % first find all the waists after the first component
            
            % find the z positions off all comps in domain, add the end of
            % zdomain to the list, loop will go through the number of
            % comps, not including the end of zdomain
            
            zmin = min([z1 z2]);
            zmax = max([z1 z2]);
            
            zComponents = [pathobj.components.z];
            zComponents = zComponents(zmin < zComponents & zComponents <=zmax);
            
            nLoop = length(zComponents);
            
            zComponents = [zComponents,zmax];
            
            zWaists = [];
            % loop through z values
            for kk = 1:nLoop
                
                % find the waist position relative to the current component
                
                beamAfterComponent = pathobj.qPropagate(zComponents(kk));
                
                waistPos = zComponents(kk) - beamAfterComponent.waistZ;
                
                % make sure it's between this component and the next, if yes,
                % keep track of that one. (> current, <= next)
                
                if waistPos > zComponents(kk) && waistPos <= zComponents(kk+1)
                    zWaists = [zWaists;waistPos]; %#ok<AGROW>
                end
            end
            
            % finally as a special case, get the waist before the first
            % component, make sure it's between the beginning of the
            % zdomain and the z of first component.
            leftStep = 1;

            if length(zComponents)<2 && zComponents(1)==Inf
                zComponents(1)=leftStep;
            end
            
            beamBeforeFirstComponent = pathobj.qPropagate(zComponents(1)-leftStep);
            
            waistPos = zComponents(1) - leftStep - beamBeforeFirstComponent.waistZ;
            
            if waistPos > zmin && waistPos <= zComponents(1)
                zWaists = [waistPos;zWaists];
            end
            
            zList = zWaists;
            
            if nargout > 1
                %qList = []; % can't do this because we want to fill with objects
                for jj = 1:length(zList)
                    qList(jj) = pathobj.qPropagate(zList(jj)); %#ok<AGROW>
                end
                if isempty(zList)
                    qList = [];
                end
            end
            
        end
        % these are for making graphics
        function plothandle = plotBeamWidth(pathobj,zdomain,varargin)
            % -- beamPath.plotBeamWidth --
            % Creates a plot of the beam width over some z domain, the seed
            % beam is used as the input. Although the seed beam may be defined anywhere
            % on the beam path.
            % Example:
            % z = 0:.01:5;
            % path1.plotBeamWidth(z,{arguments passed to MATLAB plot function})
            % This will make a plot of the beam width over the domain
            % defined by z. Extra arguments are passed to MATLAB's plotting
            % function and can be used to change the line style. One could
            % use 'r--' to make the plot appear as a red dashed line.
            %
            % If an output argument is used, the function returns the plot handle
            % of the positive side of the beam, this useful for making legends.
            % Example:
            % plot1 = path1.plotBeamWidth(z,'b')
            % plot2 = path2.plotBeamWidth(z,'r')
            % legend([plot1 plot2],'Beam Path 1','Beam Path 2')
            qplot = pathobj.qPropagate(zdomain);
            
            ploth = qplot.plotBeamWidth(zdomain,varargin{:});
            if nargout>0
                plothandle = ploth;
            end
        end
        function plothandle = plotGouyPhase(pathobj,zdomain,dontWrap,varargin)
            % -- beamPath.plotBeamWidth --
            % Creates a plot of the accumulate guoy phase in degrees over some z domain, the seed
            % beam is used as the input. Although the seed beam may be defined anywhere
            % on the beam path.
            % Example:
            % z = 0:.01:5;
            % path1.plotGouyPhase(z,'wrap',{arguments passed to MATLAB plot function})
            % This will make a plot of the beam width over the domain
            % defined by z. Extra arguments are passed to MATLAB's plotting
            % function and can be used to change the line style. One could
            % use 'r--' to make the plot appear as a red dashed line.
            %
            % the second argument is a string, if it is 'nowrap' then the Gouy Phase
            % will not wrap inside of +-180 degrees, if it is anything else, or omitted, it will
            % be wrapped.
            %
            % If an output argument is used, the function returns the plot handle.
            if nargin<3
                dontWrap = '';
            end
            gouyPlot = pathobj.gouyPhase(zdomain);
            
            if ~strcmpi(dontWrap,'nowrap')
                gouyPlot = mod(gouyPlot+180,360)-180;
            end
            
            ploth = plot(zdomain,gouyPlot,varargin{:});
            if nargout>0
                plothandle = ploth;
            end
        end
        function plothandle = plotComponents(pathobj,zdomain,yoffset,plotString)
            % -- beamPath.plotComponents -- 
            % Make plot with the positions and labels of components.
            % Example:
            % z = 0:.01:5;
            % hold on
            % path1.plotBeamWidth(z,'r')
            % path1.plotComponents(z,yoffset,'r*')
            % hold off
            % This will plot the beam width and label the plot with
            % the component positions and labels. yoffset is used to offset
            % the y position (in the units of the vertical axis). If omitted,
            % the default value is 0.
            % The third argument chooses the graphic to plot. If omitted, 'b*' is used.
            % If no arguments are given, the function will plot all components in the path.
            compsToPlot = pathobj.components;
            
            componentZ = [compsToPlot.z];
            
            if nargin<3
                yoffset = 0;
            end
            if nargin<4
                plotString = 'b*';
                if ischar(yoffset) % allow the second argument to be the plotstring
                    plotString = yoffset;
                    yoffset = 0;
                end
            end
            if nargin>1
                zmin = min(zdomain);
                zmax = max(zdomain);
                
                inRangeIndex = zmin < componentZ & componentZ <=zmax;
                compsToPlot = compsToPlot(inRangeIndex);
            end
            
            componentZcell = {compsToPlot.z};
            componentLabels = {compsToPlot.label};
            compCount = length(componentLabels);
            offsetLabel = cell(1,compCount);
            for j = 1:compCount
                if mod(j,2)==1
                    offsetLabel{j} = [componentLabels(j);{''}];
                else
                    offsetLabel{j} = [{''};componentLabels(j)];
                end
            end
            
            multiplothandle = plot([componentZcell{:}],yoffset*ones(length(componentZcell)),plotString);
            if ~isempty(multiplothandle) % don't do stuff if there's nothing to plot
                ploth = multiplothandle(1);
                cellfun(@(z,lab) text(z,yoffset,lab),componentZcell,offsetLabel)
            else
                ploth = [];
            end
            if nargout>0
                plothandle = ploth;
            end
        end
        function plotBeams(pathobj,zdomain,yoffset,colorString)
            % -- beamPath.plotBeams --
            % Annotates the plot object with textarrow objects which show the location
            % of the seed and target beams, and also displays some information about them.
            % Example:
            % path1.plotBeams(zdomain,yoffset,colorString)
            % yoffset offsets the textarrows in the y direction (default is 0). colorString sets the 
            % color of the textarrow objects ('r' will make them red etc.).
            if nargin < 2
                zdomain = [-Inf Inf];
            end
            if nargin < 3
                yoffset = 0;
            end
            if nargin < 4
                colorString = 'k';
                if ischar(yoffset) % allow the second argument to be the colorstring
                    colorString = yoffset;
                    yoffset = 0;
                end
            end

            if ~isempty(pathobj.seedq.q) && min(zdomain) < pathobj.seedz && pathobj.seedz <= max(zdomain)
                seedString = {'Seed Beam;';['w = ' num2str(pathobj.seedq.beamWidth/1e-6) '\mum'];...
                    ['ROC = ' num2str(pathobj.seedq.radiusOfCurvature) 'm']};
                if pathobj.seedq.radiusOfCurvature == Inf
                    seedString{3} = 'ROC = \infty';
                end
                [seedPlotx, seedPloty] = dsxy2figxy(gca,pathobj.seedz,yoffset);
                annotation('textarrow',[seedPlotx-.01 seedPlotx],[seedPloty+.05 seedPloty],...
                    'string',seedString,'FontSize',8,'Color',colorString);
            end

            if ~isempty(pathobj.targetq.q) && min(zdomain) < pathobj.targetz && pathobj.targetz <= max(zdomain)
                targetString = {'Target Beam;';['w = ' num2str(pathobj.targetq.beamWidth/1e-6) '\mum'];...
                    ['ROC = ' num2str(pathobj.targetq.radiusOfCurvature) 'm']};
                if pathobj.targetq.radiusOfCurvature == Inf
                    targetString{3} = 'ROC = \infty';
                end

                if ~isempty(pathobj.seedq.q)
                    targetString = [targetString;['Overlap = ' num2str(pathobj.targetOverlap)]];
                end

                [targetPlotx, targetPloty] = dsxy2figxy(gca,pathobj.targetz,yoffset);
                annotation('textarrow',[targetPlotx+.01 targetPlotx],[targetPloty-.05 targetPloty],...
                    'string',targetString,'FontSize',8,'Color',colorString);
            end

        end
        function plotWaists(pathobj,zdomain,colorString)
            arrowLength = .02;
            
            if nargin < 2
                zdomain = [-Inf Inf];
            end
            if nargin < 3
                colorString = 'k';
            end
            
            [waistZs,waistQs]=pathobj.getWaists(min(zdomain),max(zdomain));
            
            for jj = 1:length(waistZs)
                waistString = {['w_0 = ' num2str(round(waistQs(jj).beamWidth/1e-6)) '\mum']};
                [waistPlotx, waistPloty] = dsxy2figxy(gca,waistZs(jj),waistQs(jj).beamWidth);
                [waistPlotx, waistPlotyneg] = dsxy2figxy(gca,waistZs(jj),-waistQs(jj).beamWidth);
                try
                annotation('textarrow',[waistPlotx waistPlotx],[waistPloty+arrowLength waistPloty],...
                    'string',waistString,'FontSize',8,'Color',colorString);
                
                annotation('textarrow',[waistPlotx waistPlotx],[waistPlotyneg-arrowLength waistPlotyneg],...
                    'string','','FontSize',8,'Color',colorString);
                catch
                    error('You may be trying to plot a waist outside of the plot domain')
                end
            end
        end
        function plotSummary(pathobj,zdomain,varargin)
            % -- beamPath.plotSummary --
            % Plots a summary of the beampath. This is a two panel plot with beam 
            % width on the top and gouy phase on the bottom.
            % The components and beams are also plotted and the axis are labeled.
            % Example:
            % path1.plotSummary(zdomain)
            % Plots the plot summary over the domain given in zdomain, 
            % if zdomain is omitted, it chooses a domain based on what is
            % in the beam path
            
            % set plot defaults
            if nargin<3
                plotComponents = 1;
                plotBeamWidth = 1;
                plotBeams = 0;
                plotWaists = 1;
                plotGouyPhase = 1;
            else % they are choosing, set all to false
                plotComponents = 0;
                plotBeamWidth = 0;
                plotBeams = 0;
                plotWaists = 0;
                plotGouyPhase = 0;
            end
            
            % parse plot options
            for kk=1:length(varargin)
                switch lower(varargin{kk})
                    case 'components'
                        plotComponents=1;
                    case 'beamwidth'
                        plotBeamWidth=1;
                    case 'beams'
                        plotBeams=1;
                    case 'waists'
                        plotWaists=1;
                    case 'gouyphase'
                        plotGouyPhase=1;
                    case 'nobeams'
                        warning('nobeams is depreciated.');
                        plotComponents = 1;
                        plotBeamWidth = 1;
                        plotGouyPhase = 1;
                    otherwise
                        error(['I don''t understand how to plot ' varargin{kk}]);
                end
            end
            
            if plotBeamWidth && plotGouyPhase
                doSubPlots = 1;
            else
                doSubPlots = 0;
            end
            
            if ~plotBeamWidth && ~plotGouyPhase
                error('You must plot either the Beam Width or the Gouy Phase')
            end
            
            if nargin<2 || isempty(zdomain)
                zlist = [pathobj.components.z];
                if ~isempty(pathobj.seedq.q) && plotBeams
                    zlist = [zlist,pathobj.seedz];
                end
                if ~isempty(pathobj.targetq.q) && plotBeams
                    zlist = [zlist,pathobj.targetz];
                end
                
                if isempty(zlist)
                    error('Nothing to plot')
                end
                
                zmin = min(zlist);
                zmax = max(zlist);
                
                if zmin == zmax %make a very basic plot if there is only one object to plot
                    zmin = zmin - 0.5;
                    zmax = zmin + 1;
                end
                
                zlength = zmax - zmin;
                
                zdomain = linspace(zmin - 0.1*zlength,zmax + 0.1*zlength,100);
            end
            
            clf
            
            if doSubPlots
                subplot(2,1,1)
            end
            if plotBeamWidth
                hold on
                pathobj.plotBeamWidth(zdomain,'r');
                axis tight
                grid on
                if plotComponents
                    pathobj.plotComponents(zdomain,'b*');
                end
                if plotWaists
                    pathobj.plotWaists(zdomain);
                end
                if plotBeams
                    pathobj.plotBeams(zdomain);
                end
                hold off
                
                ylabel('Beam Width (m)')
            end
            title(['Beam Path Summary: ' inputname(1)])
            
            if doSubPlots
                subplot(2,1,2)
            end
            if plotGouyPhase
                hold on
                pathobj.plotGouyPhase(zdomain,'wrap','r');
                axis tight
                grid on
                if plotComponents
                    pathobj.plotComponents(zdomain,'b*');
                end
                hold off
                
                ylabel('Gouy Phase (degrees)')
                
            end
            xlabel('axial dimension, z (m)')
            
            %return to top plot
            if doSubPlots
                subplot(2,1,1)
            end
        end
        % used for finding overlap optimization
        function [optimizedPathobj,optimumOverlap] = optimizePath(pathobj,varargin)
            % -- beamPath.optimizePath --
            % This function will return a new path object after having optimized the component
            % positions of a given path object in order to maximize the overlap of 
            % the beem propagated from the seed and the beam at the target.
            % The arguments should go like this:
            % newPath = oldPath.optimizePath('lens1',[-1.5,-1],'lens3',[2,2.3],'target',[4,inf])
            % (where oldPath is an existing beamPath object with at least the components 'lens1' and 'lens3'
            % this will attempt to find the optimum modematching for the components in the given range
            % for each component. If a component is not named, it will not be moved. The target beam
            % position can also be optimized. The numbers represent the valid range of positions that can 
            % be searched for that object. (Eg. 'lens1' is allowed to be placed between z=-1.5 and z=-1)
            %
            % options:
            %   '-v' Include the string '-v' (verbose) as one of your input arguments to enable some
            %        output to the command window.
            %   '-r' The option '-r' will make the placement restriction arguments
            %        relative to the current component position. If you want to allow 'lens1' to move plus or minus
            %        10 centimeters, the arguments would go ...'lens1',[-0.1,0.1]...
            %   '-c' Use this option and a function handle as the following argument, this function will be used as a
            %        custom cost function that will be minimized while the beam path is modified, instead of maximizing the
            %        target modematching. See customcostexample.m in the alm/examples directory.
            lvargin = length(varargin);
            if lvargin<1  
                error('Input arguments are required for optimizePath.')
            end
            
            % default options
            verbose = 0;
            relative = 0;
            costFunction = @(path) 1-path.targetOverlap;
            costTransform.f = @(cost) 1-cost;
            costTransform.name = 'overlap';
            
            flagIndex = find(strncmp(varargin,'-',1)); % find the option flag
            removeInds = [];
            
            for ind=flagIndex
                if any(varargin{ind}=='v') || any(varargin{ind}=='V')
                    verbose = 1; % identify verbose flag
                end
                if any(varargin{ind}=='r') || any(varargin{ind}=='R')
                    relative = 1; % identify relative flag
                end
                if any(varargin{ind}=='c') || any(varargin{ind}=='C') % identify cost function flag
                    costFunction = varargin{ind+1};
                    costTransform.f = @(cost) cost;
                    costTransform.name = 'cost';
                    removeInds = [removeInds,ind+1];
                end;
                removeInds = [removeInds,ind];
            end
            
            for ind=sort(removeInds,'descend')
                varargin = {varargin{1:end~=(ind)}};
            end
            lvargin = length(varargin);
            
            if mod(lvargin,2)==1
                error('Number of args not even, each component needs a range, you can use inf to be unconstrained')
            end
            % this block sets the initial values of the optimization to the current component locations,
            % it also sets the upper and lower bounds equal to the current location, if the component is not 
            % named in the input arguments, he is constrained to his initial point. (The optimization algorithm 
            % treats this case in a smart way).
            z0 = [pathobj.components.z,pathobj.targetz];
            UB = z0;
            LB = z0;
            
            % this parses the arguments
            argNames = {varargin{1:2:lvargin}};
            for j = 1:length(argNames)
                originShift = 0;
                if strcmpi('target',argNames(j))
                    if relative
                        originShift = z0(end);
                    end
                    UB(end) = originShift + varargin{2*j}(2);
                    LB(end) = originShift + varargin{2*j}(1);
                    continue
                end
                namedIndex=pathobj.findComponentIndex(argNames{j});
                if relative
                    originShift = z0(namedIndex);
                end
                UB(namedIndex) = originShift + varargin{2*j}(2);
                LB(namedIndex) = originShift + varargin{2*j}(1);
            end
            

            cost = @(zVec,vargs) pathobj.applyCostFunc(costFunction,zVec,vargs);
            [zOptimized,optimumCost] = fminsearchbnd(cost,z0,LB,UB); 
            
            % obfuscated output
            if verbose
                disp([costTransform.name ' is ' num2str(costTransform.f(optimumCost)) '. for z='])
                disp(zOptimized);
            end
            
            % return optimized beampath
            optimizedPathobj = pathobj.duplicate;
            optimizedPathobj.batchMove(zOptimized);
            
            optimumOverlap = costTransform.f(optimumCost);
        end
        function [pathList,costList] = chooseComponents(pathin,varargin)
            % -- beamPath.chooseComponents --
            % This function will return an array of beamPath objects which are created 
            % by selecting components from an array of component objects and optimizing the positions using
            % the optimizePath method. Use a second output argument and the function will return an array of the mode
            % overlap calcualted for each path.
            %
            % The beamPath object which is used to call the function will provide the initial conditions
            % and should include objects in the desired starting positions. The components in the initial beam path
            % will not be part of the list that the components will be chosen from. If an empty array is used
            % instead of a list, the component will not be changed, but the placement will be changed to optimize
            % the mode overlap. 
            %
            % Each component which is being changed should be passed its own component list to choose from.
            % The algorithm recognizes when a component object is in more than one list and won't try solutions which
            % use the same component more than once. If the number of times a component is used doesn't matter, one
            % should make a new component list with the component.duplicate method and pass a duplicate list for each
            % component being replaced.
            %
            % See the example file for more information
            %
            % Syntax: path.chooseComponents(componentLabel,componentList,[lower_bound upper_bound],...)
            %       path is the beampath being optimized, 'componentLabel' is the string label of a component already
            %       defined in path which will be replaced by components in the 'componentList' array. The bounds
            %       are numbers which restrict the z coordinate of the component when optimizing.
            %
            % example: [pathlist,overlaplist] = path1.chooseComponents('lens1',lensList,[.5 .75],'lens2',lensList,[1.2 1.3])
            %
            % options:
            %   '-v' Include the string '-v' for verbose output.
            %   '-t' Use the option '-t' followed by a number as the next argument to set
            %          a threshold below which solutions will not attempt to be optimized after the initial
            %          placement of the components. This can be used to reduce the search time significantly when the number
            %          of combinations is large.
            %   '-r' will make the upper and lower Z bounds of components act relative to the current z position.
            %   '-c' Use this option and a function handle as the following argument, this function will be used as a
            %        custom cost function that will be minimized while the beam path is modified, instead of maximizing the
            %        target modematching. See customcostexample.m in the alm/examples directory.
            %

            lvargin = length(varargin);
            if lvargin<1
                error('Input arguments are required for chooseComponents.')
            end

            % check for options
            verbose = 0;
            relative = 0;
            costFunction = @(path) 1-path.targetOverlap;
            costTransform.f = @(cost) 1-cost;
            costTransform.name = 'overlap';
            costTransform.skip = @(transCost,thresh) transCost <= thresh;
            transCostThresh = 0;
            useDefaultThresh = 1;
            
            flagIndex = find(strncmp(varargin,'-',1)); % find the option flag
            removeInds = [];
            
            for ind = flagIndex
                if any(varargin{ind}=='v') || any(varargin{ind}=='V') % identify verbose flag
                    verbose = 1; 
                end
                if any(varargin{ind}=='r') || any(varargin{ind}=='R') % identify verbose flag
                    relative = 1; 
                end
                if any(varargin{ind}=='t') || any(varargin{ind}=='T') % identify threshold flag
                    transCostThresh = varargin{ind+1};
                    useDefaultThresh = 0;
                    removeInds = [removeInds,ind+1];
                end
                if any(varargin{ind}=='c') || any(varargin{ind}=='C') % identify threshold flag
                    %customCost = 1;
                    costFunction = varargin{ind+1};
                    costTransform.f = @(cost) cost;
                    costTransform.name = 'cost';
                    costTransform.skip = @(transCost,thresh) transCost >= thresh;
                    if useDefaultThresh
                        transCostThresh = Inf;
                    end
                    removeInds = [removeInds,ind+1];
                end
                removeInds = [removeInds,ind];
            end
            
            for ind=sort(removeInds,'descend')
                varargin = {varargin{1:end~=(ind)}};
            end
            lvargin = length(varargin);
            
            % This block will look in the arguments for a reference to the target
            indexTarget = find(strcmpi(varargin,'target'),2);
            if length(indexTarget)>1
                error('Please only reference the target once in arguments.')
            end

            % If the target is mentioned, make note of the range and remove from varargin
            if indexTarget
                targetRange = varargin{indexTarget+1};
                ii = 1:lvargin;
                varargin = {varargin{ii~=indexTarget & ii~=indexTarget+1}};
                lvargin = length(varargin);
            end

            % now parse the rest of the arguments
            if mod(lvargin,3)~=0
                error('Number of arguments is invalid.')
            end

            compLabels = {varargin{1:3:lvargin}};
            compLists = {varargin{2:3:lvargin}};
            compRanges = {varargin{3:3:lvargin}};
            numlists = length(compLists);
            
            % build the cell of arguments for optimizePath, also replace empty lists with the current component,
            %      also make sure that all the named components exist
            argCell = {};
            for jj = 1:numlists
                argCell = [argCell,{compLabels{jj},compRanges{jj}}]; %#ok<AGROW>
                if isempty(compLists{jj})
                    compLists{jj} = pathin.component(compLabels{jj}).duplicate;
                else
                    pathin.findComponentIndex(compLabels{jj}); %this will return an error if a component doesn't exist
                end
            end
            
            % add target to argument cell if it was named
            if indexTarget
                argCell = [argCell,{'target',targetRange}];
            end
            
            % pass relative flag if used
            if relative
                argCell = [argCell,{'-r'}];
            end
            
            % pass custom cost
            argCell = [argCell,{'-c',costFunction}];
            
            % duplicate the user's path so we don't screw up the original
            pathobj = pathin.duplicate;
            
            % set up the loop
            indicies = ones(1,numlists);
            indexlimit = cellfun(@length,compLists);
            
            step = 1;
            totalSteps = prod(indexlimit);
            costList = Inf*ones(totalSteps,1);
            bestCost = Inf;
            pathList(totalSteps,1) = beamPath;
            
            if verbose
                dispString = ['Searching through ' num2str(totalSteps) ' combinations.'];
                if ~useDefaultThresh
                    dispString = [dispString  ' Initial ' costTransform.name ' threshhold: ' num2str(transCostThresh)];
                end
                disp(dispString)
            end
            
            while indicies(end)<=indexlimit(end) % main loop
                % do things here
                
                % check that a component is not getting double-used
                doubleFlag = 0;
                for jj = 1:numlists
                    for kk = jj+1:numlists
                        if compLists{jj}(indicies(jj)) == compLists{kk}(indicies(kk))
                            doubleFlag = 1;
                            costList(step) = NaN;
                            break
                        end
                    end
                    if doubleFlag
                        break
                    end
                end

                if ~doubleFlag
                    % replace the components with components from the list
                    for jj = 1:numlists
                        % grab a component from the list
                        compCopy = compLists{jj}(indicies(jj)).duplicate;
                        pathobj.replaceComponent(compLabels{jj},compCopy);
                    end

                    % check the unoptimized cost against the threshold, if it passes, try to optimize
                    if  ~costTransform.skip(costTransform.f(costFunction(pathobj)),transCostThresh)
                        [pathList(step),costList(step)] = pathobj.optimizePath(argCell{:});
                        if verbose
                            if costList(step) < bestCost
                                bestCost = costList(step);
                            end
                            disp(['current ' costTransform.name ': ' num2str(costTransform.f(costList(step))) '. best so far: '...
                                num2str(costTransform.f(bestCost)) '. ' num2str(totalSteps-step) ' more to try.'])
                        end
                    end
                end
                
                % overall step counter
                step = step + 1;
                % increment index vector
                for kk = 1:numlists;
                    indicies(kk) = indicies(kk)+1;
                    if indicies(kk)>indexlimit(kk) && kk ~= numlists
                        indicies(kk)=1;
                        continue
                    end
                    break
                end
            end
            
            % remove illegal or skipped combinations
            ix = ~isnan(costList);
            costList = costList(ix);
            pathList = pathList(ix);
            
            if isempty(costList)
                error(['Could not find any solutions, try changing initial conditions or changing ' costTransform.name ' threshold.'])
            end
            
            % now sort path list according to cost
            [costList,ix] = sort(costList);
            costList = arrayfun(costTransform.f,costList);
            pathList = pathList(ix);
        end
        % used for beam width fitting
        function fittedPath = fitBeamWidth(pathobj,zPred,widthPred)
            % -- beamPath.fitBeamWidth --
            % Fits a given set of beam width and position data and returns a beamPath
            % object with a seed beam to match to the data.
            % Example:
            % path2 = path1.fitBeamWidth(zPred,widthPred)
            % zPred is the position data and widthPred is the width data to fit to. 
            % Both should be array vectors of the same length.
            if nargin<3
                error('Not enough input arguments for fitBeamWidth')
            end
            
            if length(zPred)~=length(widthPred)
                error('Position and beam width vectors must be same length in fitBeamWidth.')
            end
            
            waistParamsStart = [pathobj.seedq.waistSize pathobj.seedq.waistZ];
            
            lambda = pathobj.seedq.lambda;
            
            fitWaistParams = fminsearch(@(waistParams)pathobj.beamFitError(waistParams,zPred,widthPred,lambda),waistParamsStart);
            
            fittedPath = pathobj.duplicate;
            fittedPath.seedq = beamq.beamWaistAndZ(fitWaistParams(1),fitWaistParams(2),lambda);
        end
    end % methods
end % classdef