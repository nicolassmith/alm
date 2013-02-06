classdef component < handle
    % -- component --
    %
    %    The component class is used to represent optical components.
    %    The properties associated with a component object are: the 
    %    position (z), the ray transfer matrix, or ABCD matrix (M) and the
    %    label. Component objects used with the special constructors also
    %    have a 'parameters' property.
    %
    %    Methods:
    %        component - The component constructor, for making a component
    %            object. The arguments are (M,z,label)
    %        <a href="matlab:help component.duplicate">duplicate</a> - creates and identical component object.
    %        <a href="matlab:help component.lens">lens</a> - creates a lens component.
    %        <a href="matlab:help component.dielectric">dielectric</a> - creates a dielectric (thick lens) component.
    %        <a href="matlab:help component.curvedMirror">curvedMirror</a> - creates a curved mirror object.
    %        <a href="matlab:help component.flatMirror">flatMirror</a> - creates flat mirror object.
    %        <a href="matlab:help component.propagator">propagator</a> - creates free space propagator object.
    %        <a href="matlab:help component.combine">combine</a> - creates a composit component object with a 
    %            tranfer matrix which is the product of the components in the calling array
    %     
    properties (SetAccess = private) %only z and label can be changed at the moment
                                     %eventually i need to make it possible for components to have
                                     %a transfer matrix as a function of Z
        M
        type = 'other';
    end
    properties
        label = 'none';
        z
    end 
    properties (SetAccess = private)
        parameters=struct;
    end
    methods (Static)
        % these methods are to construct different types of components
        function objout = lens(focalLength,Z,label)
            % -- component.lens --
            % Create a lens component object.
            % Example:
            % mylens = component.lens(f,z,label);
            % This creates a lens component with focal length f at position
            % z. label is a string which is used to identify the component.
            if nargin < 2
                Z = 0;
            end
            
            numcomps = length(focalLength);
            if numcomps>1
                zlength=length(Z);
                if zlength~=numcomps 
                    if zlength~=1
                        error('List of focal lengths must be the same length as list of z positions')
                    else
                        Z(1:numcomps,1)=Z;
                    end
                end
                
                if nargin>2
                    lablength=length(label);
                    if ischar(label)
                        lablength = 1;
                    end
                    if lablength~=numcomps
                        if lablength~=1
                            error('List of focal lengths must be the same length as list of labels')
                        else
                            singlabel = label;
                            label = cell(numcomps,1);
                            for jj = 1:numcomps
                                label{jj}=singlabel;
                            end
                        end
                    end
                end
                
                objout(numcomps,1) = component;
                for jj = 1:numcomps
                    if nargin>2
                        objout(jj) = component.lens(focalLength(jj),Z(jj),label{jj});
                    else
                        objout(jj) = component.lens(focalLength(jj),Z(jj));
                    end
                end
                return
            end
            
            M = [ 1, 0; -1/focalLength, 1];
            objout = component(M,Z);
            objout.type='lens';
            objout.parameters = rmfield(objout.parameters,'none');
            objout.parameters.focalLength = focalLength;
            if nargin > 2
                objout.label = label;
            end
        end
        function objout = curvedMirror(radiusOfCurvature,Z,label)
            % -- component.curvedMirror --
            % Create a curved mirror component object.
            % Example:
            % mylens = component.curvedMirror(ROC,z,label);
            % This creates a lens component with radius of curvature ROC at position
            % z. label is a string which is used to identify the component.
            if nargin < 2
                Z = 0;
            end
            
            numcomps = length(radiusOfCurvature);
            if numcomps>1
                zlength=length(Z);
                if zlength~=numcomps 
                    if zlength~=1
                        error('List of radii must be the same length as list of z positions')
                    else
                        Z(1:numcomps,1)=Z;
                    end
                end
                
                if nargin>2
                    lablength=length(label);
                    if ischar(label)
                        lablength = 1;
                    end
                    if lablength~=numcomps
                        if lablength~=1
                            error('List of radii must be the same length as list of labels')
                        else
                            singlabel = label;
                            label = cell(numcomps,1);
                            for jj = 1:numcomps
                                label{jj}=singlabel;
                            end
                        end
                    end
                end
                
                objout(numcomps,1) = component;
                for jj = 1:numcomps
                    if nargin>2
                        objout(jj) = component.curvedMirror(radiusOfCurvature(jj),Z(jj),label{jj});
                    else
                        objout(jj) = component.curvedMirror(radiusOfCurvature(jj),Z(jj));
                    end
                end
                return
            end
            
            M = [ 1, 0; -2/radiusOfCurvature, 1];
            objout = component(M,Z);
            objout.type='curved mirror';
            objout.parameters=rmfield(objout.parameters,'none');
            objout.parameters.ROC = radiusOfCurvature;
            if nargin > 2
                objout.label = label;
            end
        end
        function objout = flatMirror(Z,label)
            % -- component.flatMirror --
            % Create a flat Mirror component object.
            % Example:
            % mylens = component.flatMirror(z,label);
            % This creates a flat mirror component at position
            % z. label is a string which is used to identify the component.
            if nargin < 1
                Z = 0;
            end
            M = [ 1, 0; 0, 1];
            objout = component(M,Z);
            objout.type='flat mirror';
            if nargin > 1
                objout.label = label;
            end
        end
        function objout = dielectric(R1, R2, thickness, n, Z, label)
            % -- component.dielectric --
            % Create a dielectric component object.
            % Example:
            % mylens = component.dielectric(R1, R2, th, n, Z, label);
            % This creates a dielectric (thick lens) component at position
            % z. label is a string which is used to identify the component.
            if nargin < 1
                Z = 0;
            end
            dist = @(x) [ 1, x; 0, 1];
            refract = @(R,n1,n2) [1 0 ; (n1-n2)/(R*n2) n1/n2];
            M = refract(R2, n, 1)*dist(thickness)*refract(R1,1,n);
            objout = component(M,Z);
            objout.type='dielectric';
            objout.parameters=rmfield(objout.parameters,'none');
            objout.parameters.length = thickness;
            if nargin > 1
                objout.label = label;
            end
        end
        function objout = propagator(DZ,Z,label)
            % -- component.propagator --
            % Create a free-space propagator component object.
            % Example:
            % mylens = component.propagator(dz,z,label);
            % This creates a propagator component with length dz at position
            % z. label is a string which is used to identify the component.
            if nargin < 2
                Z = 0;
            end
            M = [ 1, DZ; 0, 1];
            objout = component(M,Z);
            objout.type='propagator';
            objout.parameters=rmfield(objout.parameters,'none');
            objout.parameters.length = DZ;
            if nargin > 2
                objout.label = label;
            end
        end
    end
    methods
        function objout = component(M,Z,label)
            if nargin > 0
                objout.M = M;
                objout.z = Z;
                if nargin > 2
                    objout.label = label;
                end
            else
                objout.M = [ 1, 0; 0, 1];
                objout.z = 0;
            end
            if isempty(fieldnames(objout.parameters))
                objout.parameters.none = [];
            end
        end
        function objout = duplicate(objin)
            % -- component.duplicate --
            % make a new component (or array of components) with the
            % same properties as the original.
            % Example:
            % lens1copy = lens1.duplicate;
            [nn mm] = size(objin);
            if prod([nn mm])>0
                objout(nn,mm) = component;
                for ii = 1:nn
                    for jj = 1:mm
                        objout(ii,jj) = component(objin(ii,jj).M,objin(ii,jj).z,objin(ii,jj).label);
                        objout(ii,jj).type = objin(ii,jj).type;
                        objout(ii,jj).parameters = objin(ii,jj).parameters;
                    end
                end
            else
                objout = component;
                objout(1) = [];
            end
        end
        function objout = combine(componentTrain)
            % -- component.combine --
            % Squashes a list of components together to make a single component
            % with transfer matrix equal to the product of the list, multiplied
            % in order of index array.
            Mtrain = [ 1, 0; 0, 1];
            for j = 1:length(componentTrain)
                Mtrain = componentTrain(j).M * Mtrain;
            end
            objout = component(Mtrain,0);
            objout.type = 'composite';
        end
        function obj = set.z(obj,zin)
            if length(zin)>1 || ~isnumeric(zin)
                error('Sorry, axial position Z must be a 1x1 numeric array')
            else
                obj.z = zin;
            end
        end
        function obj = set.M(obj,Min)
            if all(size(Min)==[2 2]) && isnumeric(Min)
                obj.M = Min;
            else
                error('Sorry, transformation matrix M must be a 2x2 numeric matrix')
            end
        end
        function setLength(obj,L)
            % -- component.setLength --
            % set the real space length of an object.
            
            if isfield(obj.parameters,'none')
                obj.parameters=rmfield(obj.parameters,'none');
            end
            
            obj.parameters.length = L;
        end
        function display(obj)            
            topRow =    {'label','z (m)','type'};
            lineBreak = {'-----','-----','----'};
            parameterColumn{length(obj)+2,1} = [];
            parameterColumn{1} = 'parameters';
            parameterColumn{2} = '----------';
            for j = 1:length(obj)
                x = evalc('disp(obj(j).parameters)');
                
                parameterColumn{j+2} = x(5:min([end-2,23]));
            end            
            output = [topRow;lineBreak;[{obj.label}.',{obj.z}.',{obj.type}.']];
            output = [output,parameterColumn];           %#ok<NASGU>
            outstring = evalc('disp(output)');
            outstring = outstring(outstring~='''' & outstring~='[' & outstring~=']');
            
            disp(' ');
            disp([inputname(1) ' = '])
            disp(' ');
            disp(outstring);
            disp(' ');
        end
    end
end