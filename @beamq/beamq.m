classdef beamq
    % -- beamq --
    %
    %    The beamq class is used for calculating the properties of a 
    %    gaussian beam. It stores the complex beam q parameter used to 
    %    chatacterize a gaussian beam. By setting the 'q' and 'lambda' 
    %    properties of a  beamq object, it can return various properties
    %    of the beam.
    %
    %    Constructor Methods:
    %         Note: The default value of lambda is 1064nm.
    %         beamq(q,lambda) - returns a beamq object with the defined q
    %            value for wavelength lambda (in meters).
    %         beamq.beamWaistAndZ(w0,Z,lambda) - returns a beamq object
    %            with a waist of w0 (in meters) at position Z (in meters)
    %            with wavelength lambda (in meters).
    %         beamq.beamWaistAndR(w0,R,lambda) - returns a beamq object
    %            with a waist of w0 (in meters) and a radius of R 
    %            (in meters) at Z=0 with wavelength lambda (in meters). 
    %         beamq.beamWidthAndR(w,R,lambda) - returns a beamq object
    %            with a beam width of w (in meters) at Z=0 and a radius of 
    %            R (in meters) Z=0 with wavelength lambda (in meters).
    %
    %    Properties:
    %         beamWidth - the 1/e electric field amplitude radius.
    %         waistSize - the beam width at the waist of the beam.
    %         waistZ - the relative distance from the beam waist, this is
    %             also the real part of the q parameter.
    %         divergenceAngle - the angle between the propagation axis and
    %             the diverging radius of the beam in the far field.
    %         radiusOfCurvature - The radius of curvature of the constant 
    %             phase front of the beam.
    %         rayleighRange - The axial length scale of the beam focus,
    %             this is also the imaginary part of the q parameter.
    %
    properties
        q;
        lambda = 1064e-9;
    end
    properties (Dependent, SetAccess = private)
        waistSize;
        waistZ;
        divergenceAngle;
        radiusOfCurvature;
        beamWidth;
        rayleighRange;
    end
    methods (Static)
        % alternative variable parameterization to the constructor
        function qobjout = beamWaistAndZ(w0,Z,lambda)
            if nargin<3
                lambda = 1064e-9;
            end
            
            ZR = pi*w0^2/lambda;
            q = Z + i*ZR;
            
            qobjout = beamq(q);
            qobjout.lambda = lambda;
        end
        function qobjout = beamWaistAndR(w0,R,lambda)
            if nargin<3
                lambda = 1064e-9;
            end
            
            ZR = pi*w0^2/lambda;
            q = (1/R - i/ZR)^-1;

            qobjout = beamq(q);
            qobjout.lambda = lambda;
        end
        function qobjout = beamWidthAndR(w, R, lambda)
            if nargin<3
                lambda = 1064e-9;
            end
            Z = R/(1+(R*lambda/pi/w^2)^2);
            ZR = sqrt(Z*(R-Z));
            q = Z + 1i*ZR;
            
            qobjout = beamq(q);
            qobjout.lambda = lambda;
        end
        
        function qvalout = transformValue(qvalin,M)
            qvalout=(M(1,1).*qvalin+M(1,2))./(M(2,1).*qvalin+M(2,2));     
        end
    end
    methods
        % constructor and data access methods
        function qobj = beamq(qvalue,lambda)
            if nargin>0
                qobj.q = qvalue;
                if nargin>1
                    qobj.lambda = lambda;
                end
            end
        end
        function qobj = set.q(qobj,qvalue)
            if imag(qvalue)<0
                error('imaginary part of q parameter must be positron')
            end
            qobj.q=qvalue;
        end
        function qobj = set.lambda(qobj,newlambda)
            if newlambda<=0
                error('lambda must be positron')
            end
            qobj.lambda = newlambda;
        end
        function qnew = duplicate(qold)
            % -- beamq.duplicate --
            % Make a copy of a beamq object with the same properties 
            % as the original
            % Example:
            % beamcopy = beam1.duplicate;
            qnew = beamq(qold.q,qold.lambda);
        end
        % methods for dependent properties
        function valOut = get.waistSize(qin)
            lambda = qin.lambda;
            
            valOut = sqrt(imag(qin.q).*lambda/pi);
        end
        function valOut = get.rayleighRange(qin)
            w0 = qin.waistSize;
            lambda = qin.lambda;
            
            valOut = pi*w0.^2./lambda;
        end
        function valOut = get.divergenceAngle(qin)
            w0 = qin.waistSize;
            zR = qin.rayleighRange;
            
            valOut = w0./zR;
        end
        function valOut = get.waistZ(qin)
            valOut = real(qin.q);
        end
        function valOut = get.beamWidth(qin)
            z = qin.waistZ;
            zR = qin.rayleighRange;
            w0 = qin.waistSize;

            valOut = w0 .* sqrt( 1 + (z./zR).^2 );
        end
        function valOut = get.radiusOfCurvature(qin)
            z = qin.waistZ;
            zR = qin.rayleighRange;
            
            if z ~= 0
                valOut = z.*(1+(zR./z).^2);
            else
                valOut = Inf;
            end
        end
        % methods for making useful calculations
        function fraction = overlap(beam1,beam2)
            % -- beamq.overlap --
            % Find the overlap fraction of 2 beams (assumes axial symmetry).
            % Example:
            % fraction = overlap(beam1,beam2)
            
            q1 = beam1.q;
            q2 = beam2.q;

            w1 = beam1.waistSize;
            w2 = beam2.waistSize;
            
            lambda = beam1.lambda;
            if lambda~=beam2.lambda
                error('Cannot overlap beams of different wavelength')
            end
            
            fraction = (2*pi/lambda*w1*w2 * 1/abs(conj(q2)-q1))^2; %square for 2D modematching
        end
        function beamout = transform(beamin,M)
            % -- beamq.transform --
            % Creates a new beamq object after being transformed by an ABCD matrix.
            % Example:
            % newbeam = oldbeam.transform(M)
            % This transforms the oldbeam object and placed the new object into
            % newbeam, using the ABCD matrix M.
            qin = beamin.q;
            
            qout=beamq.transformValue(qin,M);
            
            beamout = beamin.duplicate;
            beamout.q = qout;
        end
        % plotting
        function plothandle = plotBeamWidth(qarray,zdomain,varargin)
            % -- beamq.plotBeamWidth --
            % Given an array of beamq objects, this function will plot 
            % the beam width.
            washold = ishold;
            if ~washold
                hold on
            end
            ploth = plot(zdomain,[qarray.beamWidth],varargin{:});
            plot(zdomain,-[qarray.beamWidth],varargin{:});
            if ~washold
                hold off
            end
            if nargout>0
                plothandle = ploth;
            end
        end
    end
end