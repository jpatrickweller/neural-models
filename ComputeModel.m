function [response] = ComputeModel(params,stimuli,surftype)
% Compute an LN or LNLN model (linear or nonlienar combinations of 
% cone signals followed by a Naka-Rushton function) on passed in stimuli.
% Several different models may be computed depending on specified surftype.

% Initialize some conserved parameters
response =  nan(size(stimuli,1),1);

% Fit requested surface
if strcmp(surftype, 'symmetric_LN')
    % symmetric_LN: LN model, only accepts vector as input. Allows for 2 symmetric 
    % Naka-Rushtons tail to tail. Stimuli may be positive and/or negative.
    
    % Assign parameter values
    sig = params(2);
    exp = params(3);
    bl = params(4);
    A= params(1) - bl;
    
    % Generate responses
    response = A * abs(stimuli).^exp ./ (abs(stimuli).^exp + sig.^exp) + bl;

elseif strcmp(surftype, 'asymmetric_LN')
    % asymmetric: LN model, allows for 2 NakaRushtons tail to tail. Baseline
    % and upper asymptote is shared, but other params are independent.
    % Stimuli should have both positive and negative values. For
    % half-rectified responses, send sig2 as Inf.
    
    % Assign parameter values
    sig1 = params(3);
    sig2 = params(4);
    exp1 = params(5);
    exp2 = params(6);
    bl = params(7);
    A = params(1) - bl;
    
    % Generate responses
    L = stimuli > 0; 
    response(L) = A * abs(stimuli(L)).^exp1 ./ (abs(stimuli(L)).^exp1 + sig1.^exp1) + bl;
    response(~L) = A * abs(stimuli(~L)).^exp2 ./ (abs(stimuli(~L)).^exp2 + sig2.^exp2) + bl;

elseif strcmp(surftype, 'asymmetric_2D_LN')
    % Same as asymmetric, but stimuli are 2-dimensional, sent in as a matrix.
    % The sign of sigma designates excitation or suppression.
    
    % Assign parameters
    sig1 = params(2);
    sig2 = params(3);
    exp = params(4);
    bl = params(5);
    A = params(1) - bl;
    rot = params(6);
    
    % Unpack stimuli
    rotMat = [cos(rot) sin(rot); -sin(rot) cos(rot)];
    tempRotPts = (rotMat * [stimuli(:,1) stimuli(:,2)]')';
    contrast = tempRotPts(:,1);
    posIdx = contrast >= 0;
    poscon = contrast(posIdx);
    negcon = abs(contrast(~posIdx));

    % Generate repsonses
    response(posIdx) = (A * poscon.^exp) ./ (poscon.^exp + sig1.^exp) + bl;
    if sig2 > 0
        response(~posIdx) = (A * negcon.^exp) ./ (negcon.^exp + sig2.^exp) + bl;
    else
        response(~posIdx) = -(bl * negcon.^exp) ./ (negcon.^exp + abs(sig2).^exp) + bl;
    end
    
elseif strcmp(surftype, 'asymmetric_2D_LNLN') 
    % asymmetric_2D_LNLN: All parameters are shared except for
    % c50s along each axis.  c50s change smoothly with angle, as the major
    % and minor axes of an ellipse.
    % r(theta) = sqrt((a*sin(theta)).^2 + (b*cos(theta)).^2))
    
    % Assign parameter values
    sig1 = params(2);
    sig2 = params(3);
    sig3 = params(4);
    sig4 = params(5);
    exp = params(6);
    bl = params(7);
    A = params(1) - bl;
    rot = params(8);
    
    % Unpack stimuli
    rotMat = [cos(rot) sin(rot); -sin(rot) cos(rot)];
    tempRotPts = (rotMat * [stimuli(:,1) stimuli(:,2)]')';
    [theta,rho] = cart2pol(tempRotPts(:,1),tempRotPts(:,2));
    
    % Calclulate all of the possible sigma values for each quadrant 
    sigmasQ1 = (sig1*sig2) ./ (abs((sig1*sin(theta)).^2) + abs((sig2*cos(theta)).^2)).^(1/2);
    sigmasQ2 = (sig2*sig3) ./ (abs((sig2*cos(theta)).^2) + abs((sig3*sin(theta)).^2)).^(1/2);
    sigmasQ3 = (sig3*sig4) ./ (abs((sig3*sin(theta)).^2) + abs((sig4*cos(theta)).^2)).^(1/2);
    sigmasQ4 = (sig4*sig1) ./ (abs((sig4*cos(theta)).^2) + abs((sig1*sin(theta)).^2)).^(1/2);
    
    % Assign sigma values based on quadrant
    sig = nan(size(response));
    for n = 1:numel(theta)
        if theta(n) > pi/2
            sig(n) = sigmasQ3(n);
        elseif theta(n) > 0 && theta(n) <= pi/2
            sig(n) = sigmasQ4(n);
        elseif theta(n) <= 0 && theta(n) > -pi/2
            sig(n) = sigmasQ1(n);
        elseif theta(n) <= -pi/2
            sig(n) = sigmasQ2(n);
        end        
    end
    
    % Generate responses
    response = A * rho.^exp ./ (rho.^exp + sig.^exp) + bl;
    
elseif strcmp(surftype,'conicsection_rt') || strcmp(surftype,'conicsection')
    % Building a surface whose c50 changes as a conic section.
    % *sig1 controls the preferred direction c50.
    % *sig2 controls the anti-preferred c50. Negative value is suppression.
    % c50 is sent in as reciprocal so that the variable can be contrinuous. 
    % *orthosig is the orthogonal axis. Constrained to be symmetrical.
    % Positive values are the minor axis of ellipse. Negative values are
    % width of the curved hyperbola at the focus. Reciprocal of c50 is sent
    % in so that values can be continuous. If orth = Inf, surface is 1D.
    
    % Assign parameter values
    sig1 = 1/params(2);
    sig2 = 1/params(3);
    orthosig = params(4); % converted to reciprocal below
    exp = params(5);
    bl = params(6);
    A = params(1)-bl;
    rot = params(7);
    
    % Unpack stimuli
    rotMat = [cos(rot) sin(rot); -sin(rot) cos(rot)];
    tempRotPts = (rotMat * stimuli')';
    x = tempRotPts(:,1);
    y = tempRotPts(:,2);
    posIdx = x >= 0;
    [theta,rho] = cart2pol(x,y);
    theta = abs(theta);
    
    % Calclulate all of the possible sigma values for each quadrant 
    % r(theta) = (a*b) ./ sqrt((a*sin(theta)).^2 +/- (b*cos(theata)).^2))
    sig = nan(size(theta));
    
    % Positive contrast
    if orthosig > 0 % ellipse
        nom = sig1 * 1/orthosig;
        denom = (sig1*sin(theta(posIdx))).^2 + (1/orthosig*cos(theta(posIdx))).^2;
        sig(posIdx) = nom ./ sqrt(denom);
    elseif orthosig < 0 % hypoerbola
        % there is a sign switch and half rectification here.
        nom = sig1 * abs(1/orthosig);
        denom = abs(min(0,(sig1*sin(theta(posIdx))).^2 - (1/orthosig*cos(theta(posIdx))).^2));
        sig(posIdx) = nom ./ sqrt(denom);
    else % 1D naka rushton
        %sig(posIdx) = sqrt(sig1^2 ./ (cos(theta(posIdx))).^2);
        sig(posIdx) = sig1 ./ cos(theta(posIdx));
    end
    
    % Negative contrast
    if orthosig > 0 % ellipse
        if sig2 == Inf
            sig(~posIdx) = (1/orthosig) ./ sin(theta(~posIdx));
        else
            nom = sig2 * 1/orthosig;
            denom = (sig2*sin(theta(~posIdx))).^2 + (1/orthosig*cos(theta(~posIdx))).^2;
            sig(~posIdx) = nom ./ sqrt(denom);
        end
    elseif orthosig < 0 % hypoerbola
        % there is a sign switch and half rectification here.
        nom = sig2 * abs(1/orthosig);
        denom = abs(min(0,(sig2*sin(theta(~posIdx))).^2 - (1/orthosig*cos(theta(~posIdx))).^2)); 
        sig(~posIdx) = nom ./ sqrt(denom);
    else % 1D naka rushton        
        sig(~posIdx) = abs(sig2 ./ cos(theta(~posIdx)));
    end
    
    response = A * rho.^exp ./ (rho.^exp + sig.^exp) + bl;

    
elseif strcmp(surftype,'conicsection_xy') || strcmp(surftype,'conicsection2')
    % An alternative parameterization of 'conicsection' in which cone
    % contrast is speified in cartesian coordinates.
    
    % Assign parameter values
    sig1 = (1/params(2));
    sig2 = (1/params(3));
    orthosig = (params(4)); % converted to reciprocal below
    exp = params(5);
    bl = params(6);
    A = params(1)-bl;
    rot = params(7);
        
    % Rotate axes such that the x-axis is the principle axis.
    rotMat = [cos(rot) sin(rot); -sin(rot) cos(rot)];
    tempRotPts = (rotMat * stimuli')';
    x = tempRotPts(:,1);
    y = tempRotPts(:,2); 
    y = abs(y);
    posIdx = x >= 0;
    
    % Convert 2D stim (rho,theta) into 1D stim (effective contrast)
    effcont = nan(size(x));
    
    if orthosig > 0 % ellipse
        effcont(posIdx) = sqrt(x(posIdx).^2 + (y(posIdx)./(1/orthosig/sig1)).^2);
        if sig2 == Inf
            effcont(~posIdx) = abs(y(~posIdx)./(1/orthosig/sig1));
        else
            effcont(~posIdx) = sqrt((x(~posIdx)/(sig2/sig1)).^2 + (y(~posIdx)./(1/orthosig/sig1)).^2);
        end
    elseif orthosig < 0 % hypoerbola
        effcont(posIdx) = sqrt(max(x(posIdx).^2 - (y(posIdx)./(1/abs(orthosig)/sig1)).^2,0));
        if sig2 == Inf
            effcont(~posIdx) = x(~posIdx)./(sig2/sig1);
            %effcont(~posIdx) = -y(~posIdx)./(1/abs(orthosig)/sig1);
        else
            effcont(~posIdx) = sqrt(max((x(~posIdx)./(sig2/sig1)).^2 - (y(~posIdx)./(1/abs(orthosig)/sig1)).^2,0));
        end
    else % 1D naka rushton
        effcont(posIdx) = x(posIdx);
        effcont(~posIdx) = abs(x(~posIdx)./(sig2/sig1));
    end
    if any(effcont<0)
        disp('negative effcont on non-hyperbola')
        keyboard
        effcont = max(effcont,0);
    end
    response = A * effcont.^exp ./ (effcont.^exp + sig1.^exp) + bl;

    
elseif strcmp(surftype,'conicsection_sym')
    
    % A symmetric version of 'conicsection'.
    % Building a surface whose c50 changes as a conic section.
    % *sig1 controls the preferred direction c50.
    % *sig2 controls rectification (1=full, 0=half)
    % *orthosig is the c50 along the orthogonal axis. Constrained to be symmetrical.
    % Positive values are the minor axis of ellipse. Negative values are
    % width of the curved hyperbola at the focus. Reciprocal of c50 is sent
    % in so that values can be continuous. If orth = Inf, surface is 1D.
    
    % Assign parameter values
    sig1 = params(2); % asymmetric axis c50
    sig2 = params(3); % binary (0 or 1) that makes PD half or full rectified
    orthosig = (1/params(4)); % symmetric axis c50 (input is reciprocal of c50)
    exp = params(5);
    bl = params(6);
    A = params(1)-bl;
    rot = params(7);
        
    % Rotate axes such that the x-axis is the principle axis.
    rotMat = [cos(rot) sin(rot); -sin(rot) cos(rot)];
    tempRotPts = (rotMat * stimuli')';
    x = tempRotPts(:,1);
    y = tempRotPts(:,2); 
    
    effcont = max(x,0).^2 + sig2 * x.^2 + 1/orthosig * sig1 * y.^2;

    if any(effcont<0)
        disp('negative effcont values...')
        effcont = max(effcont,0);
    end
    response = A * effcont.^exp ./ (effcont.^exp + sig1.^exp) + bl;
    
else
    error('Must specify desired fit.');
end


