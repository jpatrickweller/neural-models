 function [f] = FitModel(params,contrast,response,FITSTR,error,showplot)
% Evaluate model fit  based on the passed in parameters (params) and 
% returns an error metric of the data (specified by 'error') given the
% parameters. Contrast and response are the stimuli and responses to be
% fitted, respectively.

% Display the fit of the guesses in the fitting proceedure
if nargin < 5
    disp('Must specify type of error to be used...')
    return
elseif (nargin < 6)
    showplot = 0;
end

% Calculate responses based on given parameters, stimuli, and model-type
prediction = ComputeModel(params,contrast,FITSTR);

% Determine which type of error to calculate
if strcmp(error,'Poisson') || strcmp(error,'poisson')
    f = sum(prediction)-(log(prediction')*response);  % -1 * log-likelihood (Poisson)
elseif strcmp(error,'Gaussian') || strcmp(error,'gaussian')
    f = sum((prediction-response).^2); % Sum of squarred error (Gaussian)
elseif strcmp(error,'Bernoulli') || strcmp(error,'bernoulli')
    f = -1 * (log(prediction')*response + (log(1-prediction)'*(1-response))); % Bernoulli error
elseif strcmp(error,'NegativeBinomial') || strcmp(error,'negativebinomial')
    kappa = params(end);
    mu = prediction;
    sigsq = mu + kappa * mu.^2;
    p = (sigsq - mu) ./ sigsq;
    r = mu.^2 ./ (sigsq - mu);
    if kappa == 0
        f = sum(prediction)-(log(prediction')*response);  % -1 * log-likelihood (Poisson in special case of negative binomial)
    elseif any(r <= 0)
        disp('r <= 0')
        f = 5000000000;
        keyboard
    elseif any(p>1) || any(p<0)
        disp('p > 1 or p < 0')
        f = 5000000000;
        keyboard
    elseif any(~isreal(mu))
        disp('mu has immaginary component')
        f = 5000000000;
        keyboard
    else
        f = -sum(gammaln(r+response)) + sum(gammaln(r))...
            - (log(p)'*response) - (log(1-p)'*r); % -1 * log-likelihood (negative binomial)
    end
end

if showplot == 1
    if strcmp(FITSTR,'asymmetric') 
        figure(500); clf; hold on; grid on;
        plot(contrast,response,'ko')
        plot(contrast,prediction,'m*');
        set(gca,'ylim',[0 max(prediction)],...
            'xlim',[min(contrast(:)) max(contrast(:))],...
            'ylim',[min(contrast(:)) max(contrast(:))]);
        drawnow
    end
    if strcmp(FITSTR,'conicsection') || strcmp(FITSTR,'conicsection_xy')
        figure(500); 
        campos = get(gca,'cameraposition');
        clf; hold on; grid on;
        plot3(contrast(:,1),contrast(:,2),response,'ko')
        plot3(contrast(:,1),contrast(:,2),prediction,'m*')
        set(gca,'ylim',[0 max(prediction)],...
            'xlim',[min(contrast(:)) max(contrast(:))],...
            'ylim',[min(contrast(:)) max(contrast(:))],...
            'cameraposition',campos);
        drawnow
    end
end

% Handle bizarre response prediction cases
if (isnan(f))
    f = 500000000;
    disp('LL is nan...');
end

if (all(prediction == response))
    f = 500000000;
    disp('EXACT MATCH!');
end

if ~isreal(f)
    f = 500000000;
    disp('Irrational -LL...')
end

if any(~isreal(prediction))
    f = 500000000;
    disp('Irrational prediction...')
end

end


