
function PD_Estimate_Variance()
global data

% This script models the responses of LN neurons tuned to various 
% directions in the LM plane. By fitting the modeled responses, we can 
% estimate bias and confidence intervals in the estimated preferred
% directions.

% Set up figure
figure(1000); clf;
set(gcf,'units','normalized','pos',[.25 .2 .5 .6],'NumberTitle','off',...
    'Name','Variance of Tuning Model');
modelfig = get(gcf,'UserData');
modelfig.conpanel = uipanel('pos',[.525 .525 .45 .45],'parent',gcf,'title','Control Panel');
modelfig.surfpanel = uipanel('Pos',[.025 .525 .45 .45],'Parent',gcf);
modelfig.fitspanel = uipanel('pos',[.025 .025 .95 .45],'parent',gcf);
conpanel = get(modelfig.conpanel,'UserData');
surfpanel = get(modelfig.surfpanel,'UserData');
fitspanel = get(modelfig.fitspanel,'UserData');

% Set up controls for analysis
conpanel.uicontrols.nsamps = uicontrol('style','edit','parent',modelfig.conpanel,...
    'units','normalized','pos',[.65 .7 .2 .07],'string',10,'fontsize',10);
conpanel.labels.nsamps = uicontrol('Parent',modelfig.conpanel,'Units','normalized',...
    'pos',[.55 .8 .4 .15],'HorizontalAlignment','center',...
    'style','text','string','Samples/preferred direction','FontSize',10);
conpanel.uicontrols.nDirs = uicontrol('parent',modelfig.conpanel,'units','normalized',...
    'style','edit','pos',[.65 .35 .2 .07],'string',30,'fontsize',10);
conpanel.labels.nDirs = uicontrol('parent',modelfig.conpanel,'units','normalized',...
    'pos',[.55 .45 .4 .15],'style','text','HorizontalAlignment','center',...
    'string','# of preferred color directions','fontsize',10);

% Set up controls for parameter values
conpanel.uicontrols.upperA = uicontrol('parent',modelfig.conpanel,'units','normalized',...
    'style','edit','pos',[.15 .8 .2 .07],'string',50,'fontsize',10);
conpanel.labels.upperA = uicontrol('parent',modelfig.conpanel,'units','normalized',...
    'pos',[.05 .88 .4 .1],'style','text','HorizontalAlignment','center',...
    'string','Upper Asymptote','fontsize',10);
conpanel.uicontrols.baseline = uicontrol('parent',modelfig.conpanel,'units','normalized',...
    'style','edit','pos',[.15 .55 .2 .07],'string',5,'fontsize',10);
conpanel.labels.baseline = uicontrol('parent',modelfig.conpanel,'units','normalized',...
    'pos',[.05 .63 .4 .1],'style','text','HorizontalAlignment','center',...
    'string','Baseline','fontsize',10);
conpanel.uicontrols.exp = uicontrol('parent',modelfig.conpanel,'units','normalized',...
    'style','edit','pos',[.15 .3 .2 .07],'string',3,'fontsize',10);
conpanel.labels.exp = uicontrol('parent',modelfig.conpanel,'units','normalized',...
    'pos',[.05 .38 .4 .1],'style','text','HorizontalAlignment','center',...
    'string','Exponent','fontsize',10);
conpanel.uicontrols.kappa = uicontrol('parent',modelfig.conpanel,'units','normalized',...
    'style','edit','pos',[.15 .05 .2 .07],'string',2,'fontsize',10);
conpanel.labels.kappa = uicontrol('parent',modelfig.conpanel,'units','normalized',...
    'pos',[.05 .12 .4 .1],'style','text','HorizontalAlignment','center',...
    'string','Kappa (variance)','fontsize',10);

% Start button
conpanel.uicontrols.startanalysis = uicontrol('parent',modelfig.conpanel,'style','pushbutton',...
    'units','normalized','pos',[.6 .05 .3 .15],'string','Start Analysis','fontsize',10,...
    'backgroundColor',[.2 1 .2],'callback',@StartAnalysis);

% Set up surfpanel
surfpanel.axes = axes('parent',modelfig.surfpanel,'units','normalized','pos',[.15 .2 .7 .7],...
    'tickdir','out');
axis square; box on;
xlabel('L-cone contrast')
ylabel('M-cone contrast')
surfpanel.surftype = 'conicsection_xy';
surfpanel.errortype = 'NegativeBinomial';
upperA = str2double(conpanel.uicontrols.upperA.String);
baseline = str2double(conpanel.uicontrols.baseline.String);
exp = str2double(conpanel.uicontrols.exp.String);
kappa = str2double(conpanel.uicontrols.kappa.String);
surfpanel.realparams = [upperA nan 0 0 exp baseline nan kappa];

% Set up fitspanel
fitspanel.axes.LL = axes('parent',modelfig.fitspanel,'units','normalized',...
    'pos',[.1 .2 .8 .7],'box','on','xlim',[-pi/2 pi/2],'tickdir','out');
xlabel('Preferred Direction (deg)')
ylabel('PD Estimate Error (deg)')

try
    
    % Try to load previously saved data
    load('PD Est Var Model data','data')

    % Organize and express in degrees
    angs = data.angs./pi*180;
    angdiffs = data.angdiff./pi*180;
    fitmean = mean(angdiffs,2);
    fitstd = std(angdiffs,[],2);
    
    % Plot the previously fit data
    axes(fitspanel.axes.LL); cla; hold on; grid on;
    h = shadedErrorBar(angs,fitmean,fitstd,'r-');
    alpha(.5)
    h.edge(1).LineStyle = 'none';
    h.edge(2).LineStyle = 'none';
    plot(angs,angdiffs,'ko','ButtonDownFcn',@DispDatasets);
    xlim([min(angs) max(angs)])
    xlabel('Preferred Direction (deg)')
    ylabel('PD Estimate Error (deg)')

    
    % Fill in parameters from previous fits
    set(conpanel.uicontrols.upperA,'string',data.realparams(1))
    set(conpanel.uicontrols.baseline,'string',data.realparams(6))
    set(conpanel.uicontrols.exp,'string',data.realparams(5))
    set(conpanel.uicontrols.nDirs,'string',size(data.LL,1))
    set(conpanel.uicontrols.nsamps,'string',size(data.LL,2))
    surfpanel.surftype = data.surftype;
    surfpanel.errortype = data.errortype;
    surfpanel.realparams = data.realparams;

catch
    
    % If no previous data exists, propt the user to check the path.
    disp('No previously modeled data found.')
end

% Save variables
set(modelfig.conpanel,'UserData',conpanel);
set(modelfig.surfpanel,'UserData',surfpanel);
set(modelfig.fitspanel,'UserData',fitspanel);
set(gcf,'UserData',modelfig);

end

function StartAnalysis(~,~)
global gl data

% Load variables
modelfig = get(gcf,'UserData');
conpanel = get(modelfig.conpanel,'UserData');
fitspanel = get(modelfig.fitspanel,'UserData');
surfpanel = get(modelfig.surfpanel,'UserData');

% GL struct stores experimental variables and tracks the progress of the model.
gl.nPres = 5;
gl.nSamps = str2double(get(conpanel.uicontrols.nsamps,'string'));
gl.nRnds = 3;
npd = str2double(get(conpanel.uicontrols.nDirs,'string'));
gl.allAngs = linspace(-pi/2,pi/2,npd);

% Fill in specified parameter values
surfpanel.realparams = nan(1,8);
surfpanel.realparams(1) = str2double(conpanel.uicontrols.upperA.String);
surfpanel.realparams(6) = str2double(conpanel.uicontrols.baseline.String);
surfpanel.realparams(5) = str2double(conpanel.uicontrols.exp.String);
surfpanel.realparams(8) = str2double(conpanel.uicontrols.kappa.String);

% Initialize structure 'data' for saving model fit data.
data.angs = gl.allAngs;
data.realparams = surfpanel.realparams;
data.surftype = surfpanel.surftype;
data.errortype = surfpanel.errortype;
data.angdiff = nan(npd,gl.nSamps);
data.params = cell(npd,gl.nSamps);
data.LL = nan(npd,gl.nSamps);
data.stim = cell(npd,gl.nSamps);
data.resp = cell(npd,gl.nSamps);

% Iterate through each angle and number of samples
for rot = 1:numel(gl.allAngs)
    gl.currentAng = rot;
    for sampn = 1:gl.nSamps
        
        gl.currentSamp = sampn;
        
        disp(['Angle ' num2str(gl.currentAng) ' of ' num2str(numel(gl.allAngs))])
        disp(['Sample ' num2str(gl.currentSamp) ' of ' num2str(gl.nSamps)])
        
        disp('Choosing Stim...')
        ChooseLMStimuli;
        
        disp('Modeling Responses...')
        CreateModelSurface;
        
        disp('Fitting Model Data...')
        FitModelData()
        
        % Load in figure variables
        fitspanel = get(modelfig.fitspanel,'UserData');
        surfpanel = get(modelfig.surfpanel,'UserData');
        angdiff = fitspanel.params(end-1) - gl.allAngs(gl.currentAng);
        if angdiff > pi/2
            angdiff = angdiff - pi;
        elseif angdiff < -pi/2
            angdiff = angdiff + pi;
        end
        data.angdiff(gl.currentAng,gl.currentSamp) = angdiff;
        data.params{gl.currentAng,gl.currentSamp} = fitspanel.params;
        data.LL(gl.currentAng,gl.currentSamp) = fitspanel.LL;
        data.stim{gl.currentAng,gl.currentSamp} = [surfpanel.Lcc surfpanel.Mcc];
        data.resp{gl.currentAng,gl.currentSamp} = surfpanel.nsp;
        
    end
end

disp('Saving data...')
save('PD Est Var Model data','data')

% Organize and express in degrees
angs = data.angs./pi*180;
angdiffs = data.angdiff./pi*180;
fitmean = mean(angdiffs,2);
fitstd = std(angdiffs,[],2);

% Plot
axes(fitspanel.axes.LL); cla; hold on; grid on;
h = shadedErrorBar(angs,fitmean,fitstd,'r-');
alpha(.5)
h.edge(1).LineStyle = 'none';
h.edge(2).LineStyle = 'none';
plot(angs,angdiffs,'ko','ButtonDownFcn',@DispDatasets);
xlim([min(angs) max(angs)])
xlabel('Preferred Direction (deg)')
ylabel('PD Estimate Error (deg)')

disp('Analysis finished.')

end

function ChooseLMStimuli()
global gl

% Load variables
modelfig = get(gcf,'UserData');
conpanel = get(modelfig.conpanel,'UserData');
surfpanel = get(modelfig.surfpanel,'UserData');
fitspanel = get(modelfig.fitspanel,'UserData');

% Pre-Defined Variables
thetaspace = pi/4;
rhospace = .5;
colmaxcc = .09;
lummaxcc = .7;

% Construct Polar Grid
if gl.nRnds == 2
    rhospace = rhospace/2;
elseif gl.nRnds == 3
    rhospace = rhospace/2;
    thetaspace = thetaspace/2;
else
    disp('Hardcoded rounds 2 and 3, but not others yet...')
    keyboard
end
thetas = shiftdim(0:thetaspace:2*pi-thetaspace);
rhos = shiftdim(rhospace:rhospace:1);

% Ennumerate all conditions
PolRhoIdx = fullfact([numel(thetas) numel(rhos)]);
rhothetalist = repmat([thetas(PolRhoIdx(:,1)) rhos(PolRhoIdx(:,2))],gl.nPres,1);
thetas = rhothetalist(:,1);
rhos = rhothetalist(:,2);
[tempLcc,tempMcc] = pol2cart(thetas,rhos);

% Scale Cone Contrast Units for Monitor
scale = lummaxcc*colmaxcc./sqrt((colmaxcc.*cos(thetas-pi/4)).^2 ...
    +(lummaxcc.*sin(thetas-pi/4)).^2);
Lcc = tempLcc .* scale;
Mcc = tempMcc .* scale;
surfpanel.Lcc = Lcc;
surfpanel.Mcc = Mcc;

% Save variables
set(modelfig.conpanel,'UserData',conpanel);
set(modelfig.surfpanel,'UserData',surfpanel);
set(modelfig.fitspanel,'UserData',fitspanel);
set(gcf,'UserData',modelfig);

end

function CreateModelSurface()
global gl

% Load Figure Variables
modelfig = get(gcf,'UserData');
%conpanel = get(modelfig.conpanel,'UserData');
surfpanel = get(modelfig.surfpanel,'UserData');
%fitspanel = get(modelfig.fitspanel,'UserData');

% Specify parameters for new neuron
params = surfpanel.realparams;
params(end-1) = gl.allAngs(gl.currentAng);
[pdx,pdy] = pol2cart(params(end-1),1);
proj = [pdx pdy] * [surfpanel.Lcc surfpanel.Mcc]';
params(2) = 1/(max(proj)/2);
gensig = ComputeNakaRushtonJPW(params(1:7),[surfpanel.Lcc surfpanel.Mcc],surfpanel.surftype);

% Use generator signal as mean of a draw from negative binomial distribution
kappa = surfpanel.realparams(end);
mu = gensig;
sigsq = mu + kappa * mu.^2;
p = (sigsq - mu) ./ sigsq;
r = mu.^2 ./ (sigsq - mu);
surfpanel.nsp = nbinrnd(r,1-p);

% Save Figure Variables
%set(modelfig.conpanel,'UserData',conpanel);
set(modelfig.surfpanel,'UserData',surfpanel);
%set(modelfig.fitspanel,'UserData',fitspanel);
set(gcf,'UserData',modelfig);

end

function FitModelData()

% Load Figure Variables
modelfig = get(gcf,'UserData');
surfpanel = get(modelfig.surfpanel,'UserData');
fitspanel = get(modelfig.fitspanel,'UserData');

% Load data parameters
Lcc = surfpanel.Lcc;
Mcc = surfpanel.Mcc;
nsp = surfpanel.nsp;
angs = linspace(-pi/2,pi/2,9);
GOF = nan(numel(angs),1);

% Set up some variables
[~,rho] = cart2pol(Lcc,Mcc);
ub = [max(nsp)        300 0 0 10 max(nsp)  pi 5];
lb = [min(nsp) 1/max(rho) 0 0  0       0  -pi 0];
a(1,:) = [0 -1 1 0 0 0 0 0]; % make sure sig1 is greater than sig2
a(2,:) = [-1 0 0 0 0 1 0 0]; % make sure upperA is above baseline
b = zeros(size(a,1),1);
options = optimset('Algorithm','interior-point','MaxFunEvals',5000,...
    'MaxIter',5000,'FinDiffType','central','Hessian','bfgs','display','off',...
    'TolFun',0,'TolCon',0,'FunValCheck','on','AlwaysHonorConstraints','bounds');
Aguess = max(nsp)*.8;
sigguess = max(rho)/2;
expguess = 2;
blguess = min(nsp);
kappaguess = 1;
params = nan(numel(angs),numel(lb));
guessIdx = fullfact([numel(sigguess) numel(angs)]);

% Fit using a variety of initial guesses
for rot = 1:size(guessIdx,1)
    paramsGuess = [Aguess 1/sigguess(guessIdx(rot,1)) 1/sigguess(guessIdx(rot,1)) 0 expguess blguess angs(guessIdx(rot,2)) kappaguess];
    [f1,fval] = fmincon('FitModel',paramsGuess,a,b,[],[],lb,ub,[],options,[Lcc Mcc],nsp,surfpanel.surftype,surfpanel.errortype);
    params(rot,:) = f1;
    GOF(rot) = fval;
end

% Find the best fitting surface
[~,bestIdx] = min(GOF);
params1 = params(bestIdx,:);

% Calculate surface
x = linspace(-max(Lcc),max(Lcc),50);
[xx,yy] = meshgrid(x,x);
surface = ComputeNakaRushtonJPW(params1,[xx(:) yy(:)],surfpanel.surftype);
surface = reshape(surface,size(xx));

% Set colormap and record current view
cmap = repmat(linspace(.7,0,32)',1,3);
colormap(cmap)
prevview = get(gca,'view');

% Plot current dataset
axes(surfpanel.axes); cla; hold on; grid on;
uniquestim = unique([Lcc Mcc],'rows');
maxnsp = max(nsp);
for i = 1:size(uniquestim,1)
    L = Lcc == uniquestim(i,1) &  Mcc == uniquestim(i,2);
    mn = mean(nsp(L))/maxnsp*10;
    h = plot3(uniquestim(i,1),uniquestim(i,2),mean(nsp(L)),'ko');
    set(h,'MarkerFaceColor','black','MarkerSize',max([2*mn 0]+4),'MarkerEdgeColor','none')
end
h = surfc(xx,yy,surface);
set(h(1),'edgecolor','none')
alpha(.5);
set(gca,'view',prevview);
xlabel('L-cone contrast');
ylabel('M-cone contrast');
zlabel('# of spikes');
drawnow

% Save results
fitspanel.params = params1;
fitspanel.LL = -GOF(bestIdx);

% Save Figure Variables
set(modelfig.surfpanel,'UserData',surfpanel);
set(modelfig.fitspanel,'UserData',fitspanel);
set(gcf,'UserData',modelfig);

end

function DispDatasets(~,~)
global data

% Grab current point before anything else
whichpt = get(gca,'CurrentPoint');
whichpt = whichpt(1,[1 2]);

% Load Figure Variables
modelfig = get(gcf,'UserData');
surfpanel = get(modelfig.surfpanel,'UserData');
fitspanel = get(modelfig.fitspanel,'UserData');

% Organize and express in degrees
angs = data.angs./pi*180;
angdiffs = data.angdiff./pi*180;
fitmean = mean(angdiffs,2);
fitstd = std(angdiffs,[],2);

% Which datafile corresponds to selected point?
[~,PDidx] = min(abs(whichpt(1)-angs));
[~,sampIdx] = min(abs(whichpt(2)-angdiffs(PDidx,:)));

% Plot the previously fit data
axes(fitspanel.axes.LL); cla; hold on; grid on;
h = shadedErrorBar(angs,fitmean,fitstd,'r-');
alpha(.5)
h.edge(1).LineStyle = 'none';
h.edge(2).LineStyle = 'none';
plot(angs,angdiffs,'ko','ButtonDownFcn',@DispDatasets);
plot(angs(PDidx),angdiffs(PDidx,sampIdx),'r*');
xlim([min(angs) max(angs)])
xlabel('Preferred Direction (deg)')
ylabel('PD Estimate Error (deg)')


% Load saved data for chosen stimulus set
nsp = data.resp{PDidx,sampIdx};
params = data.params{PDidx,sampIdx};
stim = data.stim{PDidx,sampIdx};
Lcc = stim(:,1);
Mcc = stim(:,2);

% Calculate surface
x = linspace(-max(Lcc),max(Lcc),50);
[xx,yy] = meshgrid(x,x);
surface = ComputeNakaRushtonJPW(params,[xx(:) yy(:)],surfpanel.surftype);
surface = reshape(surface,size(xx));

% Plot current dataset
axes(surfpanel.axes); 
prevview = get(gca,'view');
cla; hold on; grid on;
Lcc = stim(:,1);
Mcc = stim(:,2);
uniquestim = unique([Lcc Mcc],'rows');
maxnsp = max(nsp);
for i = 1:size(uniquestim,1)
    L = Lcc == uniquestim(i,1) &  Mcc == uniquestim(i,2);
    mn = mean(nsp(L))/maxnsp*10;
    h = plot3(uniquestim(i,1),uniquestim(i,2),mean(nsp(L)),'ko');
    set(h,'MarkerFaceColor','black','MarkerSize',max([2*mn 0]+4),'MarkerEdgeColor','none')
end
h = surfc(xx,yy,surface);
set(h(1),'edgecolor','none')
alpha(.5);
set(gca,'view',prevview);
xlabel('L-cone contrast');
ylabel('M-cone contrast');
zlabel('# of spikes');



end
