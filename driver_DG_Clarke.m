function results = driver_DG_Clarke(probID)
% driver_DG_Clarke
%
% Driver for running DG_Clarke on one problem from TE/TE.mat.
% It runs both direction modes:
%
%   dir = 1 : steepest discrete-gradient direction
%   dir = 2 : MatCSG direction
%
% and plots qf versus nf for both variants.
%
% Folder layout expected:
%
%   DG_Clarke/
%       DG_Clarke.m
%       driver_DG_Clarke.m
%       TE/TE.mat
%       HIT/RMLMMAX001.mat, ..., RMLMMAX080.mat
%
% Example:
%
%   driver_DG_Clarke(1)
%
% or simply:
%
%   driver_DG_Clarke
%
% which runs RMLMMAX001.

clc

if nargin < 1 || isempty(probID)
    probID = 1;
end

% ============================================================
% Root folder of this driver.
% This makes the driver independent of the current MATLAB path.
% ============================================================
root = fileparts(mfilename('fullpath'));

if isempty(root)
    root = pwd;
end

solverFile = fullfile(root,'DG_Clarke.m');
TEfile     = fullfile(root,'TE','TE.mat');
HITdir     = fullfile(root,'HIT');

if ~exist(solverFile,'file')
    error('DG_Clarke.m was not found at: %s', solverFile);
end

if ~exist(TEfile,'file')
    error('TE.mat was not found at: %s', TEfile);
end

if ~exist(HITdir,'dir')
    error('HIT folder was not found at: %s', HITdir);
end

addpath(root)
addpath(fullfile(root,'TE'))
addpath(HITdir)

fprintf('\n==============================================================\n');
fprintf('DG_Clarke driver\n');
fprintf('root   : %s\n', root);
fprintf('TEfile : %s\n', TEfile);
fprintf('HITdir : %s\n', HITdir);
fprintf('==============================================================\n');

% ============================================================
% Load TE
% ============================================================
S = load(TEfile,'TE');

if ~isfield(S,'TE')
    error('The file %s does not contain variable TE.', TEfile);
end

TE = S.TE;

pname = sprintf('RMLMMAX%03d',probID);

if ~isfield(TE.problem,pname)
    error('Problem %s was not found in TE.problem.', pname);
end

P = TE.problem.(pname);

fprintf('\nSelected problem\n');
fprintf('  name       : %s\n', pname);

if isfield(P,'dataset')
    fprintf('  dataset    : %s\n', P.dataset);
end

if isfield(P,'loss')
    fprintf('  loss       : %s\n', P.loss);
end

if isfield(P,'dim')
    fprintf('  dim        : %d\n', P.dim);
end

if isfield(P,'m')
    fprintf('  max terms  : %d\n', P.m);
end

% ============================================================
% Load hitlist
% ============================================================
hitfile = fullfile(HITdir,[pname,'.mat']);

if ~exist(hitfile,'file')
    error('Hitlist file was not found: %s', hitfile);
end

H = load(hitfile);

if isfield(H,'hitlist')
    hitlist = H.hitlist;
else
    hitlist = H;
end

fprintf('\nLoaded hitlist\n');
fprintf('  file       : %s\n', hitfile);

% ============================================================
% Objective function
% ============================================================
if isfield(P,'funf') && isa(P.funf,'function_handle')
    baseFun = @(x) P.funf(x(:));
elseif isfield(P,'data') && exist('TE_realML_minmax_eval','file') == 2
    baseFun = @(x) TE_realML_minmax_eval(x(:),P.data);
else
    error(['No valid objective handle found. Expected P.funf or ', ...
           'TE_realML_minmax_eval(x,P.data).']);
end

% ============================================================
% Initial point from TE
% ============================================================
initialPointName = 'xr';

n = getProblemDimension(P,hitlist);

if isfield(P,'points') && isstruct(P.points) && ...
        isfield(P.points,initialPointName)

    x0 = P.points.(initialPointName);

else

    error('Initial point %s not found in TE.problem.%s.points.', ...
          initialPointName,pname);

end

x0 = x0(:);

if numel(x0) ~= n
    error('Initial point %s for %s has length %d, but n = %d.', ...
          initialPointName,pname,numel(x0),n);
end

% Clip to bounds if available.
[low,upp] = getBounds(P,hitlist,n);

if ~isempty(low) && ~isempty(upp)
    x0 = min(upp,max(low,x0));
end

f0 = baseFun(x0);

if ~isfinite(f0)
    error('Initial objective value is not finite: f0 = %g.', f0);
end

fprintf('\nInitial point\n');
fprintf('  point name : %s\n', initialPointName);
fprintf('  n          : %d\n', n);
fprintf('  f0         : %.16e\n', f0);
fprintf('  ||x0||     : %.6e\n', norm(x0));
% ============================================================
% Target value from hitlist
% ============================================================
[ftarget,hasTarget] = extractTargetFromHitlist(hitlist);

if hasTarget && isfinite(ftarget) && f0 > ftarget
    accf = 1e-4;
    fprintf('\nBenchmark target\n');
    fprintf('  ftarget    : %.16e\n', ftarget);
    fprintf('  accf       : %.1e\n', accf);
else
    hasTarget = false;
    ftarget = f0 - max(1,abs(f0));
    accf = -inf;

    fprintf('\nBenchmark target\n');
    fprintf('  usable target not found or f0 <= ftarget\n');
    fprintf('  accuracy stopping disabled\n');
end

% ============================================================
% Tuning parameters
% ============================================================
tuneBase = struct();

% ============================================================
% Stopping and printing structure
% ============================================================
STbase = struct();

STbase.finit    = f0;
STbase.ftarget  = ftarget;
STbase.accf     = accf;

STbase.nfmax    = 50000;
STbase.secmax   = 300;
STbase.prt      = -1;

STbase.nf       = 0;
STbase.qf       = inf;
STbase.sec      = 0;
STbase.done     = false;
STbase.initTime = cputime;

% ============================================================
% Run both direction modes
% ============================================================
dirs = [1,2];

% Store outputs in a cell array.
results = cell(numel(dirs),1);

for r = 1:numel(dirs)

    dirValue = dirs(r);

    fprintf('\n==============================================================\n');
    fprintf('Running DG_Clarke on %s with dir = %d\n', pname, dirValue);
    fprintf('==============================================================\n');

    R = runOneDirection(baseFun,x0,tuneBase,STbase, ...
                        dirValue,ftarget,hasTarget, ...
                        f0,pname,probID);

    results{r} = R;

    fprintf('\nFinished dir = %d\n', dirValue);
    fprintf('  fsol        = %.16e\n', R.fsol);
    fprintf('  qsol        = %.6e\n', R.qsol);
    fprintf('  logged nf   = %d\n', R.nfFinal);

    fprintf('  driver flag = %s\n', R.driverStatus);
    
    if isstruct(R.outAll) && isfield(R.outAll,'flag')
        fprintf('  solver flag = %s\n', R.outAll.flag);
    end

end

% ============================================================
% Plot log10(qf) versus nf
% ============================================================
fig = figure('Color','w');
set(fig, 'Renderer', 'painters');  

ax = axes('Parent',fig);

hold(ax,'on')
grid(ax,'on')
set(ax,'Box','on')

lineWidth = 1.8;

for r = 1:numel(results)

    R = results{r};

    if hasTarget
        y = R.qfPlot;
        y = max(y, realmin);
        ylog = log10(y);
        plot(ax, R.nfPlot, ylog, 'LineWidth', lineWidth);
    else
        plot(ax, R.nfPlot, R.fbestPlot, 'LineWidth', lineWidth);
    end

end

xlabel(ax,'Number of function evaluations');

if hasTarget
    ylabel(ax,'log_{10}(q_f)');
else
    ylabel(ax,'Best objective value');
end

title(ax,sprintf('%s: DG-Clarke comparison',pname), ...
      'Interpreter','none');

legend(ax,{'dir=1','dir=2'}, 'Location','best');

plotPNG = ...
   fullfile(root,sprintf('plot_%s_logqf_vs_nf_dir1_dir2.png',pname));
plotFIG = ...
    fullfile(root,sprintf('plot_%s_logqf_vs_nf_dir1_dir2.fig',pname));

saveas(fig,plotPNG);
savefig(fig,plotFIG);

fprintf('\nSaved plot PNG : %s\n', plotPNG);
fprintf('Saved plot FIG : %s\n', plotFIG);
% ============================================================
% Save comparison result
% ============================================================
resultFile = ...
    fullfile(root,sprintf('result_%s_DG_Clarke_dir1_dir2.mat',pname));

save(resultFile, ...
    'results','f0','ftarget','hasTarget','pname','probID', ...
    'initialPointName','x0', ...
    'tuneBase','STbase','TEfile','hitfile','plotPNG','plotFIG');

fprintf('Saved result   : %s\n', resultFile);
end
% =====================================================================
function R = runOneDirection(baseFun,x0,tuneBase,STbase,dirValue, ...
                             ftarget,hasTarget,f0,pname,probID)

% Use the same random seed for both direction modes.
rng(1,'twister');

fLog = [];

% Track best evaluated point seen by the driver.
xBestLog  = x0(:);
fBestLog  = f0;
nfBestLog = 0;

tune = tuneBase;
tune.dir = dirValue;

ST = STbase;
ST.nf = 0;
ST.qf = inf;
ST.sec = 0;
ST.done = false;
ST.initTime = cputime;

[xReturned,outAll] = DG_Clarke(@loggedFun,x0,tune,ST);

elapsed = cputime - ST.initTime;

% Returned iterate value.
fReturned = baseFun(xReturned);

% Decide the reported/best solution.
% Use the best evaluated objective value.
if isfinite(fBestLog) && fBestLog < fReturned
    xsol = xBestLog;
    fsol = fBestLog;
    solutionSource = 'best_logged_evaluation';
else
    xsol = xReturned;
    fsol = fReturned;
    solutionSource = 'returned_iterate';
end

nfEval = (1:numel(fLog)).';

% Include nf=0 value for plotting.
fAll = [f0; fLog(:)];
nfPlot = [0; nfEval];

fbestPlot = cummin(fAll);

if hasTarget
    denom = f0 - ftarget;

    if denom <= 0 || ~isfinite(denom)
        denom = max(1,abs(f0-ftarget));
    end

    qfPlot = (fbestPlot - ftarget) ./ denom;
    qsol = (fsol - ftarget) ./ denom;
    qReturned = (fReturned - ftarget) ./ denom;
    qBestLogged = (fBestLog - ftarget) ./ denom;
else
    qfPlot = NaN(size(fbestPlot));
    qsol = NaN;
    qReturned = NaN;
    qBestLogged = NaN;
end

% Driver-side status.
if hasTarget && qsol <= ST.accf
    driverStatus = 'accuracy reached';
elseif numel(fLog) >= ST.nfmax
    driverStatus = 'nfmax reached';
elseif elapsed >= ST.secmax
    driverStatus = 'secmax reached';
elseif isstruct(outAll) && isfield(outAll,'flag')
    driverStatus = outAll.flag;
else
    driverStatus = 'unknown';
end

R = struct();

R.dir = dirValue;
R.problem = pname;
R.probID = probID;

R.xsol = xsol;
R.fsol = fsol;
R.solutionSource = solutionSource;

R.xReturned = xReturned;
R.fReturned = fReturned;

R.xBestLogged = xBestLog;
R.fBestLogged = fBestLog;
R.nfBestLogged = nfBestLog;

R.f0 = f0;
R.ftarget = ftarget;
R.hasTarget = hasTarget;

R.qsol = qsol;
R.qReturned = qReturned;
R.qBestLogged = qBestLogged;

R.fLog = fLog(:);
R.nfEval = nfEval;
R.nfFinal = numel(fLog);
R.elapsed = elapsed;

R.nfPlot = nfPlot;
R.fbestPlot = fbestPlot;
R.qfPlot = qfPlot;

R.driverStatus = driverStatus;
R.outAll = outAll;
R.tune = tune;
R.ST = ST;

    function f = loggedFun(x)
        x = x(:);
        f = baseFun(x);

        fLog(end+1,1) = f;

        if isfinite(f) && f < fBestLog
            fBestLog = f;
            xBestLog = x;
            nfBestLog = numel(fLog);
        end
    end

end
% =====================================================================
function n = getProblemDimension(P,hitlist)

if isfield(P,'dim') && isnumeric(P.dim) && isscalar(P.dim)
    n = P.dim;
elseif isfield(P,'n') && isnumeric(P.n) && isscalar(P.n)
    n = P.n;
elseif isfield(P,'low')
    n = numel(P.low);
elseif isfield(hitlist,'low')
    n = numel(hitlist.low);
elseif isfield(hitlist,'xopt')
    n = numel(hitlist.xopt);
else
    error('Could not determine problem dimension.');
end

end

% =====================================================================
function x = findVectorOfLength(S,n)

x = [];

if isnumeric(S)
    if isvector(S) && numel(S) == n
        x = S(:);
        return
    end

    if ismatrix(S)
        if size(S,1) == n && size(S,2) >= 1
            x = S(:,1);
            return
        elseif size(S,2) == n && size(S,1) >= 1
            x = S(1,:).';
            return
        end
    end
end

if isstruct(S)
    preferred = {'x0','xinit','xstart','x','point','start','startingPoint'};

    fields = fieldnames(S);
    ordered = {};

    for i = 1:numel(preferred)
        if ismember(preferred{i},fields)
            ordered{end+1} = preferred{i}; %#ok<AGROW>
        end
    end

    for i = 1:numel(fields)
        if ~ismember(fields{i},ordered)
            ordered{end+1} = fields{i}; %#ok<AGROW>
        end
    end

    for i = 1:numel(ordered)
        x = findVectorOfLength(S.(ordered{i}),n);
        if ~isempty(x)
            return
        end
    end
end

end

% =====================================================================
function [low,upp] = getBounds(P,hitlist,n)

low = [];
upp = [];

if isfield(P,'low') && isfield(P,'upp') && ...
        numel(P.low) == n && numel(P.upp) == n

    low = P.low(:);
    upp = P.upp(:);
    return
end

if isfield(hitlist,'low') && isfield(hitlist,'upp') && ...
        numel(hitlist.low) == n && numel(hitlist.upp) == n

    low = hitlist.low(:);
    upp = hitlist.upp(:);
    return
end

end

% =====================================================================
function [ftarget,hasTarget] = extractTargetFromHitlist(hitlist)

ftarget = [];
hasTarget = false;

preferred = {'fopt','ftar','ftarget','fmin','fbest','bestf','minf'};

for i = 1:numel(preferred)
    name = preferred{i};

    if isfield(hitlist,name)
        val = extractScalarTarget(hitlist.(name));

        if ~isempty(val) && isfinite(val)
            ftarget = val;
            hasTarget = true;
            return
        end
    end
end

end

% =====================================================================
function val = extractScalarTarget(S)

val = [];

if isnumeric(S) && isscalar(S) && isfinite(S)
    val = S;
    return
end

if ~isstruct(S)
    return
end

preferred = {'fopt','ftar','ftarget','fmin','fbest','bestf','f','val','value'};

fields = fieldnames(S);

for i = 1:numel(preferred)
    name = preferred{i};

    if ismember(name,fields)
        val = extractScalarTarget(S.(name));
        if ~isempty(val)
            return
        end
    end
end

end