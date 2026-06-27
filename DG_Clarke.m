%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DG_Clarke.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DG_Clarke
%
% Discrete-gradient Clarke-stationarity wrapper with an inner
% DG-MatCSG solver.
%
% This routine applies an outer sequence of stationarity,
% line-search, and discrete-gradient parameters. At each outer
% iteration, it calls DG_MatCSG to compute an approximate stationary
% point using a bundle of discrete gradients, a minimal-norm
% convex-combination residual, and one of two direction modes:
%
%   tune.dir = 1 : normalized steepest discrete-gradient direction,
%   tune.dir = 2 : matrix-CG / MatCSG direction with beta-theta memory.
%
% Usage:
%
%   [x,outAll] = DG_Clarke(fun,x0,tune,ST)
%
% Inputs:
%
%   fun  : function handle. Given x in R^n, fun(x) returns f(x).
%
%   x0   : initial point; numeric vector in R^n.
%
%   tune : algorithmic-parameter structure. Missing fields are
%          completed by initTune. Main fields are listed below.
%
%          Outer parameters:
%             delta0, epsLiS0, lambda0, pwt0
%             outerDecay, warmOuter
%             deltaMin, epsLiSMin, lambdaMin
%             deltaSeq, epsSeq, lambdaSeq
%
%          Inner iteration control:
%             dir, maxIter, maxLiS
%             disableInnerStationarity, storeOuter
%
%          Discrete-gradient and line-search parameters:
%             lambda, pwt, mu1, mu2
%             epsLiS, tLowerLiS, zetaLiS, pLiS, tInitLiS
%
%          MatCSG / matrix-CG parameters:
%             tDL, epsDelta, varrho
%             angleDegTol, corrNormTol
%             clipNegativeBeta, restartAfterNull
%             matEigMinTol, matCondMax
%             betaMax, thetaMax, betaMin, thetaMin
%             useCoeffDecay, coeffDecayPower, coeffDecay
%
%          Bundle and convex-hull QP parameters:
%             qpReg, simplexIter, simplexTol, maxBundle, weightTol
%
%          Diagonal scaling for dir = 2:
%             useBestScaleDir2, bestScaleMemory, bestScaleMethod
%             bestScaleMin, bestScaleMax
%             betaDomFactor, thetaDomFactor
%
%   ST   : stopping and printing structure. Missing fields are
%          completed by initST. Main fields are listed below.
%
%             finit    : initial objective value used in q_f,
%             ftarget  : benchmark target value,
%             accf     : stopping tolerance for q_f,
%             nfmax    : maximum number of function evaluations,
%             secmax   : maximum CPU time,
%             prt      : printing level,
%             nf       : cumulative function-evaluation counter,
%             sec      : elapsed CPU time,
%             qf       : current benchmark ratio,
%             done     : logical stopping flag.
%
% Outputs:
%
%   x      : final point returned by the solver. If a better
%            internally evaluated point is available, that point is
%            returned.
%
%   outAll : output structure, or cell array if tune.storeOuter = true.
%            Main fields are listed below.
%
%             flag, iter, fval
%             xcurrent, fcurrent
%             xbestSolver, fbestSolver
%             V, hist, tune, lastLiS, vhat
%             outerNu, outerDelta, outerEpsLiS
%             outerLambda, outerPwt, outerDir
%
% Local subfunctions:
%
%   checkInput              : input validation.
%   initST                  : stopping and printing defaults.
%   initTune                : algorithmic-parameter defaults.
%   getStopFlagFromST       : stopping-flag conversion.
%
%   DG_MatCSG               : inner discrete-gradient solver.
%   DG_TPLiS                : two-point line search with DG enrichment.
%   DG_DiscreteGradient     : Bagirov-style discrete gradient.
%
%   MatCSG_direction        : direction computation.
%   enforceAngle            : bounded-angle correction.
%   minNormConv             : minimal-norm element of conv(V).
%   projectSimplex          : simplex projection.
%
%   updateBestPointArchive  : archive of best evaluated points.
%   computeBestPointScale   : scaling vector for dir = 2.
%   betaMatrixStats         : tridiagonal matrix diagnostics.
%   randUnit                : random unit-vector generator.
%
%   printIterDynamic        : iteration printing.
%   printMemoryDiagnostics  : memory diagnostics.
%   printLiSTrace           : line-search trace printing.
%   getScalarField          : safe scalar-field extraction.
%
% Nested subfunctions inside DG_TPLiS:
%
%   buildLiSInfo            : line-search diagnostic structure.
%   stopReturn              : early return when ST.done is triggered.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [x,outAll] = DG_Clarke(fun,x,tune,ST)

% check input 
[x,outAll,tune,ST] = checkInput(fun,x,tune,ST);

% Complete missing stopping test parameters and print level
ST  = initST(ST);
prt = ST.prt;


if prt >= 0
    disp(' ')
    disp('==============================================================')
    disp('start DG-MatCSG')
    disp('==============================================================')
end
x = x(:); n = length(x);


% Complete missing tuning parameters
tune = initTune(tune,n);
nu = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% infinite outer Clarke loop %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

while true

    nu = nu + 1;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % choose delta
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if tune.disableInnerStationarity
        tune.delta = -inf;
    elseif ~isempty(tune.deltaSeq)

        if nu <= length(tune.deltaSeq), tune.delta = tune.deltaSeq(nu);
        else, tune.delta = tune.deltaSeq(end);
        end
    else
        if nu <= tune.warmOuter, tune.delta = tune.delta0;
        else
            j           = nu - tune.warmOuter;
            decayFactor = exp(-tune.outerDecay*j);
            tune.delta  = max(tune.deltaMin,tune.delta0*decayFactor);
        end

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % choose epsLiS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ~isempty(tune.epsSeq)
        if nu <= length(tune.epsSeq), tune.epsLiS = tune.epsSeq(nu);
        else, tune.epsLiS = tune.epsSeq(end);
        end
    else
        if nu <= tune.warmOuter, tune.epsLiS = tune.epsLiS0;
        else
            j           = nu - tune.warmOuter;
            decayFactor = exp(-tune.outerDecay * j);
            tune.epsLiS = max(tune.epsLiSMin, tune.epsLiS0 * decayFactor);
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % choose lambda
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ~isempty(tune.lambdaSeq)
        if nu<=length(tune.lambdaSeq), tune.lambda = tune.lambdaSeq(nu);
        else, tune.lambda = tune.lambdaSeq(end);
        end
    else
        if nu <= tune.warmOuter, tune.lambda = tune.lambda0;
        else
            j           = nu - tune.warmOuter;
            decayFactor = exp(-tune.outerDecay * j);
            tune.lambda = max(tune.lambdaMin,tune.lambda0*decayFactor);
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % coordinate perturbation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tune.pwt = tune.pwt0;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % print outer information
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if prt >= 1
        fprintf('\n==================================================\n');
        fprintf('DG_Clarke outer iteration nu = %d\n',nu);

        if isinf(tune.delta) && tune.delta < 0
            fprintf('delta  = -inf  (inner stationarity disabled)\n');
        else
            fprintf('delta  = %.3e\n',tune.delta);
        end

        fprintf('epsLiS = %.3e\n',tune.epsLiS);
        fprintf('lambda = %.3e\n',tune.lambda);
        fprintf('pwt    = %.3e\n',tune.pwt);
        fprintf('dir    = %d\n',tune.dir);
        fprintf('maxIter inner = %d\n',tune.maxIter);
        fprintf('==================================================\n');
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % call inner DG-MatCSG
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    [x,out,ST] = DG_MatCSG(fun,x,tune,ST);
    
    % Set the stopping flag.
    if ST.done
        out.flag = getStopFlagFromST(ST);
    end

    out.outerNu     = nu;
    out.outerDelta  = tune.delta;
    out.outerEpsLiS = tune.epsLiS;
    out.outerLambda = tune.lambda;
    out.outerPwt    = tune.pwt;
    out.outerDir    = tune.dir;
    

    if tune.storeOuter
        outAll{nu,1} = out;
    else
        outAll = out;
    end
    
    if ST.done, break; end
    
end

ST.status_of_converge = getStopFlagFromST(ST);

end

% ============================================================
function flag = getStopFlagFromST(ST)
% getStopFlagFromST
% Converts the external stopping structure into a solver flag.

flag = 'unknown';

if isfield(ST,'qf') && isfield(ST,'accf') && ...
        isfinite(ST.qf) && isfinite(ST.accf) && ST.qf <= ST.accf

    flag = 'accuracy reached';

elseif isfield(ST,'nf') && isfield(ST,'nfmax') && ST.nf >= ST.nfmax

    flag = 'nfmax reached';

elseif isfield(ST,'sec') && isfield(ST,'secmax') && ST.sec >= ST.secmax

    flag = 'secmax reached';

elseif isfield(ST,'done') && ST.done

    flag = 'stopped';

end

end

% ============================================================
function [x,out,ST] = DG_MatCSG(fun,x,tune,ST)
% DG_MatCSG
% Discrete-gradient derivative-free nonsmooth solver.
%
% Direction modes:
%   tune.dir = 1 : steepest descent direction only
%   tune.dir = 2 : matrix CG / MatCSG direction with beta and theta
%

out = struct(); x = x(:); n = length(x); prt = ST.prt;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% printing-level info  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if prt >= 0
    disp(' ')
    disp('printing level:')
    disp('  prt = -1 : silent')
    disp('  prt =  0 : start/finish only')
    disp('  prt =  1 : compact dynamic iteration line')
    disp('  prt =  2 : detailed dynamic diagnostics')
    disp(['  prt =  3 :',... 
          'detailed diagnostics + memory/vector norms + LiS trace'])
end

%%%%%%%%%%%%%%%%%%%%%%%%
%%%% main algorithm %%%%
%%%%%%%%%%%%%%%%%%%%%%%%

% Initial objective value and base discrete-gradient sample.
% IMPORTANT:
%   fcur is the cached value f(x). Do not recompute fun(x) inside the
%   main loop or inside DG_TPLiS.
if ~isfield(ST,'nf') || isempty(ST.nf)
    ST.nf = 0;
end

fcur  = fun(x);
ST.nf = ST.nf+1;
% check stopping test
sec       = (cputime-ST.initTime);
ST.done   = (sec>ST.secmax)|(ST.nf>=ST.nfmax);
qf        = (fcur-ST.ftarget)/(ST.finit-ST.ftarget);
ST.qf     = qf;
ST.done   = (ST.done|ST.qf<=ST.accf);
ST.sec    = sec;
if ST.done, return; end

[vhat,~,bestDG,ST] = ...
   DG_DiscreteGradient(fun,x,tune.lambda,randUnit(n),fcur,tune,ST);
 
if ST.done, return; end

% Opportunistic acceptance of an improving point found during the
% initial discrete-gradient construction.
if bestDG.improved && bestDG.f < fcur
    if prt >= 1
        fprintf(['DG opportunistic initial improvement |', ...
        ' source=%s | f_old=%+.16e | f_new=%+.16e | dec=%+.3e\n'], ...
        bestDG.source, fcur, bestDG.f, bestDG.f - fcur);
    end
  
    x = bestDG.x; fcur = bestDG.f;
    
    % Recompute the base DG at the new accepted point so that vhat and x
    % are consistent. Ignore further opportunistic improvements here to
    % avoid turning initialization into an uncontrolled coordinate search.
    [vhat, ~, ~,ST] = ...
    DG_DiscreteGradient(fun, x, tune.lambda, randUnit(n), fcur, tune,ST);
    if ST.done, return; end
end

V = vhat;

% Memory for CG / MatCSG.
mem       = struct();
mem.dprev = [];
mem.sprev = [];
mem.yprev = [];
mem.sold  = [];
mem.yold  = [];
mem.lastStepWasNull = false;

% Pending accepted step.
pending = struct();
pending.active = false;

hist.f     = [];
hist.normv = [];
hist.dir   = [];
hist.stepType = {};
hist.betaCG   = [];
hist.thetaCG  = [];
hist.Delta = [];
hist.denDL = [];
hist.normd = [];
hist.vTd   = [];
hist.cosvd = [];
hist.lis_t = [];
hist.lis_innerIter = [];
hist.bundleSize    = [];
hist.matEigMin     = [];
hist.matEigMax     = [];
hist.matCondAbs    = [];
hist.usedFallback  = [];
hist.usedRestart   = [];
hist.usedAngleCorrection = [];
hist.usedMatrixStabilityFallback = [];
hist.matrixStabilityReason = {};
hist.betaWasClipped    = [];
hist.betaWasNonfinite  = [];
hist.thetaWasNonfinite = [];
hist.nonfiniteReason   = {};

% DG internal/opportunistic improvement diagnostics
hist.dgInternalImprovement       = [];
hist.dgInternalImprovementValue  = [];
hist.dgInternalImprovementSource = {};

% Extra line-search / safeguard diagnostics
hist.lis_barT = [];
hist.lis_q    = [];
hist.lis_qBeforeSafeguard = [];
hist.lis_usedDirectionalSafeguard   = [];
hist.lis_usedComponentFlipSafeguard = [];
hist.lis_usedFullRestartSafeguard   = [];
hist.lis_numFlippedComponents       = [];

hist.lis_f_t_minus_fx      = [];
hist.lis_f_bar_minus_fx    = [];
hist.lis_armijo_t_margin   = [];
hist.lis_armijo_bar_margin = [];
hist.lis_actualDecBar      = [];
hist.lis_requiredDecBar    = [];

hist.lis_gTd           = [];
hist.lis_enrich_rhs    = [];
hist.lis_enrich_margin = [];

hist.lis_tl = [];
hist.lis_tu = [];
hist.lis_t0Paper = [];
hist.lis_epsLiS  = [];
hist.lis_tLower  = [];
hist.lis_zeta    = [];
hist.lis_pLiS    = [];

hist.lis_localFunEvals = [];
hist.lis_localDGEvals  = [];
hist.lis_localBarEvals = [];
hist.lis_localApproxTotalEvals = [];
hist.lis_forcedByMaxLiS = [];

% New line-search acceptance / forcing diagnostics
hist.lis_acceptedStepSource = {};
hist.lis_qIsZero    = [];
hist.lis_qIsDescent = [];
hist.lis_qZeroTol   = [];
hist.lis_rhoForce   = [];
hist.lis_tTestType  = {};
hist.lis_barTestType = {};


k = 0; out.flag = 'unknown'; lisInfo = [];

xbestSolver = x; fbestSolver = fcur;
% ============================================================
% Best-point archive for dir = 2 scaling.
% Stores up to tune.bestScaleMemory best points and function values.
% ============================================================
bestHistX = xbestSolver;
bestHistF = fbestSolver;
scDir2 = ones(n,1);
while k < tune.maxIter
    
    
    % ============================================================
    % Default diagnostics for this outer iteration.
    % These are stored ONCE near the end of the loop.
    % ============================================================
    dgInternalImprovementFlag = 0;
    dgInternalImprovementValue = 0;
    dgInternalImprovementSource = 'none';
    
    % ============================================================
    % Compute minimal-norm element of conv(V).
    % This vector is used ONLY for the stationarity test.
    % ============================================================
    vstat = minNormConv(V, tune);
    nv    = norm(vstat);
    
    % ============================================================
    % If a serious step was accepted in the previous iteration,
    % update MatCSG memory.
    % ============================================================
    if pending.active
        mem.yold = mem.yprev;
        mem.sold = mem.sprev;
        
        % Update memory consistently with the minimal-norm bundle vector
        % used in the MatCSG direction.
        mem.yprev = vstat - pending.vold;
        mem.sprev = pending.step;
        mem.dprev = pending.d;
        
        pending.active = false;
    end
    
    hist.f(end+1,1)     = fcur;
    hist.normv(end+1,1) = nv;
    hist.dir(end+1,1)   = tune.dir;
    
    % ============================================================
    % Stationarity test based on approximate minimal-norm vector.
    % ============================================================
    if nv <= tune.delta
        out.flag = 'stationary';
        
        if prt >= 1
            fprintf(['k=%5d | f=%+.5e | ||vstat||=%8.2e | ',...
                     'STOP: stationarity\n'], k, fcur, nv);
        end
        
        break;
    end
    
    % ============================================================
    % Reset persistent MatCSG coefficient counter at the first iteration
    % of each solver call.
    % ============================================================
    if k == 0
        tune.resetMatCSGCounter = true;
    else
        tune.resetMatCSGCounter = false;
    end
    
    % ------------------------------------------------------------
    % Build componentwise scaling vector for dir = 2 from 
    % previous best points.
    % sc is based on the median/max absolute difference
    % between the current best point and the other saved best points.
    % ------------------------------------------------------------
    if tune.dir == 2 && tune.useBestScaleDir2
        scDir2 = computeBestPointScale(bestHistX, bestHistF, tune, n);
    else
        scDir2 = ones(n,1);
    end
    
    tune.scDir2 = scDir2;
    
    [d, dinfo] = MatCSG_direction(vstat, mem, tune);
    
    % ============================================================
    % Two-point line search.
    % ============================================================
    xold = x;
    
    % Use the same minimal-norm bundle vector in the line-search model
    % as the one used to generate the search direction.
    vold = vstat;
    
    [s, I, lisInfo,ST] = DG_TPLiS(fun, x, d, vstat, fcur, tune,ST);
    
    if ST.done, return; end
    % ============================================================
    % Capture best point evaluated inside DG_TPLiS.
    %
    % Store the best point evaluated during the line search.
    % ============================================================
    if isfield(lisInfo,'bestPoint') && ...
            isfield(lisInfo.bestPoint,'improved') && ...
            lisInfo.bestPoint.improved && ...
            isfinite(lisInfo.bestPoint.f) && ...
            lisInfo.bestPoint.f < fbestSolver
        
        xbestSolver = lisInfo.bestPoint.x;
        fbestSolver = lisInfo.bestPoint.f;
        
        [bestHistX, bestHistF] = updateBestPointArchive( ...
        bestHistX,bestHistF,xbestSolver,fbestSolver,tune.bestScaleMemory);
        
        if prt >= 1
            fprintf(['DG-TPLiS internal best update | source=%s | ', ...
                'fbestSolver=%+.16e | innerIter=%d\n'], ...
                lisInfo.bestPoint.source, ...
                fbestSolver, ...
                lisInfo.bestPoint.innerIter);
        end
    end
    
    % ============================================================
    % DG_TPLiS may safeguard the direction if vhat'*d is nonnegative
    % or nonfinite. If so, use the actual direction in memory updates
    % and diagnostics.
    % ============================================================
    if isfield(lisInfo,'d')
        d = lisInfo.d;
    end
    
    if isfield(lisInfo,'usedDirectionalSafeguard') && ...
            lisInfo.usedDirectionalSafeguard
        
        if ~isfield(dinfo,'wasLineSearchSafeguarded') ...
                || ~dinfo.wasLineSearchSafeguarded
            dinfo.directionType = ...
                   ['line_search_safeguarded: ', dinfo.directionType];
        end
        
        dinfo.wasLineSearchSafeguarded = true;
        dinfo.usedRestart = true;
        
    else
        
        dinfo.wasLineSearchSafeguarded = false;
        
    end
        
    % ============================================================
    % Prepare print values.
    % ============================================================
    if I == 1
        stepKindPrint = 'serious';
        stepNormPrint = norm(s - x);
        ftrialPrint = lisInfo.ftrial;
    elseif I == 0
        stepKindPrint = 'null';
        stepNormPrint = 0;
        ftrialPrint = lisInfo.ftrial;
    else
        stepKindPrint = 'line_search_failed';
        stepNormPrint = 0;
        ftrialPrint = lisInfo.ftrial;
    end
    
    if prt >= 1
        printIterDynamic(k, fcur, vstat, vhat, d, V, dinfo, lisInfo, ...
            stepKindPrint, stepNormPrint, ftrialPrint, prt);
    end
    
    if prt >= 3
        printMemoryDiagnostics(k, x, V, vhat, d, mem, pending);
    end
    
    % ============================================================
    % Store direction / CG / matrix diagnostics.
    % ============================================================
    hist.normd(end+1,1) = norm(d);
    
    % Direction is based on vstat, so diagnostics should also use vstat.
    hist.vTd(end+1,1) = vstat' * d;
    
    if norm(vstat) > 0 && norm(d) > 0
        hist.cosvd(end+1,1) = (vstat' * d) / (norm(vstat) * norm(d));
    else
        hist.cosvd(end+1,1) = NaN;
    end
    
    hist.betaCG(end+1,1)  = dinfo.betaCG;
    hist.thetaCG(end+1,1) = dinfo.thetaCG;
    hist.Delta(end+1,1)   = dinfo.Delta;
    hist.denDL(end+1,1)   = dinfo.denDL;
    hist.bundleSize(end+1,1) = size(V,2);
    
    hist.matEigMin(end+1,1) = dinfo.matEigMin;
    hist.matEigMax(end+1,1) = dinfo.matEigMax;
    hist.matCondAbs(end+1,1) = dinfo.matCondAbs;
    
    hist.usedFallback(end+1,1) = dinfo.usedFallback;
    hist.usedRestart(end+1,1) = dinfo.usedRestart;
    hist.usedAngleCorrection(end+1,1) = dinfo.usedAngleCorrection;
    hist.usedMatrixStabilityFallback(end+1,1) = ...
                                 dinfo.usedMatrixStabilityFallback;
    
    hist.matrixStabilityReason{end+1,1} = dinfo.matrixStabilityReason;
    
    hist.betaWasClipped(end+1,1) = dinfo.betaWasClipped;
    hist.betaWasNonfinite(end+1,1) = dinfo.betaWasNonfinite;
    hist.thetaWasNonfinite(end+1,1) = dinfo.thetaWasNonfinite;
    hist.nonfiniteReason{end+1,1} = dinfo.nonfiniteReason;
    
    % ============================================================
    % Store line-search diagnostics.
    % ============================================================
    if isfield(lisInfo,'t')
        hist.lis_t(end+1,1) = lisInfo.t;
    else
        hist.lis_t(end+1,1) = NaN;
    end
    
    if isfield(lisInfo,'innerIter')
        hist.lis_innerIter(end+1,1) = lisInfo.innerIter;
    else
        hist.lis_innerIter(end+1,1) = NaN;
    end
    
    hist.lis_barT(end+1,1) = ...
        getScalarField(lisInfo,'barT',NaN);
    hist.lis_q(end+1,1) = ...
        getScalarField(lisInfo,'q',NaN);
    hist.lis_qBeforeSafeguard(end+1,1) = ...
        getScalarField(lisInfo,'qBeforeSafeguard',NaN);
    
    hist.lis_usedDirectionalSafeguard(end+1,1) = ...
        getScalarField(lisInfo,'usedDirectionalSafeguard',NaN);
    hist.lis_usedComponentFlipSafeguard(end+1,1) = ...
        getScalarField(lisInfo,'usedComponentFlipSafeguard',NaN);
    hist.lis_usedFullRestartSafeguard(end+1,1) = ...
        getScalarField(lisInfo,'usedFullRestartSafeguard',NaN);
    hist.lis_numFlippedComponents(end+1,1) = ...
        getScalarField(lisInfo,'numFlippedComponents',NaN);
    
    hist.lis_f_t_minus_fx(end+1,1) = ...
        getScalarField(lisInfo,'f_t_minus_fx',NaN);
    hist.lis_f_bar_minus_fx(end+1,1) = ...
        getScalarField(lisInfo,'f_bar_minus_fx',NaN);
    hist.lis_armijo_t_margin(end+1,1) = ...
        getScalarField(lisInfo,'armijo_t_margin',NaN);
    hist.lis_armijo_bar_margin(end+1,1) = ...
        getScalarField(lisInfo,'armijo_bar_margin',NaN);
    hist.lis_actualDecBar(end+1,1) = ...
        getScalarField(lisInfo,'actualDecBar',NaN);
    hist.lis_requiredDecBar(end+1,1) = ...
        getScalarField(lisInfo,'requiredDecBar',NaN);
    
    hist.lis_gTd(end+1,1) = ...
        getScalarField(lisInfo,'gTd',NaN);
    hist.lis_enrich_rhs(end+1,1) = ...
        getScalarField(lisInfo,'enrich_rhs',NaN);
    hist.lis_enrich_margin(end+1,1) = ...
        getScalarField(lisInfo,'enrich_margin',NaN);
    
    hist.lis_tl(end+1,1) = getScalarField(lisInfo,'tl',NaN);
    hist.lis_tu(end+1,1) = getScalarField(lisInfo,'tu',NaN);
    hist.lis_t0Paper(end+1,1) = getScalarField(lisInfo,'t0Paper',NaN);
    hist.lis_epsLiS(end+1,1) = getScalarField(lisInfo,'epsLiS',NaN);
    hist.lis_tLower(end+1,1) = getScalarField(lisInfo,'tLower',NaN);
    hist.lis_zeta(end+1,1) = getScalarField(lisInfo,'zeta',NaN);
    hist.lis_pLiS(end+1,1) = getScalarField(lisInfo,'pLiS',NaN);
    
    hist.lis_localFunEvals(end+1,1) = ...
        getScalarField(lisInfo,'localFunEvals',NaN);
    hist.lis_localDGEvals(end+1,1) = ...
        getScalarField(lisInfo,'localDGEvals',NaN);
    hist.lis_localBarEvals(end+1,1) = ...
        getScalarField(lisInfo,'localBarEvals',NaN);
    hist.lis_localApproxTotalEvals(end+1,1) =...
        getScalarField(lisInfo,'localApproxTotalEvals',NaN);
    hist.lis_forcedByMaxLiS(end+1,1) = ...
        getScalarField(lisInfo,'forcedByMaxLiS',NaN);
    
    % Line-search acceptance and forcing diagnostics
    if isfield(lisInfo,'acceptedStepSource')
        hist.lis_acceptedStepSource{end+1,1} = ...
                                      lisInfo.acceptedStepSource;
    else
        hist.lis_acceptedStepSource{end+1,1} = 'unknown';
    end
    
    hist.lis_qIsZero(end+1,1) = ...
        getScalarField(lisInfo,'qIsZero',NaN);
    hist.lis_qIsDescent(end+1,1) = ...
        getScalarField(lisInfo,'qIsDescent',NaN);
    hist.lis_qZeroTol(end+1,1) = ...
        getScalarField(lisInfo,'qZeroTol',NaN);
    hist.lis_rhoForce(end+1,1) = ...
        getScalarField(lisInfo,'rhoForce',NaN);
    
    if isfield(lisInfo,'tTestType')
        hist.lis_tTestType{end+1,1} = lisInfo.tTestType;
    else
        hist.lis_tTestType{end+1,1} = 'unknown';
    end
    
    if isfield(lisInfo,'barTestType')
        hist.lis_barTestType{end+1,1} = lisInfo.barTestType;
    else
        hist.lis_barTestType{end+1,1} = 'unknown';
    end
    
    % ============================================================
    % Serious/null/failed update.
    % ============================================================
    if I == 1
        
        % --------------------------------------------------------
        % Serious step.
        % --------------------------------------------------------
        x = s;
        
        % Reuse accepted trial value from DG_TPLiS.
        fcur = lisInfo.ftrial;
        
        % --------------------------------------------------------
        % Reset bundle at the new point.
        % --------------------------------------------------------
        [gnew, ~, bestDG,ST] = ...
            DG_DiscreteGradient(fun, x, tune.lambda, d, fcur, tune, ST);
        
        if ST.done, return; end

        
        if prt >= 3
          fprintf(['post-serious DG check | fcur=%+.16e | ',...
                  'bestDG.f=%+.16e | ','improved=%d | source=%s\n'],...
                   fcur, bestDG.f, bestDG.improved, bestDG.source);
        end
        
        if bestDG.improved && bestDG.f < fcur
            
            fOldDG = fcur;
            fNewDG = bestDG.f;
            
            dgInternalImprovementFlag = 1;
            dgInternalImprovementValue = fNewDG - fOldDG;
            dgInternalImprovementSource = bestDG.source;
            
            if prt >= 1
                fprintf(['DG opportunistic post-serious improvement ',...
                    '| source=%s | f_old=%+.16e | f_new=%+.16e',...
                    ' | dec=%+.3e\n'],bestDG.source, fOldDG, ...
                    fNewDG, fNewDG - fOldDG);
            end
            
            % Accept the best point found internally
            % during DG construction.
            x = bestDG.x;
            fcur = bestDG.f;
            
            % Recompute base DG at the accepted improved point.
            [gnew, ~, ~,ST] = ...
                 DG_DiscreteGradient(fun,x,tune.lambda,d,fcur,tune,ST);
            if ST.done, return; end
        end
        
        vhat = gnew; V = vhat;
        
        % Update solver-level best.
        if fcur < fbestSolver
            xbestSolver = x;
            fbestSolver = fcur;
            
            [bestHistX, bestHistF] = updateBestPointArchive( ...
                bestHistX, bestHistF, xbestSolver, fbestSolver, ...
                tune.bestScaleMemory);
        end
        
        % Store pending memory update for next iteration.
        %
        % Use the previous minimal-norm bundle vector
        % in the memory update
        
        pending.active = true;
        pending.vold   = vold;
        actualStep     = x - xold;
        pending.step   = actualStep;
        
        if norm(actualStep) > 0 && all(isfinite(actualStep))
            pending.d = actualStep / norm(actualStep);
        else
            pending.d = d;
        end
        
        mem.lastStepWasNull = false;
        
        hist.stepType{end+1,1} = 'serious';
        
    elseif I == 0
        
        % --------------------------------------------------------
        % Null / enrichment step.
        % --------------------------------------------------------
        
        % Current x and fcur remain unchanged.
        % The returned discrete-gradient vector becomes 
        % the new base vector.
        vhat = s;
        V = [V, vhat];
        
        if size(V,2) > tune.maxBundle
            V = V(:, end-tune.maxBundle+1:end);
        end
        
        mem.lastStepWasNull = true;
        
        hist.stepType{end+1,1} = 'null';
        
    else
        
        % --------------------------------------------------------
        % Failed line-search step.
        %
        % Keep the current point, reset the bundle and memory, and do not
        % add an invalid enrichment vector.
        % --------------------------------------------------------
        
        if prt >= 1
            fprintf(['DG_TPLiS failed: reset bundle/memory',...
                     ' without adding enrichment.\n']);
        end
        
        [vhat, ~, ~,ST] = DG_DiscreteGradient(fun,x,tune.lambda, ...
                       randUnit(length(x)), fcur, tune,ST);
        if ST.done, return; end

        V = vhat;
        
        mem.dprev = [];
        mem.sprev = [];
        mem.yprev = [];
        mem.sold  = [];
        mem.yold  = [];
        mem.lastStepWasNull = false;
        
        pending.active = false;
        
        hist.stepType{end+1,1} = 'line_search_failed';
        
    end
    % ============================================================
    % Store DG internal/opportunistic improvement diagnostics.
    %
    % Store diagnostics after the iterate update, since serious steps
    % may modify the internal-improvement flag.
    % ============================================================
    hist.dgInternalImprovement(end+1,1) = ...
        dgInternalImprovementFlag;
    hist.dgInternalImprovementValue(end+1,1)= ...
        dgInternalImprovementValue;
    hist.dgInternalImprovementSource{end+1,1}= ...
         dgInternalImprovementSource;
    
    % ============================================================
    % Advance outer loop.
    % ============================================================
    k = k + 1;
    
end
if k >= tune.maxIter && strcmp(out.flag,'unknown')
    out.flag = 'maxIter';
end

out.iter = k;

% Preserve current iterate information.
out.xcurrent = x;
out.fcurrent = fcur;

% Preserve best internally evaluated point.
out.xbestSolver = xbestSolver;
out.fbestSolver = fbestSolver;
out.bestHistX = bestHistX;
out.bestHistF = bestHistF;
out.lastScDir2 = scDir2;

% Return the best internally evaluated point if it is better than current.
if isfinite(fbestSolver) && fbestSolver < fcur
    
    if prt >= 1
        fprintf(['Returning best internally evaluated point | ', ...
            'fcurrent=%+.16e | fbestSolver=%+.16e | dec=%+.3e\n'], ...
            fcur, fbestSolver, fbestSolver - fcur);
    end
    
    x = xbestSolver;
    out.fval = fbestSolver;
    
else
    
    out.fval = fcur;
    
end

out.V = V;
out.hist = hist;
out.tune = tune;
out.lastLiS = lisInfo;
out.vhat = vhat;

if prt >= 0
    disp('==============================================================')
    disp('finish DG-MatCSG')
    disp('==============================================================')
end

end

% ============================================================
function [s, I, info,ST] = DG_TPLiS(fun, x, d, vbase, fx, tune,ST)
% DG_TPLiS
% Paper-style two-point line-search with discrete gradients.
%
% Line-search logic:
%
%   1) If Armijo/forcing succeeds at bar_t_i, accept bar_t_i.
%
%   2) If bar_t_i fails but Armijo/forcing succeeds at t_i,
%      accept x + t_i*d as a serious step BEFORE allowing enrichment.
%
%   3) If q = vbase'*d is numerically zero, use a forcing function:
%
%          f(x+t*d) - f(x) <= -rho*t^2
%
%      instead of the degenerate condition
%
%          f(x+t*d) - f(x) <= mu*t*q = 0.
%
% If q < 0, the usual Armijo term mu*t*q is used.
% If q > 0, component-wise sign repair is attempted first,
%           then full restart.

x = x(:);
d = d(:);
vbase = vbase(:);

n = length(x);

% ------------------------------------------------------------
% Parameters from the paper
% ------------------------------------------------------------
epsLiS = tune.epsLiS;
tLower = tune.tLowerLiS;
zeta   = tune.zetaLiS;
pLiS   = tune.pLiS;

% Optional forcing coefficient for q approximately zero:
% condition is f(x+t*d)-f(x) <= -rhoForce*t^2.
if isfield(tune,'forceCoeffLiS')
    rhoForce = tune.forceCoeffLiS;
else
    rhoForce = 1.0;
end

rhoForce = max(rhoForce, eps);

% Optional zero tolerance for q.
if isfield(tune,'qZeroTolLiS')
    qZeroTolUser = tune.qZeroTolLiS;
else
    qZeroTolUser = 1e-12;
end

% Defensive projection into valid ranges.
epsLiS = min(max(epsLiS, 10*eps), 1 - 10*eps);
tLower = min(max(tLower, 10*eps), epsLiS - 10*eps);
zeta   = min(max(zeta, 10*eps), 0.5 - 10*eps);
pLiS   = max(1, round(pLiS));

% Choose t_0 in (underline t, epsilon).
if isfield(tune,'tInitLiS') && ~isempty(tune.tInitLiS)
    t0 = tune.tInitLiS;
else
    t0 = sqrt(tLower * epsLiS);
end

t0 = min(max(t0, tLower + 10*eps), epsLiS - 10*eps);

% ------------------------------------------------------------
% Directional model derivative q = vbase' * d.
%
% Important:
%   - If q > 0, repair/restart.
%   - If q is approximately zero, use the forcing function.
%     -rho*t^2 in the sufficient-decrease test.
% ------------------------------------------------------------
q = vbase' * d;

usedDirectionalSafeguard = false;
usedComponentFlipSafeguard = false;
usedFullRestartSafeguard = false;

qBeforeSafeguard = q;
flipIdx = false(size(d));

qScale = max(1, norm(vbase) * norm(d));
qZeroTol = qZeroTolUser * qScale;

% ------------------------------------------------------------
% Safeguard 1:
% If q > 0, first try component-wise sign repair:
%
%   if vbase_i*d_i > 0, set d_i <- -d_i.
% ------------------------------------------------------------
if isfinite(q) && q > qZeroTol

    flipIdx = (vbase .* d) > 0;

    if any(flipIdx)
        d(flipIdx) = -d(flipIdx);
        q = vbase' * d;

        usedDirectionalSafeguard = true;
        usedComponentFlipSafeguard = true;
    end
end

% ------------------------------------------------------------
% Safeguard 2:
% If q is still positive or nonfinite, restart fully to -vbase.
%
% Note:
%   q approximately zero is NOT restarted. It is handled by the forcing
%   function in the line-search tests.
% ------------------------------------------------------------
if ~isfinite(q) || q > qZeroTol

    nvbase = norm(vbase);

    if nvbase > 0 && isfinite(nvbase)
        d = -vbase / nvbase;
    else
        d = -vbase / max(nvbase, eps);
    end

    q = vbase' * d;

    usedDirectionalSafeguard = true;
    usedFullRestartSafeguard = true;
end

% Classify directional derivative after safeguard.
qIsZero = isfinite(q) && abs(q) <= qZeroTol;
qIsDescent = isfinite(q) && q < -qZeroTol;

% ------------------------------------------------------------
% Local evaluation counters.
% These do not affect the external fginfo counter; they only diagnose cost.
% One DG call costs approximately n function evaluations because f(base)
% is already passed into DG_DiscreteGradient.
% ------------------------------------------------------------
localFunEvals = 0;
localDGEvals = 0;
localBarEvals = 0;

% ------------------------------------------------------------
% Best point evaluated inside this line search.
% This does not change the line-search decision. It only protects the
% solver from losing a good sampled point if the external nfmax stops
% during the current line-search call.
% ------------------------------------------------------------
bestPoint = struct();
bestPoint.x = x;
bestPoint.f = fx;
bestPoint.improved = false;
bestPoint.source = 'base';
bestPoint.innerIter = -1;

% ------------------------------------------------------------
% Default line-search outputs.
% ------------------------------------------------------------
s = x;
I = -1;

info = struct();
info.type = 'default_line_search_output';
info.innerIter = -1;
info.t = NaN;
info.barT = NaN;
info.ftrial = fx;
info.f_t = NaN;
info.f_bar = NaN;
info.acceptedStepSource = 'none';
info.forcedByMaxLiS = false;
info.bestPoint = bestPoint;

% ------------------------------------------------------------
% Inner trace arrays
% ------------------------------------------------------------
trace_i = [];
trace_t = [];
trace_barT = [];
trace_tl = [];
trace_tu = [];
trace_tl_next = [];
trace_tu_next = [];

trace_f_t_minus_fx = [];
trace_f_bar_minus_fx = [];
trace_armijo_t_rhs = [];
trace_armijo_bar_rhs = [];
trace_armijo_t_margin = [];
trace_armijo_bar_margin = [];

trace_actualDecBar = [];
trace_requiredDecBar = [];

trace_gTd = [];
trace_enrich_rhs = [];
trace_enrich_margin = [];

trace_armijo_t_pass = [];
trace_armijo_bar_pass = [];
trace_enrich_pass = [];
trace_bar_eval_reused = [];

trace_t_test_type = {};
trace_bar_test_type = {};

% ------------------------------------------------------------
% Safe default scalar diagnostics.
% These are overwritten inside the while-loop before buildLiSInfo()
% is used in ordinary serious/null/failed decisions.
% ------------------------------------------------------------
armijo_t_lhs = NaN;
armijo_t_rhs = NaN;
armijo_t_margin = NaN;
armijo_t_pass = false;

armijo_bar_lhs = NaN;
armijo_bar_rhs = NaN;
armijo_bar_margin = NaN;
armijo_bar_pass = false;

actualDecBar = NaN;
requiredDecBar = NaN;

gTd = NaN;
enrich_rhs = NaN;
enrich_margin = NaN;
enrich_pass = false;

tTestType = 'undefined';
barTestType = 'undefined';

f_t = NaN;
f_bar = NaN;
barEvalReused = false;

% ------------------------------------------------------------
% Initialize paper variables:
%   t_0^l = 0, t_0^u = epsilon, bar_t_0 = 1, i = 0.
% ------------------------------------------------------------
tl = 0;
tu = epsLiS;

t = t0;
barT = 1.0;
i = 0;

% Initial enrichment discrete gradient g_0 at x + t_0*d.
y_t = x + t*d;
f_t = fun(y_t);
localFunEvals = localFunEvals + 1;
ST.nf = ST.nf + 1;

if isfinite(f_t) && f_t < bestPoint.f
    bestPoint.x = y_t;
    bestPoint.f = f_t;
    bestPoint.improved = true;
    bestPoint.source = 't_point';
    bestPoint.innerIter = i;
end

% check stopping test
sec     = (cputime - ST.initTime);
ST.done = (sec > ST.secmax) | (ST.nf >= ST.nfmax);
qf      = (f_t - ST.ftarget) / (ST.finit - ST.ftarget);
ST.qf   = qf;
ST.done = (ST.done | ST.qf <= ST.accf);
ST.sec  = sec;

if ST.done
    [s,I,info] = stopReturn('stopped_after_initial_t_eval', ...
                            f_t, NaN, i, t, barT);
    return
end

dDG = randUnit(n);
[g,~,bestDG0,ST] = ...
    DG_DiscreteGradient(fun,y_t,tune.lambda,dDG,f_t,tune,ST);
localDGEvals = localDGEvals + 1;

% Capture best point found inside the DG construction at y_t.
if exist('bestDG0','var') && isstruct(bestDG0) && ...
        isfield(bestDG0,'improved') && bestDG0.improved && ...
        isfield(bestDG0,'f') && isfinite(bestDG0.f) && ...
        bestDG0.f < bestPoint.f

    bestPoint.x = bestDG0.x;
    bestPoint.f = bestDG0.f;
    bestPoint.improved = true;
    bestPoint.source = ['DG_initial_', bestDG0.source];
    bestPoint.innerIter = i;
end

if ST.done
    [s,I,info] = stopReturn('stopped_after_initial_DG', ...
                            f_t, NaN, i, t, barT);
    return
end

lastG = g;
lastFt = f_t;
lastT = t;
lastBarT = barT;
lastI = i;

while true

    % --------------------------------------------------------
    % Sufficient-decrease model at t_i.
    %
    % If q < 0:
    %   f(x+t*d)-f(x) <= mu1*t*q.
    %
    % If q approximately zero:
    %   f(x+t*d)-f(x) <= -rhoForce*t^2.
    % --------------------------------------------------------
    armijo_t_lhs = f_t - fx;

    if qIsDescent
        armijo_t_rhs = tune.mu1 * t * q;
        tTestType = 'Armijo-q';
    else
        armijo_t_rhs = -rhoForce * t^2;
        tTestType = 'forcing-t2';
    end

    armijo_t_margin = armijo_t_lhs - armijo_t_rhs;
    armijo_t_pass = armijo_t_lhs <= armijo_t_rhs;

    if armijo_t_pass
        tl_next = t;
        tu_next = tu;
    else
        tl_next = tl;
        tu_next = t;
    end

    % --------------------------------------------------------
    % Serious-step test at bar_t_i.
    % Avoid duplicate evaluation if barT == t.
    % --------------------------------------------------------
    if abs(barT - t) <= 10*eps*max(1, abs(t))
        y_bar = y_t;
        f_bar = f_t;
        barEvalReused = true;
    else
        y_bar = x + barT*d;
        f_bar = fun(y_bar);

        ST.nf = ST.nf + 1;
        localFunEvals = localFunEvals + 1;
        localBarEvals = localBarEvals + 1;
        barEvalReused = false;

        if isfinite(f_bar) && f_bar < bestPoint.f
            bestPoint.x = y_bar;
            bestPoint.f = f_bar;
            bestPoint.improved = true;
            bestPoint.source = 'barT_point';
            bestPoint.innerIter = i;
        end

        % check stopping test
        sec     = (cputime - ST.initTime);
        ST.done = (sec > ST.secmax) | (ST.nf >= ST.nfmax);
        qf      = (f_bar - ST.ftarget) / (ST.finit - ST.ftarget);
        ST.qf   = qf;
        ST.done = (ST.done | ST.qf <= ST.accf);
        ST.sec  = sec;

        if ST.done
            [s,I,info] = stopReturn('stopped_after_barT_eval', ...
                                    f_t, f_bar, i, t, barT);
            return
        end
    end

    if isfinite(f_bar) && f_bar < bestPoint.f
        bestPoint.x = y_bar;
        bestPoint.f = f_bar;
        bestPoint.improved = true;
        bestPoint.source = 'barT_point';
        bestPoint.innerIter = i;
    end

    armijo_bar_lhs = f_bar - fx;

    if qIsDescent
        armijo_bar_rhs = tune.mu1 * barT * q;
        barTestType = 'Armijo-q';
    else
        armijo_bar_rhs = -rhoForce * barT^2;
        barTestType = 'forcing-t2';
    end

    armijo_bar_margin = armijo_bar_lhs - armijo_bar_rhs;

    actualDecBar = fx - f_bar;
    requiredDecBar = -armijo_bar_rhs;

    armijo_bar_pass = ...
          (armijo_bar_lhs <= armijo_bar_rhs) && (barT >= tLower);

    % --------------------------------------------------------
    % Enrichment test using g_i computed at x + t_i*d.
    %
    % For q approximately zero, keep the enrichment threshold at zero:
    %   g_i^T d >= 0.
    %
    % This avoids mixing derivative units with function-decrease units.
    % The forcing function is used only for accepting serious steps.
    % --------------------------------------------------------
    gTd = g' * d;

    if qIsDescent
        enrich_rhs = tune.mu2 * q;
    else
        enrich_rhs = 0;
    end

    enrich_margin = gTd - enrich_rhs;
    enrich_pass = gTd >= enrich_rhs;

    % --------------------------------------------------------
    % Store trace for this inner iteration.
    % --------------------------------------------------------
    trace_i(end+1,1) = i;
    trace_t(end+1,1) = t;
    trace_barT(end+1,1) = barT;
    trace_tl(end+1,1) = tl;
    trace_tu(end+1,1) = tu;
    trace_tl_next(end+1,1) = tl_next;
    trace_tu_next(end+1,1) = tu_next;

    trace_f_t_minus_fx(end+1,1) = armijo_t_lhs;
    trace_f_bar_minus_fx(end+1,1) = armijo_bar_lhs;
    trace_armijo_t_rhs(end+1,1) = armijo_t_rhs;
    trace_armijo_bar_rhs(end+1,1) = armijo_bar_rhs;
    trace_armijo_t_margin(end+1,1) = armijo_t_margin;
    trace_armijo_bar_margin(end+1,1) = armijo_bar_margin;

    trace_actualDecBar(end+1,1) = actualDecBar;
    trace_requiredDecBar(end+1,1) = requiredDecBar;

    trace_gTd(end+1,1) = gTd;
    trace_enrich_rhs(end+1,1) = enrich_rhs;
    trace_enrich_margin(end+1,1) = enrich_margin;

    trace_armijo_t_pass(end+1,1) = armijo_t_pass;
    trace_armijo_bar_pass(end+1,1) = armijo_bar_pass;
    trace_enrich_pass(end+1,1) = enrich_pass;
    trace_bar_eval_reused(end+1,1) = barEvalReused;

    trace_t_test_type{end+1,1} = tTestType;
    trace_bar_test_type{end+1,1} = barTestType;

    % --------------------------------------------------------
    % Serious-step acceptance priority.
    %
    % 1) First accept bar_t_i if it succeeds.
    % 2) If bar_t_i fails but t_i succeeds, accept t_i.
    % 3) Only then allow enrichment/null step.
    %
    % Serious-step candidates are tested before enrichment so that an
    % already evaluated decreasing point is not discarded.
    % --------------------------------------------------------

    % --------------------------------------------------------
    % Serious-step acceptance.
    %
    % If both barT and t_i satisfy sufficient decrease, accept the one
    % with smaller objective value. This uses only points already evaluated.
    % --------------------------------------------------------
    t_pass_serious = armijo_t_pass && (t >= tLower);
    bar_pass_serious = armijo_bar_pass;

    if bar_pass_serious || t_pass_serious

        if bar_pass_serious && t_pass_serious

            % Both are valid serious-step candidates.
            % Choose the better objective value.
            if f_t < f_bar
                acceptSource = 't_i';
                acceptPoint = x + t*d;
                acceptF = f_t;
            else
                acceptSource = 'barT';
                acceptPoint = x + barT*d;
                acceptF = f_bar;
            end

        elseif bar_pass_serious

            acceptSource = 'barT';
            acceptPoint = x + barT*d;
            acceptF = f_bar;

        else

            acceptSource = 't_i';
            acceptPoint = x + t*d;
            acceptF = f_t;

        end

        s = acceptPoint;
        I = 1;

        info = buildLiSInfo();
        info.innerIter = i;
        info.type = ['successful step at ', acceptSource];
        info.t = t;
        info.barT = barT;
        info.ftrial = acceptF;
        info.f_t = f_t;
        info.f_bar = f_bar;
        info.acceptedStepSource = acceptSource;
        info.forcedByMaxLiS = false;

        return
    end

    % --------------------------------------------------------
    % Enrichment/null step.
    % Check enrichment after serious-step tests.
    % --------------------------------------------------------
    if enrich_pass

        s = g;
        I = 0;

        info = buildLiSInfo();
        info.innerIter = i;
        info.type = 'enrichment';
        info.t = t;
        info.barT = barT;
        info.ftrial = f_t;
        info.f_t = f_t;
        info.f_bar = f_bar;
        info.acceptedStepSource = 'none';
        info.forcedByMaxLiS = false;

        return
    end

    % --------------------------------------------------------
    % Forced enrichment cap for implementation.
    % Cap the number of inner line-search iterations.
    % --------------------------------------------------------
    if i >= tune.maxLiS

        % At the line-search cap, accept the best improved point
        % if one was found; otherwise report failure without
        % adding an enrichment vector.

        if bestPoint.improved && isfinite(bestPoint.f) ...
                              && bestPoint.f < fx

            s = bestPoint.x;
            I = 1;

            info = buildLiSInfo();
            info.innerIter = lastI;
            info.type = 'bestPoint fallback serious';
            info.t = lastT;
            info.barT = lastBarT;
            info.ftrial = bestPoint.f;
            info.f_t = f_t;
            info.f_bar = f_bar;
            info.acceptedStepSource = 'bestPoint';
            info.forcedByMaxLiS = true;

            return

        else

            s = x;
            I = -1;

            info = buildLiSInfo();
            info.innerIter = lastI;
            info.type = 'line search failed';
            info.t = lastT;
            info.barT = lastBarT;
            info.ftrial = fx;
            info.f_t = f_t;
            info.f_bar = f_bar;
            info.acceptedStepSource = 'none';
            info.forcedByMaxLiS = true;

            return

        end

    end

    gap = tu_next - tl_next;

    a = tl_next + zeta * gap;
    b = tu_next - zeta * gap;

    t_next = sqrt(max(a*b, 0));

    % Numerical safety.
    t_next = min(max(t_next, 10*eps), epsLiS - 10*eps);

    i_next = i + 1;
    barT_next = exp((i_next / pLiS) * log(t0));

    % --------------------------------------------------------
    % Select new DG direction d^{(i+1)} in D and compute g_{i+1}.
    % --------------------------------------------------------
    dDG = randUnit(n);

    y_t_next = x + t_next*d;
    f_t_next = fun(y_t_next);
    localFunEvals = localFunEvals + 1;

    ST.nf = ST.nf + 1;

    if isfinite(f_t_next) && f_t_next < bestPoint.f
        bestPoint.x = y_t_next;
        bestPoint.f = f_t_next;
        bestPoint.improved = true;
        bestPoint.source = 't_next_point';
        bestPoint.innerIter = i_next;
    end

    % check stopping test
    sec     = (cputime - ST.initTime);
    ST.done = (sec > ST.secmax) | (ST.nf >= ST.nfmax);
    qf      = (f_t_next - ST.ftarget) / (ST.finit - ST.ftarget);
    ST.qf   = qf;
    ST.done = (ST.done | ST.qf <= ST.accf);
    ST.sec  = sec;

    if ST.done
        [s,I,info] = stopReturn('stopped_after_t_next_eval', ...
                                f_t_next, NaN, i_next, t_next, barT_next);
        return
    end

    [g_next, ~, bestDGNext,ST] = ...
    DG_DiscreteGradient(fun, y_t_next, tune.lambda, dDG, f_t_next, tune,ST);
    localDGEvals = localDGEvals + 1;

    % Capture best point found inside the DG construction at y_t_next.
    if exist('bestDGNext','var') && isstruct(bestDGNext) && ...
            isfield(bestDGNext,'improved') && bestDGNext.improved && ...
            isfield(bestDGNext,'f') && isfinite(bestDGNext.f) && ...
            bestDGNext.f < bestPoint.f

        bestPoint.x = bestDGNext.x;
        bestPoint.f = bestDGNext.f;
        bestPoint.improved = true;
        bestPoint.source = ['DG_next_', bestDGNext.source];
        bestPoint.innerIter = i_next;
    end

    if ST.done
        [s,I,info] = stopReturn('stopped_after_next_DG', ...
                                f_t_next, NaN, i_next, t_next, barT_next);
        return
    end

    % Save candidate for forced enrichment.
    lastG = g_next;
    lastFt = f_t_next;
    lastT = t_next;
    lastBarT = barT_next;
    lastI = i_next;

    % Advance loop.
    tl   = tl_next;
    tu   = tu_next;
    t    = t_next;
    barT = barT_next;
    f_t  = f_t_next;
    g    = g_next;
    i    = i_next;
    y_t  = y_t_next;

end

% ------------------------------------------------------------
% Nested helper: collect all diagnostic fields consistently.
% ------------------------------------------------------------
function info = buildLiSInfo()

info = struct();

info.q = q;
info.qBeforeSafeguard = qBeforeSafeguard;
info.qZeroTol = qZeroTol;
info.qIsZero = qIsZero;
info.qIsDescent = qIsDescent;
info.rhoForce = rhoForce;
info.d = d;

info.bestPoint = bestPoint;

info.usedDirectionalSafeguard = usedDirectionalSafeguard;
info.usedComponentFlipSafeguard = usedComponentFlipSafeguard;
info.usedFullRestartSafeguard = usedFullRestartSafeguard;
info.numFlippedComponents = nnz(flipIdx);
info.flipIdx = flipIdx;

info.tl = tl;
info.tu = tu;
info.t0Paper = t0;
info.epsLiS = epsLiS;
info.tLower = tLower;
info.zeta = zeta;
info.pLiS = pLiS;

info.f_t_minus_fx = armijo_t_lhs;
info.f_bar_minus_fx = armijo_bar_lhs;

info.armijo_t_rhs = armijo_t_rhs;
info.armijo_bar_rhs = armijo_bar_rhs;
info.armijo_t_margin = armijo_t_margin;
info.armijo_bar_margin = armijo_bar_margin;

info.actualDecBar = actualDecBar;
info.requiredDecBar = requiredDecBar;

info.gTd = gTd;
info.enrich_rhs = enrich_rhs;
info.enrich_margin = enrich_margin;

info.armijo_t_pass = armijo_t_pass;
info.armijo_bar_pass = armijo_bar_pass;
info.enrich_pass = enrich_pass;

info.tTestType = tTestType;
info.barTestType = barTestType;

info.localFunEvals = localFunEvals;
info.localDGEvals = localDGEvals;
info.localBarEvals = localBarEvals;
info.localApproxTotalEvals = localFunEvals + n*localDGEvals;

info.trace.i = trace_i;
info.trace.t = trace_t;
info.trace.barT = trace_barT;
info.trace.tl = trace_tl;
info.trace.tu = trace_tu;
info.trace.tl_next = trace_tl_next;
info.trace.tu_next = trace_tu_next;

info.trace.f_t_minus_fx = trace_f_t_minus_fx;
info.trace.f_bar_minus_fx = trace_f_bar_minus_fx;
info.trace.armijo_t_rhs = trace_armijo_t_rhs;
info.trace.armijo_bar_rhs = trace_armijo_bar_rhs;
info.trace.armijo_t_margin = trace_armijo_t_margin;
info.trace.armijo_bar_margin = trace_armijo_bar_margin;

info.trace.actualDecBar = trace_actualDecBar;
info.trace.requiredDecBar = trace_requiredDecBar;

info.trace.gTd = trace_gTd;
info.trace.enrich_rhs = trace_enrich_rhs;
info.trace.enrich_margin = trace_enrich_margin;

info.trace.armijo_t_pass = trace_armijo_t_pass;
info.trace.armijo_bar_pass = trace_armijo_bar_pass;
info.trace.enrich_pass = trace_enrich_pass;
info.trace.bar_eval_reused = trace_bar_eval_reused;

info.trace.t_test_type = trace_t_test_type;
info.trace.bar_test_type = trace_bar_test_type;

end

% ------------------------------------------------------------
% Nested helper: safe early return when ST.done is triggered
% before a normal line-search decision is reached.
% ------------------------------------------------------------
function [sStop,IStop,infoStop] = stopReturn(stopType, fTlocal, fBarLocal, ...
                                             iterLocal, tLocal, barTLocal)

if bestPoint.improved && isfinite(bestPoint.f) && bestPoint.f < fx
    sStop = bestPoint.x;
    IStop = 1;
    ftrialStop = bestPoint.f;
    acceptedSourceStop = bestPoint.source;
else
    sStop = x;
    IStop = -1;
    ftrialStop = fx;
    acceptedSourceStop = 'none';
end

infoStop = struct();

infoStop.q = q;
infoStop.qBeforeSafeguard = qBeforeSafeguard;
infoStop.qZeroTol = qZeroTol;
infoStop.qIsZero = qIsZero;
infoStop.qIsDescent = qIsDescent;
infoStop.rhoForce = rhoForce;
infoStop.d = d;

infoStop.bestPoint = bestPoint;

infoStop.usedDirectionalSafeguard = usedDirectionalSafeguard;
infoStop.usedComponentFlipSafeguard = usedComponentFlipSafeguard;
infoStop.usedFullRestartSafeguard = usedFullRestartSafeguard;
infoStop.numFlippedComponents = nnz(flipIdx);
infoStop.flipIdx = flipIdx;

infoStop.tl = tl;
infoStop.tu = tu;
infoStop.t0Paper = t0;
infoStop.epsLiS = epsLiS;
infoStop.tLower = tLower;
infoStop.zeta = zeta;
infoStop.pLiS = pLiS;

infoStop.type = stopType;
infoStop.innerIter = iterLocal;
infoStop.t = tLocal;
infoStop.barT = barTLocal;
infoStop.ftrial = ftrialStop;
infoStop.f_t = fTlocal;
infoStop.f_bar = fBarLocal;
infoStop.acceptedStepSource = acceptedSourceStop;
infoStop.forcedByMaxLiS = false;

infoStop.f_t_minus_fx = NaN;
infoStop.f_bar_minus_fx = NaN;

if isfinite(fTlocal)
    infoStop.f_t_minus_fx = fTlocal - fx;
end

if isfinite(fBarLocal)
    infoStop.f_bar_minus_fx = fBarLocal - fx;
end

infoStop.armijo_t_rhs = NaN;
infoStop.armijo_bar_rhs = NaN;
infoStop.armijo_t_margin = NaN;
infoStop.armijo_bar_margin = NaN;

infoStop.actualDecBar = NaN;
infoStop.requiredDecBar = NaN;

infoStop.gTd = NaN;
infoStop.enrich_rhs = NaN;
infoStop.enrich_margin = NaN;

infoStop.armijo_t_pass = false;
infoStop.armijo_bar_pass = false;
infoStop.enrich_pass = false;

infoStop.tTestType = 'stopped';
infoStop.barTestType = 'stopped';

infoStop.localFunEvals = localFunEvals;
infoStop.localDGEvals = localDGEvals;
infoStop.localBarEvals = localBarEvals;
infoStop.localApproxTotalEvals = localFunEvals + n*localDGEvals;

infoStop.trace.i = trace_i;
infoStop.trace.t = trace_t;
infoStop.trace.barT = trace_barT;
infoStop.trace.tl = trace_tl;
infoStop.trace.tu = trace_tu;
infoStop.trace.tl_next = trace_tl_next;
infoStop.trace.tu_next = trace_tu_next;

infoStop.trace.f_t_minus_fx = trace_f_t_minus_fx;
infoStop.trace.f_bar_minus_fx = trace_f_bar_minus_fx;
infoStop.trace.armijo_t_rhs = trace_armijo_t_rhs;
infoStop.trace.armijo_bar_rhs = trace_armijo_bar_rhs;
infoStop.trace.armijo_t_margin = trace_armijo_t_margin;
infoStop.trace.armijo_bar_margin = trace_armijo_bar_margin;

infoStop.trace.actualDecBar = trace_actualDecBar;
infoStop.trace.requiredDecBar = trace_requiredDecBar;

infoStop.trace.gTd = trace_gTd;
infoStop.trace.enrich_rhs = trace_enrich_rhs;
infoStop.trace.enrich_margin = trace_enrich_margin;

infoStop.trace.armijo_t_pass = trace_armijo_t_pass;
infoStop.trace.armijo_bar_pass = trace_armijo_bar_pass;
infoStop.trace.enrich_pass = trace_enrich_pass;
infoStop.trace.bar_eval_reused = trace_bar_eval_reused;

infoStop.trace.t_test_type = trace_t_test_type;
infoStop.trace.bar_test_type = trace_bar_test_type;

end

end

% ============================================================
function [dg, f2, best,ST] = DG_DiscreteGradient(fun,x,sl,g,f1,tune,ST)
% DG_DiscreteGradient
% Bagirov-style discrete-gradient approximation with nonfinite handling.
%
% Nonfinite policy:
%   +Inf  -> +largeMag
%   -Inf  -> -largeMag
%   NaN   -> median(finite(dg)) * 1e-3 (fallback 0 if no finite entries)

x = x(:);
g = g(:);
m = length(x);

% Default outputs.
% These must be assigned before any possible early return.
dg = zeros(m,1);
f2 = f1;

% Best sampled point tracker (opportunistic)
best = struct();
best.x = x;
best.f = f1;
best.improved = false;
best.source = 'base';
best.index = NaN;

% Normalize direction
ng = norm(g);
if ng == 0 || ~isfinite(ng)
    g = randUnit(m);
else
    g = g / ng;
end

pwt = tune.pwt;

[~, imax] = max(abs(g));
if abs(g(imax)) <= eps
    error('DG_DiscreteGradient: direction has no usable nonzero component.');
end

% First directional point
x1    = x + sl*g;
f2    = fun(x1);
ST.nf = ST.nf+1;

if isfinite(f2) && f2 < best.f
    best.x = x1;
    best.f = f2;
    best.improved = true;
    best.source = 'directional';
    best.index = imax;
end

% check stopping test
sec       = (cputime-ST.initTime);
ST.done  = (sec>ST.secmax)|(ST.nf>=ST.nfmax);
qf        = (f2-ST.ftarget)/(ST.finit-ST.ftarget);
ST.qf     = qf;
ST.done   = (ST.done|ST.qf<=ST.accf);
ST.sec    = sec;
if ST.done
    return
end

dsum = 0;
r4 = f2;

% Coordinate-perturbation sequence
for k = 1:m
    if k ~= imax
        r3    = r4;
        x1(k) = x1(k) + pwt;
        r4    = fun(x1);
        ST.nf = ST.nf+1;
        
        if isfinite(r4) && r4 < best.f
            best.x = x1;
            best.f = r4;
            best.improved = true;
            best.source = sprintf('coordinate_%d', k);
            best.index = k;
        end
        
        dg(k) = (r4 - r3) / pwt;
        dsum = dsum + dg(k)*g(k);
        
        % check stopping test
        sec       = (cputime-ST.initTime);
        ST.done  = (sec>ST.secmax)|(ST.nf>=ST.nfmax);
        qf        = (r4-ST.ftarget)/(ST.finit-ST.ftarget);
        ST.qf     = qf;
        ST.done   = (ST.done|ST.qf<=ST.accf);
        ST.sec    = sec;
        if ST.done
            % ---------------- Nonfinite handling before early return ----------------
            largeMag = 10;
            
            infPos = isinf(dg) & dg > 0;
            infNeg = isinf(dg) & dg < 0;
            dg(infPos) = +largeMag;
            dg(infNeg) = -largeMag;
            
            nanIdx = isnan(dg);
            if any(nanIdx)
                finiteIdx = isfinite(dg);
                if any(finiteIdx)
                    medVal = median(dg(finiteIdx));
                else
                    medVal = 0;
                end
                dg(nanIdx) = medVal * 1e-3;
            end
            
            return
        end
    end
    
    
end

% imax component from secant relation
den_imax = sl*g(imax);
if den_imax == 0 || ~isfinite(den_imax)
    dg(imax) = (f2 - f1)/max(den_imax, eps) - dsum/max(g(imax), eps);
else
    dg(imax) = (f2 - f1)/den_imax - dsum/g(imax);
end

% ---------------- Nonfinite handling ----------------
largeMag = 10;

% Clamp Inf to +/- largeMag
infPos = isinf(dg) & dg > 0;
infNeg = isinf(dg) & dg < 0;
dg(infPos) = +largeMag;
dg(infNeg) = -largeMag;

% Replace NaNs with median(finite(dg))*1e-3 (fallback 0 if none finite)
nanIdx = isnan(dg);
if any(nanIdx)
    finiteIdx = isfinite(dg);
    if any(finiteIdx)
        medVal = median(dg(finiteIdx));
    else
        medVal = 0;
    end
    dg(nanIdx) = medVal * 1e-3;
end

end

% ============================================================
function d = randUnit(n)
d = randn(n,1);
d = d / max(norm(d), eps);
end

% ============================================================
function [dnew, info] = MatCSG_direction(v, mem, tune)
% MatCSG_direction
%
% tune.dir = 1:
%   d = -v
%
%   tune.dir = 2 : matrix CG / MatCSG direction with beta and theta
%
% Matrix-CG safeguards:
%
%   1) Matrix-CG stability is checked on raw beta/theta before decay.
%
%   2) If the raw matrix correction is unstable, fallback to DL is used,
%      and the same coefficient decay is then applied to the fallback.
%
%   3) Dominance adjustment for dir = 2:
%      beta and theta are clipped so that beta*dprev and theta*dtilde
%      cannot dominate the robust component -v.

v = v(:);

persistent matcsgCoeffCounter

if isempty(matcsgCoeffCounter), matcsgCoeffCounter = 0; end

if isfield(tune,'resetMatCSGCounter') && tune.resetMatCSGCounter
    matcsgCoeffCounter = 0;
end


info = struct();

info.mode = tune.dir;
info.directionType = 'undefined';

info.betaCG = 0;
info.thetaCG = 0;
info.betaRawCG = 0;
info.thetaRawCG = 0;

info.Delta = NaN;
info.denDL = NaN;

info.coeffCounter = matcsgCoeffCounter;
info.coeffScale = 1;

info.usedFallback = false;
info.usedRestart = false;
info.usedAngleCorrection = false;
info.usedMatrixStabilityFallback = false;
info.matrixStabilityReason = 'none';

info.betaWasClipped = false;
info.thetaWasClipped = false;
info.betaWasNonfinite = false;
info.thetaWasNonfinite = false;
info.nonfiniteReason = 'none';

info.matrixSystemWasRepaired = false;
info.nonfiniteACount = 0;
info.nonfiniteBCount = 0;

info.matEntryMin = NaN;
info.matEntryMax = NaN;
info.matEigMin = NaN;
info.matEigMax = NaN;
info.matCondAbs = NaN;

info.rawMatEntryMin = NaN;
info.rawMatEntryMax = NaN;
info.rawMatEigMin = NaN;
info.rawMatEigMax = NaN;
info.rawMatCondAbs = NaN;

% ------------------------------------------------------------
% dir = 1: steepest descent only
% ------------------------------------------------------------
if tune.dir == 1
    
    dnew = -v;
    dnew = dnew / max(norm(dnew), eps);
    
    info.directionType = 'steepest';
    
    [info.matEntryMin, info.matEntryMax, ...
        info.matEigMin, info.matEigMax, info.matCondAbs] = ...
        betaMatrixStats(0, 0, length(v));
    
    return
    
end

% ------------------------------------------------------------
% Optional safeguard: restart after null step
% ------------------------------------------------------------
if tune.restartAfterNull
    
    if isfield(mem,'lastStepWasNull') && mem.lastStepWasNull
        
        dnew = -v;
        dnew = dnew / max(norm(dnew), eps);
        
        info.directionType = 'restart_after_null_step';
        info.usedRestart = true;
        
        [info.matEntryMin, info.matEntryMax, ...
            info.matEigMin, info.matEigMax, info.matCondAbs] = ...
            betaMatrixStats(0, 0, length(v));
        
        return
        
    end
    
end

% ------------------------------------------------------------
% Missing-memory restart
% ------------------------------------------------------------
if isempty(mem.dprev) || isempty(mem.yprev) || isempty(mem.sprev)
    
    dnew = -v;
    dnew = dnew / max(norm(dnew), eps);
    
    info.directionType = 'restart_missing_memory';
    info.usedRestart = true;
    
    [info.matEntryMin, info.matEntryMax, ...
        info.matEigMin, info.matEigMax, info.matCondAbs] = ...
        betaMatrixStats(0, 0, length(v));
    
    return
    
end

dprev = mem.dprev(:);
yprev = mem.yprev(:);
sprev = mem.sprev(:);

n = length(dprev);

% ------------------------------------------------------------
% Dai--Liao scalar beta
% ------------------------------------------------------------
denDL = dprev' * yprev;
info.denDL = denDL;

if ~isfinite(denDL) || abs(denDL) < eps
    
    dnew = -v;
    dnew = dnew / max(norm(dnew), eps);
    
    info.directionType = 'restart_bad_denDL';
    info.usedRestart = true;
    info.betaWasNonfinite = ~isfinite(denDL);
    info.nonfiniteReason = 'denDL_nan_or_inf';
    
    [info.matEntryMin, info.matEntryMax, ...
        info.matEigMin, info.matEigMax, info.matCondAbs] = ...
        betaMatrixStats(0, 0, n);
    
    return
    
end

zprev = yprev - tune.tDL * sprev;
betaDL = (v' * zprev) / denDL;

% Mandatory NaN/Inf beta safeguard.
if ~isfinite(betaDL)
    
    dnew = -v;
    dnew = dnew / max(norm(dnew), eps);
    
    info.directionType = 'restart_nonfinite_betaDL';
    info.usedRestart = true;
    info.betaWasNonfinite = true;
    info.nonfiniteReason = 'betaDL_nan_or_inf';
    
    info.betaCG = 0;
    info.thetaCG = 0;
    
    [info.matEntryMin, info.matEntryMax, ...
        info.matEigMin, info.matEigMax, info.matCondAbs] = ...
        betaMatrixStats(0, 0, n);
    
    return
    
end

% Optional practical safeguard for nonsmooth/nonconvex cases.
if tune.clipNegativeBeta && betaDL < 0
    betaDL = 0;
    info.betaWasClipped = true;
end

% Enforce signed lower bound for betaDL only if nonzero.
if betaDL > 0
    betaDL = max(betaDL, tune.betaMin);
elseif betaDL < 0
    betaDL = min(betaDL, -tune.betaMin);
end

info.betaCG = betaDL;
info.thetaCG = 0;
info.betaRawCG = betaDL;
info.thetaRawCG = 0;

% ============================================================
% dir = 2: matrix beta/theta with symmetric tridiagonal action
% ============================================================

dhat   = [0; dprev(1:n-1)];
dcheck = [dprev(2:n); 0];
dtilde = dhat + dcheck;

% ------------------------------------------------------------
% If only one displacement is available, use scalar DL fallback.
% ------------------------------------------------------------
if isempty(mem.yold) || isempty(mem.sold)
    
    beta = betaDL;
    theta = 0;
    
    info.usedFallback = true;
    info.directionType = 'matrix_CG_fallback_to_DL';
    
else
    
    yold = mem.yold(:);
    sold = mem.sold(:);
    zold = yold - tune.tDL * sold;
    
    A = [dprev' * yprev,  dtilde' * yprev;
         dprev' * yold,   dtilde' * yold];
    
    b = [v' * zprev;
         v' * zold];
    
    % --------------------------------------------------------
    % Correct nonfinite matrix-system safeguard.
    % Do not repair A or b using medians.
    % --------------------------------------------------------
    if any(~isfinite(A(:))) || any(~isfinite(b(:)))
        
        beta = betaDL;
        theta = 0;
        
        info.usedFallback = true;
        info.directionType = ...
            'matrix_CG_fallback_nonfinite_system';
        info.betaWasNonfinite = true;
        info.thetaWasNonfinite = true;
        info.nonfiniteReason = ...
            'matrix_system_nan_or_inf_fallback_to_DL';
        info.nonfiniteACount = nnz(~isfinite(A));
        info.nonfiniteBCount = nnz(~isfinite(b));
        
    else
        
        Delta = det(A);
        info.Delta = Delta;
        
        if ~isfinite(Delta)
            
            beta = betaDL;
            theta = 0;
            
            info.usedFallback = true;
            info.directionType = ...
                'matrix_CG_fallback_nonfinite_Delta';
            info.betaWasNonfinite = true;
            info.thetaWasNonfinite = true;
            info.nonfiniteReason = ...
                'matrix_delta_nan_or_inf_fallback_to_DL';
            
        elseif abs(Delta) >= tune.epsDelta
            
            pars = A\b;
            
            if any(~isfinite(pars))
                
                beta = betaDL;
                theta = 0;
                
                info.usedFallback = true;
                info.directionType = ...
                    'matrix_CG_fallback_nonfinite_parameters';
                info.betaWasNonfinite = true;
                info.thetaWasNonfinite = true;
                info.nonfiniteReason = ...
                        'matrix_beta_theta_nan_or_inf_fallback_to_DL';
                
            else
                
                beta = pars(1);
                theta = pars(2);
                
                info.directionType = 'matrix_CG';
                
            end
            
        else
            
            beta = betaDL;
            theta = 0;
            
            info.usedFallback = true;
            info.directionType = 'matrix_CG_fallback_small_Delta';
            
        end
        
    end
    
end

% ------------------------------------------------------------
% Final raw beta/theta finite check.
% ------------------------------------------------------------
if ~isfinite(beta) || ~isfinite(theta)
    
    dnew = -v;
    dnew = dnew / max(norm(dnew), eps);
    
    info.directionType = 'restart_nonfinite_beta_or_theta';
    info.usedRestart = true;
    
    info.betaWasNonfinite = ~isfinite(beta);
    info.thetaWasNonfinite = ~isfinite(theta);
    info.nonfiniteReason = 'beta_or_theta_nan_or_inf';
    
    info.betaCG = 0;
    info.thetaCG = 0;
    info.betaRawCG = 0;
    info.thetaRawCG = 0;
    
    [info.matEntryMin, info.matEntryMax, ...
        info.matEigMin, info.matEigMax, info.matCondAbs] = ...
        betaMatrixStats(0, 0, n);
    
    return
    
end

% ------------------------------------------------------------
% Enforce signed lower bound for beta/theta only if nonzero.
% ------------------------------------------------------------
if beta > 0
    beta = max(beta, tune.betaMin);
elseif beta < 0
    beta = min(beta, -tune.betaMin);
end

if theta > 0
    theta = max(theta, tune.thetaMin);
elseif theta < 0
    theta = min(theta, -tune.thetaMin);
end

% ------------------------------------------------------------
% Matrix-CG stability safeguard on RAW beta/theta.
% ------------------------------------------------------------
betaRaw = beta;
thetaRaw = theta;

info.betaRawCG = betaRaw;
info.thetaRawCG = thetaRaw;

[info.rawMatEntryMin, info.rawMatEntryMax, ...
    info.rawMatEigMin, info.rawMatEigMax, info.rawMatCondAbs] = ...
    betaMatrixStats(betaRaw, thetaRaw, n);

unstableMatrix = false;
reason = 'none';

if info.rawMatEigMin < tune.matEigMinTol
    unstableMatrix = true;
    reason = 'negative_or_too_small_matrix_eigenvalue';
elseif info.rawMatCondAbs > tune.matCondMax
    unstableMatrix = true;
    reason = 'large_matrix_condition_number';
elseif abs(betaRaw) > tune.betaMax
    unstableMatrix = true;
    reason = 'large_beta';
elseif abs(thetaRaw) > tune.thetaMax
    unstableMatrix = true;
    reason = 'large_theta';
end

if unstableMatrix
    
    betaRaw = betaDL;
    thetaRaw = 0;
    
    info.usedFallback = true;
    info.usedMatrixStabilityFallback = true;
    info.matrixStabilityReason = reason;
    info.directionType = ['matrix_CG_stability_fallback_to_DL: ', reason];
    
    info.betaRawCG = betaRaw;
    info.thetaRawCG = thetaRaw;
    
end

% ------------------------------------------------------------
% Dominance adjustment for the matrix-CG direction.
%
% The direction is
%
%   dnew = -v + beta*dprev + theta*dtilde.
%
% We clip beta and theta so that neither beta*dprev nor theta*dtilde
% dominates the robust component -v.
% ------------------------------------------------------------
nv = norm(v);
ndprev = norm(dprev);
ndtilde = norm(dtilde);

betaDomMax  = tune.betaDomFactor  * nv / max(ndprev, eps);
thetaDomMax = tune.thetaDomFactor * nv / max(ndtilde, eps);

betaMaxEff  = min(tune.betaMax,  betaDomMax);
thetaMaxEff = min(tune.thetaMax, thetaDomMax);

if betaRaw > betaMaxEff
    betaRaw = betaMaxEff;
    info.betaWasClipped = true;
elseif betaRaw < -betaMaxEff
    betaRaw = -betaMaxEff;
    info.betaWasClipped = true;
end

if thetaRaw > thetaMaxEff
    thetaRaw = thetaMaxEff;
    info.thetaWasClipped = true;
elseif thetaRaw < -thetaMaxEff
    thetaRaw = -thetaMaxEff;
    info.thetaWasClipped = true;
end

info.betaRawCG = betaRaw;
info.thetaRawCG = thetaRaw;

% ------------------------------------------------------------
% Persistent coefficient decay for matrix CG.
% ------------------------------------------------------------
if tune.useCoeffDecay
    matcsgCoeffCounter = matcsgCoeffCounter + 1;
    coeffScale = ...
       tune.coeffDecay/ (1 + matcsgCoeffCounter)^tune.coeffDecayPower;
else
    coeffScale = 1;
end

beta  = coeffScale * betaRaw;
theta = coeffScale * thetaRaw;

info.coeffCounter = matcsgCoeffCounter;
info.coeffScale = coeffScale;

info.betaCG = beta;
info.thetaCG = theta;

[info.matEntryMin, info.matEntryMax, ...
    info.matEigMin, info.matEigMax, info.matCondAbs] = ...
    betaMatrixStats(beta, theta, n);

% ------------------------------------------------------------
% Final matrix direction.
% ------------------------------------------------------------
dnew = -v + beta*dprev + theta*dtilde;

% ------------------------------------------------------------
% Best-point diagonal scaling for dir = 2.
%
% Apply sc BEFORE enforceAngle.
% Otherwise sc can destroy the bounded-angle condition.
% ------------------------------------------------------------
if isfield(tune,'scDir2') && numel(tune.scDir2) == n
    
    sc = tune.scDir2(:);
    sc = abs(sc);
    
    badSc = ~isfinite(sc) | (sc == 0);
    sc(badSc) = 1;
    
    sc = min(tune.bestScaleMax, max(tune.bestScaleMin, sc));
    
    dnew = sc .* dnew;
    
end

% ------------------------------------------------------------
% Final direction safety check after scaling.
% ------------------------------------------------------------
if any(~isfinite(dnew))
    
    dnew = -v;
    dnew = dnew / max(norm(dnew), eps);
    
    info.directionType = 'restart_nonfinite_final_direction';
    info.usedRestart = true;
    info.betaWasNonfinite = true;
    info.thetaWasNonfinite = true;
    info.nonfiniteReason = 'final_direction_nan_or_inf';
    
    info.betaCG = 0;
    info.thetaCG = 0;
    
    [info.matEntryMin, info.matEntryMax, ...
        info.matEigMin, info.matEigMax, info.matCondAbs] = ...
        betaMatrixStats(0, 0, n);
    
    return
    
end

% ------------------------------------------------------------
% Final bounded-angle safeguard.
%
% This must be applied AFTER diagonal scaling.
% enforceAngle normalizes the final direction.
% ------------------------------------------------------------
[dnew, angleInfo] = enforceAngle(dnew, v, tune);

info.usedRestart = angleInfo.usedRestart;
info.usedAngleCorrection = angleInfo.usedAngleCorrection;

end

% ============================================================
function [d, info] = enforceAngle(d, v, tune)
% enforceAngle
% Theorem-based bounded-angle correction:
%
%   d <- d - varpi*v
%
% where varpi is chosen so that the corrected direction satisfies
%
%   v' d / (||v|| ||d||) = -varrho
%
% in the nondegenerate correction branch.
%
% Degenerate positive-collinearity and near-zero corrected direction
% cases are protected by restart to -v.

info = struct();
info.usedRestart = false;
info.usedAngleCorrection = false;

v = v(:);
d = d(:);

nv = norm(v);
nd = norm(d);

if nv == 0
    d = zeros(size(v));
    return
end

if nd == 0 || ~isfinite(nd)
    d = -v / nv;
    info.usedRestart = true;
    return
end

varrho = tune.varrho;
varrho = min(max(varrho, 1e-12), 1 - 1e-12);

omega1 = nv^2;
omega2 = nd^2;
omega  = v' * d;

if ~isfinite(omega1) || ~isfinite(omega2) || ~isfinite(omega)
    d = -v / nv;
    info.usedRestart = true;
    return
end

% Already satisfies bounded-angle condition.
if omega <= -varrho * sqrt(omega1 * omega2)
    d = d / nd;
    return
end

info.usedAngleCorrection = true;

c = omega / sqrt(omega1 * omega2);
c = min(max(c, -1), 1);

% Degenerate positive-collinearity case:
% d is almost a positive multiple of v.
% The theoretical correction may produce nearly zero.
if abs(1 - c) <= tune.angleDegTol
    d = -v / nv;
    info.usedRestart = true;
    return
end

w = omega1 * omega2 * (1 - c^2) / (1 - varrho^2);
w = max(w, 0);

if w <= eps || ~isfinite(w)
    d = -v / nv;
    info.usedRestart = true;
    return
end

varpi = (omega + varrho * sqrt(w)) / omega1;

if ~isfinite(varpi)
    d = -v / nv;
    info.usedRestart = true;
    return
end

d = d - varpi*v;

if norm(d) <= tune.corrNormTol || v' * d >= 0 || any(~isfinite(d))
    d = -v / nv;
    info.usedRestart = true;
    return
end

d = d / norm(d);

end

% ============================================================
function v = minNormConv(V, tune)
% minNormConv (safe stabilized version)
%
% Solves:
%
%     min ||V p||  subject to  p >= 0,  sum(p)=1.
%
% Cached objective value at the current point.

[n, m] = size(V);

if m == 0, v = zeros(n,1); return; end

if m == 1, v = V(:,1); return; end

% ------------------------------------------------------------
% Remove only invalid columns.
% Do NOT remove rank-dependent columns.
% ------------------------------------------------------------
finiteCols = all(isfinite(V),1);

if ~all(finiteCols)
    
    V = V(:,finiteCols);
    m = size(V,2);
    
    if m == 0, v = zeros(n,1); return; end
    
    if m == 1, v = V(:,1); return; end
    
end

% ------------------------------------------------------------
% Gram matrix and adaptive regularization.
%
% Original objective:
%
%     min ||Vp||^2 = p'*(V'V)*p.
%
% quadprog solves:
%
%     min 0.5*p'*H*p + f'*p.
%
% Therefore H = 2*(V'V + reg*I).
% The regularization is only for numerical stabilization.
% Final v is still computed as V*p.
% ------------------------------------------------------------
G = V' * V;
G = 0.5*(G + G');

Gscale = max(1, norm(G,'fro'));
reg = tune.qpReg * Gscale;

H = 2*(G + reg*eye(m));
H = 0.5*(H + H');

f = zeros(m,1);
Aeq = ones(1,m);
beq = 1;
lb = zeros(m,1);
ub = [];

% ------------------------------------------------------------
% Try quadprog with adaptive regularization retry.
% ------------------------------------------------------------
p = [];

if exist('quadprog','file') == 2
    
    regTry = reg;
    
    for attempt = 1:2
        
        try
            
            Htry = 2*(G + regTry*eye(m));
            Htry = 0.5*(Htry + Htry');
            
            qopts = optimoptions( ...
                'quadprog', ...
                'Display','off', ...
                'Algorithm','interior-point-convex');
            
            p = quadprog(Htry, f, [], [], Aeq, beq, lb, ub, [], qopts);
            
            if ~isempty(p) && all(isfinite(p))
                break
            end
            
        catch
            p = [];
        end
        
        regTry = 10*regTry;
        
    end
    
end

% ------------------------------------------------------------
% Fallback: projected gradient on simplex.
% ------------------------------------------------------------
if isempty(p) || any(~isfinite(p))
    
    p = ones(m,1)/m;
    
    Hnorm = norm(H,2);
    eta = 1 / max(Hnorm, 1);
    
    oldObj = Inf;
    
    for it = 1:tune.simplexIter
        
        grad = H*p;
        
        pnew = projectSimplex(p - eta*grad);
        
        stepNorm = norm(pnew - p);
        obj = 0.5*pnew'*H*pnew;
        
        if stepNorm <= tune.simplexTol && ...
                abs(oldObj - obj) <= tune.simplexTol*max(1,abs(oldObj))
            
            p = pnew;
            break
            
        end
        
        p = pnew;
        oldObj = obj;
        
    end
    
end

% ------------------------------------------------------------
% Repair simplex feasibility.
% ------------------------------------------------------------
p(~isfinite(p)) = 0;
p = max(p,0);

s = sum(p);

if s <= 0
    p = ones(m,1)/m;
else
    p = p/s;
end

% ------------------------------------------------------------
% Optional tiny-weight cleanup.
%
% Keep tune.weightTol = 0 by default.
% Tiny weights can be essential for cancellation.
% ------------------------------------------------------------
if tune.weightTol > 0
    
    p(p < tune.weightTol) = 0;
    
    s = sum(p);
    
    if s <= 0
        p = ones(m,1)/m;
    else
        p = p/s;
    end
    
end

% ------------------------------------------------------------
% Return the original unregularized convex combination.
% ------------------------------------------------------------
v = V*p;

end

% ============================================================
function p = projectSimplex(y)
% projectSimplex
% Projection onto {p : p >= 0, sum(p)=1}.

y = y(:);
m = length(y);

u = sort(y,'descend');
cssv = cumsum(u);

rho = find(u > (cssv - 1)./(1:m)', 1, 'last');

if isempty(rho)
    p = ones(m,1)/m;
    return
end

theta = (cssv(rho) - 1)/rho;

p = max(y - theta, 0);

s = sum(p);

if s <= 0
    p = ones(m,1)/m;
else
    p = p/s;
end

end

% ============================================================
function printIterDynamic(k, fk, vstat, vhat, d, V, dinfo, lisInfo, ...
    stepKind, stepNorm, ftrial, prt)
% printIterDynamic
% Diagnostic-heavy printing.
% Does not change algorithmic behavior.

dn = norm(d);
nvstat = norm(vstat);
nvhat = norm(vhat);

if dn > 0 && nvstat > 0
    vd = vstat' * d;
    cosvd = vd / (nvstat * dn);
else
    vd = NaN;
    cosvd = NaN;
end

bundleSize = size(V,2);

t = getScalarField(lisInfo,'t',NaN);
barT = getScalarField(lisInfo,'barT',NaN);
innerIter = getScalarField(lisInfo,'innerIter',NaN);
q = getScalarField(lisInfo,'q',NaN);
qBefore = getScalarField(lisInfo,'qBeforeSafeguard',NaN);

f_t_minus_fx = getScalarField(lisInfo,'f_t_minus_fx',NaN);
f_bar_minus_fx = getScalarField(lisInfo,'f_bar_minus_fx',NaN);

armijo_t_margin = getScalarField(lisInfo,'armijo_t_margin',NaN);
armijo_bar_margin = getScalarField(lisInfo,'armijo_bar_margin',NaN);

actualDecBar = getScalarField(lisInfo,'actualDecBar',NaN);
requiredDecBar = getScalarField(lisInfo,'requiredDecBar',NaN);

gTd = getScalarField(lisInfo,'gTd',NaN);
enrich_rhs = getScalarField(lisInfo,'enrich_rhs',NaN);
enrich_margin = getScalarField(lisInfo,'enrich_margin',NaN);

tl = getScalarField(lisInfo,'tl',NaN);
tu = getScalarField(lisInfo,'tu',NaN);
t0Paper = getScalarField(lisInfo,'t0Paper',NaN);

usedDirectionalSafeguard = ...
              getScalarField(lisInfo,'usedDirectionalSafeguard',NaN);
usedComponentFlipSafeguard = ...
               getScalarField(lisInfo,'usedComponentFlipSafeguard',NaN);
usedFullRestartSafeguard =  ...
              getScalarField(lisInfo,'usedFullRestartSafeguard',NaN);
numFlippedComponents = ...
               getScalarField(lisInfo,'numFlippedComponents',NaN);

localFunEvals = getScalarField(lisInfo,'localFunEvals',NaN);
localDGEvals = getScalarField(lisInfo,'localDGEvals',NaN);
localBarEvals = getScalarField(lisInfo,'localBarEvals',NaN);
localApproxTotalEvals = ...
                 getScalarField(lisInfo,'localApproxTotalEvals',NaN);

forcedByMaxLiS = getScalarField(lisInfo,'forcedByMaxLiS',NaN);

if isfield(lisInfo,'type')
    lisType = lisInfo.type;
else
    lisType = 'unknown';
end

if isnan(ftrial)
    df = NaN;
else
    df = ftrial - fk;
end

% ------------------------------------------------------------
% Compact one-line print.
% Sign convention:
%   Armijo margin <= 0 means pass.
%   Enrichment margin >= 0 means pass.
% ------------------------------------------------------------
if prt >= 1
    fprintf(['k=%5d | f=%+.5e | df=%+.2e | ||vstat||=%8.2e |',...
        ' ||vhat||=%8.2e | ', ...
        '||d||=%8.2e | vstatTd=%+.2e | cos=%+.2e | ', ...
        'q=%+.2e | t=%8.2e | barT=%8.2e | LiS=%2d | B=%2d | %s'], ...
        k, fk, df, nvstat, nvhat, dn, vd, cosvd, ...
        q, t, barT, innerIter, bundleSize, stepKind);
    
    if usedComponentFlipSafeguard == 1
        fprintf(' | flip=%d', numFlippedComponents);
    end
    
    if usedFullRestartSafeguard == 1
        fprintf(' | fullRestart');
    end
    
    if dinfo.usedMatrixStabilityFallback
        fprintf(' | MFB: %s', dinfo.matrixStabilityReason);
    end
    
    if forcedByMaxLiS == 1
        fprintf(' | forcedMaxLiS');
    end
    
    fprintf('\n');
end

% ------------------------------------------------------------
% Detailed diagnostics.
% ------------------------------------------------------------
if prt >= 2
    disp(' ')
    disp('---------------- dynamic diagnostics ----------------')
    
    fprintf('iteration k                     = %d\n',k)
    fprintf('direction mode dir              = %d\n',dinfo.mode)
    fprintf('direction type                  = %s\n',dinfo.directionType)
    
    fprintf('step type                       = %s\n', stepKind)
    fprintf('line-search type                = %s\n', lisType)
    fprintf('line-search inner iter          = %d\n', innerIter)
    
    fprintf('f(x_k)                          = %+.16e\n', fk)
    fprintf('f(trial)-f(x_k)                 = %+.16e\n', df)
    fprintf('accepted step norm              = %.16e\n', stepNorm)
    
    fprintf('||vstat_k||                     = %.16e\n', nvstat)
    fprintf('||vhat_k||                      = %.16e\n', nvhat)
    fprintf('||d_k||                         = %.16e\n', dn)
    fprintf('vstat_k^T d_k                   = %+.16e\n', vd)
    fprintf('cos(vstat_k,d_k)                = %.16e\n', cosvd)
    
    disp(' ')
    disp('---- line-search scalar state ----')
    fprintf('q before safeguard              = %+.16e\n', qBefore)
    fprintf('q after safeguard               = %+.16e\n', q)
    fprintf('t_i                              = %.16e\n', t)
    fprintf('bar_t_i                          = %.16e\n', barT)
    fprintf('paper t0                         = %.16e\n', t0Paper)
    fprintf('bracket tl                       = %.16e\n', tl)
    fprintf('bracket tu                       = %.16e\n', tu)
    
    disp(' ')
    disp('---- Armijo diagnostics ----')
    fprintf('f(x+t_i d)-f(x)                 = %+.16e\n', f_t_minus_fx)
    fprintf('Armijo rhs at t_i               = %+.16e\n', ...
            getScalarField(lisInfo,'armijo_t_rhs',NaN))
    fprintf('Armijo margin at t_i            = %+.16e   (<=0 pass)\n',...
           armijo_t_margin)
    
    fprintf('f(x+bar_t_i d)-f(x)             = %+.16e\n',f_bar_minus_fx)
    fprintf('Armijo rhs at bar_t_i           = %+.16e\n',...
           getScalarField(lisInfo,'armijo_bar_rhs',NaN))
    fprintf('Armijo margin at bar_t_i        = %+.16e   (<=0 pass)\n',...
          armijo_bar_margin)
    fprintf('actual decrease at bar_t_i      = %+.16e\n', actualDecBar)
    fprintf('required decrease at bar_t_i    = %+.16e\n', requiredDecBar)
    
    disp(' ')
    disp('---- enrichment diagnostics ----')
    fprintf('g_i^T d                         = %+.16e\n', gTd)
    fprintf('mu2*q                           = %+.16e\n', enrich_rhs)
    fprintf('enrichment margin               = %+.16e   (>=0 pass)\n',...
            enrich_margin)
    
    disp(' ')
    disp('---- safeguard diagnostics ----')
    fprintf('used directional safeguard       = %d\n', ...
        usedDirectionalSafeguard)
    fprintf('used component flip safeguard    = %d\n', ...
        usedComponentFlipSafeguard)
    fprintf('used full restart safeguard      = %d\n', ...
        usedFullRestartSafeguard)
    fprintf('number flipped components        = %d\n', ...
        numFlippedComponents)
    
    disp(' ')
    disp('---- CG / matrix diagnostics ----')
    fprintf('betaCG                          = %+.16e\n', ...
        dinfo.betaCG)
    fprintf('thetaCG                         = %+.16e\n', ...
        dinfo.thetaCG)
    fprintf('Delta                           = %+.16e\n', ...
        dinfo.Delta)
    fprintf('denDL                           = %+.16e\n', ...
        dinfo.denDL)
    
    fprintf('matrix beta entry min           = %+.16e\n', ...
        dinfo.matEntryMin)
    fprintf('matrix beta entry max           = %+.16e\n', ...
        dinfo.matEntryMax)
    fprintf('matrix beta eig min             = %+.16e\n', ...
        dinfo.matEigMin)
    fprintf('matrix beta eig max             = %+.16e\n', ...
        dinfo.matEigMax)
    fprintf('matrix beta abs cond            = %+.16e\n', ...
        dinfo.matCondAbs)
    
    fprintf('used fallback                   = %d\n', ...
        dinfo.usedFallback)
    fprintf('used restart                    = %d\n', ...
        dinfo.usedRestart)
    fprintf('used angle correction           = %d\n', ...
        dinfo.usedAngleCorrection)
    fprintf('matrix stability fallback       = %d\n', ...
        dinfo.usedMatrixStabilityFallback)
    fprintf('matrix stability reason         = %s\n', ...
        dinfo.matrixStabilityReason)
    fprintf('beta was clipped                = %d\n', ...
        dinfo.betaWasClipped)
    fprintf('beta was nonfinite              = %d\n', ...
        dinfo.betaWasNonfinite)
    fprintf('theta was nonfinite             = %d\n', ...
        dinfo.thetaWasNonfinite)
    fprintf('nonfinite reason                = %s\n', ...
        dinfo.nonfiniteReason)
    
    disp(' ')
    disp('---- local evaluation-cost diagnostics ----')
    fprintf('local direct fun evals in LiS    = %d\n', ...
        localFunEvals)
    fprintf('local barT fun evals             = %d\n', ...
        localBarEvals)
    fprintf('local DG calls in LiS            = %d\n', ...
        localDGEvals)
    fprintf('approx total local evals         = %d\n', ...
        localApproxTotalEvals)
    
    disp('-----------------------------------------------------')
end

% ------------------------------------------------------------
% Full inner line-search trace.
% This is the most useful part for debugging stagnation.
% ------------------------------------------------------------
if prt >= 3 && isfield(lisInfo,'trace')
    printLiSTrace(lisInfo.trace);
end

end

% ============================================================
function printMemoryDiagnostics(k, x, V, v, d, mem, pending)
% printMemoryDiagnostics
% Extra dynamic diagnostics for debugging memory updates.

disp(' ')
disp('*************** memory/vector diagnostics ********************')
fprintf('iteration k              = %d\n', k)
fprintf('dimension n              = %d\n', length(x))
fprintf('bundle size              = %d\n', size(V,2))
fprintf('||x||                    = %.16e\n', norm(x))
fprintf('||v||                    = %.16e\n', norm(v))
fprintf('||d||                    = %.16e\n', norm(d))
fprintf('pending.active           = %d\n', pending.active)
fprintf('mem.lastStepWasNull      = %d\n', mem.lastStepWasNull)

if isempty(mem.dprev)
    fprintf('mem.dprev                = []\n')
else
    fprintf('||mem.dprev||            = %.16e\n', norm(mem.dprev))
end

if isempty(mem.sprev)
    fprintf('mem.sprev                = []\n')
else
    fprintf('||mem.sprev||            = %.16e\n', norm(mem.sprev))
end

if isempty(mem.yprev)
    fprintf('mem.yprev                = []\n')
else
    fprintf('||mem.yprev||            = %.16e\n', norm(mem.yprev))
end

if isempty(mem.sold)
    fprintf('mem.sold                 = []\n')
else
    fprintf('||mem.sold||             = %.16e\n', norm(mem.sold))
end

if isempty(mem.yold)
    fprintf('mem.yold                 = []\n')
else
    fprintf('||mem.yold||             = %.16e\n', norm(mem.yold))
end

disp('**************************************************************')

end

% ============================================================
function val = getScalarField(S, fieldName, defaultVal)
% getScalarField
% Safe scalar field extractor for diagnostics.

if isstruct(S) && isfield(S, fieldName)
    tmp = S.(fieldName);
    
    if isempty(tmp)
        val = defaultVal;
    elseif isnumeric(tmp) || islogical(tmp)
        val = tmp(1);
    else
        val = defaultVal;
    end
else
    val = defaultVal;
end

end

% ============================================================
function printLiSTrace(tr)
% printLiSTrace
% Prints the full inner DG-TPLiS trace.
%
% Sign convention:
%   armijo_*_margin <= 0 means Armijo passed.
%   enrich_margin   >= 0 means enrichment passed.

disp(' ')
disp('++++++++++++++++ DG-TPLiS inner trace ++++++++++++++++')
disp(['    i', ...
    '          t_i', ...
    '        barT_i', ...
    '            tl', ...
    '            tu', ...
    '      f_t-fx', ...
    '    f_bar-fx', ...
    '   ArmTmargin', ...
    ' ArmBarMargin', ...
    '          gTd', ...
    '      EnMargin', ...
    '  AT AB EN'])

m = length(tr.i);

for r = 1:m
    fprintf(['%5d  %12.4e  %12.4e  %12.4e  %12.4e  ', ...
        '%+12.4e  %+12.4e  %+12.4e  %+12.4e  ', ...
        '%+12.4e  %+12.4e   %1d  %1d  %1d\n'], ...
        tr.i(r), tr.t(r), tr.barT(r), tr.tl(r), tr.tu(r), ...
        tr.f_t_minus_fx(r), tr.f_bar_minus_fx(r), ...
        tr.armijo_t_margin(r), tr.armijo_bar_margin(r), ...
        tr.gTd(r), tr.enrich_margin(r), ...
        tr.armijo_t_pass(r), tr.armijo_bar_pass(r), tr.enrich_pass(r));
end

disp('AT = Armijo at t_i, AB = Armijo at barT_i, EN = enrichment')
disp('Armijo margin <= 0 passes; enrichment margin >= 0 passes.')
disp('+++++++++++++++++++++++++++++++++++++++++++++++++++++++')
end

% ============================================================
function [entryMin, entryMax, eigMin, eigMax, condAbs] = ...
    betaMatrixStats(beta, theta, n)
% betaMatrixStats
% Matrix beta is the symmetric tridiagonal matrix
%
%   B = beta*I + theta*(subdiag + superdiag).
%
% We report both entry range and eigenvalue range.

if n <= 0
    entryMin = NaN;
    entryMax = NaN;
    eigMin = NaN;
    eigMax = NaN;
    condAbs = NaN;
    return
end

entryVals = [0; beta; theta];
entryMin = min(entryVals);
entryMax = max(entryVals);

if n == 1
    eigsB = beta;
else
    j = (1:n)';
    eigsB = beta + 2*theta*cos(j*pi/(n+1));
end

eigMin = min(eigsB);
eigMax = max(eigsB);

absEig = abs(eigsB);
if min(absEig) <= eps
    condAbs = Inf;
else
    condAbs = max(absEig) / min(absEig);
end

end


% ============================================================
function [bestHistX, bestHistF] = ...
   updateBestPointArchive(bestHistX, bestHistF, xnew, fnew, maxMem)
% updateBestPointArchive
%
% Keeps at most maxMem best points, sorted by function value.

xnew = xnew(:);

if isempty(bestHistX)
    bestHistX = xnew;
    bestHistF = fnew;
    return
end

% Append new candidate.
bestHistX = [bestHistX, xnew];
bestHistF = [bestHistF(:); fnew];

% Remove nonfinite function values.
finiteIdx = isfinite(bestHistF);
bestHistX = bestHistX(:, finiteIdx);
bestHistF = bestHistF(finiteIdx);

% Sort by function value.
[bestHistF, ord] = sort(bestHistF, 'ascend');
bestHistX = bestHistX(:, ord);

% Remove near-duplicate points, keeping the best occurrence.
keep = true(length(bestHistF),1);

for i = 1:length(bestHistF)
    if ~keep(i)
        continue
    end
    
    for j = i+1:length(bestHistF)
        if norm(bestHistX(:,i) - bestHistX(:,j), inf) ...
                          <= 1e-14*max(1,norm(bestHistX(:,i),inf))
            keep(j) = false;
        end
    end
end

bestHistX = bestHistX(:, keep);
bestHistF = bestHistF(keep);

% Keep only maxMem best points.
m = min(maxMem, length(bestHistF));
bestHistX = bestHistX(:, 1:m);
bestHistF = bestHistF(1:m);

end

% ============================================================
function sc = computeBestPointScale(bestHistX, bestHistF, tune, n)
% computeBestPointScale
%
% Builds a componentwise scaling vector from the saved best points.
%
% Let x_best be the best saved point. The vector sc is computed from
% abs(x_j - x_best), where x_j are the other saved best points.
%
% If fewer than two best points are available, sc = ones(n,1).

sc = ones(n,1);

if isempty(bestHistX) || isempty(bestHistF)
    return
end

if size(bestHistX,2) < 2
    return
end

% Sort archive by function value.
[~, ord] = sort(bestHistF(:), 'ascend');
bestHistX = bestHistX(:, ord);

m = min(tune.bestScaleMemory, size(bestHistX,2));

if m < 2
    return
end

xbest = bestHistX(:,1);
Xother = bestHistX(:,2:m);

D = abs(Xother - xbest);

switch lower(tune.bestScaleMethod)
    case 'max'
        sc = max(D, 2);
    case  'min'
        sc = min(D, 2);
    case  'mean'
        sc = mean(D, 2);
    case 'median'
         sc = median(D, 2);
    otherwise
        error('bestScaleMethod must be max, min, mean, median')
       
end

sc = sc(:);

% Replace zero and nonfinite entries by one.
bad = ~isfinite(sc) | (sc == 0);
sc(bad) = 1;

% Apply uniform componentwise bounds to the scaling vector.
sc = min(tune.bestScaleMax, max(tune.bestScaleMin, sc));

end

% ============================================================
function tune = initTune(tune, n)
% initTune
% Initializes default tuning parameters for DG-Clarke / DG-MatCSG.
%
% Input:
%   tune : user-defined tuning structure
%   n    : problem dimension
%
% Output:
%   tune : tuning structure with missing fields completed

if nargin < 1 || isempty(tune), tune = struct(); end
if nargin < 2 || isempty(n), n = 1; end

if ~isstruct(tune)
    error('initTune: tune must be a structure.');
end

% Controls whether DG_MatCSG may stop by its internal stationarity test.
if ~isfield(tune,'disableInnerStationarity'),
    tune.disableInnerStationarity = true; 
end

% Initial stationarity tolerance for the outer Clarke loop.
if ~isfield(tune,'delta0'), tune.delta0 = 1e-1; end

% Initial line-search radius for the outer Clarke loop.
if ~isfield(tune,'epsLiS0'), tune.epsLiS0 = 1e-1; end

% Initial discrete-gradient displacement parameter.
if ~isfield(tune,'lambda0'), tune.lambda0 = 1e-2; end

% Exponential decay rate for outer parameters.
if ~isfield(tune,'outerDecay'), tune.outerDecay = 0.25; end

% Number of outer iterations before decay starts.
if ~isfield(tune,'warmOuter'), tune.warmOuter = 3; end

% Lower bound for the outer stationarity tolerance.
if ~isfield(tune,'deltaMin'), tune.deltaMin = 1e-20; end

% Lower bound for the line-search radius.
if ~isfield(tune,'epsLiSMin'), tune.epsLiSMin = 1e-8; end

% Lower bound for the discrete-gradient displacement parameter.
if ~isfield(tune,'lambdaMin'), tune.lambdaMin = 1e-6; end

% Preserve a user-supplied coordinate perturbation as the initial value.
if isfield(tune,'pwt') && ~isfield(tune,'pwt0'), tune.pwt0 = tune.pwt; end

% Initial coordinate perturbation used 
% in the discrete-gradient construction.
if ~isfield(tune,'pwt0'), tune.pwt0 = 1e-8; end

% Optional user-defined sequence for delta.
if ~isfield(tune,'deltaSeq'), tune.deltaSeq = []; end

% Optional user-defined sequence for epsLiS.
if ~isfield(tune,'epsSeq'), tune.epsSeq = []; end

% Optional user-defined sequence for lambda.
if ~isfield(tune,'lambdaSeq'), tune.lambdaSeq = []; end

% Controls whether all outer outputs are stored.
if ~isfield(tune,'storeOuter'), tune.storeOuter = false; end

% Direction mode: 1 = steepest discrete-gradient direction,
%                 2 = MatCSG direction.
if ~isfield(tune,'dir'), tune.dir = 1; end

% Validate the direction mode.
if ~ismember(tune.dir,[1,2]),
    error('initTune: tune.dir must be 1 or 2.');
end

% Internal stationarity tolerance used by DG_MatCSG.
if ~isfield(tune,'delta'), tune.delta = -inf; end

% Maximum number of DG_MatCSG iterations per call.
if ~isfield(tune,'maxIter'), tune.maxIter = 100; end

% Maximum number of inner DG-TPLiS iterations.
if ~isfield(tune,'maxLiS'), tune.maxLiS = 50; end

% Discrete-gradient displacement parameter.
if ~isfield(tune,'lambda'), tune.lambda = 0.01; end

% Coordinate perturbation used in the discrete-gradient construction.
if ~isfield(tune,'pwt'), tune.pwt = tune.lambda^2; end

% Armijo sufficient-decrease parameter.
if ~isfield(tune,'mu1'), tune.mu1 = 1e-4; end

% Enrichment parameter in DG-TPLiS.
if ~isfield(tune,'mu2'), tune.mu2 = 0.1; end

% Line-search radius used in DG-TPLiS.
if ~isfield(tune,'epsLiS'), tune.epsLiS = 0.5; end

% Lower accepted-step threshold in DG-TPLiS.
if ~isfield(tune,'tLowerLiS'), tune.tLowerLiS = 1e-2; end

% Interval-contraction parameter in DG-TPLiS.
if ~isfield(tune,'zetaLiS'), tune.zetaLiS = 0.25; end

% Exponent parameter used in the auxiliary serious-step sequence.
if ~isfield(tune,'pLiS'), tune.pLiS = 5; end

% Optional initial trial step for DG-TPLiS.
if ~isfield(tune,'tInitLiS'), tune.tInitLiS = []; end

% Dai--Liao parameter used in the scalar fallback.
if ~isfield(tune,'tDL'), tune.tDL = 1; end

% Determinant threshold for the two-parameter matrix system.
if ~isfield(tune,'epsDelta'), tune.epsDelta = 1e-12; end

% Bounded-angle parameter.
if ~isfield(tune,'varrho'), tune.varrho = 1e-12; end

% Degeneracy tolerance for angle correction.
if ~isfield(tune,'angleDegTol'), tune.angleDegTol = 1e-12; end

% Minimum acceptable norm after angle correction.
if ~isfield(tune,'corrNormTol'), tune.corrNormTol = 1e-14; end

% Controls whether negative beta values are clipped.
if ~isfield(tune,'clipNegativeBeta'), tune.clipNegativeBeta = false; end

% Controls whether MatCSG restarts after null steps.
if ~isfield(tune,'restartAfterNull'), tune.restartAfterNull = false; end

% Minimum admissible eigenvalue for the matrix-CG parameter matrix.
if ~isfield(tune,'matEigMinTol'), tune.matEigMinTol = -1e-12; end

% Maximum admissible absolute condition number for 
% the matrix-CG parameter matrix.
if ~isfield(tune,'matCondMax'), tune.matCondMax = 1e3; end

% Upper bound for the beta coefficient.
if ~isfield(tune,'betaMax'), tune.betaMax = 5; end

% Upper bound for the theta coefficient.
if ~isfield(tune,'thetaMax'), tune.thetaMax = 5; end

% Regularization parameter for the convex-hull minimal-norm QP.
if ~isfield(tune,'qpReg'), tune.qpReg = 1e-12; end

% Maximum number of projected-gradient iterations
% for the simplex fallback.
if ~isfield(tune,'simplexIter'), tune.simplexIter = 500; end

% Maximum number of stored discrete gradients in the bundle.
if ~isfield(tune,'maxBundle'), tune.maxBundle = max(n+3,25); end

% Stopping tolerance for the simplex projection fallback.
if ~isfield(tune,'simplexTol'), tune.simplexTol = 1e-12; end

% Threshold for optional removal of tiny simplex weights.
if ~isfield(tune,'weightTol'), tune.weightTol = 0; end

% Lower bound for nonzero beta values.
if ~isfield(tune,'betaMin'), tune.betaMin = 0; end

% Lower bound for nonzero theta values.
if ~isfield(tune,'thetaMin'), tune.thetaMin = 0; end

% Controls whether beta and theta are damped across iterations.
if ~isfield(tune,'useCoeffDecay'), tune.useCoeffDecay = true; end

% Decay exponent for beta and theta.
if ~isfield(tune,'coeffDecayPower'), tune.coeffDecayPower = 0.85; end

% Multiplicative coefficient-decay factor.
if ~isfield(tune,'coeffDecay'), tune.coeffDecay = 1; end

% Controls whether best-point diagonal scaling is used for MatCSG.
if ~isfield(tune,'useBestScaleDir2'), tune.useBestScaleDir2 = true; end

% Number of best points stored for diagonal scaling.
if ~isfield(tune,'bestScaleMemory'), tune.bestScaleMemory = 3; end

% Componentwise statistic used to compute the scaling vector.
if ~isfield(tune,'bestScaleMethod'),tune.bestScaleMethod='median';end

% Lower bound for the diagonal scaling vector.
if ~isfield(tune,'bestScaleMin'), tune.bestScaleMin = 1; end

% Upper bound for the diagonal scaling vector.
if ~isfield(tune,'bestScaleMax'), tune.bestScaleMax = 1e8; end

% Dominance factor for the beta memory component.
if ~isfield(tune,'betaDomFactor'), tune.betaDomFactor = 0.5; end

% Dominance factor for the theta memory component.
if ~isfield(tune,'thetaDomFactor'), tune.thetaDomFactor = 0.25; end

end

% ============================================================
function ST=initST(ST)


% ST           % structure with stop and print criteria
%              %   (indefinite run if no stopping criterion is given)
%  .secmax     %   stop if sec>=secmax (default: inf)
%  .nfmax      %   stop if nf>=nfmax   (default: inf)
%  .acc        %   stop if qf<=acc     (default: 1e-4)
%  .ftarget    %   function value accepted as optimal (default: 0)
%  .prt        %   Printing levels:
%                    -1 : silent
%                     0 : start/finish only
%                     1 : compact dynamic iteration line
%                     2 : detailed dynamic diagnostics
%                     3 : detailed diagnostics + memory/
%                         vector norms + LiS trace


if ~isfield(ST,'secmax'), ST.secmax=inf; end
if ~isfield(ST,'nfmax'), ST.nfmax=inf; end
if ~isfield(ST,'ftarget'), ST.ftarget=0; end
if ~isfield(ST,'accf'), ST.accf=-inf; end
if ~isfield(ST,'prt'), ST.prt=-1; end

if ~isfield(ST,'nf'), ST.nf = 0; end
if ~isfield(ST,'sec'), ST.sec = 0; end
if ~isfield(ST,'qf'), ST.qf = inf; end
if ~isfield(ST,'done'), ST.done = false; end
if ~isfield(ST,'initTime'), ST.initTime = cputime; end
end

% ============================================================
function [x,outAll,tune,ST]=checkInput(fun,x,tune,ST)

outAll = {};

if nargin < 1
    message = 'fun and x must be as input!';
    disp(message)
    x = [];
    outAll.flag = 'input_error';
    outAll.message = message;
end

if nargin < 2
    message = 'x must be as input!';
    disp(message)
    x = [];
    outAll.flag = 'input_error';
    outAll.message = message;
end
if nargin < 3 || isempty(tune)
    tune = [];
end

if nargin < 4 || isempty(ST)
    ST = [];
end

if ~isstruct(tune)
    message = 'tune should be a structure';
    disp(message)
    x = [];
    outAll.flag = 'input_error';
    outAll.message = message;
    return
end


% check function handle
if isempty(fun)
    message = 'DG-MatCSG needs the function handle fun to be defined';
    disp(message)
    x = [];
    outAll.flag = 'input_error';
    outAll.message = message;
    return
elseif ~isa(fun,'function_handle')
    message = 'fun should be a function handle';
    disp(message)
    x = [];
    outAll.flag = 'input_error';
    outAll.message = message;
    return
end

% starting point
if isempty(x)
    message = 'starting point must be defined';
    disp(message)
    x = [];
    outAll.flag = 'input_error';
    outAll.message = message;
    return
elseif ~isa(x,'numeric') || ~isvector(x)
    message = 'x should be a numeric vector';
    disp(message)
    x = [];
    outAll.flag = 'input_error';
    outAll.message = message;
    return
end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
