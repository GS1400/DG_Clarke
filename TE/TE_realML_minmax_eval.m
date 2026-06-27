function f = TE_realML_minmax_eval(x, data)
% TE_realML_minmax_eval
%
% Objective evaluator for real-data ML finite-max problems.
%
% Each objective has the form
%
%     f(x) = max_i ell_i(x).
%
% Inputs:
%   x    : decision vector
%   data : structure with fields A, b, y, loss, param
%
% Output:
%   f    : scalar objective value
%
% No subgradient is returned.

x = x(:);

A = data.A;
b = data.b;
y = data.y;

lossName = data.loss;
param = data.param;

Ax = A*x;
r  = Ax - b;

switch lossName

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 'max_abs_residual'

        phi = abs(r);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 'max_squared_residual'

        phi = 0.5*r.^2;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 'max_huber_residual'

        delta = param.delta;
        ar = abs(r);

        phi = zeros(size(r));

        idx1 = ar <= delta;
        idx2 = ar >  delta;

        phi(idx1) = 0.5*r(idx1).^2;
        phi(idx2) = delta*(ar(idx2) - 0.5*delta);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 'max_quantile'

        tau = param.tau;
        phi = max(tau*r, (tau-1)*r);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 'max_epsilon_insensitive'

        epsilon = param.epsilon;
        phi = max(0, abs(r) - epsilon);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 'max_hinge'

        margin = 1 - y.*Ax;
        phi = max(0, margin);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 'max_squared_hinge'

        margin = 1 - y.*Ax;
        phi = max(0, margin).^2;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 'max_logistic'

        z = y.*Ax;

        % Stable log(1+exp(-z)).
        phi = max(0, -z) + log(1 + exp(-abs(z)));

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 'max_relu_residual'

        phi = max(0, r);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 'max_fractional_residual'

        epsilon = param.epsilon;
        phi = (r.^2 + epsilon).^(1/4);

    otherwise

        error('Unknown loss name: %s', lossName);

end

f = max(phi);

if ~isfinite(f)
    error('Objective returned nonfinite value.');
end

end