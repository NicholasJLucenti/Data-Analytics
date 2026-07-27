function [rational, library_names] = run_implicit_sindy(X, dXdt, lambda, poly_order)
%RUN_IMPLICIT_SINDY Discover rational governing equations (e.g.
%Michaelis-Menten kinetics Vmax*x/(Km+x)) via an implicit/SINDy-PI style
%formulation, reusing the same stlsq_solve.m core as every other variant.
%
%   [rational, library_names] = RUN_IMPLICIT_SINDY(X, dXdt, lambda, poly_order)
%
%   THE TRICK: standard/weak SINDy assume dx/dt = Theta(x)*xi -- linear in
%   the unknowns, explicit in dx/dt. Rational dynamics like
%   dx/dt = a*x/(b+x) are NOT linear in (a,b). But multiplying through by
%   the denominator makes it linear again:
%       dx/dt * (b + x) = a*x
%       dx/dt*b + dx/dt*x - a*x = 0
%   Generalizing: build an extended library containing both the ordinary
%   polynomial terms phi_i(x) AND every product phi_i(x)*dx/dt. Every
%   candidate implicit equation is then linear in its coefficients. To
%   avoid the trivial zero solution, fix the coefficient on the "1*dxdt"
%   term to 1 and move it to the RHS as the regression target:
%       sum_i c_i*phi_i(x)  +  sum_{j>=2} d_j*phi_j(x)*dxdt  =  -dxdt
%   Solving this with STLSQ (exactly like standard SINDy, just a
%   different design matrix) and rearranging gives:
%       dxdt = [ -sum_i c_i*phi_i(x) ] / [ 1 + sum_{j>=2} d_j*phi_j(x) ]
%   i.e. a genuine rational function falls out directly -- the numerator
%   is the fitted c's, the denominator is 1 plus the fitted d's. This is
%   why rational terms don't need to be pre-specified: any (Km+x)-style
%   denominator up to poly_order shows up automatically if the data
%   supports it.
%
%   Inputs:
%     X          - state matrix (N x D)
%     dXdt       - derivative matrix (N x D), same preprocessing
%                  requirements as run_standard_sindy.m (finite-difference
%                  derivatives on smoothed data). A weak-form version of
%                  this trick is possible but not implemented here -- see
%                  note at bottom.
%     lambda     - sparsity threshold, applied to both numerator and
%                  denominator coefficients (default 0.1)
%     poly_order - polynomial library order underlying both numerator and
%                  denominator candidate terms (default 2)
%
%   Output:
%     rational(d) - struct array, one entry per state dimension d:
%       .numerator_Xi     - (M x 1) coefficients on library_names, the
%                            fitted numerator polynomial
%       .denominator_Xi   - (M x 1) coefficients on library_names, the
%                            fitted denominator polynomial. Entry 1 (the
%                            constant term) is always exactly 1 by
%                            construction, NOT a fitted/thresholded value.
%       .min_abs_denominator - smallest |denominator(X)| observed across
%                            the data. If this is close to zero, the
%                            fitted model has a near-singularity within
%                            the observed range and forward simulation
%                            will be unreliable -- check this before
%                            trusting the model.
%     library_names - cell array of term names shared by every
%                     numerator_Xi/denominator_Xi (same polynomial basis
%                     used for both, from build_polynomial_library.m)
%
%   INTEGRATION NOTE: this variant returns a rational model, not a single
%   coefficient matrix like run_standard_sindy.m/run_weak_sindy.m --
%   benchmarking/simulate_trajectory.m's RHS builder currently assumes
%   dxdt = Theta(x)*Xi and does not yet know how to evaluate a
%   numerator/denominator pair. Forward-simulating an implicit model
%   needs a small dedicated RHS (evaluate build_polynomial_library at the
%   current state for both numerator and denominator, divide) before this
%   variant can be scored by run_grid_search.m/select_best_model.m
%   alongside the other two. That RHS is the natural next piece to add
%   once this variant needs benchmarking, not before.

if nargin < 4 || isempty(poly_order)
    poly_order = 2;
end
if nargin < 3 || isempty(lambda)
    lambda = 0.1;
end

[N, D] = size(X);
if ~isequal(size(dXdt), [N, D])
    error('run_implicit_sindy:sizeMismatch', ...
        'dXdt must be the same size as X (%d x %d).', N, D);
end

[Theta_poly, library_names] = build_library(X, poly_order);
M = size(Theta_poly, 2);

rational = struct('numerator_Xi', {}, 'denominator_Xi', {}, 'min_abs_denominator', {});

for d = 1:D
    dxdt_d = dXdt(:, d);

    % Every poly term times dxdt_d, EXCLUDING the trivial "1*dxdt_d"
    % column (index 1), which is the term we fixed to 1 and moved to the
    % RHS to avoid the zero solution.
    interaction = Theta_poly(:, 2:end) .* dxdt_d;

    L_d = [Theta_poly, interaction];  % (N x (M + (M-1)))
    Y_d = -dxdt_d;

    Xi_d = stlsq_solve(L_d, Y_d, lambda, 10);

    numerator_Xi = Xi_d(1:M);
    denominator_Xi = [1; Xi_d(M+1:end)]; % constant term fixed by construction, not thresholded

    denom_values = Theta_poly * denominator_Xi;
    min_abs_denominator = min(abs(denom_values));

    rational(d).numerator_Xi = numerator_Xi;
    rational(d).denominator_Xi = denominator_Xi;
    rational(d).min_abs_denominator = min_abs_denominator;

    if min_abs_denominator < 1e-3
        warning('run_implicit_sindy:nearSingularDenominator', ...
            ['Dimension %d: fitted denominator gets as small as %.4g over the observed data -- ' ...
             'this model has a near-singularity within the training range and forward simulation ' ...
             'will likely be unreliable.'], d, min_abs_denominator);
    end
end

end