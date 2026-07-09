function opts = get_resolution_preset(level)
%GET_RESOLUTION_PRESET Named grid-density presets for run_grid_search.m/
%run_full_search.m, so sweep resolution can be dialed with one argument
%instead of hand-tuning every grid separately.
%
%   opts = GET_RESOLUTION_PRESET(level)
%
%   level: 'fast' | 'balanced' | 'thorough'
%     'fast'     - coarsest grids, fewest combinations, quickest
%                  turnaround. Good for iterating on the pipeline itself
%                  or a first pass on new data.
%     'balanced' - the previous hard-coded defaults. Reasonable middle
%                  ground for routine use. (Default if unspecified.)
%     'thorough' - denser grids, more combinations, best chance of
%                  finding a genuinely optimal model at real cost in
%                  runtime. Use once you already trust the pipeline and
%                  want the best answer, not while iterating.
%
%   Approximate total combinations across all 3 variants (standard + weak
%   + implicit) with the default LibraryFlavorGrid, for a 2-channel
%   system like Hes1:
%     fast     ~ 40-50 combinations
%     balanced ~ 300-350 combinations
%     thorough ~ 1200-1300 combinations
%   Weak-SINDy's window_points x test_function_order dimensions are what
%   drive most of the growth between levels -- if 'thorough' is still too
%   slow, dropping WindowPointsGrid/TestFunctionOrderGrid back down is the
%   single highest-leverage manual override.
%
%   Output: opts struct with fields matching run_grid_search.m's
%   name-value options: LambdaGrid, PolyOrderGrid, WindowPointsGrid,
%   TestFunctionOrderGrid, HillCoeffGrid, HillKCandidates, NumDensePoints,
%   LibraryFlavorGrid.

switch level
    case 'fast'
        opts.LambdaGrid = logspace(-2, 0, 3);
        opts.PolyOrderGrid = [1 2];
        opts.WindowPointsGrid = 20;
        opts.TestFunctionOrderGrid = 4;
        opts.HillCoeffGrid = [1 2];
        opts.HillKCandidates = 2;
        opts.NumDensePoints = 200;
        opts.LibraryFlavorGrid = {'poly_only', 'hill_mixed'};

    case 'balanced'
        opts.LambdaGrid = logspace(-2, 0, 5);
        opts.PolyOrderGrid = [1 2 3];
        opts.WindowPointsGrid = [15 25];
        opts.TestFunctionOrderGrid = [2 4];
        opts.HillCoeffGrid = [1 2 4];
        opts.HillKCandidates = 3;
        opts.NumDensePoints = 300;
        opts.LibraryFlavorGrid = {'poly_only', 'hill_mixed'};

    case 'thorough'
        opts.LambdaGrid = logspace(-2, 0, 8);
        opts.PolyOrderGrid = [1 2 3];
        opts.WindowPointsGrid = [15 22 30];
        opts.TestFunctionOrderGrid = [2 4 6];
        opts.HillCoeffGrid = [1 2 4];
        opts.HillKCandidates = 4;
        opts.NumDensePoints = 400;
        opts.LibraryFlavorGrid = {'poly_only', 'hill_mixed'};

    otherwise
        error('get_resolution_preset:badLevel', ...
            'level must be ''fast'', ''balanced'', or ''thorough'' (got ''%s'').', level);
end

end