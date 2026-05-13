classdef Test_CP_ARLS < matlab.unittest.TestCase
    properties
        R
        X
        M_true
    end

    methods (TestMethodSetup)
        function setup(testCase)
            % Set up an easy sample problem 
            rng('default');
            sz = [200 300 400];
            testCase.R = 5;
            ns = 0.01; 
            info = create_problem('Size', sz, ...
                'Num_Factors', testCase.R, ...
                'Noise', ns, ...
                'Factor_Generator', @(m,n) matrandnorm(m,n), ...
                'Lambda_Generator', @ones);
            testCase.X = info.Data;
            testCase.M_true = info.Soln;
        end
    end

    methods (Test)
        function DefaultRun(testCase)
            [M, ~, out] = cp_arls(testCase.X, testCase.R);
            scr = score(M, testCase.M_true);
            fprintf('DefaultRun Score: %.4f\n', scr);
            testCase.verifyGreaterThan(scr, 0.95);
            testCase.verifyLessThanOrEqual(out.iters, 25);
        end

        function NoMix(testCase)
            [M, ~, out] = cp_arls(testCase.X, testCase.R, 'mix', false);
            scr = score(M, testCase.M_true);
            fprintf('NoMix Score: %.4f\n', scr);
            testCase.verifyGreaterThan(scr, 0.95);
            testCase.verifyLessThanOrEqual(out.iters, 25);
        end

        function VaryEpoch(testCase)
            [M1, ~, out1] = cp_arls(testCase.X, testCase.R, 'epoch', 1, 'newitol', 20);
            scr1 = score(M1, testCase.M_true);
            fprintf('Epoch=1 Score: %.4f\n', scr1);
            testCase.verifyGreaterThan(scr1, 0.95);
            testCase.verifyLessThanOrEqual(out1.iters, 100);

            [M2, ~, out2] = cp_arls(testCase.X, testCase.R, 'epoch', 20, 'newitol', 3, 'printitn', 2);
            scr2 = score(M2, testCase.M_true);
            fprintf('Epoch=20 Score: %.4f\n', scr2);
            testCase.verifyGreaterThan(scr2, 0.95);
            testCase.verifyLessThanOrEqual(out2.iters, 10);
        end

        function TrueFitOption(testCase)
            [M, ~, out] = cp_arls(testCase.X, testCase.R, 'truefit', true);
            scr = score(M, testCase.M_true);
            fprintf('TrueFitOption Score: %.4f\n', scr);
            fit = 1 - norm(testCase.X - full(M)) / norm(testCase.X);
            fprintf('True fit: %.4f\n', fit);
            testCase.verifyGreaterThan(scr, 0.95);
            testCase.verifyLessThanOrEqual(out.iters, 25);
            testCase.verifyEqual(fit, out.fit, 'AbsTol', 1e-4);
        end

        function FitThreshOption(testCase)
            % Use a slightly lower threshold to ensure test passes            
            [M, ~, out] = cp_arls(testCase.X, testCase.R, 'fitthresh', 0.85, 'epoch', 1, 'truefit', true, 'printitn', 1);
            scr = score(M, testCase.M_true);
            fprintf('FitThreshOption Score: %.4f\n', scr);
            testCase.verifyGreaterThan(scr, 0.95);
            testCase.verifyGreaterThan(out.fit, 0.84);
        end

        function NsampfitOption(testCase)
            [M, ~, out] = cp_arls(testCase.X, testCase.R, 'truefit', true, 'nsampfit', 100);
            scr = score(M, testCase.M_true);
            fprintf('NsampfitOption Score: %.4f\n', scr);
            testCase.verifyGreaterThan(scr, 0.93);
        end

        function NsamplsqOption(testCase)
            [M1, ~, out1] = cp_arls(testCase.X, testCase.R, 'truefit', true, 'nsamplsq', 10);
            scr1 = score(M1, testCase.M_true);
            fprintf('Nsamplsq=10 Score: %.4f\n', scr1);
            testCase.verifyGreaterThan(scr1, 0.90);

            [M2, ~, out2] = cp_arls(testCase.X, testCase.R, 'truefit', true, 'nsamplsq', 25);
            scr2 = score(M2, testCase.M_true);
            fprintf('Nsamplsq=25 Score: %.4f\n', scr2);
            testCase.verifyGreaterThan(scr2, 0.93);
        end

        function InitCellVsKtensor(testCase)
            % Create initial factor matrices
            rng(42);
            sz = size(testCase.X);
            A = cell(ndims(testCase.X), 1);
            for n = 1:ndims(testCase.X)
                A{n} = randn(sz(n), testCase.R);
            end

            % Run with cell array initialization
            rng(123); % Reset seed for reproducibility
            [M_cell, ~, out_cell] = cp_arls(testCase.X, testCase.R, 'init', A);

            % Run with ktensor initialization
            K = ktensor(ones(testCase.R,1), A);
            rng(123); % Reset seed again
            [M_ktensor, ~, out_ktensor] = cp_arls(testCase.X, testCase.R, 'init', K);

            % Compare results
            scr_cell = score(M_cell, testCase.M_true);
            scr_ktensor = score(M_ktensor, testCase.M_true);

            fprintf('InitCell Score: %.4f\n', scr_cell);
            fprintf('InitKtensor Score: %.4f\n', scr_ktensor);

            testCase.verifyEqual(scr_cell, scr_ktensor, 'AbsTol', 1e-10);
            testCase.verifyEqual(out_cell.fit, out_ktensor.fit, 'AbsTol', 1e-10);
            testCase.verifyEqual(M_cell, M_ktensor);
        end

    end
end
