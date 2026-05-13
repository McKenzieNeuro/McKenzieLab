classdef Test_CP_OPT < matlab.unittest.TestCase
    % Test_CP_OPT - Unit tests for cp_opt.m in Tensor Toolbox
    % Covers scenarios from cp_opt_doc.m and cp_opt.m help

    properties
        R;
        X;
        M_true;
        M_init;
    end

    methods (TestMethodSetup)
        function setup(testCase)
            % Create a reproducible test problem for CP-OPT
            rng(1);
            testCase.R = 5;
            problem = create_problem('Size', [50 40 30], 'Num_Factors', testCase.R, 'Noise', 0.10);
            testCase.X = problem.Data;
            testCase.M_true = problem.Soln;
            testCase.M_init = create_guess('Data', testCase.X, 'Num_Factors', testCase.R, 'Factor_Generator', 'nvecs');
        end
    end

    methods (Test)
        function testBasicSolve(testCase)
            % Basic call with nvecs initialization
            [M, M0, info] = cp_opt(testCase.X, testCase.R, 'init', testCase.M_init, 'printitn', 0);
            testCase.verifyClass(M, 'ktensor');
            testCase.verifyEqual(ndims(M), 3);
            testCase.verifyEqual(ncomponents(M), testCase.R);
            fit = 1 - sqrt(info.f);
            testCase.verifyGreaterThan(fit, 0.85); % Should be close to 0.9 due to 10% noise
        end

        function testReproducibility(testCase)
            % Test that using info.params reproduces the result
            [M1, M0, info] = cp_opt(testCase.X, testCase.R, 'init', testCase.M_init, 'printitn', 0);
            params = info.params;
            params.init = M0;
            M2 = cp_opt(testCase.X, testCase.R, params);
            testCase.verifyEqual(M1, M2);
        end

        function testScoreAgainstTruth(testCase)
            % Score the solution against the true model
            [M, ~, ~] = cp_opt(testCase.X, testCase.R, 'init', testCase.M_init, 'printitn', 0);
            scr = score(M, testCase.M_true);
            testCase.verifyGreaterThan(scr, 0.85); % Should be close to 1, but not perfect due to noise
        end

		function testDifferentMethodFminunc(testCase)
			% Only run if Optimization Toolbox is available
			if ~license('test', 'Optimization_Toolbox')
				testCase.assumeFail('Optimization Toolbox is not available.');
			end
			% Test with fminunc method
			[M, ~, info] = cp_opt(testCase.X, testCase.R, 'init', testCase.M_init, 'printitn', 0, 'method', 'fminunc');
			testCase.verifyClass(M, 'ktensor');
			fit = 1 - sqrt(info.f);
			testCase.verifyGreaterThan(fit, 0.85);
		end

        function testHigherRankGuess(testCase)
            % Test with R+1 factors (overestimate rank)
            rng(2);
            M_plus_init = create_guess('Data', testCase.X, 'Num_Factors', testCase.R+1, 'Factor_Generator', 'nvecs');
            [M_plus, ~, info] = cp_opt(testCase.X, testCase.R+1, 'init', M_plus_init, 'printitn', 0);
            fit = 1 - sqrt(info.f);
            testCase.verifyGreaterThan(fit, 0.85);
            scr = score(M_plus, testCase.M_true);
            testCase.verifyGreaterThan(scr, 0.7); % Should still be reasonable
        end

        function testNonnegativeOption(testCase)
            % Test with nonnegative factors using lower bound
            rng(3);
            problem2 = create_problem('Size', [50 40 30], 'Num_Factors', testCase.R, 'Noise', 0.10, ...
                'Factor_Generator', 'rand', 'Lambda_Generator', 'rand');
            X2 = problem2.Data .* (problem2.Data > 0); % Force nonnegative
            M_true2 = problem2.Soln;
            [M, ~, info] = cp_opt(X2, testCase.R, 'init', 'rand', 'lower', 0, 'printitn', 0);
            fit = 1 - sqrt(info.f);
            testCase.verifyGreaterThan(fit, 0.85);
            scr = score(M, M_true2);
            testCase.verifyGreaterThan(scr, 0.7);
        end

        function testDifferentMethodLbfgs(testCase)
            % Test with lbfgs method (POBLANO, pure MATLAB)
            [M, ~, info] = cp_opt(testCase.X, testCase.R, 'init', testCase.M_init, 'printitn', 0, 'method', 'lbfgs');
            testCase.verifyClass(M, 'ktensor');
            fit = 1 - sqrt(info.f);
            testCase.verifyGreaterThan(fit, 0.85);
        end

        function testDifferentMethodCompact(testCase)
            % Test with compact method (Limited-Memory Compact Representation)
            [M, ~, info] = cp_opt(testCase.X, testCase.R, 'init', testCase.M_init, 'printitn', 0, 'method', 'compact');
            testCase.verifyClass(M, 'ktensor');
            fit = 1 - sqrt(info.f);
            testCase.verifyGreaterThan(fit, 0.85);
        end
    end
end
