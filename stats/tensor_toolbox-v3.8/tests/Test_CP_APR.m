classdef Test_CP_APR < matlab.unittest.TestCase
    properties
			R; % Rank of the tensor
			X; % Data tensor
			M_true; % True solution tensor
	end
	
	methods (TestMethodSetup)
		function setup(testCase)
			% Create a problem for testing CP-APR that's easy to solve.
            rng('default')
			A = cell(3,1);
			A{1} = rand(10,2);
			A{1}(1:4:end,1) = 10;
			A{1}(3:4:end,2) = 10;
			A{2} = rand(8,2);
			A{2}(1:3:end,1) = 10;
			A{2}(2:3:end,2) = 10;
			A{3} = rand(6,2);
			A{3}(1:3:end,1) = 10;
			A{3}(2:3:end,2) = 10;
			lambda = [1; 1];
			sz = [10 8 6];
			S = ktensor(lambda, A);
			S = normalize(S,'sort',1);
			% Create sparse test problem based on provided solution.
			nz = 10000;
			info = create_problem('Soln', S, 'Sparse_Generation', nz);

			% Extract data and solution
			testCase.R = 2;
			testCase.X = info.Data;
			testCase.M_true = info.Soln;
        end
    end

    methods (Test)
        function BasicSolve(testCase)
			rng('default'); 
            M = cp_apr(testCase.X, testCase.R);
			testCase.verifyClass(M, 'ktensor');			
			scr = score(M, testCase.M_true, 'greedy', true);
			fprintf('Score: %.4f\n', scr);
			testCase.verifyGreaterThan(scr, 0.95, 'Score is below threshold');
        end

        function SolvePQNR(testCase)
            rng('default');
            M = cp_apr(testCase.X, testCase.R, 'alg', 'pqnr');
            testCase.verifyClass(M, 'ktensor');
            scr = score(M, testCase.M_true, 'greedy', true);
			fprintf('Score: %.4f\n', scr);
            testCase.verifyGreaterThan(scr, 0.95, 'PQNR: Score is below threshold');
        end

        function SolvePDNR(testCase)
            rng('default');
            M = cp_apr(testCase.X, testCase.R, 'alg', 'pdnr', 'maxiters', 50);
            testCase.verifyClass(M, 'ktensor');
            scr = score(M, testCase.M_true, 'greedy', true);
			fprintf('Score: %.4f\n', scr);
            testCase.verifyGreaterThan(scr, 0.95, 'PDNR: Score is below threshold');
        end

        function SolveMU(testCase)
            rng('default');
            M = cp_apr(testCase.X, testCase.R, 'alg', 'mu', 'maxiters', 200, 'kappa', 50);
            testCase.verifyClass(M, 'ktensor');
            scr = score(M, testCase.M_true, 'greedy', true);
            fprintf('Score: %.4f\n', scr);
            testCase.verifyGreaterThan(scr, 0.95, 'MU: Score is below threshold');
        end

        function SolveWithKtensorInit(testCase)
            rng('default');
            % Create a ktensor initial guess with correct size and nonnegative values
            sz = size(testCase.X);
            F = cell(ndims(testCase.X),1);
            for n = 1:ndims(testCase.X)
                F{n} = rand(sz(n), testCase.R);
            end
            Minit = ktensor(F);
            M = cp_apr(testCase.X, testCase.R, 'init', Minit);
            scr = score(M, testCase.M_true, 'greedy', true);
			fprintf('Score: %.4f\n', scr);
            testCase.verifyGreaterThan(scr, 0.95, 'ktensor init: Score is below threshold');
        end

        function SolveWithCellInit(testCase)
            rng('default');
            % Create a cell array initial guess with correct size and nonnegative values
            sz = size(testCase.X);
            F = cell(ndims(testCase.X),1);
            for n = 1:ndims(testCase.X)
                F{n} = rand(sz(n), testCase.R);
            end
            M = cp_apr(testCase.X, testCase.R, 'init', F);
            scr = score(M, testCase.M_true, 'greedy', true);
			fprintf('Score: %.4f\n', scr);
            testCase.verifyGreaterThan(scr, 0.95, 'cell array init: Score is below threshold');
        end
    end
end