classdef Test_CP_ALS < matlab.unittest.TestCase
    % Test_CP_ALS - Unit tests for cp_als.m in Tensor Toolbox
    % Covers scenarios from cp_als_doc.m and cp_als.m help

    methods (Test)
        function testRandomInitBasic(testCase)
            % Basic call with random initialization
			rng(0);
            K = ktensor(@rand, [4 5 3], 2);
			X = full(K);
            M = cp_als(X, 2, 'printitn', 0);
            testCase.verifyClass(M, 'ktensor');
            testCase.verifyEqual(ndims(M), 3);
            testCase.verifyEqual(ncomponents(M), 2);
			relerr = norm(full(X) - full(M)) / norm(full(X));
			testCase.verifyLessThanOrEqual(relerr, 0.02);
        end

        function testDimorderOption(testCase)
			rng(0);
            K = ktensor(@rand, [4 5 3], 2);
			X = full(K);
			dimorders = perms(1:3);
			for i = 1:size(dimorders,1)
				M = cp_als(X, 2, 'dimorder', dimorders(i,:), 'printitn', 0);
				testCase.verifyClass(M, 'ktensor');
				relerr = norm(full(X) - full(M)) / norm(full(X));
				testCase.verifyLessThanOrEqual(relerr, 0.02);
			end
        end

        function testNvecsInit(testCase)
            rng(0);
			K = ktensor(@rand, [5 4 3], 2);
			X = full(K);
            M = cp_als(X, 2, 'init', 'nvecs', 'printitn', 0);
            relerr = norm(full(X) - full(M)) / norm(full(X));
            testCase.verifyLessThanOrEqual(relerr, 0.02);
        end

        function testCellArrayInit(testCase)
			rng(0);
            K = ktensor(@rand, [5 4 3], 2);
			X = full(K);
            U0 = {rand(5,2), rand(4,2), []};
            M = cp_als(X, 2, 'dimorder', [3 2 1], 'init', U0, 'printitn', 0);
            relerr = norm(full(X) - full(M)) / norm(full(X));
            testCase.verifyLessThanOrEqual(relerr, 0.02);
        end

        function testKtensorInit(testCase)
            rng(0);
            K = ktensor(@rand, [5 4 3], 2);
			X = full(K);            
            M = cp_als(X, 2, 'init', K, 'printitn', 0);
			relerr = norm(full(X) - full(M)) / norm(full(X));
			testCase.verifyLessThanOrEqual(relerr, 1e-10);
        end

        function testMaxitersAndTol(testCase)
            rng(0);
            K = ktensor(@rand, [5 4 3], 2);
			X = full(K);
            M = cp_als(X, 2, 'maxiters', 100, 'tol', 1e-8, 'printitn', 0);
			relerr = norm(full(X) - full(M)) / norm(full(X));
			testCase.verifyLessThanOrEqual(relerr, 0.02);
        end
        
        function testFixsignsOption(testCase)
            rng(0);
            K = ktensor(@rand, [5 4 3], 2);
			X = full(K);
            M = cp_als(X, 2, 'fixsigns', false, 'printitn', 0);
			relerr = norm(full(X) - full(M)) / norm(full(X));
			testCase.verifyLessThanOrEqual(relerr, 0.02);
		end

        function testReproducibility(testCase)
            rng(0);
            K = ktensor(@rand, [5 4 3], 2);
			X = full(K);
            [M1, U0, out] = cp_als(X, 2, 'dimorder', [3 2 1], 'init', 'random', 'printitn', 0);
			params = out.params;
			params.init = U0;
            M2 = cp_als(X, 2, params); % same params as previous run
            testCase.verifyEqual(M1,M2);
        end

        function testSparseTensorInput(testCase)
            rng(0);
            X = sptenrand([5 4 3], 10); 
            M = cp_als(X, 2, 'printitn', 0);
            testCase.verifyClass(M, 'ktensor');
            testCase.verifyEqual(ndims(M), 3);
            testCase.verifyEqual(ncomponents(M), 2);
            relerr = norm(full(X) - full(M)) / norm(full(X));
            testCase.verifyGreaterThanOrEqual(relerr, 0);
            testCase.verifyLessThanOrEqual(relerr, 1);
        end

        function testKtensorInput(testCase)
            rng(0);
            K = ktensor(@rand, [5 4 3], 2);
            M = cp_als(K, 2, 'printitn', 0);
            testCase.verifyClass(M, 'ktensor');
            testCase.verifyEqual(ndims(M), 3);
            testCase.verifyEqual(ncomponents(M), 2);
            relerr = norm(full(K) - full(M)) / norm(full(K));
            testCase.verifyLessThanOrEqual(relerr, 0.02);
        end

    end
end
