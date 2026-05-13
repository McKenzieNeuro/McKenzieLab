classdef Test_CreateProblem < matlab.unittest.TestCase
    % Test_CreateProblem Test cases for create_problem and create_problem_binary
    %   Tests for the standard and binary problem creation functions
    
    methods (Test)
        
        function testDefaultCreateProblem(testCase)
            % Test default parameters for create_problem
            rng('default');
            info = create_problem();
            
            % Check that the solution is a ktensor by default
            testCase.verifyClass(info.Soln, 'ktensor');
            
            % Data should be a tensor
            testCase.verifyClass(info.Data, 'tensor');
            
            % Verify sizes match
            testCase.verifyEqual(size(info.Data), size(info.Soln));
            
            % Default size should be [100 100 100]
            testCase.verifyEqual(size(info.Data), [100 100 100]);
            
            % Default number of factors should be 2
            testCase.verifyEqual(size(info.Soln.lambda, 1), 2);

            % Check error is within specified default
            M = info.Soln;
            Mfull = full(M);
            X = info.Data;
            relerr = norm(X - Mfull) / norm(Mfull);
            testCase.verifyLessThanOrEqual(relerr, 0.1+eps);
        end
        
        function testTuckerProblem(testCase)
            % Test creating a Tucker problem
            info = create_problem('Type', 'Tucker', 'Size', [50 40 30], 'Num_Factors', [5 4 3]);
            
            % Verify that the solution is a ttensor
            testCase.verifyClass(info.Soln, 'ttensor');
            
            % Check correct sizes
            testCase.verifyEqual(size(info.Data), [50 40 30]);
            testCase.verifyEqual(size(info.Soln.core), [5 4 3]);
            
            % Check factor matrix sizes
            testCase.verifyEqual(size(info.Soln.U{1}), [50 5]);
            testCase.verifyEqual(size(info.Soln.U{2}), [40 4]);
            testCase.verifyEqual(size(info.Soln.U{3}), [30 3]);

            % Check error is within specified default
            M = info.Soln;
            Mfull = full(M);
            X = info.Data;
            relerr = norm(X - Mfull) / norm(Mfull);
            testCase.verifyLessThanOrEqual(relerr, 0.1+eps);

        end
        
        function testExplicitSoln(testCase)
            % Test providing an explicit solution
            rng('default');
            sz = [10 15 20];
            lambda = [1; 0.5];
            U = {rand(sz(1), 2), rand(sz(2), 2), rand(sz(3), 2)};
            K = ktensor(lambda, U);
            info = create_problem('Soln', K, 'Noise', 0);
            
            % Verify solution is preserved
            testCase.verifyEqual(info.Soln, K);
            
            % Verify data tensor matches solution (no noise)
            testCase.verifyEqual(double(info.Data), double(full(K)), 'AbsTol', 1e-12);
        end
        
        function testMissingData(testCase)
            % Test problem with missing data
            rng('default');
            info = create_problem('Size', [30 40 50], 'M', 0.3, 'Noise', 0);
            
            % Verify pattern is a tensor or sptensor
            testCase.verifyTrue(isa(info.Pattern, 'tensor') || isa(info.Pattern, 'sptensor'));
            
            % Verify pattern has correct size
            testCase.verifyEqual(size(info.Pattern), [30 40 50]);
            
            % Verify approximately 30% of data is missing
            pattern_values = double(info.Pattern);
            missing_ratio = 1 - nnz(pattern_values) / numel(pattern_values);
            testCase.verifyEqual(missing_ratio, 0.3, 'RelTol', 0.05);
            
            % Verify data tensor has zeros at missing entries
            data_values = double(info.Data);
            testCase.verifyEqual(data_values(pattern_values == 0), zeros(nnz(pattern_values == 0), 1));

            % Verify non-missing entries match the solution
            M = info.Soln;
            Mfull = full(M);
            non_missing_indices = find(pattern_values > 0);
            testCase.verifyEqual(data_values(non_missing_indices), Mfull(non_missing_indices), 'AbsTol', 1e-12);
        end
        
        function testSparseGeneration(testCase)
            % Test sparse generation option
            rng('default');
            
            % Factors must be nonnegative for sparse generation
            info = create_problem('Size', [100 100 100], 'Sparse_Generation', 1000, 'Factor_Generator', 'rand'); 
            
            % Verify data is sparse
            testCase.verifyClass(info.Data, 'sptensor');
            
            % Verify number of nonzeros is approximately the requested number
            testCase.verifyEqual(nnz(info.Data), 1000);
        end
        
        function testSymmetric(testCase)
            % Test symmetric problem
            rng('default');
            info = create_problem('Size', [40 40 30], 'Symmetric', [1 2]);
            
            % Verify factor matrices for symmetric modes are identical
            testCase.verifyEqual(info.Soln.U{1}, info.Soln.U{2});
            
            % Verify data tensor is symmetric in the specified modes
            X = info.Data;
            testCase.verifyTrue(issymmetric(X,[1 2]))
        end
        
        function testFactorGenerators(testCase)
            % Test different factor generators
            generators = {'rand', 'randn', 'orthogonal', 'stochastic'};
            
            for i = 1:length(generators)
                rng('default');
                info = create_problem('Size', [10 10 10], 'Factor_Generator', generators{i});
                                
                % Check specific properties based on generator type
                if strcmp(generators{i}, 'stochastic')
                    % Column sums should be approximately 1
                    for j = 1:3
                        col_sums = sum(info.Soln.U{j}, 1);
                        testCase.verifyEqual(col_sums, ones(1, 2), 'AbsTol', 1e-10);
                    end
                elseif strcmp(generators{i}, 'orthogonal')
                    % Columns should be orthogonal
                    for j = 1:3
                        U = info.Soln.U{j};
                        testCase.verifyEqual(U'*U, eye(2), 'AbsTol', 1e-10);
                    end
                end
            end
        end
        
        function testDefaultCreateProblemBinary(testCase)
            % Test default parameters for create_problem_binary
            rng('default');
            sz = [10 15 20];
            r = 3;
            [X, M] = create_problem_binary(sz, r);
            
            % Verify outputs
            testCase.verifyClass(X, 'sptensor');
            testCase.verifyClass(M, 'ktensor');
            
            % Verify sizes
            testCase.verifyEqual(size(X), sz);
            testCase.verifyEqual(size(M), sz);
            
            % Verify rank
            testCase.verifyEqual(length(M.lambda), r);
            
            % Verify binary values
            testCase.verifyEqual(unique(X.vals), 1);
        end
        
        function testCreateProblemBinaryParams(testCase)
            % Test various parameters for create_problem_binary
            rng('default');
            sz = [10 15 20];
            r = 3;
            opts = {'loprob', 0.05, 'hiprob', 0.8, 'density', 0.2, 'verbosity', 0};
            [X, M, info] = create_problem_binary(sz, r, opts{:});
            
            % Verify correct params were stored
            testCase.verifyEqual(info.params.loprob, 0.05);
            testCase.verifyEqual(info.params.hiprob, 0.8);
            testCase.verifyEqual(info.params.density, 0.2);
            
            % Verify outputs
            testCase.verifyClass(X, 'sptensor');
            testCase.verifyClass(M, 'ktensor');
        end
        
        function testCreateProblemBinaryUserSpecifiedM(testCase)
            % Test create_problem_binary with user-specified M
            rng('default');
            sz = [10 15 20];
            r = 3;
            A = cell(3,1);
            for k = 1:3
                A{k} = rand(sz(k), r);
            end
            userM = ktensor(ones(r,1), A);
            
            [X, M, info] = create_problem_binary(sz, r, 'Mtrue', userM);
            
            % Verify M matches user-specified M
            testCase.verifyEqual(M, userM);
        end
        
    end
end
