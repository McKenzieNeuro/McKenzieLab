classdef Test_CP_NMU < matlab.unittest.TestCase
    properties
        R
        X
        M_true
    end

    methods (TestMethodSetup)
        function setup(testCase)
            % Create a simple problem for testing CP-NMU
            rng('default');
            sz = [10 8 6];
            testCase.R = 2;
            ns = 0.01; 
            info = create_problem('Size', sz, ...
                'Num_Factors', testCase.R, ...
                'Noise', ns, ...
                'Factor_Generator', @rand, ...
                'Lambda_Generator', @rand);
            testCase.X = info.Data;
            testCase.M_true = info.Soln;
        end
    end

    methods (Test)
        function testDefaultRun(testCase)
            [M, ~, out] = cp_nmu(testCase.X, testCase.R,'printitn',0);
            scr = score(M, testCase.M_true);
            % fprintf('DefaultRun Score: %.4f\n', scr);
            testCase.verifyGreaterThan(scr, 0.95);
            testCase.verifyGreaterThan(out.final_fit, 0.95);
        end

        function testCustomInit(testCase)
            sz = size(testCase.X); 
            R_local = testCase.R;
            Uinit = cell(3,1);
            for n = 1:3
                Uinit{n} = rand(sz(n), R_local);
            end
            [M, U0, out] = cp_nmu(testCase.X, R_local, 'init', Uinit,'printitn',0);
            scr = score(M, testCase.M_true);
            % fprintf('DefaultRun Score: %.4f\n', scr);
            testCase.verifyGreaterThan(scr, 0.95);
            testCase.verifyGreaterThan(out.final_fit, 0.95);
            testCase.verifyEqual(U0,Uinit);

            % Check that it produces the same resutl with a ktensor version of U0.
            kt_init = ktensor(Uinit);
            [M2, U0_kt, out2] = cp_nmu(testCase.X, R_local, 'init', kt_init,'printitn',0);
            scr2 = score(M2, testCase.M_true);
            % fprintf('Ktensor Init Score: %.4f\n', scr2);
            testCase.verifyGreaterThan(scr2, 0.95);
            testCase.verifyGreaterThan(out2.final_fit, 0.95);
            testCase.verifyEqual(M, M2);

        end

    end

    
end
