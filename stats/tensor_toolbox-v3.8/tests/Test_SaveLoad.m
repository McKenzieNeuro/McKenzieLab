classdef Test_SaveLoad < matlab.unittest.TestCase
    methods (Test)
        function testSaveLoadTensor(testCase)
            % Create a tensor object
            data = rand(3, 3, 3);
            t = tensor(data);

            % Save the tensor object
            save('tensor_test.mat', 't');

            % Load the tensor object
            loadedData = load('tensor_test.mat');
            t_loaded = loadedData.t;

            % Verify the loaded object is of class tensor
            testCase.verifyClass(t_loaded, 'tensor');

            % Verify the loaded tensor object
            testCase.verifyEqual(t.data, t_loaded.data);
            testCase.verifyEqual(t.size, t_loaded.size);

            % Clean up
            delete('tensor_test.mat');
        end

        function testSaveLoadSptensor(testCase)
            % Create an sptensor object
            subs = [1 1 1; 2 2 2];
            vals = [1; 2];
            sz = [3, 3, 3];
            t = sptensor(subs, vals, sz);

            % Save the sptensor object
            save('sptensor_test.mat', 't');

            % Load the sptensor object
            loadedData = load('sptensor_test.mat');
            t_loaded = loadedData.t;

            % Verify the loaded object is of class sptensor
            testCase.verifyClass(t_loaded, 'sptensor');

            % Verify the loaded sptensor object
            testCase.verifyEqual(t.subs, t_loaded.subs);
            testCase.verifyEqual(t.vals, t_loaded.vals);
            testCase.verifyEqual(t.size, t_loaded.size);

            % Clean up
            delete('sptensor_test.mat');
        end

        function testSaveLoadKtensor(testCase)
            % Create a ktensor object
            lambda = [1; 2];
            U = {rand(3, 2), rand(3, 2), rand(3, 2)};
            t = ktensor(lambda, U);

            % Save the ktensor object
            save('ktensor_test.mat', 't');

            % Load the ktensor object
            loadedData = load('ktensor_test.mat');
            t_loaded = loadedData.t;

            % Verify the loaded object is of class ktensor
            testCase.verifyClass(t_loaded, 'ktensor');

            % Verify the loaded ktensor object
            testCase.verifyEqual(t.lambda, t_loaded.lambda);
            testCase.verifyEqual(t.u, t_loaded.u);

            % Clean up
            delete('ktensor_test.mat');
        end

        function testSaveLoadTenmat(testCase)
            % Create a tenmat object
            T = tensor(rand(3,3,3));
            A = tenmat(T, [1 2]);
            save('tenmat_test.mat', 'A');
            loaded = load('tenmat_test.mat');
            A_loaded = loaded.A;
            testCase.verifyClass(A_loaded, 'tenmat');
            testCase.verifyEqual(A.tsize, A_loaded.tsize);
            testCase.verifyEqual(A.rindices, A_loaded.rindices);
            testCase.verifyEqual(A.cindices, A_loaded.cindices);
            testCase.verifyEqual(A.data, A_loaded.data);
            delete('tenmat_test.mat');
        end

        function testSaveLoadSptenmat(testCase)
            % Create a sptenmat object
            T = sptensor([1 1 1; 2 2 2], [1; 2], [3 3 3]);
            A = sptenmat(T, [1], [2 3]);
            save('sptenmat_test.mat', 'A');
            loaded = load('sptenmat_test.mat');
            A_loaded = loaded.A;
            testCase.verifyClass(A_loaded, 'sptenmat');
            testCase.verifyEqual(A.tsize, A_loaded.tsize);
            testCase.verifyEqual(A.rdims, A_loaded.rdims);
            testCase.verifyEqual(A.cdims, A_loaded.cdims);
            testCase.verifyEqual(A.subs, A_loaded.subs);
            testCase.verifyEqual(A.vals, A_loaded.vals);
            delete('sptenmat_test.mat');
        end

        function testSaveLoadTtensor(testCase)
            % Create a ttensor object
            core = tensor(rand(2,2,2));
            U = {rand(3,2), rand(3,2), rand(3,2)};
            t = ttensor(core, U);
            save('ttensor_test.mat', 't');
            loaded = load('ttensor_test.mat');
            t_loaded = loaded.t;
            testCase.verifyClass(t_loaded, 'ttensor');
            testCase.verifyEqual(t.core, t_loaded.core);
            testCase.verifyEqual(t.u, t_loaded.u);
            delete('ttensor_test.mat');
        end

        function testSaveLoadSumtensor(testCase)
            % Create a sumtensor object
            T1 = tensor(rand(3,3,3));
            T2 = tensor(rand(3,3,3));
            t = sumtensor(T1, T2);
            save('sumtensor_test.mat', 't');
            loaded = load('sumtensor_test.mat');
            t_loaded = loaded.t;
            testCase.verifyClass(t_loaded, 'sumtensor');
            testCase.verifyEqual(numel(t.part), numel(t_loaded.part));
            for i = 1:numel(t.part)
                testCase.verifyEqual(t.part{i}, t_loaded.part{i});
            end
            delete('sumtensor_test.mat');
        end

        function testSaveLoadSymktensor(testCase)
            % Create a symktensor object
            lambda = [1; 2];
            U = rand(3,2);
            m = 3;
            t = symktensor(lambda, U, m);
            save('symktensor_test.mat', 't');
            loaded = load('symktensor_test.mat');
            t_loaded = loaded.t;
            testCase.verifyClass(t_loaded, 'symktensor');
            testCase.verifyEqual(t.lambda, t_loaded.lambda);
            testCase.verifyEqual(t.u, t_loaded.u);
            testCase.verifyEqual(t.m, t_loaded.m);
            delete('symktensor_test.mat');
        end

        function testSaveLoadSymtensor(testCase)
            % Create a symtensor object
            vals = (1:10)';
            m = 2;
            n = 4;
            t = symtensor(vals, m, n);
            save('symtensor_test.mat', 't');
            loaded = load('symtensor_test.mat');
            t_loaded = loaded.t;
            testCase.verifyClass(t_loaded, 'symtensor');
            testCase.verifyEqual(t.val, t_loaded.val);
            testCase.verifyEqual(t.m, t_loaded.m);
            testCase.verifyEqual(t.n, t_loaded.n);
            delete('symtensor_test.mat');
        end
    end
end