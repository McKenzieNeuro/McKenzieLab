classdef test_sptensor < matlab.unittest.TestCase
    methods (Test)
        function testEmptyConstructor(testCase)
            % Test the empty constructor
            t = sptensor();
            testCase.verifyEmpty(t.subs);
            testCase.verifyEmpty(t.vals);
            testCase.verifyEmpty(t.size);
            testCase.verifyEqual(t.type, 'sparse');
        end
        
        function testSizeConstructor(testCase)
            % Test constructor with size
            sz = [3, 3, 3];
            t = sptensor(sz);
            testCase.verifyEmpty(t.subs);
            testCase.verifyEmpty(t.vals);
            testCase.verifyEqual(t.size, sz);
            testCase.verifyEqual(t.type, 'sparse');
        end
        
        function testSubsValsConstructor(testCase)
            % Test constructor with subscripts and values
            subs = [1 1 1; 2 2 2];
            vals = [1; 2];
            sz = [3, 3, 3];
            t = sptensor(subs, vals, sz);
            testCase.verifyEqual(t.subs, subs);
            testCase.verifyEqual(t.vals, vals);
            testCase.verifyEqual(t.size, sz);
            testCase.verifyEqual(t.type, 'sparse');
        end
        
        function testCopyConstructor(testCase)
            % Test copy constructor
            subs = [1 1 1; 2 2 2];
            vals = [1; 2];
            sz = [3, 3, 3];
            t1 = sptensor(subs, vals, sz);
            t2 = sptensor(t1);
            testCase.verifyEqual(t2.subs, subs);
            testCase.verifyEqual(t2.vals, vals);
            testCase.verifyEqual(t2.size, sz);
            testCase.verifyEqual(t2.type, 'sparse');
        end
        
        function testSparseMatrixConstructor(testCase)
            % Test constructor with sparse matrix
            A = sparse([1 0 0; 0 2 0; 0 0 3]);
            t = sptensor(A);
            testCase.verifyEqual(t.subs, [1 1; 2 2; 3 3]);
            testCase.verifyEqual(t.vals, [1; 2; 3]);
            testCase.verifyEqual(t.size, [3 3]);
            testCase.verifyEqual(t.type, 'sparse');
        end
        
        function testTensorConstructor(testCase)
            % Test constructor with tensor
            A = rand(3, 3, 3);
            t = sptensor(tensor(A));
            [subs, vals] = find(tensor(A));
            testCase.verifyEqual(t.subs, subs);
            testCase.verifyEqual(t.vals, vals);
            testCase.verifyEqual(t.size, size(A));
            testCase.verifyEqual(t.type, 'sparse');
        end
        
        function testRandomTensorConstructor(testCase)
            % Test constructor with random tensor
            fh = @rand;
            sz = [3, 3, 3];
            nv = 5;
            t = sptensor(fh, nv, sz);
            testCase.verifyEqual(t.size, sz);
            testCase.verifyEqual(t.type, 'sparse');
            testCase.verifySize(t.subs, [nv, 3]);
            testCase.verifySize(t.vals, [nv, 1]);
        end
        
        function testAccumulateConstructor(testCase)
            % Test constructor with accumulation function
            subs = [1 1 1; 1 1 1; 2 2 2];
            vals = [1; 2; 3];
            sz = [3, 3, 3];
            t = sptensor(subs, vals, sz, @sum);
            testCase.verifyEqual(t.subs, [1 1 1; 2 2 2]);
            testCase.verifyEqual(t.vals, [3; 3]);
            testCase.verifyEqual(t.size, sz);
            testCase.verifyEqual(t.type, 'sparse');
        end
        
        function testIncompleteTensorConstructor(testCase)
            % Test constructor with incomplete tensor
            subs = [1 1 1; 2 2 2];
            vals = [1; 0];
            sz = [3, 3, 3];
            t = sptensor(subs, vals, sz, [], 'incomplete');
            testCase.verifyEqual(t.subs, subs);
            testCase.verifyEqual(t.vals, vals);
            testCase.verifyEqual(t.size, sz);
            testCase.verifyEqual(t.type, 'incomplete');
        end
        
        function testIncompleteTensorWithAccumulation(testCase)
            % Test constructor with incomplete tensor and accumulation function
            subs = [1 1 1; 1 1 1; 2 2 2];
            vals = [1; 2; 0];
            sz = [3, 3, 3];
            t = sptensor(subs, vals, sz, @sum, 'incomplete');
            testCase.verifyEqual(t.subs, [1 1 1; 2 2 2]);
            testCase.verifyEqual(t.vals, [3; 0]);
            testCase.verifyEqual(t.size, sz);
            testCase.verifyEqual(t.type, 'incomplete');
        end
        
        function testZeroValsSparse(testCase)
            % Test that zeros are removed for 'sparse' type
            subs = [1 1 1; 2 2 2; 3 3 3];
            vals = [1; 0; 2];
            sz = [3, 3, 3];
            t = sptensor(subs, vals, sz);
            testCase.verifyEqual(t.subs, [1 1 1; 3 3 3]);
            testCase.verifyEqual(t.vals, [1; 2]);
            testCase.verifyEqual(t.size, sz);
            testCase.verifyEqual(t.type, 'sparse');
        end

        function testZeroValsIncomplete(testCase)
            % Test that zeros are retained for 'incomplete' type
            subs = [1 1 1; 2 2 2; 3 3 3];
            vals = [1; 0; 2];
            sz = [3, 3, 3];
            t = sptensor(subs, vals, sz, [], 'incomplete');
            testCase.verifyEqual(t.subs, subs);
            testCase.verifyEqual(t.vals, vals);
            testCase.verifyEqual(t.size, sz);
            testCase.verifyEqual(t.type, 'incomplete');
        end

        function testFullFunctionIncomplete(testCase)
            % Test full function for incomplete tensor
            subs = [1 1 1; 2 2 2];
            vals = [5; 7];
            sz = [3, 3, 3];
            t = sptensor(subs, vals, sz, [], 'incomplete');
            % The full function should fill missing entries with NaN for incomplete tensors
            A = full(t);
            expected = nan(sz);
            expected(1,1,1) = 5;
            expected(2,2,2) = 7;
            expected = tensor(expected);
            testCase.verifyEqual(A, expected);
        end

        function testFullFunctionSparse(testCase)
            % Test full function for sparse tensor
            subs = [1 1 1; 2 2 2];
            vals = [5; 7];
            sz = [3, 3, 3];
            t = sptensor(subs, vals, sz);
            % The full function should fill missing entries with 0 for sparse tensors
            A = full(t);
            expected = zeros(sz);
            expected(1,1,1) = 5;
            expected(2,2,2) = 7;
            expected = tensor(expected);
            testCase.verifyEqual(A, expected);
        end

        function testDoubleSparse(testCase)
            % Test double conversion for sparse tensor
            subs = [1 1 1; 2 2 2];
            vals = [5; 7];
            sz = [3, 3, 3];
            t = sptensor(subs, vals, sz);
            A = double(t);
            expected = zeros(sz);
            expected(1,1,1) = 5;
            expected(2,2,2) = 7;
            testCase.verifyEqual(A, expected);
        end

        function testDoubleIncomplete(testCase)
            % Test double conversion for incomplete tensor
            subs = [1 1 1; 2 2 2];
            vals = [5; 7];
            sz = [3, 3, 3];
            t = sptensor(subs, vals, sz, [], 'incomplete');
            A = double(t);
            expected = nan(sz);
            expected(1,1,1) = 5;
            expected(2,2,2) = 7;
            testCase.verifyEqual(A, expected);
        end

        function testDoubleEmpty(testCase)
            % Test double conversion for empty tensor
            sz = [2, 2, 2];
            t = sptensor(sz);
            A = double(t);
            expected = zeros(sz);
            testCase.verifyEqual(A, expected);
        end

        function testDoubleScalar(testCase)
            % Test double conversion for scalar tensor
            subs = [1 1 1];
            vals = 42;
            sz = [1 1 1];
            t = sptensor(subs, vals, sz);
            A = double(t);
            expected = zeros(sz);
            expected(1,1,1) = 42;
            testCase.verifyEqual(A, expected);
        end

        function testDoubleAllZeroSparse(testCase)
            % Test double conversion for all-zero sparse tensor
            sz = [2, 2];
            t = sptensor([], [], sz);
            A = double(t);
            expected = zeros(sz);
            testCase.verifyEqual(A, expected);
        end

        function testDoubleAllZeroIncomplete(testCase)
            % Test double conversion for all-zero incomplete tensor
            sz = [2, 2];
            t = sptensor([], [], sz, [], 'incomplete');
            A = double(t);
            expected = nan(sz);
            testCase.verifyEqual(A, expected);
        end

        function testSaveLoad(testCase)
            % Test saveobj and loadobj methods
            subs = [1 2; 2 1];
            vals = [10; 20];
            sz = [2 2];
            t = sptensor(subs, vals, sz, [], 'incomplete');
            s = saveobj(t);
            t2 = sptensor.loadobj(s);
            testCase.verifyEqual(t2.subs, t.subs);
            testCase.verifyEqual(t2.vals, t.vals);
            testCase.verifyEqual(t2.size, t.size);
            testCase.verifyEqual(t2.type, t.type);
        end

        function testScalarValsExpansion(testCase)
            % Test scalar value expansion
            subs = [1 1 1; 2 2 2; 3 3 3];
            vals = 4;
            sz = [3, 3, 3];
            t = sptensor(subs, vals, sz);
            testCase.verifyEqual(t.subs, subs);
            testCase.verifyEqual(t.vals, [4; 4; 4]);
            testCase.verifyEqual(t.size, sz);
        end

    end
end