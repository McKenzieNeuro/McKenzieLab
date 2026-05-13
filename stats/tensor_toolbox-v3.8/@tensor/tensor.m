%TENSOR Class for dense tensors.
%
%TENSOR Methods:
%   and         - Logical AND (&) for tensors.
%   collapse    - Collapse tensor along specified dimensions.
%   contract    - Contract tensor along two dimensions (array trace).
%   disp        - Command window display of a tensor.
%   display     - Command window display of a tensor.
%   double      - Convert tensor to double array.
%   end         - Last index of indexing expression for tensor.
%   eq          - Equal (==) for tensors.
%   exp         - Exponential for tensors.
%   fibers      - Extracts specified mode-k fibers and creates matrix.
%   find        - Find subscripts of nonzero elements in a tensor.
%   full        - Convert to a (dense) tensor.
%   ge          - Greater than or equal (>=) for tensors.
%   gt          - Greater than (>) for tensors.
%   innerprod   - Efficient inner product with a tensor.
%   isequal     - for tensors.
%   isscalar    - False for tensors.
%   issymmetric - Verify that a tensor X is symmetric in specified modes.
%   ldivide     - Left array divide for tensor.
%   le          - Less than or equal (<=) for tensor.
%   lt          - Less than (<) for tensor.
%   mask        - Extract values as specified by a mask tensor.
%   minus       - Binary subtraction (-) for tensors.
%   mldivide    - Slash left division for tensors.
%   mrdivide    - Slash right division for tensors.
%   mtimes      - tensor-scalar multiplication.
%   mttkrp      - Matricized tensor times Khatri-Rao product for tensor.
%   mttkrps     - Sequence of MTTKRP calculations for a tensor.
%   ndims       - Return the number of dimensions of a tensor.
%   ne          - Not equal (~=) for tensors.
%   nnz         - Number of nonzeros for tensors.
%   norm        - Frobenius norm of a tensor.
%   not         - Logical NOT (~) for tensors.
%   nvecs       - Compute the leading mode-n vectors for a tensor.
%   or          - Logical OR (|) for tensors.
%   permute     - Permute tensor dimensions.
%   plus        - Binary addition (+) for tensors.
%   power       - Elementwise power (.^) operator for a tensor.
%   rdivide     - Right array divide for tensors.
%   reshape     - Change tensor size.
%   scale       - Scale along specified dimensions of tensor.
%   size        - Tensor dimensions.
%   squeeze     - Remove singleton dimensions from a tensor.
%   subsasgn    - Subscripted assignment for a tensor.
%   subsref     - Subscripted reference for tensors.
%   symmetrize  - Symmetrize a tensor X in specified modes.
%   tenfun      - Apply a function to each element in a tensor.
%   tensor      - Create tensor.
%   times       - Array multiplication for tensors.
%   transpose   - is not defined on tensors.
%   ttm         - Tensor times matrix.
%   ttsv        - Tensor times same vector in multiple modes.
%   ttt         - Tensor mulitplication (tensor times tensor).
%   ttv         - Tensor times vector.
%   uminus      - Unary minus (-) for tensors.
%   unfold      - Unfold a tensor into a matrix.
%   uplus       - Unary plus (+) for tensors.
%   vec         - Vectorize a tensor.
%   xor         - Logical EXCLUSIVE OR for tensors.
%
%   <a href="matlab:web(strcat('file://',fullfile(getfield(what('tensor_toolbox'),'path'),'doc','html','tensor_doc.html')))">Documentation page for tensor class</a>
%
%   See also TENSOR_TOOLBOX
%
%   Reference:
%   * B.W. Bader and T.G. Kolda. Algorithm 862: MATLAB Tensor Classes for
%     Fast Algorithm Prototyping, ACM Trans. Mathematical Software,
%     32:635-653, 2006, http://dx.doi.org/10.1145/1186785.1186794.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>
classdef tensor
    properties
        data
        size
    end
    methods
        function t = tensor(input, sz)            
            %TENSOR Create tensor.
            %
            %   X = TENSOR(A,SZ) creates a tensor from the multidimensional array A.
            %   The SZ argument is a size vector specifying the desired shape
            %   of A.
            %
            %   X = TENSOR(F,SZ) createa a tensor of size SZ using the function
            %   handle F to create the data. The function F must take a size vector as
            %   input.
            %
            %   X = TENSOR(A) creates a tensor from the multidimensional array A, using
            %   SZ = size(A).
            %
            %   X = TENSOR(S) copies a tensor S.
            %
            %   X = TENSOR(A) converts an sptensor, ktensor, ttensor, or tenmat object
            %   to a tensor.
            %
            %   X = TENSOR creates an empty, dense tensor object.
            %
            %   Examples:
            %       X = tensor(rand(3, 4, 2)); % Tensor of size 3 x 4 x 2
            %       X = tensor(@rand, [3 4 2]); % Equivalent
            %       Y = tensor(zeros(3, 1), 3); % Tensor of size 3
            %       Y = tensor(@zeros, [3 1]);
            %       Z = tensor(ones(12, 1), [3 4 1]); % Tensor of size 3 x 4 x 1
            %       Z = tensor(@ones, [3 4 1]); % Equivalent
            %
            %   See also TENSOR, TENSOR/NDIMS.
            %
            %Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>
            arguments
                input = [] % Input data or function handle
                sz = [] % Size of tensor
            end

            % EMPTY/DEFAULT CONSTRUCTOR
            if isempty(input) && isempty(sz)
                t.data = [];
                t.size = [];
                return;
            end

            % Adapt to work with derived classes
            tmp = superclasses(input);
            if isempty(tmp)
                class_input = class(input);
            else
                class_input = tmp{1};
            end

            switch class_input
                
                case 'function_handle'
                    % FUNCTION HANDLE AND SIZE
                    fh = input;

                    % Check size
                    if ~isvector(sz)
                        error('Argument ''sz'' must be a vector');
                    end
                    if ~isrow(sz)
                        sz = sz';
                    end

                    % Generate data
                    t.data = fh([sz 1]);
                    t.size = sz;

                case 'tensor'
                    % COPY CONSTRUCTOR
                    t.data = input.data;
                    t.size = input.size;

                case {'ktensor','ttensor','sptensor','sumtensor','symtensor','symktensor'}
                    % CONVERSION
                    t = full(input);

                case 'tenmat'
                    % RESHAPE TENSOR-AS-MATRIX
                    % Here we just reverse what was done in the tenmat constructor.
                    % First we reshape the data to be an MDA, then we un-permute
                    % it using ipermute.
                    isz = tsize(input);
                    order = [input.rdims,input.cdims];
                    data = reshape(input.data, [isz(order) 1 1]);
                    if numel(order) >= 2
                        t.data = ipermute(data,order);
                    else
                        t.data = data;
                    end
                    t.size = isz;

                otherwise

                    if isnumeric(input) || islogical(input)
                        % CONVERT A MULTIDIMENSIONAL ARRAY

                        if isempty(sz)
                            sz = size(input);
                        elseif ~isvector(sz)
                            error('Argument ''sz'' must be a vector.');
                        elseif ~isrow(sz)
                            sz = sz';
                        end

                        t.data = input;
                        t.size = sz;

                    else

                        error('Unsupported use of function TENSOR.');
                    end
            end

            if isempty(t.size)
                if ~isempty(t.data)
                    error('Size is empty but data is not');
                end

            else
                if prod(t.size) ~= numel(t.data)
                    error('Number of elements in ''input'' does not match ''sz''.');
                end
                
                % Make sure the data is indeed the right shape
                t.data = reshape(t.data,[t.size 1]);

            end

            return;

        end % constructor

        function s = saveobj(obj)
            %SAVEOBJ Save tensor for MAT-file.
            %   S = SAVEOBJ(OBJ) saves the tensor object OBJ for use with
            %   the LOAD function. The output S is a structure with fields
            %   'data' and 'size', which contain the data and size of the
            %   tensor, respectively. This is used to ensure compatibility
            %   with future versions of the Tensor Toolbox.
            s.data = obj.data;
            s.size = obj.size;
        end

    end

    methods (Static)
        function obj = loadobj(s)
            %LOADOBJ Load tensor from MAT-file.
            %   OBJ = LOADOBJ(S) loads a tensor object from the structure S
            %   created by the SAVEOBJ method. The structure S must contain
            %   fields 'data' and 'size', which are used to reconstruct the
            %   tensor object.
            if isstruct(s)
                obj = tensor(s.data,s.size);
            else
                obj = s;
            end
        end
    end
end
