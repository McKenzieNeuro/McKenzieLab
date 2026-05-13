%SPTENSOR Class for sparse and incomplete tensors.
%
%SPTENSOR Methods:
%   and          - Logical AND (&) for sptensors.
%   collapse     - Collapse sptensor values along specified dimensions.
%   contract     - Contract sparse tensor along two dimensions (array trace).
%   disp         - Command window display of a sptensor.
%   display      - Command window display of a sptensor.
%   divide       - Divide an sptensor by a nonnegative KTENSOR.
%   double       - Converts a sptensor to a dense multidimensional array.
%   elemfun      - Manipulate the elements of a sptensor.
%   end          - Last index of indexing expression for sptensor.
%   eq           - Equal (==) for sptensors.
%   fibers       - Extracts specified mode-k fibers and creates matrix.
%   find         - Find subscripts of nonzero elements in a sptensor.
%   findices     - Compute mode-k unfolding column index for every value.
%   full         - Convert a sptensor to a (dense) tensor.
%   ge           - Greater than or equal for sptensors.
%   gt           - Greater than for sptensors.
%   innerprod    - Efficient inner product with a sparse tensor.
%   isequal      - Compare sptensors for equality.
%   isincomplete - True if X is an incomplete tensor (versus sparse).
%   isscalar     - False for sptensors.
%   issparse     - True if X is a sparse tensor (versus incomplete).
%   ldivide      - Array right division for sparse tensors.
%   le           - Less than or equal for sptensors.
%   lt           - Less than for sptensors.
%   mask         - Extract values as specified by a mask tensor.
%   minus        - Binary subtraction for sparse tensors. 
%   mldivide     - Slash left division for sparse tensors.
%   mrdivide     - Slash right division for sparse tensors.
%   mtimes       - sptensor-scalar multiplication.
%   mttkrp       - Matricized tensor times Khatri-Rao product for sparse tensor.
%   ndims        - Number of dimensions of a sptensor.
%   ne           - Not equal (~=) for sptensors.
%   nnz          - Number of values in sptensor.
%   norm         - Frobenius norm of a sparse tensor.
%   not          - Logical NOT (~) for sptensors.
%   nvecs        - Compute the leading mode-n vectors for a sparse tensor.
%   ones         - Replace values of sptensor with ones.
%   or           - Logical OR (|) for sptensors.
%   permute      - Rearrange the dimensions of a sptensor.
%   plus         - Binary addition for sparse tensors. 
%   rdivide      - Array right division for sparse tensors.
%   reshape      - Reshape sptensor.
%   rrf          - Produce matrix via sparse randomized range finder in mode-k.
%   scale        - Scale along specified dimensions for sptensors.
%   size         - Dimensions of sptensor.
%   spmatrix     - Converts a two-way sparse tensor to sparse matrix.
%   spones       - Replace sptensor elements with ones.
%   sptensor     - Create a sparse or incomplete tensor.
%   squash       - Remove empty slices from a sptensor.
%   squeeze      - Remove singleton dimensions from a sptensor.
%   subsasgn     - Subscripted assignment for sptensor.
%   subsref      - Subscripted reference for a sptensor.
%   times        - Array multiplication for sptensors.
%   ttm          - Sparse tensor times matrix.
%   ttt          - Sparse tensor times sparse tensor.
%   ttv          - Sparse tensor times vector.
%   uminus       - Unary minus (-) for sptensor.
%   uplus        - Unary plus (+) for sptensor.
%   xor          - Logical XOR for sptensors.
%
%   <a href="matlab:web(strcat('file://',fullfile(getfield(what('tensor_toolbox'),'path'),'doc','html','sptensor_doc.html')))">Documentation page for sptensor Class</a>
%
%   See also TENSOR_TOOLBOX
%
%   How to cite the sptensor class:
%   * B.W. Bader and T.G. Kolda. Efficient MATLAB Computations with Sparse
%     and Factored Tensors, SIAM J. Scientific Computing, 30:205-231, 2007,
%     http://dx.doi.org/10.1137/060676489.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>

classdef sptensor

    properties
        subs,
        vals,
        size,
        type
    end

    methods

        function t = sptensor(varargin)
            %SPTENSOR Create a sparse or incomplete tensor.
            %
            %   X = SPTENSOR(SUBS, VALS, SZ, ACC, TYPE) creates a sparse or incomplete
            %   tensor as follows (using [] yields the default value):
            %
            %    o SUBS - p x n array specifying the subscripts of the values.
            %    o VALS - p x 1 array of values or a scalar for all values.
            %    o SZ   - 1 x n array specifying the size. Default: max(SUBS,[],1).
            %    o ACC  - function to accumulate repeats. Default: @sum.
            %    o TYPE - 'sparse' or 'incomplete'. Default: 'sparse'.
            %
            %   X = SPTENSOR(Y) or X = SPTENSOR(Y,TYPE) copies/converts Y if it is
            %   another compatible object. Note that a MATLAB row vector will be
            %   interpreted as a size (see previous constructor).
            %
            %   X = SPTENSOR(FUN,P,SZ,TYPE) uses FUN to create the values in a
            %   sptensor with randomly generated unique subscripts. The value P can be
            %   an integer or a proportion of entries.
            %
            %   Examples
            %
            %   % Setup
            %   subs = [1 1 1; 1 2 3; 2 1 1; 1 1 1; 3 1 1; 4 2 1; 3 1 1];
            %   vals = [1 0 2 4 0 2 1]';
            %   siz = [4 4 4];
            %
            %   % Sparse 4 x 4 x 4 tensor: zeros ignored, repeats summed
            %   X = sptensor(subs,vals,siz)
            %
            %   % Incomplete 4 x 4 x 4 tensor: repeats summed
            %   X = sptensor(subs,vals,siz,[],'incomplete')
            %
            %   % Set every value to be 1 (repeats summed)
            %   X = sptensor(subs,1,siz)
            %
            %   % Use max accumulation instead
            %   X = sptensor(subs,1,[],@max)
            %
            %   % Sparse 4 x 4 x 4 tensor: zeros ignored, min of repeats
            %   X = sptensor(subs,vals,siz,@min)
            %
            %   % Sparse and incomplete tensors with random generation
            %   X = sptensor(@rand,[10 10 10],0.01) %<- Proportion
            %   X = sptensor(@randn,[10 10 10],10,[],'incomplete') %<- Integer
            %
            %   See also SPTENSOR, SPTENRAND.
            %
            %Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>

            % Other help..
            %
            %   X = SPTENSOR(SZ) abbreviates X = SPTENSOR([],[],SZ).
            %
            %   X = SPTENSOR(FUN,SZ,P,...) would also work via backwards
            %   compatibitlity.
            %
            %   X = SPTENSOR(SUBS, VALS, SZ, 0) or
            %   X = SPTENSOR(SUBS, VALS, SZ, TYPE, 0) skips the checks for repeats and
            %   just trusts that everything is correct. Use with caution.


            % Defaults!
            t.subs = [];
            t.vals = [];
            t.size = [];
            t.type = 'sparse';

            % EMPTY Constructor
            if (nargin == 0) || ((nargin == 1) && isempty(varargin{1}))
                %t = class(t,'sptensor');
                return;
            end

            % SINGLE ARGUMENT
            if (nargin == 1)

                source = varargin{1};

                switch(class(source))

                    % COPY CONSTRUCTOR
                    case 'sptensor'
                        t.subs = source.subs;
                        t.vals = source.vals;
                        t.size = source.size;
                        t.type = source.type;
                        %t = class(t, 'sptensor');
                        return;

                        % CONVERT SPTENMAT
                    case 'sptenmat'

                        % Extract the tensor size and order
                        siz = source.tsize;

                        if isempty(source.subs) %There are no nonzero terms
                            subs = [];
                        else % Convert the 2d-subscipts into nd-subscripts
                            if ~isempty(source.rdims)
                                subs(:,source.rdims) = ...
                                    tt_ind2sub(siz(source.rdims),source.subs(:,1));
                            end
                            if ~isempty(source.cdims)
                                subs(:,source.cdims) = ...
                                    tt_ind2sub(siz(source.cdims),source.subs(:,2));
                            end
                        end
                        % Copy the values (which do not need to be modified)
                        vals = source.vals;

                        % Store everything
                        t.subs = subs;
                        t.vals = vals;
                        t.size = siz;
                        %t = class(t, 'sptensor');
                        return;

                        % CONVERT TENSOR
                    case 'tensor'
                        [subs,vals] = find(source);
                        t.subs = subs;
                        t.vals = vals;
                        t.size = size(source);
                        %t = class(t, 'sptensor');
                        return;

                        % SPARSE MATRIX, SIZE, or MDA
                    case {'numeric','logical','double'}

                        % Case 1: SPARSE MATRIX
                        if issparse(source)
                            [i,j,s] = find(source);
                            siz = size(source);
                            t.subs = [i,j];
                            t.vals = s;
                            t.size = siz;
                            %t = class(t,'sptensor');
                            return;
                        end

                        % Case 2: SPECIFYING THE SIZE
                        if tt_sizecheck(source)
                            t.size = source;
                            %t = class(t, 'sptensor');
                            return;
                        end

                        % Case 3: An MDA
                        t = sptensor(tensor(source));
                        return;

                end % switch

            end % nargin == 1

            % SPECIAL CASE (explicit fifth argument as zero) for INTERACTION WITH MEX
            % FILES OR DIRECT CREATION OF SPTENSOR WITHOUT ANY SORTING OR OTHER
            % STANDARD CHECKS
            if (nargin == 5) && (isnumeric(varargin{5})) && (varargin{5} == 0)

                % Store everything
                t.subs = varargin{1};
                t.vals = varargin{2};
                t.size = varargin{3};
                t.type = varargin{4};

                return;

            end

            % SAVED FOR BACKWARDS COMPATIBILITY
            if (nargin == 4) && (isnumeric(varargin{4})) && (varargin{4} == 0)

                % Store everything
                t.subs = varargin{1};
                t.vals = varargin{2};
                t.size = varargin{3};
                t.type = 'sparse';

                return;

            end

            % Create with function handle and size
            if (nargin >= 3) && isa(varargin{1},'function_handle')
                fh = varargin{1};

                % backwards compatibility - it used to be that size and nz were swapped
                if isscalar(varargin{3})
                    nv = varargin{3};
                    sz = varargin{2};
                else
                    nv = varargin{2};
                    sz = varargin{3};
                end

                if (nargin >= 4) && ~isempty(varargin{4})
                    type = varargin{4};
                    typecheck(type);
                else
                    type = 'sparse';
                end

                if (nv <= 0) || (nv >= prod(sz))
                    error('Requested number of values must be positive and less than the total size')
                elseif (nv < 1)
                    nv = ceil(prod(sz) * nv);
                else
                    nv = floor(nv);
                end

                % Keep iterating until we find enough unique nonzeros or we give up
                subs = [];
                cnt = 0;
                while (size(subs,1) < nv) && (cnt < 10)
                    newsubs = ceil( rand(nv, size(sz,2)) * diag(sz) );
                    subs = unique([subs; newsubs], 'rows');
                    cnt = cnt + 1;
                end

                if size(subs,1) < nv
                    warning(['Could only generate %d unique subscripts rather' ...
                        ' than the desired %d'], size(subs,1), nv);
                    nv = min(nv, size(subs,1));
                end

                % eliminate extra subscripts, if any
                subs = subs(1:nv,:);

                % Generate the values
                vals = fh(nv,1);

                if strcmp(type,'sparse') && any(vals==0)

                        zidx = find(vals == 0);
                    zcnt = length(zidx);
                    warning(['Generation function created %d zeros, reducing the ' ...
                        'number of final values to %d'], zcnt, nv-zcnt);

                    % Eliminate any zeros
                    idx = find(vals);
                    subs = subs(idx,:);
                    vals = vals(idx);

                end

                

                % Store everything
                t.subs = subs;
                t.vals = vals;
                t.size = sz;
                t.type = type;

                % Create the tensor
                %t = class(t, 'sptensor');
                return;
            end

            % CONVERT A SET OF INPUTS
            if nargin >= 2

                % Extract the subscripts and values
                subs = varargin{1};
                vals = varargin{2};

                tt_subscheck(subs);
                tt_valscheck(vals);
                if isscalar(vals)
                    vals = vals * ones(size(subs,1),1);
                end
                if (size(vals,1) ~= size(subs,1))
                    error('Number of subscripts and values must be equal');
                end

                % Extract the size
                if nargin < 3 || isempty(varargin{3})
                    siz = max(subs,[],1);
                else
                    siz = varargin{3};
                    tt_sizecheck(siz);
                end

                % Check for wrong input
                if ~isempty(subs) && size(subs,2) ~= size(siz,2)
                    error('Number of subscripts does not match size')
                end

                % Check for subscripts out of range
                for j = 1:numel(siz)
                    if ~isempty(subs) && max(subs(:,j)) > siz(j)
                        error('Subscript exceeds sptensor size')
                    end
                end

                % Extract the 'combiner' function handle
                if (nargin< 4) || isempty(varargin{4})
                    fun = @sum;
                else
                    fun = varargin{4};
                end

                % Extract the type
                if (nargin < 5) || isempty(varargin{5})
                    type = 'sparse';
                else
                    type = varargin{5};
                    typecheck(type);
                end

                if isempty(subs)
                    newsubs = [];
                    newvals = [];
                else
                    % Identify only the unique indices
                    [newsubs,~,loc] = unique(subs,'rows');

                    % Accumulate repeated values
                    newvals = accumarray(loc,vals,[size(newsubs,1) 1],fun);
                end

                % Remove any zero if sptensor
                if strcmp(type,'sparse')
                    nzidx = find(newvals);
                    newsubs = newsubs(nzidx,:);
                    newvals = newvals(nzidx);
                end

                % Store everything
                t.subs = newsubs;
                t.vals = newvals;
                t.size = siz;
                t.type = type;

                % Create the tensor
                %t = class(t, 'sptensor');

                return;
            end

            error('Unsupported use of function SPTENSOR.');

        end

        function s = saveobj(obj)
            %SAVEOBJ Save the sptensor object for MAT-file.
            %   S = SAVEOBJ(OBJ) saves the sptensor object OBJ for use with
            %   the load function. The saved object S is a structure with the
            %   fields 'subs', 'vals', 'size', and 'type', which contain the
            %   subscripts, values, size, and type of the sptensor, respectively.
            s.subs = obj.subs;
            s.vals = obj.vals;
            s.size = obj.size;
            s.type = obj.type;
        end
    end

    methods (Static)

        function obj = loadobj(s)
            %LOADOBJ Load the sptensor object from MAT-file.
            %   OBJ = LOADOBJ(S) loads a sptensor object from the structure S
            %   created by the SAVEOBJ method. The structure S must contain
            %   fields 'subs', 'vals', 'size', and 'type', which are used to
            %   reconstruct the sptensor object. If S is not a structure,
            %   it is returned directly as the sptensor object.
            if isstruct(s)
                obj = sptensor(s.subs,s.vals,s.size,s.type,0);
            else
                obj = s;
            end

        end
    end

end
