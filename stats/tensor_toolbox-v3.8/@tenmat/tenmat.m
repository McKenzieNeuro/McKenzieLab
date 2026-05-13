%TENMAT Store tensor as a matrix.
%
%TENMAT Methods:
%   ctranspose - Complex conjugate transpose for tenmat.
%   disp       - Command window display of a matricized tensor (tenmat).
%   display    - Command window display of a tenmat.
%   double     - Convert tenmat to double array.
%   end        - Last index of indexing expression for tenmat.
%   minus      - Binary subtraction (-) for tenmat.
%   mtimes     - Multiplies two tenmat objects.
%   norm       - Frobenius norm of a tenmat.
%   plus       - Binary addition (+) for tenmat. 
%   size       - Size of tenmat.
%   subsasgn   - Subscripted assignment for tenmat.  
%   subsref    - Subscripted reference for tenmat.
%   tenmat     - Create a matricized tensor.
%   tsize      - Tensor size of tenmat.
%   uminus     - Unary minus (-) for tenmat.
%   uplus      - Unary plus (+) for tenmat.
%
%   <a href="matlab:web(strcat('file://',fullfile(getfield(what('tensor_toolbox'),'path'),'doc','html','tenmat_doc.html')))">Documentation page for tensor-as-matrix class</a>
%
%   See also TENSOR_TOOLBOX
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>

classdef tenmat

    properties
        tsize    %< Size of the original tensor.
        rindices %< Row indices in the original tensor.
        cindices %< Column indices in the original tensor.
        data     %< Matrix data.
    end

    methods
        function A = tenmat(varargin)
            %TENMAT Create a matricized tensor.
            %
            %   A = TENMAT(T, RDIMS) creates a matrix representation of a tensor
            %   T.  The dimensions (or modes) specified in RDIMS map to the rows
            %   of the matrix, and the remaining dimensions (in ascending order)
            %   map to the columns.
            %
            %   A = TENMAT(T, CDIMS, 't') does the same as above, but instead the
            %   column dimensions are specified, and the remaining dimensions (in
            %   ascending order) map to the rows.
            %
            %   A = TENMAT(T, RDIMS, CDIMS) creates a matrix representation of
            %   tensor T.  The dimensions specified in RDIMS map to the rows of
            %   the matrix, and the dimensions specified in CDIMS map to the
            %   columns, in the order given.
            %
            %   A = TENMAT(T, RDIM, STR) creates the same matrix representation as
            %   above, except only one dimension in RDIM maps to the rows of the
            %   matrix, and the remaining dimensions span the columns in an order
            %   specified by the string argument STR as follows:
            %
            %     'fc' - Forward cyclic.  Order the remaining dimensions in the
            %            columns by [RDIM+1:ndims(T), 1:RDIM-1].  This is the
            %            ordering defined by Kiers.
            %
            %     'bc' - Backward cyclic.  Order the remaining dimensions in the
            %            columns by [RDIM-1:-1:1, ndims(T):-1:RDIM+1].  This is the
            %            ordering defined by De Lathauwer, De Moor, and Vandewalle.
            %
            %   A = TENMAT(A, RDIMS, CDIMS, TSIZE) creates a tenmat from a matrix
            %   A along with the mappings of the row (RDIMS) and column indices
            %   (CDIMS) and the size of the original tensor (TSIZE).
            %
            %   A = TENMAT(B) is the copy constructor for B also a tenmat.
            %
            %   A = TENMAT is the empty constructor.
            %
            %   See also TENMAT.
            %
            %Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>


            %
            % Includes improvements offered by Marcus Brubaker.


            % Case 0a: Empty Contructor
            if (nargin == 0)
                A.tsize = [];
                A.rindices = [];
                A.cindices = [];
                A.data = [];
                return;
            end

            % Case 0b: Copy Constructor
            if (nargin == 1)
                B = varargin{1};
                if isa(B, 'tenmat')
                    A.tsize = B.tsize;
                    A.rindices = B.rindices;
                    A.cindices = B.cindices;
                    A.data = B.data;
                else
                    error('Single argument constructor expects a tenmat object.');
                end
                return;
            end

            % Case I: Matrix to tenmat
            if (nargin == 4)
                data_in = varargin{1};
                if ~isnumeric(data_in) || (ndims(data_in) ~= 2)
                    error('A must be a matrix.');
                end
                rdims = varargin{2};
                cdims = varargin{3};
                tsize_in = varargin{4};
                n = numel(tsize_in);
                if ~isequal(1:n, sort([rdims cdims]))
                    error('Incorrect specification of dimensions');
                elseif (prod(tsize_in(rdims)) ~= size(data_in,1))
                    error('SIZE(A,1) does not match size specified by RDIMS and SIZE.');
                elseif (prod(tsize_in(cdims)) ~= size(data_in,2))
                    error('SIZE(A,2) does not match size specified by CDIMS and SIZE.');
                end
                A.tsize = tsize_in;
                A.rindices = rdims;
                A.cindices = cdims;
                A.data = data_in;
                return;
            end

            % Case II: MDA to tenmat (via tensor)
            if isa(varargin{1},'double')
                temp_A = tenmat(tensor(varargin{1}),varargin{2:nargin});
                A.tsize = temp_A.tsize;
                A.rindices = temp_A.rindices;
                A.cindices = temp_A.cindices;
                A.data = temp_A.data;
                return;
            end

            % Case III: tensor to tenmat
            if (nargin < 2)  ||  (nargin > 3)
                error('Incorrect number of arguments.');
            end
            T = varargin{1};
            if ~isa(T,'tensor')
                error('First argument must be a tensor, a tenmat, a matrix, or empty.');
            end
            tsize_val = size(T);
            n = ndims(T);
            if (nargin == 2)
                rdims = varargin{2};
                tmp = true(1,n); 
                tmp(rdims) = false; 
                cdims = find(tmp);
            elseif isa(varargin{3},'char')
                switch varargin{3}
                    case 't' % Transpose
                        cdims = varargin{2};
                        tmp = true(1,n); 
                        tmp(cdims) = false; 
                        rdims = find(tmp);
                    case 'fc' % Forward cyclic
                        rdims = varargin{2};
                        if (numel(rdims) ~= 1)
                            error('Only one row dimension if third argument is ''fc''.');
                        end
                        cdims = [rdims+1:n, 1:rdims-1];
                    case 'bc' % Backward cyclic
                        rdims = varargin{2};
                        if (numel(rdims) ~= 1)
                            error('Only one row dimension if third argument is ''bc''.');
                        end
                        cdims = [rdims-1:-1:1, n:-1:rdims+1];
                    otherwise
                        error('Unrecognized option');
                end
            else
                rdims = varargin{2};
                cdims = varargin{3};
            end
            if ~isequal(1:n, sort([rdims cdims]))
                error('Incorrect specification of dimensions');
            end
            data_val = reshape(double(permute(T,[rdims cdims])), prod(tsize_val(rdims)), prod(tsize_val(cdims)));
            A.tsize = tsize_val;
            A.rindices = rdims;
            A.cindices = cdims;
            A.data = data_val;
        end

        function b = saveobj(a)
            %SAVEOBJ Save a tenmat object.
            b.tsize = a.tsize;
            b.rindices = a.rindices;
            b.cindices = a.cindices;
            b.data = a.data;
        end
    end

    methods (Static)
        function t = loadobj(s)
            %LOADOBJ Load a tenmat object.
            if isa(s,'tenmat')
                t = s;
            else
                t = tenmat(s.data, s.rindices, s.cindices, s.tsize);
            end
        end
    end
end
