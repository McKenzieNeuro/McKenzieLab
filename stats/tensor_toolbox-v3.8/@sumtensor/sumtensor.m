%SUMTENSOR Class for implicit sum of other tensors.
%
%SUMTENSOR Methods:
%   disp      - Command window display of a sumtensor.
%   display   - Command window display of a sumtensor.
%   double    - Convert sumtensor to double array.
%   full      - Convert a sumtensor to a (dense) tensor.
%   innerprod - Efficient inner product with a sumtensor.
%   isscalar  - False for sumtensors.
%   mttkrp    - Matricized tensor times Khatri-Rao product for sumtensor.
%   ndims     - Return the number of dimensions for a sumtensor.
%   norm      - Frobenius norm of a sumtensor.
%   plus      - Plus for sumtensor.
%   size      - Size of a sumtensor.
%   subsref   - Subscript reference for sumtensor.
%   sumtensor - Tensor stored as sum of tensors.
%   ttv       - Tensor times vector for sumtensor.
%   uminus    - Unary minus for sumtensor.
%   uplus     - Unary plus for sumtensor.
%
%   <a href="matlab:web(strcat('file://',fullfile(getfield(what('tensor_toolbox'),'path'),'doc','html','sumtensor_doc.html')))">Documentation page for sum of tensors class</a>
%
%   See also TENSOR_TOOLBOX
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>

classdef sumtensor

    properties
        part %< Cell array of tensors to be summed.
    end

    methods
        function t = sumtensor(varargin)
            %SUMTENSOR Tensor stored as sum of tensors.
            %
            %   The SUMTENSOR class is limited to certain operations that easily
            %   decompose as sums: INNERPROD, MTTKRP, TTV. Note that the NORM function
            %   is not easily computed for a SUMTENSOR.
            %
            %   T = SUMTENSOR(T1,T2,...) creates a tensor that is the sum of its
            %   constituent parts. The tensor is stored implicitly, i.e., each
            %   component is retained. This may lead to storage and computation
            %   efficiency. All input tensors must be the same size, but they can be
            %   any type of tensor. 
            %
            %   T = SUMTENSOR(S) creates a SUMTENSOR by copying an existing
            %   SUMTENSOR.
            %
            %   T = SUMTENSOR(C) where C is a cell array of tensors creates a
            %   SUMTENSOR from the tensors in C. 
            %
            %   T = SUMTENSOR is the empty constructor.
            %
            %   Examples
            %   T1 = tensor(rand(4,3,3));
            %   T2 = sptensor([1 1 1; 3 1 2; 4 3 3], 1, [4,3,3]);
            %   T = sumtensor(T1,T2); %<--A sumtensor with parts T1 and T2
            %
            %   See also TENSOR, SPTENSOR, TTENSOR, KTENSOR
            %
            %Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>

            % Empty constructor
            if (nargin == 0)
                t.part = cell(0);
                return;
            end

            % Copy constructor
            if (nargin == 1) && isa(varargin{1}, 'sumtensor')
                t.part = varargin{1}.part;
                return;
            end

            % Determine the actual parts to process
            parts_to_process = {};
            if (nargin == 1) && iscell(varargin{1})
                % Case: sumtensor(cell_array_of_tensors) - e.g., from loadobj
                parts_to_process = varargin{1};
            else
                % Case: sumtensor(T1, T2, ...) or sumtensor(T1)
                parts_to_process = varargin;
            end
            
            if isempty(parts_to_process)
                t.part = cell(0); % Handle sumtensor({}) or sumtensor(empty_cell_from_loadobj)
                return;
            end

            % Process the parts
            t.part = cell(length(parts_to_process),1);
            first_tensor_size = [];
            for i = 1:length(parts_to_process)
                current_item = parts_to_process{i};
                cl = class(current_item);
                if ismember(cl,'double') % Convert an MDA
                    current_item = tensor(current_item);
                elseif ~ismember(cl, {'tensor','sptensor','ktensor','ttensor'})
                    error('Inputs must be tensors. Symtensors are not supported.');
                end
                
                if (i == 1) % In MATLAB, length(parts_to_process) >= 1 here
                    first_tensor_size = size(current_item);
                else
                    if ~isequal(size(current_item), first_tensor_size)
                        error('All inputs must be the same size.');
                    end
                end
                t.part{i} = current_item;  
            end
        end

        function b = saveobj(a)
            %SAVEOBJ Save a sumtensor object.
            %
            %   B = SAVEOBJ(A) is called by SAVE when a sumtensor object is
            %   saved. The result B is a struct containing the sumtensor data.
            %
            %   See also SUMTENSOR/LOADOBJ, SAVE, LOAD.
            %
            %Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>
            
            b.part = a.part;
        end

    end % methods

    methods (Static)
        function t = loadobj(s)
            %LOADOBJ Load a sumtensor object.
            %
            %   T = LOADOBJ(S) is called by LOAD when a sumtensor object is
            %   loaded. S is a structure containing the sumtensor data and T
            %   is the restored sumtensor object. If S is an object, then the
            %   object is simply returned.
            %
            %   See also SUMTENSOR/SAVEOBJ, SAVE, LOAD.
            %
            %Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>
            
            if isa(s,'sumtensor')
                t = s;
            else
                % s is a struct, create a new sumtensor object
                % The constructor is designed to handle a cell array as a single argument
                t = sumtensor(s.part);
            end
        end
    end % static methods
end
