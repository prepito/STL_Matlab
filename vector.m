classdef vector < handle % matlab.mixin.Copyable % 
% vector class
% Imitates STL vector container.
% 
% Note: use .size_ instead of size to get the STL-vector-size. Creation of 
% matrix of vectors is supported by e.g.
%   vecArr = copy( repmat(vector(zeros(5,1),2), 5,2) );
% so the size(vecArr) will correctly return [5 2].
% 
% Implemented methods:
%  obj = vector(oneElem,n)        %constructor
%  // outObjs = copy(inObjs)      %copy method
%  [obj2, obj] = swap(obj, obj2)
%  cpt = capacity(obj)
%  obj = reserve(obj,n)
%  obj = resize(obj,n [,val])
%  obj = pop_back(obj)
%  obj = push_back(obj,elems)
%  elems = getElems(obj, whereIdc) or 
%  elems = getElems(obj, whereIdc1, whereIdc2, ...)
%  obj = setElems(obj, whereIdc, what) or
%  obj = setElems(obj, whereIdc1, whereIdc2, ..., what)
%  elem = back(obj)
%  elem = front(obj)
%  obj = clear(obj)
%  isEmpty = empty(obj)
%  obj = insert(obj, whereIdc, what)
%  obj = erase(obj,idc)
%  obj = shrink_to_fit(obj)
%
% getElems and setElems support ':' (in apostrophes) as index
% specification. But does not support end as index; e.g. 1:2:end is not
% supported. Use 1:2:vec.size_ instead.
% 
% The data_ is a matrix of elements.  Elements are concatenated along the
% first dimension of length 1 or a new dimension is added if necessary to
% concatenate along.
%
% Developed in the R2016b Matlab version.
% 29-th Dec 2020, Jozef Lukac, CTU in Prague, FEE.


% handle class together with copy-method enables deep-copy (needed e.g. 
% when vector of vectors is created) whereas matlab.mixin.Copyable and its
% copy method does that not.
% 
% But vector derived from matlab.mixin.Copyable is a bit faster. Just
% comment-out the copy-method to use it.
%
% It does NOT fully work in Octave -- e.g. the copy method does not
% work.
%
% matlab.mixin.Copyable allows the handle-object to be (deep-)copied (but
% only to the first level).
% using:  obj2Arr = copy(obj1Arr); 
% Note, obj2Arr = obj1Arr copies just the handle/handle-array.
% 
% Test functionality by running the prepared tests in test_vector.m file.
% Run by calling
%
% runtests
%
% in the command window. In Octave, only manual 'copy-paste' tests work. :(

    properties (SetAccess = private)
       data_ = []; %data storage
       Ndims_ = 0; %number of dimensions of data
       
       %index of the dimension to concatenate elemenets along
       idxCatDim_ = 0;
       capacity_ = 0; %capacity of vector-container
       size_ = 0; %size of vector-container
       oneElemExample_ = []; %example of one element in the vector
       oneElem_size_ = []; %output of size(oneElemExample_)
       repMat_dimsCell_ = {}; %auxiliary cell array for the use in repmat
       %auxiliary cell array for the assignment of the new data
       assign_dimsCell_ = {};
    end
       

    
    methods
    function obj = vector(oneElem,n)
    % constructor VECTOR
    % obj = vector([oneElem,n])
    % Constructor of vector.
    % 
    % Usage (examples):
    % oneElemExample = [1 2 0 -1.1];
    % n = 3;
    % vec1 = vector(oneElemExample,n); %vector of 3 elements
    % vec2 = vector(oneElemExample); %default n is 0
    % vec3 = vector(); %default oneElemExample is [0]
    %
    % Inputs:
    %  <oneElem> -- an example of element (the default element used for
    %  resize method).
    %  <n> -- number of elements to make in the vector.
    % Outputs:
    %  <obj> -- vector handle.
    %
        %check input
        if nargin < 2
           n = 0; 
        end
        if nargin < 1
            oneElem = 0; 
        end
        sz = size(oneElem);
        
        if any(sz == 0)
           error(['vector::vector. An example of the element have',...
               ' to have all dimensions nonzero.']); 
        end
        
        obj.idxCatDim_ = find(sz == 1, 1); %find first dimension of the length 1
        NelemDims = numel(sz); %number of dimensions of the vector-element
        if isempty(obj.idxCatDim_)
            obj.idxCatDim_ = NelemDims +1; %add dimension to concatenate along
            obj.Ndims_ = NelemDims + 1;
        else
            obj.Ndims_ = NelemDims; %use first 1-dimension to concatenate along
        end
        
        obj.repMat_dimsCell_= repmat({1},obj.Ndims_,1);
        obj.repMat_dimsCell_{ obj.idxCatDim_ } = n;
        obj.data_ = repmat_general(oneElem,obj.repMat_dimsCell_{:});
        obj.capacity_ = n;
        obj.size_ = n;
        obj.oneElemExample_ = repmat_general(oneElem,1,1);
        obj.oneElem_size_ = sz;
        
        obj.assign_dimsCell_ = obj.repMat_dimsCell_;
        obj.assign_dimsCell_(1:numel(obj.oneElem_size_)) = {':'};
    end
    
    function outObjs = copy(inObjs)
    % method COPY
    % outObjs = copy(inObjs)
    % Makes the deep-copy of vector-array given as input.
    % Inputs:
    %  <inObjs> -- array/matrix of (handles of) vector objects.
    % Outputs:
    %  <outObjs> -- array/matrix of (handles of) copied vector objects.
    %
        sz = size(inObjs);
        szCell = num2cell(sz);
        NumEl = numel(inObjs);
        outObjs = repmat(vector(), szCell{:});
        for iObj = 1:NumEl
            outObjs(iObj) = vector(inObjs(iObj).oneElemExample_);
            restrCell = inObjs(iObj).assign_dimsCell_;
            restrCell{inObjs(iObj).idxCatDim_} = 1:inObjs(iObj).size_;
            outObjs(iObj).data_ = repmat_general(inObjs(iObj).data_(restrCell{:}),1,1);
            outObjs(iObj).size_ = inObjs(iObj).size_;
            outObjs(iObj).capacity_ = inObjs(iObj).size_;
        end
    end
    
    function [obj2, obj] = swap(obj, obj2)
    % method SWAP
    % [obj2, obj] = swap(obj, obj2)
    % Swaps the vectors (their handles).
    %
    % usage: [vec1,vec2] = vec1.swap(vec2), or
    % [vec1,vec2] = swap(vec1,vec2)
    %
    % Inputs:
    %  <vec2> -- the second vector (handle).
    % Outputs:
    %  <obj2>, <obj> -- handles of vectors. The correct assignment of
    %  output is necessary.
        
        %no need to copy whole data, just swap the handles/pointers
    end
    
%     function sz = size(obj)
%         sz = obj.size_;
%     end
    
    function cpt = capacity(obj)
    % method CAPACITY
    % Returns current vector capacity.
    % Inputs:
    % 
    % Outputs:
    %  <cpt> -- vector capacity.
    % 
        cpt = obj.capacity_;
    end
    
    function obj = reserve(obj,n)
    % method RESERVE
    % obj = reserve(obj,n)
    % Reserves space for <n> objects -- capacity is <n>, size_ remains.
    % Inputs:
    %  <n> -- number of elements to make space for. (so after the vector-
    %  enlargement the new memory-reallocation is not necessary).
    % Outputs:
    % 
       obj.reserve_priv(n);
    end
    
    function obj = resize(obj,n,val)
        if n < 0 || round(n) ~= n
           error('vector::resize Invalid n.'); 
        end
        if n <= obj.size_
            obj.size_ = n;
        else
            if nargin > 2
               sz = size(val);
               if any(sz ~= obj.oneElem_size_)
                  error('vector::resize Invalid dimensions of val.');
               end
            else
               val = obj.oneElemExample_;
            end
            Nadded = obj.reserve_priv(n,val);
            %if the non-default init-value was used, update new elements
            if nargin > 2
                rep_Cell = obj.repMat_dimsCell_;
                rep_Cell{ obj.idxCatDim_ } = n - obj.size_ - Nadded;
                asgn_Cell = obj.assign_dimsCell_;
                asgn_Cell{obj.idxCatDim_} = obj.size_+1 : n-Nadded;
                obj.data_( asgn_Cell{:} ) = repmat_general(val, rep_Cell{:});
            end
            obj.size_ = n;
        end
    end
    
    
    function obj = pop_back(obj)
    % method POP_BACK
    % obj = pop_back(obj)
    % Retudes the size_ of the vector (from its end).
    % Inpputs:
    % 
    % Outputs:
    % 
       if obj.size_ > 0
           obj.size_ = obj.size_ - 1;
       else
          error('vector::pop_back Already no element to pop.'); 
       end
    end    
    
    function obj = push_back(obj,elems)
    % method PUSH_BACK
    % obj = push_back(obj,elems)
    % Pushes back the element(s) given as input.
    % Inputs:
    %  <elems> -- an element (or correctly concatenated elements) to be
    %  pushed back to the vector.
    % Outputs:
    % 
        [Nadd, isOK] = checkElems_priv(obj,elems);
        if ~isOK
           error('vector::push_back Incorrect elems'); 
        end
        
        if obj.capacity_ < obj.size_ + Nadd
           %allocate additional space
           Nalloc = obj.size_ + max(obj.size_, Nadd);
           obj.reserve(Nalloc);
        end
        %assign new data
        auxCellRange = obj.assign_dimsCell_;
        auxCellRange{ obj.idxCatDim_ } = obj.size_+1:obj.size_+Nadd;
        obj.size_ = obj.size_ + Nadd;
        obj.data_(auxCellRange{:}) = repmat_general(elems,1,1);
    end
    
    
    function elems = getElems(obj, varargin)
    % method GETELEMS
    % elems = getElems(obj, whereIdc) or
    % elems = getElems(obj, whereIdc1, whereIdc2, ...)
    % Returns the elements specified by their indices.
    % Inputs:
    %  <whereIdc> -- array of indices of elements to be returned
    %  (concatenated). You can also use ':' (in apostrophes) to specify
    %  all the data.
    % Or you can index directly in the data_ matrix by specifying all the
    % necessary indices.
    % It does not support 'end'-index, e.g. 1:2:end is not allowed as an
    % input. Use 1:2:vec.size_ instead.
    %  <whereIdc1>, <whereIdc2>, ... -- arrays of indices of (parts of) 
    %  elements to be returned.
    % Outputs:
    %  <elems> -- (parts of) elements indexed by indices.
    %
       auxCell = obj.assign_dimsCell_;
       if numel(varargin) == 1
           auxCell{ obj.idxCatDim_ } = varargin{1};
       else
           auxCell = varargin;
       end
       %restrict data to the valid range
       restrCell = obj.assign_dimsCell_;
       restrCell{obj.idxCatDim_} = 1:obj.size_;
       S = struct('type','()','subs',{auxCell});
       %index in the restricted data
       elems = subsref(obj.data_(restrCell{:}), S);
    end
    
    % (indices, vals) or (idc1, idc2, ..., vals)
    function obj = setElems(obj,varargin)
    % method SETELEMS
    % obj = setElems(obj, whereIdc, what) or
    % obj = setElems(obj, whereIdc1, whereIdc2, ..., what)
    % Assigns the new values to the data.
    % Inputs:
    %  <whereIdc> -- indices of elements to update
    %  <what> -- correctly concatenated elements to be set at the 
    %  <whereIdc>-positions.
    % Or you can directly index in the underlying data_ by specifying
    % indices in the whole data_ matrix.
    % It does not support 'end'-index, e.g. 1:2:end is not allowed as an
    % input. Use 1:2:vec.size_ instead.
    %   <whereIdc1>, <whereIdc2>, ... -- indices in data_ matrix
    %   <what> -- elements to be used as the update.
    % Outputs:
    % 
       auxCell = obj.assign_dimsCell_;
       if numel(varargin) == 2
           auxCell{ obj.idxCatDim_ } = varargin{1};
       else
           auxCell = varargin(1:end-1);
       end
       vals = repmat_general(varargin{end},1,1);
       %restrict data to the valid range
       restrCell = obj.assign_dimsCell_;
       restrCell{obj.idxCatDim_} = 1:obj.size_;
       S = struct('type','()','subs',{auxCell});
       %update in the restricted data
       obj.data_(restrCell{:}) = subsasgn(obj.data_(restrCell{:}),S,vals);
    end
    
    function elem = back(obj)
    % method BACK
    % elem = back(obj)
    % Returns the last element in the vector.
    % Intputs:
    % 
    % Outputs:
    %  <elem> -- the last element in the vector.
    %
       elem = obj.getElems(obj.size_); 
    end
    
    function elem = front(obj)
    % method FRONT
    % elem = front(obj)
    % Returns the first element in the vector.
    % Intputs:
    % 
    % Outputs:
    %  <elem> -- the first element in the vector.
    %
        elem = obj.getElems(1);
    end
    
    function obj = clear(obj)
    % method CLEAR
    % obj = clear(obj)
    % Clears the vector. The capacity remains.
    % Intpus:
    % 
    % Outputs:
    %
       obj.size_ = 0; 
    end
    
    function isEmpty = empty(obj)
    % method EMPTY
    % isEmpty = empty(obj)
    % Checks whether the vector-container is empty.
    % Inputs:
    % 
    % Outputs:
    %  <isEmpty> -- flag whether the vector is empty.
    %
       isEmpty = obj.size_ == 0; 
    end
    
    function obj = insert(obj, whereIdc, what)
    % method INSERT
    % obj = insert(obj, whereIdc, what)
    % Inserts elmenets in <what> to the positions given in <whereIdc>.
    % Inputs:
    %  <whereIdc> -- array of indices to insert elements to. It does not
    %  need to be sorted. To insert more elements before the given element
    %  just duplicate the index, e.g. whereIdc = [2 2].
    %  <what> -- correctly concatenated matrix of elements to be inserted
    %  in front of elements given by their indices in <whereIdc>.
    % Outputs:
    %
       [what_numel, isOK] = checkElems_priv(obj,what);
       assert(isOK && what_numel > 0, ...
           'vector::insert Invalid data type or dimension.');
       whereIdc_numel = numel(whereIdc);
       assert(whereIdc_numel > 0, ...
           'vector::insert Zero number of indices.');
       NInsert = max(whereIdc_numel,what_numel);
       caseLogical = [whereIdc_numel == 1, what_numel == 1, ...
           whereIdc_numel == what_numel];
       assert(any(caseLogical), ...
           'vector::insert Incompatible sizes of whereIdc and what.');
       assert(isnumeric(whereIdc),...
           'vector::insert whereIdc is not numeric.');
       if caseLogical(1)
           whereIdc = repmat(whereIdc, NInsert,1);
       elseif caseLogical(2)
           rep_cell = obj.repMat_dimsCell_;
           rep_cell{obj.idxCatDim_} = NInsert;
           what = repmat_general(what, rep_cell{:});
       end
       NNewSize = obj.size_ + NInsert;
       NNewCap = obj.size_ + max(obj.size_,NInsert);
       
       [idc_oldData, idc_newData] = insert_getIdc(whereIdc,...
           obj.size_);
       
       %everything should be ok, so do the insertion
       obj.reserve(NNewCap);
       %make spaces in old data for the new data
       idc_oldData_Cell = obj.assign_dimsCell_;
       idc_oldData_noSpaces_Cell = obj.assign_dimsCell_;
       idc_newData_Cell = obj.assign_dimsCell_;
       
       idc_oldData_Cell{obj.idxCatDim_} = idc_oldData;
       idc_oldData_noSpaces_Cell{obj.idxCatDim_} = 1:obj.size_;
       idc_newData_Cell{obj.idxCatDim_} = idc_newData;
       %make spaces
       obj.data_(idc_oldData_Cell{:}) = ...
           obj.data_(idc_oldData_noSpaces_Cell{:});
       %insert new data to the prepared spaces
       obj.data_(idc_newData_Cell{:}) = repmat_general(what,1,1);
       
       %set the new size
       obj.size_ = NNewSize;
    end
    
    function obj = erase(obj,idc)
    % method ERASE
    % obj = erase(obj,idc)
    % Erases elements with indices in <idc>
    % Inputs: 
    %  <idc> -- array of element-indices to erase. Does not need to be
    %  sorted.
    % Outputs:
    %
        currVecIdc = (1:obj.size_).';
        currVecIdc(idc) = [];
        NnewVecSize = numel(currVecIdc);
        if NnewVecSize < obj.size_
           asgnCell = obj.assign_dimsCell_;
           asgnCell_subseqData = obj.assign_dimsCell_;
           
           asgnCell{obj.idxCatDim_} = currVecIdc;
           asgnCell_subseqData{obj.idxCatDim_} = (1:NnewVecSize).';
           
           obj.data_(asgnCell_subseqData{:}) = obj.data_(asgnCell{:});
           obj.size_ = NnewVecSize;
        end
    end
    
    function obj = shrink_to_fit(obj)
    % method SHRINK_TO_FIT
    % When necessary it reassignes data_ to that size_ == capacity. It does 
    % not change the data otherwise.
    % Inputs:
    %
    % Outputs:
    %
       if obj.size_ ~= obj.capacity_
           asgnCell = obj.assign_dimsCell_;           
           asgnCell{obj.idxCatDim_} = (1:obj.size_).';
           obj.data_ = obj.data_(asgnCell{:});
           obj.capacity_ = obj.size_;
       end
    end
    
    end
    
    
     
    
    methods (Access = private)
    
    function [Nadded, obj] = reserve_priv(obj,n,val)
    % method RESERVE_PRIV
    % [Nadded, obj] = reserve_priv(obj,n[,val])
    % Reserves space in the vector.
    % reserve(n) -- updates possible new elements by the default value
    % reserve(n,val) -- updates possible new elements by the val value
    % Inputs:
    %  <n> -- number of elements to make space for.
    %  <val> -- value -- an element to be used as a dummy-element
    %  for the space-reservation.
    % Outputs:
    %  <Nadded> -- number of elements the capacity was expanded by
    %
       Nadded = 0;
       if obj.capacity_ < n
          obj.repMat_dimsCell_{obj.idxCatDim_} = n;
          Nadded = n - obj.capacity_;
          if nargin < 3
             val = obj.oneElemExample_;
          end              
          newData = repmat_general(val,obj.repMat_dimsCell_{:});
          %if there are any data, copy them
          if obj.size_ > 0
              auxCellRange = obj.assign_dimsCell_;
              auxCellRange{ obj.idxCatDim_ } = 1:obj.size_;
              newData(auxCellRange{:}) = repmat_general(obj.data_(auxCellRange{:}),1,1);
          end
          obj.data_ = newData;
          obj.capacity_ = obj.repMat_dimsCell_{obj.idxCatDim_};
       end
    end
    
    function [Nelems, isOK] = checkElems_priv(obj,elems)
    % method CHECKeLEMS_PRIV
    % [Nelems, isOK] = checkElems_priv(obj,elems)
    % Checks whether the given element(s) are valid and counts their
    % number.
    % Inputs:
    %  <elems> -- vector-element or correctly concatenated matrix of 
    %  vector-elements
    % Outputs:
    %  <Nelems> -- number of elements in <elems>
    %  <isOK> -- flag whether the <elems> can be processed by the vector.
    %
        isClassOK = strcmp(class(obj.oneElemExample_), class(elems));
        sz = size(elems);
        diffSz = sz - obj.oneElem_size_;
        idcNeq0 = find(diffSz ~= 0);
        lenIdcNeq0 = length(idcNeq0);
        
        
        Nelems = 0;
        isOK = false;
        if lenIdcNeq0 <= 1 && isClassOK
            if lenIdcNeq0 == 0
                Nelems = 1;
                isOK = true;
            elseif idcNeq0 == obj.idxCatDim_
                Nelems = diffSz(idcNeq0) +1;
                isOK = true;
            end
        end
    end
    
    
    end
        
end

%--------- AUXILIARY FUNCTIONS ---------
function [idc_oldData, idc_newData] = insert_getIdc(idc, sz)
% function INSERT_GETIDC_PRIV
% [idc_oldData, idc_newData] = insert_getIdc(idc, sz)
% returns auxiliary indices for the insert vector-method
% Inputs:
%  <idc> -- indices for the insertion
%  <sz> -- current vector size
%
% Outputs:
%  <idc_oldData> -- assignment-indices of the current data
%  <idc_newData> -- assignment-indices of the old data
%
    Nidc = numel(idc);
    [srtdIdc, srtd_I] = sort(idc);
    unqIdc = unique(srtdIdc(:));
    NunqIdc = numel(unqIdc);
    freqUnqIdc = sum(unqIdc * ones(1,Nidc) == repmat(idc(:).',NunqIdc,1),2);

    ones_arr = ones(sz+1,1);
    ones_arr(unqIdc) = 1 + freqUnqIdc;

    idc_oldData = cumsum(ones_arr(1:end-1));
    idc_newData = (1:sz +Nidc).';
    idc_newData(idc_oldData) =[];
    idc_newData(srtd_I) = idc_newData;
end

function repElem = repmat_general(elem, varargin)
% function REPMAT_GENERAL
% repElem = repmat_general(elem, sz1 [, sz2, ...])
% It works the same as repmat but deep-copies/replicates also handle
% objects. It is also used as the copy for both the handle and the 
% non-handle objects (the construction repmat_genera(A,1,1) is used).
% Inputs:
%  <elem> -- element/matrix to be repicated
%  <sz1>, <sz2>, ... --- numbers to replicate <elem> (for each dimension)
%
   repElem = repmat(elem,varargin{:});
   if isa(elem,'handle')
        repElem = copy(repElem);
   end
end
