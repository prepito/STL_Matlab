function tests = test_vector
tests = functiontests(localfunctions);
end

function test_creation(testCase)
    e1 = zeros(5,1);
    e2 = cell(4,1);
    
    vec1 = vector(e1);
    vec2 = vector(5,3);
    vec3 = vector(ones(3,2),2);
    vec4 = vector(e2,3);
    vec5 = vector( vector([5 3], 3),7 );
    
    
    assert(isempty(vec1.getElems(':')), ...
        'Incorrect vector creation. vec1');
    assert(all(vec2.getElems(':') == [5 5 5].'), ...
        'Incorrect vector creation. vec2');
    auxLogical = all(all(all( vec3.getElems(':') ==...
        cat(3,ones(3,2),ones(3,2)))));
    assert(auxLogical,'Incorrect vector creation. vec3');
    assert(iscell(vec4.getElems(':')) && ...
        all(size(vec4.getElems(':')) == [4 3]),...
        'Incorrect vector creation. vec4');
    assert(...
        all(size(vec5.getElems(':')) == [7 1]),...
        'Incorrect vector creation. vec5');
    
    assert(vec1.size_() == 0, 'Incorrect vector size. vec1');
    assert(vec2.size_() == 3, 'Incorrect vector size. vec2');
    assert(vec3.size_() == 2, 'Incorrect vector size. vec3');
    assert(vec4.size_() == 3, 'Incorrect vector size. vec4');
    assert(vec5.size_() == 7, 'Incorrect vector size. vec5');
    
end

function test_vectorArray(testCase)
    vecArr = copy(repmat( vector(zeros(3,1),2), 5,1) );
    vecArr(1).push_back(2*ones(3,1));
    vecArr(2).push_back(repmat([4 5 6 10], 3,1));
    vecArr(3).resize(20, [ 5 6 7].');
    assert(vecArr(1).size_ == 3);
    assert(vecArr(2).size_ == 6);
    assert(vecArr(3).size_ == 20);
    assert(all(size(vecArr) == [5,1]));
end

function test_copy(testCase)
%nesting of vector is not fully supported -- 
% -- copy makes the deep-copy only at the first level,
% other levels are shallow-copied
    vec1 = vector( vector(zeros(3,1)),3 );
    vec1_Elem2 = vec1.getElems(2);
    vec1_Elem2.push_back([1 2 3].');
    vec2 = copy(vec1);
    vec2_Elem2 = vec2.getElems(2);
    vec2_Elem2.push_back([2 2 2; 3 3 3].');
    
    vec_toPush = vector(ones(3,1),5);
    vec1.push_back(vec_toPush);
    vec1_back = vec1.back();
    vec1_back.clear();
    
    assert(vec1.size_() == 4,'Incorrect vector size. vec1');
    assert(vec2.size_() == 3, 'Incorrect vector size. vec2');
    vec1_Elem2_aux = vec1.getElems(2);
    vec2_Elem2_aux = vec2.getElems(2);
    % copy in matlab.mixin.Copyable does not deep-copy the 'second-level' 
    % vectors, it just copies their handles
    % -- vec1_Elem2_aux.size_() == 3
    %
    % whereas the copy maunaly implemented (for the handle base-class)
    % deep-copies also the 'second-level' vectors
    %  -- vec1_Elem2_aux.size_() == 1
    assert(vec1_Elem2_aux.size_() == 1,'Incorrect vector size.');
    assert(vec2_Elem2_aux.size_() == 3,'Incorrect vector size.');
    
    
    auxVect = vec1.back(); %vec1.getElems(6);  the clear-ed vector
    assert(auxVect.size_() == 0, 'Incorrect vector size.');

    vec1 = vector(zeros(4,1),3);
    vec1.push_back(ones(4,1));
    vec2 = copy(vec1);
    vec2.erase([2 1]);
    assert(vec1.size_() == 4 && vec2.size_ == 2, 'incorrect vector sizes');
end

function test_swap(testCase)
    vec1 = vector(zeros(1,3,4),10);
    vec2 = vector(zeros(3,4,5,6),8);
    [vec1, vec2] = swap(vec1,vec2);
    assert(vec1.size_ == 8 && vec2.size_ == 10 && ...
        all(size(vec1.getElems(':')) == [3 4 5 6 8]) && ...
        all(size(vec2.getElems(':')) == [10 3 4]),...
        'Swap error.');
    
end

function test_reserve_capacity(testCase)
    vec1 = vector(zeros(5,1,7,19,7),10);
    assert(vec1.size_ == 10 && vec1.capacity() == 10);
    vec1.reserve(25);
    assert(vec1.size_ == 10 && vec1.capacity() == 25);
end

function test_resize(testCase)
    vec = vector(zeros(5,1));
    vec.resize(6);
    assert(vec.capacity() == 6);
    vec.resize(4, 2*ones(5,1));
    vec.resize(5, 3*ones(5,1));
    assert(vec.capacity() == 6);
    assert(all(all(vec.data_ == repmat([0 0 0 0 3 0], 5,1))));
    
    vec.resize(8, 4*ones(5,1));
    assert(vec.capacity() == 8 && vec.size_ == 8);
    assert(all(all(vec.data_ == repmat([0 0 0 0 3 4 4 4], 5,1))));
end

function test_popPush_back(testCase)
    vec = vector(zeros(4,1),3);
    vec.push_back(repmat([1 2 3 4 5], 4,1)); %multiple push_back
    assert(vec.size_ == 8);
    vec.pop_back();
    vec.pop_back();
    assert(vec.size_ == 6);
    vec.push_back(3 * ones(4,1));
    assert(vec.size_ == 7);
    assert(all(all(vec.getElems(':') == ...
        repmat([0 0 0 1 2 3 3],4,1))));
end

function test_getSetElems(testCase)
    vec = vector(zeros(1,4),0);
    vec.resize(5,3* ones(1,4));
    vec.setElems([2 1], repmat( [8;9], 1,4));
    assert(all(all(vec.data_ == ...
        repmat([9 8 3 3 3].',1,4))));
    vec.push_back(10*ones(1,4));
    
    assert(all(all(vec.getElems(':') == ...
        repmat([9 8 3 3 3 10].',1,4))));
    
    assert(all(all(vec.getElems(1:2,[4 3]) == ...
        repmat([9 8].', 1,2))));
    
    vec.setElems(':',[4 3],[0 0 0 0 0 0; 1 1 1 1 1 1].');
    assert(all(all(vec.getElems(2:3,':') ==...
        [8 8 1 0;3 3 1 0])));
end

function test_push_back(testCase)
    vec = vector(zeros(4,1),3);
    for i=1:4
       vec.push_back(i*ones(4,1)); 
    end
    vec.push_back(ones(4,3));
    
    assert(vec.size_ == 10,'error vector size');
    expectedData = repmat([0 0 0 1 2 3 4 1 1 1],4,1);
    assert(all(all(vec.getElems(':') == expectedData)),'Error data content');
end

function test_front_clear_empty_insert(testCase)
    vec = vector(0);
    assert(vec.empty());
    vec.resize(5,3);
    vec.insert([6,2],[10; 12]);
    assert(all(vec.getElems(':') == ...
        [3 12 3 3 3 3 10].'));
    assert(vec.front() == 3);
    vec.erase([4 3 7]);
    assert(all(vec.getElems(':') == ...
        [3 12 3 3].'));
    vec.clear();
    assert(vec.empty());
    
    vec2 = vector(zeros(3,1),3);
    vec2.insert([2 2 2 2], ones(3,1)); %should replicate <what> input
    assert(all(all(vec2.getElems(':') == repmat([0 1 1 1 1 0 0],3,1))));
    vec2.insert(5, 2 * ones(3,2)); %whould replicate <whereIdc> input
    assert(all(all(vec2.getElems(':') == repmat([0 1 1 1 2 2 1 0 0],3,1))));
    
     %insert at the end (push_back) and in front of the first element
    vec2.insert([10,1],repmat([-5 7],3,1));
    assert(all(all(vec2.getElems(':') == repmat([7 0 1 1 1 2 2 1 0 0 -5],3,1))));
end

function test_insert_erase_shrink_to_fit(testCase)
    vec = vector([10 11],5);
    vec.push_back([88 88]);
    assert(vec.size_ == 6 && vec.capacity() == 10);
    vec.shrink_to_fit();
    assert(vec.size_ == 6 && vec.capacity() == 6);
end
