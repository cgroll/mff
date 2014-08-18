function areEqual = allVectorEntriesEqual(vector)
%allVectorEntriesEqual tests whether all entries of numeric vector
% are equal
%
% Args:
%	vector           (1 x n) or (n x 1) numeric vector
%
% Output:
% 	areEqual         boolean, indicating whether elements
%                   are equal or not
%

% throw error for empty vector
if (isempty(vector))
    error('Vector:Empty', 'Testing empty vector for equal entries.');
end

% get number of entries
n = numel(vector);

% how many elements are equal to first entry
nEqual = sum(vector == vector(1));

% is number of equal elements equal to vector length 
areEqual = (n == nEqual);
