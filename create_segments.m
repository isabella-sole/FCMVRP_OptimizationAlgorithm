function [idxs, dist, tr_time, battery_tr, Max_flight_time] =  create_segments(xloc, yloc, vertex_idxs, cs_idxs, v)

% For example
% first of all, inside the for loop, a 3D matrix is build. At each layer a
% 25x2 matrix is build where there are all the combination of number from 1
% to 25 taken two by two (two column). There is a special attention
% eliminating the row where the same i-th number is repeated, for example if
% we are dealing with i = 10 the matrix stored for that iteration will be:
%     10     1
%     10     2
%     10     3
%     10     4
%     10     5
%     10     6
%     10     7
%     10     8
%     10     9
%     10    11
%     10    12
%     10    13
%     10    14
%     10    15
%     10    16
%     10    17
%     10    18
%     10    19
%     10    20
%     10    21
%     10    22
%     10    23
%     10    24
%     10    25
counter = 1; 
for i = 1:length(vertex_idxs) % for each vertex (target+cs)
    
    % C = nchoosek(v,k) returns a matrix containing all possible combinations of the elements of vector v taken k at a time. Matrix C has k columns and m!/((m–k)! k!) rows, where m is length(v)
    % idxs is a matrix containing all possible combinations of the elements of vector [i,vertex(1:end)] (from 1 to 25) taken 2 at a time
    % idxs contains all possible combinations from 1 to 25 taken 2 per time (in fact the last row is 24 25 since all the other combinations fo the number 25 are already written for the other numbers)
    % but the first rows are the combinations of the i-th number (wihch is also repeated when it is called from vertex(1:end))
    idxs = nchoosek([i,vertex_idxs(1:end)], 2); % WARNING: after some rows is is overwritten
    a = idxs(1:length(vertex_idxs), :); % extracts the submatrix 25x2: extracts the first 25 rows and both the columns of the i-th number in the iteration 
    b = a(:, 1)-a(:, 2); % 1x25, performs the subtraction for each row of the first and the second column values
    b = find(b == 0); % check where is the element == 0 and returns the index of b where there is the 0
    a1 = a(1 : b - 1, :);  % extract the submatrix of the i-th term up to the row before where there is the same number, in other words where the subtraction of the two columns is zero (b-1)
    a2 = a(b +1 : end, :); % % extract the submatrix of the i-th term from the row after there is the same number, in other words where the subtraction of the two columns is zero (b+1)
    a = [a1; a2]; % the new matrix a is composed by all the combination for the i-th term, i-th term in the first column and from 1 to 25 in the second one, without the row i-th i-th that was elimataed (the b line where the subtraction is 0)
    indexes(: ,:, i) = a;  % put this matrix in a 3D matrix that contains this type of matrix for each i-th element 
end

% Create the idxs matrix that is nothing more then all the indexes matrices
% put together into a unique 600x2 matrix (instead of havig a 3D matrix)
% idxs CONTAINS ALL THE POSSIBLE LINKS BETWEEN TWO VERTEXES (INCLUDING TARGETS, CS AND START)
idxs = []; % initialize a new matrix
idxs = zeros(size(indexes, 1)*size(indexes, 3), 2); % fill the matrix with 240x2 zeros
counter = 1;
% indexes is a possible_links x 2 matrix
for i = 1 : size(indexes, 3) % depth
    for j = 1: size(indexes, 1) % row
        for k = 1: size(indexes, 2) % column
         idxs(counter, k) = indexes(j, k, i); % put the current value of indexes in idxs
        end
        counter = counter + 1;
    end
end

% Evaluates the distance of each path 
% dist is a vector that contains the distance of each link contained in the
% idxs matrix (distanza associated to each segment composed by the element 
% of the first column and the element of the second column of each row) )
dist = hypot(xloc(idxs(:, 1)) - xloc(idxs(:, 2)), ...
             yloc(idxs(:, 1)) - yloc(idxs(:, 2))); % hypot computation of the distance
lendist1 = length(dist); 
dist = dist'.*1000; % [Km] express the distances in kilometers
tr_time = dist./v; % [sec]  vector that contains the travelling time for each link 

% Add inspection time in target node
fprintf("\nPlease add the required inspection time (in seconds):")
my_input = input('\nEnter your input:'); 
insp_time = my_input;
for t = 1:length(tr_time)
    tr_time(t) = tr_time(t) + insp_time;
end    

%Definisco peso arco in termini energetici
%U = 40 * 60; %40 min autonomia --> 2400 sec tra [20% e 80 %]
%battery_tr = tr_time./U ; %definisce percentiale batteria consumata nell 'arco (la batteria consumata è indipendente dal peso aggiuntivo in termini di tempo aggiunto negli archi da e per le cs

% Retrieve the battery required to travel each link and the maximum flight time, 
% passing the time required to travel across each link and the velocity that is used 
[battery_tr, Max_flight_time] = Energy_consumption_modelization(tr_time, v);

% Aggiungo travelling time per raggiungere le CS, è il costo aggiuntivo se
% visitate (diviso a meta tra arco per CS e arco da CS in modo da non avere
% archi vuoti da CS)
% if a cs is visited add the cost required to visit it (time)
for i = 1: length(cs_idxs)  % for all the cs
    a = find(idxs(:,2) == cs_idxs(i)); % check if the specific cs i-th is present in the second column of idxs (end of a segment) and memorize its index in a
    b = find(idxs(:,1) == cs_idxs(i));  % check if the specific cs i-th is present in the first column of idxs (start of a segment)and memorize its index in b
    tr_time(a) = tr_time(a) + 900;  % for the visit at that cs add time spent to visit it
    tr_time(b) = tr_time(b) + 900;  % for the visit at that cs add time spent to visit it
end


% travelling time vector: vector with the timining required to travel each link
% the first half of the vactor has all zeros since they are associated to the z_ij variable
tr_time = [zeros(length(battery_tr),1);tr_time];