% This program aims to find the optimal routes for multiple drones, satisfying 
% all the constraints of the case, that will be later illustrated 

close all
clear 
clc

global CS Target % global variables for charging stations, targets
              
% Initialize empty matrices for the optimization constraints
Aeq= [];
beq= [];
Aineq = [];
bineq = [];
x_temp = [];

% Constant variables
velocity = 8.3; % velocity m/s --> 30 Km/h
      
% Initialise number of cs and targets and their positions

% NOTICE that the coordinates of a cs corresponds also to the
% initial and ending position of the vehicles, since each cs has a
% vehicle that starts and ends its route there.

% Insert CS and Target x-y coords
CS = 4; % number of charging station
Target = 14; % number of target
target_x = [1, 3, 1, 5, 7, 9, 6, 7, 12, 1, 2, 10, 14, 9]; % targets x coord
target_y = [1, 2, 7, 1, 8, 3, 4, 3, 14, 10, 15, 12, 4, 11]; % targets y coord
cs_x = [5, 9, 2, 2]; % cs x coord
cs_y = [2, 9, 13, 5]; % cs y coord




%% Create map
% This function:
% -displays a map of charging stations and targets
% -creates vector contaning the vetexes, cs and target indexes
% given the number of cs, target and their x and y coordinates

% Dimensional check
if length(target_x) ~= Target
    fprintf("\nThe number of target does not match the inserted x coordinates")
end

    
if length(target_y) ~= Target
    fprintf("\nThe number of target does not match the inserted y coordinates")
end

if length(cs_x) ~= CS
    fprintf("\nThe number of cs does not match the inserted x coordinates")
end
if length(cs_y) ~= CS
    fprintf("\nThe number of cs does not match the inserted y coordinates")
end
    
xloc = [target_x, cs_x]; % builds the vector of x targets+cs locations
yloc = [target_y, cs_y]; % builds the vector of y targets+cs locations

% Build target labels
target_label = [];
for i = 1: Target
    target_label = [target_label, i]; 
end 
target_label = string(target_label);
  
% Build cs labels
cs_lalbel = [];
for j = Target+1:length(xloc)
    cs_lalbel = [cs_lalbel,j];
end  
cs_lalbel = string(cs_lalbel);

[vertex_idxs, target_idxs, cs_idxs , h_figure] = create_map(CS, Target, xloc, yloc, target_label, cs_lalbel);
CS
Target

vertex_idxs
target_idxs
cs_idxs

%% Create segments
% This function creates:
% -a matrix idxs that contains all possible links between vertexes
% -vectors containing the distance, battery, flight time of each of these segments
% given the target and cs locations together with the indexes array of vertexes and charging stations and the velocity
[idxs, dist, tr_time, battery_tr, Max_flight_time] =  create_segments(xloc, yloc, vertex_idxs, cs_idxs, velocity);
fprintf('The Maximum Flight time is: %d\n',Max_flight_time);
lendist = length(dist);  % size of the vector containing all the distances (which is also the number of links)

Vertexes = Target + CS; % total number of stops, vertexes
idxs_total = idxs; % decision variable x(i,j) --> possibile aggiungere anche segmenti tra CS e CS
A_columns = length(battery_tr) + lendist; % number of A columns --> length of decision vector
half_of_A = zeros(length(battery_tr),1); % % matrix of convenient size to perform operations

     
%% Constraint 1
% EQUALITY Aeq, beq, x_ij variable
% The sum of arcs that enter in the depot/charging station has to be
% equal to the sum of the arcs that goes out from the depot/cs: in
% other words the number of drones that arrive in a cs has to be the
% same as the number of drones that leave the cs.
% This constraint assures that each drone will return to its depot. DA CAMBIARE

outgoing_arc = []; % vector that takes into account the outgoing arcs
incoming_arc = []; % vector that takes into account the incoming arcs
xtemp = []; % temporary x
btemp = []; % temporary b

for j = 1 : length(cs_idxs)  % for all cs
    for i = 1 : length(vertex_idxs) % for all vertexes
        if (vertex_idxs(i) == cs_idxs(j)) % if the couple is composed by equal numbers pass to the next iteration and don't do the following commands
            continue;
        end
        u = [cs_idxs(j), vertex_idxs(i), vertex_idxs(i), cs_idxs(j)]; % create a vector containing the outgoing arc and the incoming arc
        outgoing_arc(:,i) = sum(idxs(:,1) == cs_idxs(j),2) & ...
            sum(idxs(:,2) == vertex_idxs(i),2); % create a matrix with 1 in correspondence of the outgoing arc
        incoming_arc(:,i) = sum(idxs(:,1) == vertex_idxs(i),2) & ...
            sum(idxs(:,2) == cs_idxs(j),2); % create a matrix with 1 in correspondence of the incoming arc
    end
    
    
    outgoing_arc = sum(outgoing_arc,2); % create a unique vector with 1s in correspondence of all the outgoing arcs
    incoming_arc = sum(incoming_arc,2); % create a unique vector with 1s in correspondence of all the incoming arcs
    incoming_arc = -1*incoming_arc; % change the incoming arcs value to -1
    x_temp = outgoing_arc + incoming_arc; % create a unique vector of the outc and inc arcs
    xtemp = [half_of_A; x_temp]; % total x temporary, with also the first half of A (still empty since I am not considering yet the z_ij variable)
    Aeq(j,:) = xtemp; % fill at each iteration the A matrix with the xtemp vector
    btemp(j,1) = 0; % fill the btemporary vector
    
end

beq = btemp;


 %% PROVA: New Constraint 1
% % INEQUALITY Aineq, bineq, x_ij variable 
% outgoing_arc = []; % vector that takes into account the outgoing arcs
% incoming_arc = []; % vector that takes into account the incoming arcs
% xtemp_out = []; % temporary x
% xtemp_in = []; % temporary x
% xtemp = []; % temporary x
% btemp = []; % temporary b
% Aineq_in = [];
% Aineq_out = [];
% s_out = size(Aineq_out,1);  % s keeps track of A dimension
% counter_out = 1;
% s_in = size(Aineq_in,1);  % s keeps track of A dimension
% counter_in = 1;
% 
% % incoming arcs
% inc_arc_matrix = zeros(length(idxs),1);
% for j = 1 : length(target_idxs)  % for all cs
%     for i = 1 : length(cs_idxs) % for all vertexes
%         if (cs_idxs(i) == target_idxs(j)) % if the couple is composed by equal numbers pass to the next iteration and don't do the following commands
%             continue;
%         end
%         u = [target_idxs(j), cs_idxs(i), cs_idxs(i), target_idxs(j)]; % create a vector containing the outgoing arc and the incoming arc
%         incoming_arc(:,i) = sum(idxs(:,1) == target_idxs(j),2) & ...
%             sum(idxs(:,2) == cs_idxs(i),2); % create a matrix with 1 in correspondence of the outgoing arc
%         inc_arc_matrix = inc_arc_matrix + incoming_arc(:,i);
%     end
%     
%     x_temp_in = -1 * inc_arc_matrix; % create a unique vector of the outc and inc arcs
%     xtemp_in = [half_of_A; x_temp_in]; % total x temporary, with also the first half of A (still empty since I am not considering yet the z_ij variable) 
%     Aineq_in(s_in + counter_in,:) = xtemp_in; % fill at each iteration the A matrix with the xtemp vector
%     counter_in = counter_in + 1; 
%     
% end
% 
% 
% % outgoing arcs
% for j = 1 : length(cs_idxs)  % for all cs
%     for i = 1 : length(vertex_idxs) % for all vertexes
%         if (vertex_idxs(i) == cs_idxs(j)) % if the couple is composed by equal numbers pass to the next iteration and don't do the following commands
%             continue;
%         end
%         u = [cs_idxs(j), vertex_idxs(i), vertex_idxs(i), cs_idxs(j)]; % create a vector containing the outgoing arc and the incoming arc
%         outgoing_arc(:,i) = sum(idxs(:,1) == cs_idxs(j),2) & ...
%             sum(idxs(:,2) == vertex_idxs(i),2); % create a matrix with 1 in correspondence of the outgoing arc
%     end
%     
%     
%     outgoing_arc = sum(outgoing_arc,2); % create a unique vector with 1s in correspondence of all the outgoing arcs
%     x_temp_out = outgoing_arc; % create a unique vector of the outc and inc arcs
%     xtemp_out = [half_of_A; x_temp_out]; % total x temporary, with also the first half of A (still empty since I am not considering yet the z_ij variable) 
%     Aineq_out(s_out + counter_out,:) = xtemp_out; % fill at each iteration the A matrix with the xtemp vector
%     counter_out = counter_out + 1; 
%     
% end
% 
% for k = 1 : length(cs_idxs)
%     Aineq(k, :) = Aineq_out(k, :) + Aineq_in(end,:) ;
% end
% 
% bineq = [bineq; zeros(length(cs_idxs),1)]; % fill beq with ones
     
%% Constraint 2
%% Constraint 2.1
% EQUALITY Aeq, beq, x_ij variable
% There has to be only one arc arriving at the target
s = size(Aeq,1);  % s keeps track of A dimension
counter = 1;
incoming_arc = []; % vector that takes into account the incoming arcs

% Link between all vertexes and targets
for i = 1: length(target_idxs) % for all the targets
    for j = 1 : length(vertex_idxs) % for all the vertexes
        u = [vertex_idxs(j), target_idxs(i)];  % couple of j-th vertex and the i-th target
        incoming_arc(:,j) = sum(idxs(:,1) == vertex_idxs(j),2) & ...
            sum(idxs(:,2) == target_idxs(i),2); % create a matrix with 1 in correspondence of the incoming arcs (arrives in target)
    end
    incoming_arc = sum(incoming_arc,2); % create a unique vector with 1s in correspondence of all the incoming arcs
    xtemp  = sparse([half_of_A;incoming_arc]); % total x temporary, with also the first half of A
    Aeq(s + counter,:) = xtemp; % fill at each iteration the A matrix with the xtemp vector
    counter = counter + 1;
end

beq = [beq; ones(length(target_idxs),1)]; % fill beq with ones

     
%% Constraint 2.2
% EQUALITY Aeq, beq, x_ij variable
% There has to be only one arc leaving the target

s = size(Aeq,1);  % s keeps track of A dimension
counter = 1;
outgoing_arc = []; % vector that takes into account the outgoing arcs

% Link between targets and all vertexes
for i = 1: length(target_idxs) % for all the targets
    for j = 1 : length(vertex_idxs) % for all the vertexes
        u = [target_idxs(i),vertex_idxs(j)];  % pair the i-th target and the j-th vertex
        outgoing_arc(:,j) = sum(idxs(:,1) == target_idxs(i),2) & ...
            sum(idxs(:,2) == vertex_idxs(j),2); % create a matrix with 1 in correspondence of the outgoing arcs (leaves the target)
    end
    outgoing_arc = sum(outgoing_arc,2); % create a unique vector with 1s in correspondence of all the outgoing arcs
    xtemp  = sparse([half_of_A;outgoing_arc]); % total x temporary, with also the first half of A
    Aeq(s + counter,:) = xtemp;  % fill at each iteration the A matrix with the xtemp vector
    counter = counter + 1;
end

beq = [beq; ones(length(target_idxs),1)]; % fill beq with ones

     
%% Constraint 3
% EQUALITY Aeq, beq, x_ij variable, z_ij variable
% This constraint states that the consumption spent up to that point is
% equal to the consumption to travel the last arc plus the one to get there.
% This constraint is also responsible for the elimination of subtours.

s = size(Aeq,1);  % s keeps track of A dimension
counter = 1;
outgoing_arc = []; % vector that takes into account the outgoing arcs
incoming_arc = []; % vector that takes into account the incoming arcs
outgoing_arc_fuel = []; % vector that takes into account the outgoing arcs for the fuel consunption
xtemp = []; % temporary x

for j = 1 : length(target_idxs)  % for all targets
    for i = 1 : length(vertex_idxs) % for all vertexes
        if (vertex_idxs(i) == target_idxs(j)) % if the couple is composed by equal numbers pass to the next iteration and don't do the following commands
            continue;
        end
        u = [target_idxs(j), vertex_idxs(i), vertex_idxs(i), target_idxs(j)]; % create a vector containing the outgoing arc and the incoming arc
        outgoing_arc(:,i) = sum(idxs(:,1) == target_idxs(j),2) & ...
            sum(idxs(:,2) == vertex_idxs(i),2); % create a matrix with 1 in correspondence of the outgoing arc
        outgoing_arc_fuel(:,i) = sum(idxs(:,1) == target_idxs(j),2) & ...
            sum(idxs(:,2) == vertex_idxs(i),2); % create a matrix with 1 in correspondence of the outgoing arc fuel consumption
        incoming_arc(:,i) = sum(idxs(:,1) == vertex_idxs(i),2) & ...
            sum(idxs(:,2) == target_idxs(j),2); % create a matrix with 1 in correspondence of the incoming arc
    end
    
    
    outgoing_arc = sum(outgoing_arc,2); % create a unique vector with 1s in correspondence of all the outgoing arcs
    outgoing_arc_fuel = sum(outgoing_arc_fuel,2); % create a unique vector with 1s in correspondence of all the outgoing arcs fuel consumption
    for k = 1 : length(outgoing_arc_fuel)
        if (outgoing_arc_fuel(k) == 1)
            outgoing_arc_fuel(k) = -battery_tr(k);  % value of the fuel/battery consumption in the selected arcs
        end
    end
    incoming_arc = sum(incoming_arc,2); % create a unique vector with 1s in correspondence of all the incoming arcs
    incoming_arc = -1*incoming_arc; % change the incoming arcs value to -1
    x_temp = outgoing_arc + incoming_arc; % create a unique vector of the outc and inc arcs
    xtemp = [x_temp; outgoing_arc_fuel]; % total x temporary, with also the second half of A (still empty since I am not considering yet the x_ij variable)
    Aeq(s + counter,:) = xtemp; % fill at each iteration the A matrix with the xtemp vector
    counter = counter + 1;
    
end

beq = [beq; zeros(length(target_idxs),1)]; % fill beq with zeros

     
%% Constraint 5
% EQUALITY Aeq, beq, x_ij variable, z_ij variable
% Together with constraint 4 impose 0<=zij<=F and they ensure that the fuel
% consumed by the vehicle to travel up to a depot does not exceed the fuel capacity F

s = size(Aeq,1);  % s keeps track of A dimension
counter = 1;
outgoing_arc = []; % vector that takes into account the outgoing arcs
outgoing_arc_fuel = [];

% Link between targets and all vertexes
for i = 1: length(cs_idxs) % for all the cs
    for j = 1 : length(target_idxs) % for all the targets
        u = [cs_idxs(i),target_idxs(j)];  % pair the i-th cs and the j-th target
        outgoing_arc(:,1) = sum(idxs(:,1) == cs_idxs(i),2) & ...
            sum(idxs(:,2) == target_idxs(j),2); % create a matrix with 1 in correspondence of the outgoing arcs (leaves the target)
        outgoing_arc_fuel(:,1) = sum(idxs(:,1) == cs_idxs(i),2) & ...  % create a matrix with 1 in correspondence of the outgoing arcs (leaves the target)
            sum(idxs(:,2) == target_idxs(j),2);
        for k = 1 : length(outgoing_arc)
            if (outgoing_arc(k) == 1)
                outgoing_arc(k) = -battery_tr(k);  % insert instead of ones the real value of the fuel consumption in the selected arcs
            end
        end
        xtemp  = sparse([outgoing_arc_fuel; outgoing_arc]); % total x temporary, first half of A with fuel, second half with ones (selected arcs)
        Aeq(s + counter,:) = xtemp;  % fill at each iteration the A matrix with the xtemp vector
        counter = counter + 1;
    end
end

beq = [beq; zeros(length(target_idxs)*length(cs_idxs),1)]; % fill beq with ones


%% Constraint 4
%% Constraint 4.1
% INEQUALITY Aineq, bineq, x_ij variable, z_ij variable
% Together with constraint 5 impose 0<=zij<=F and they ensure that the fuel
% consumed by the vehicle to travel up to a depot does not exceed the fuel capacity F

s = size(Aineq,1);  % s keeps track of A dimension
counter = 1;
arcs = []; % vector that takes into account all the arcs

% Link between all vertexes-> all possinble arcs
for i = 1: length(vertex_idxs) % for all the vertexes
    for j = 1 : length(vertex_idxs) % for all the vertexes
        if (vertex_idxs(i) == vertex_idxs(j)) % if the couple is composed by equal numbers pass to the next iteration and don't do the following commands
            continue;
        end
        u = [vertex_idxs(i),vertex_idxs(j)];  % pair the i-th vertexes and the j-th vertexes
        arcs(:,1) = sum(idxs(:,1) == vertex_idxs(i),2) & ...
            sum(idxs(:,2) == vertex_idxs(j),2); % create a matrix with 1 in correspondence of the selected arc
        xtemp  = sparse([arcs; half_of_A]); % total x temporary, with also the first half of A
        xtemp = -xtemp;
        Aineq(s + counter,:) = xtemp; % fill at each iteration the A matrix with the xtemp vector
        counter = counter + 1;
    end
end

bineq = [bineq; zeros(length(idxs),1)]; % fill beq with ones


%% Constraint 4.2
% INEQUALITY Aineq, bineq, x_ij variable, z_ij variable
% Together with constraint 5 impose 0<=zij<=F and they ensure that the fuel
% consumed by the vehicle to travel up to a depot does not exceed the fuel capacity F

F = 1; % total fuel capacity--> total battery capacity = 100%
s = size(Aineq,1);  % s keeps track of A dimension
counter = 1;
arcs = []; % vector that takes into account all the arcs
total_fuel_capacity = [];

% Link between all vertexes-> all possinble arcs
for i = 1: length(vertex_idxs) % for all the vertexes
    for j = 1 : length(vertex_idxs) % for all the vertexes
        if (vertex_idxs(i) == vertex_idxs(j)) % if the couple is composed by equal numbers pass to the next iteration and don't do the following commands
            continue;
        end
        u = [vertex_idxs(i),vertex_idxs(j)];  % pair the i-th vertexes and the j-th vertexes
        arcs(:,1) = sum(idxs(:,1) == vertex_idxs(i),2) & ...
            sum(idxs(:,2) == vertex_idxs(j),2); % create a matrix with 1 in correspondence of the selected arc
        total_fuel_capacity(:,1) = sum(idxs(:,1) == vertex_idxs(i),2) & ...
            sum(idxs(:,2) == vertex_idxs(j),2); % create a matrix with 1 in correspondence of the selected arc
        
        for k = 1 : length(total_fuel_capacity)
            if (total_fuel_capacity(k) == 1)
                total_fuel_capacity(k) = -F;
            end
        end
        
        xtemp  = sparse([arcs; total_fuel_capacity]); % total x temporary, with also the first half of A
        Aineq(s + counter,:) = xtemp; % fill at each iteration the A matrix with the xtemp vector
        counter = counter + 1;
    end
end

bineq = [bineq; zeros(length(idxs),1)]; % fill beq with ones

%% Mixed-Integer Linear Programming Algorithms
% Solve the optimisation algorithm:
%   MILP-> I want to solve a linear objective function min_x f' * x subject to:
%   - linear constraints: A_ineq * x <= b_ineq, A_eq * x = b_eq
%   - bounds: lb <= x <= up
%   - restrictions on some x to have int value: intcon

intcon = length(battery_tr):2*lendist; % restrictions on some x to have int value -> Variabile binaria x per selezione percorso
lb = zeros(2*lendist,1); % lower bound: vector of zeros
ub = ones(2*lendist,1); % upper bound : vector of ones

tic

% set optimization options 
%opts = optimoptions('intlinprog','Display','on');
%opts = optimoptions('intlinprog','CutGeneration','intermediate','CutMaxIterations', 10,'IntegerPreprocess','advanced','Display','none','PlotFcn',@optimplotmilp);%'PlotFcn','RootLPAlgorithm','primal-simplex','BranchRule','strongpscost','Heuristics','intermediate','Display','off','PlotFcn',@optimplotmilp);
%opts = optimoptions('intlinprog','Display','off','AbsoluteGapTolerance',100);
opts = optimoptions('intlinprog','PlotFcn',@optimplotmilp,'AbsoluteGapTolerance',100); % ??? CHIEDI

% Pass to intlinprog: Linear objective function f, Restrictions on some x to have int value, Aineq, bineq, Aeq, beq, lowerbound, upperbound, options
[x_tsp, costopt, exitflag,output] = intlinprog(tr_time,intcon,Aineq,bineq,Aeq,beq,lb,ub,opts);

toc

if exitflag < 1 || isempty(x_tsp) % see meanings of exitflag by typing in the command window help intlinprog
    fprintf('# No possible solution, exitflag < 1');
end
segments = find(x_tsp(length(battery_tr):end)); % Get indices of x_tsp's rows wiht 1s (optimal paths)
xopt = x_tsp(length(battery_tr) + 1:end);
lh = updateSalesmanPlot(x_tsp(length(battery_tr) + 1:end),idxs_total,xloc,yloc); % x_tsp(length(battery_tr) + 1:end) is the second half of x_tsp vector
Lat = lh(:,1);
Lon = lh(:,2);
h_figure = figure;
plot(xloc(1:Target), yloc(1:Target), 'rs', xloc(Target+1:length(xloc)), yloc(Target+1:length(yloc) ), 'bo');
text(xloc(1:Target),yloc(1:Target),target_label,'VerticalAlignment','bottom','HorizontalAlignment','right');
text(xloc(Target+1:length(xloc)), yloc(Target+1:length(yloc)),cs_lalbel,'VerticalAlignment','bottom','HorizontalAlignment','right');
lgnd = legend('Target', 'Depot/Charging Station', 'Location', 'EastOutside');
lgnd.AutoUpdate = 'off';
xlim([0 20]);ylim([0 20])
hold on
plot(Lat,Lon,'k:','LineWidth',2);
title('Solution');

ub = ones(lendist,1);

  
     

