function [drone_paths] =  dronePath_assignment(CS, Target, Vertex, idxs, vertex_idxs, target_idxs, cs_idxs, x_solution)

% Build the start_cs vector
start_cs = [];
for i = (length(idxs)-(Vertex-1)*CS +1): length(idxs)
    add_start = idxs(i,1);
    start_cs = [start_cs, add_start];
end

% Build the end_cs vector
end_cs = [];
for i = (length(idxs)-(Vertex-1)*CS +1): length(idxs)
    add_end = idxs(i,2);
    end_cs = [end_cs, add_end];
end

% Build the arc_cs vector
arc_cs = [];
for i = (length(idxs)-(Vertex-1)*CS +1): length(idxs)
    add_arc = x_solution(i);
    arc_cs = [arc_cs, add_arc];
end

% Build the start vector
start = [];
for i = 1: length(idxs)
    add_start = idxs(i,1);
    start = [start, add_start];
end

% Build the thend vector
thend = [];
for i = 1: length(idxs)
    add_end = idxs(i,2);
    thend = [thend, add_end];
end

arc = x_solution;
drone_paths = []; % all paths

% Since we know that the paths will start from a CS, the examinations starts from the CS
for j = 1: length(arc_cs) % check inside the arc_cs array where are the selected arcs
    my_loop = []; % to keep track of internal loops 
    arc_cs(j) = cast(arc_cs(j), 'int8');
    if (arc_cs(j) == 1)
        my_start = start_cs(j);
        my_end = end_cs(j);
        my_loop = [my_loop, my_start];
        my_loop = [my_loop, my_end];
%         counter = 0;
%         maxiterations = 500;
        while(my_end ~= my_start && counter < maxiterations)
            for i = ((my_end-1)*(Vertex-1) +1): (my_end)*(Vertex-1) % check inside the arc array where are the selected arcs
                arc(i) = cast(arc(i), 'int8');
                if (arc(i) == 1) % and thend[i] != my_start
                    my_end = thend(i); % memorize the end point of this link
                    my_loop = [my_loop, my_end];
                end
            end
%             counter = counter + 1;
        end
        fprintf("i")  % NON ESCE DAL WHILE
        drone_paths = [drone_paths, my_loop];
    end  
end


% Drones' list
drones = [];
for m = 1: CS
    value = m;
    drones = [drones, value];
end
fprintf("\nThe list of drones is: ")
drones

drone_paths

% Divide the drone_paths every time a new path starts
paths_cellarr = {}; % initialize a cell array
paths_cellarr{1} = drone_paths(1);
for p = 2:length(drone_paths) 
    % if drone_paths(p) ~= drone_paths(p-1)
    % if (not(ismember(drone_paths(p), cs_idxs)) & not(ismember(drone_paths(p-1), cs_idxs))) % || ismember(drone_paths(p),cs_idxs) (se i due sono successivi) o (controlla se i due successivi sono apparteneti cs_idxs) 
    if not(ismember(drone_paths(p), cs_idxs))
        size = length(paths_cellarr);
        paths_cellarr{size} = [paths_cellarr{size},drone_paths(p)];
    elseif (ismember(drone_paths(p), cs_idxs) & not(ismember(drone_paths(p-1), cs_idxs)))
        size = length(paths_cellarr);
        paths_cellarr{size} = [paths_cellarr{size},drone_paths(p)];
    elseif (ismember(drone_paths(p), cs_idxs) & ismember(drone_paths(p-1), cs_idxs))
        size = length(paths_cellarr)+1;
        paths_cellarr{size} = [drone_paths(p)];
    end
    
end
fprintf('\nThe drone paths are:')
paths_cellarr


% Create a cell array containing the path associated to a drone
assoc_paths{1} = "Path[] Drone[n]";
% Ask if the user wants to associate a drone to a specific path
for n = 1 : length(paths_cellarr) % for all the paths
    for o = 1 : length(drones) % for all the drones
        if (paths_cellarr{n}(1) == (drones(o)+Target))
            fprintf('\nWould you like to assign this path: ')
            paths_cellarr{n}
            fprintf('to drone number?')
            drones(o)
            fprintf("Type y for yes and n for no");
            my_input = input('\nEnter your input:', 's');
            % isstring(my_input)
            if (my_input == "y")
                 size = length(assoc_paths)+1;
                 % assoc_paths{size}(1) = drones(o);
                 assoc_paths{size} = paths_cellarr{n};  
                 u = length(assoc_paths{size});
                 assoc_paths{size}(u+1) = drones(o);
             end
        end
    end
end
fprintf("\nThe associations path-drone are: ")
assoc_paths

end

