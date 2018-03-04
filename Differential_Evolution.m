clear all
% f - optimal fitness
% X - optimal solution
% CONTROL PARAMETERS %
D = 3; % dimension of problem
NP =100; % size of population
F = 0.8; % differentiation constant
CR = 0.9; % crossover constant
GEN = 20000; % number of generations
L = [3 1 1.01]; % low boundary constraint
H = [6 2 3]; % high boundary constraint
% *************************** %
% ** ALGORITHM’S VARIABLES ** %
% *************************** %
X = zeros(D,1); % trial vector
Pop = zeros(D,NP); % population
Fit = zeros(1,NP); % fitness of the population
iBest = 1; % index of the best solution
r = zeros(3,1); % randomly selected indices
% *********************** %
% ** CREATE POPULATION ** %
% *********************** %
% initialize random number generator
for j = 1:NP % initialize each individual
    Pop(:,j) = L + (H-L)*rand(D,1); % within b.constraints
    Fit(1,j) = fnc4(Pop(:,j)); % and evaluate fitness
end
% ****************** %
% ** OPTIMIZATION ** %
% ****************** %
for g = 1:GEN % for each generation
    % choose three random individuals from population,
    % mutually different and different from j
    r(1) = floor(rand()* NP) + 1;
    while r(1)==j
        r(1) = floor(rand()* NP) + 1;
    end
    r(2) = floor(rand()* NP) + 1;
    while (r(2)==r(1))||(r(2)==j)
        r(2) = floor(rand()* NP) + 1;
    end
    r(3) = floor(rand()* NP) + 1;
    while (r(3)==r(2))||(r(3)==r(1))||(r(3)==j)
        r(3) = floor(rand()* NP) + 1;
    end
    % create trial individual
    % in which at least one parameter is changed
    Rnd = floor(rand()*D) + 1;
    for i = 1:D
        if ( rand()<CR ) || ( Rnd==i )
            X(i) = Pop(i,r(3))+F*(Pop(i,r(1))-Pop(i,r(2)));
        else
                X(i) = Pop(i,j);
        end
    end
    % verify boundary constraints
    for i = 1:D
        if (X(i)<L(i))||(X(i)>H(i))
            X(i) = L(i) + (H(i)-L(i))*rand();
        end
    end
    % select the best individual
    % between trial and current ones
    % calculate fitness of trial individual
    f = fnc4(X);
    % if trial is better or equal than current
    if f <= Fit(j)
        Pop(:,j) = X; % replace current by trial
        Fit(j) = f ;
        % if trial is better than the best
        if f <= Fit(iBest)
            iBest = j ; % update the best’s index
        end
    end
end
 
% ************* %
% ** RESULTS ** %
% ************* %
f = Fit(iBest)
X = Pop(:,iBest)
