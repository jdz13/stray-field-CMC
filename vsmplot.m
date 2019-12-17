%% Needs the X and Y inputs for the plot. Just feed it the data you want the DC removed from. 

% If you just want to use the data, save the ydat variable in the LHS of
% the function, and use outside the function.

function [] = vsmplot(inputx,inputy)

    av = 10; %number of averages that you want to do 
    st = 2; % start point to average at  
    gap = 5; % gap between points;
    
    % initalise the variables
    gradient = zeros(1,av); intercept = zeros(1,av);

    for n = st:st+av
        coefficients = polyfit(inputx(n:n+gap), inputy(n:n+gap), 1);
        gradient(n-st+1) = coefficients (1);
        intercept(n-st+1) = coefficients (2);
    end 
    
    avgrad = mean(gradient); % calculate the average gradient
    avint = mean(intercept); % calculate the average intercept
    
    ydat = inputy-avint-(avgrad.*inputx); % take the DC off the data
    
    plot(inputx,ydat) % will plot the renormalised data
    
end

