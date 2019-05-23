function [] = vsmplot(inputx,inputy)

    av = 10; st = 2;
    gradient = zeros(1,av); intercept = zeros(1,av);

    for n = st:st+av
        coefficients = polyfit(inputx(n:n+5), inputy(n:n+5), 1);
        gradient(n-st+1) = coefficients (1);
        intercept(n-st+1) = coefficients (2);
    end 
    
    avgrad = mean(gradient);
    avint = mean(intercept); 
    
    ydat = inputy-avint-(avgrad.*inputx);
    
    plot(inputx,ydat)
    
end

