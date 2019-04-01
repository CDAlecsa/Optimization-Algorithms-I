%% Main Function
    
    clear ;
    close all ;
    clc ;


%% Choose Functions

strOk = 0 ;

while strOk == 0
    
    strFct = input('Function : ', 's') ;
    
    if strcmp(strFct, 'Rastrigin')
        
        F = @(X) 20 + (X(1)^2 - 10*cos(2*pi*X(1))) + (X(2)^2-10*cos(2*pi*X(2))) ;
        gradF = @(X) [ 2*X(1) + 20 * pi * sin(2*pi*X(1)) , 2*X(2) + 20 * pi * sin(2*pi*X(2))] ;
        strOk = 1;
        
    elseif strcmp(strFct, 'Logistic Regression with L2 regularization')
        
        m = input('Sample batch size for reg. log. reg (default 1500) : ') ;
        n = input('Dim. of the domain of the obj. fct. (default 50) : ') ;
        
        if (m <= 0) || (n <= 0)
            m = 1500 ;
            n = 50 ;
        end
        
        C = 1/m ;            % normalizing term
 
        set = [-1, 1] ;
        y = zeros(1,m) ;
        X = zeros(n,m) ;
        for i = 1 : m
            y(i) = set(randi(length(set))) ; 
            X(:,i) = randn(1,n) ;
        end

        y = y' ;
        vectOnes = ones(m, 1) ;

        tildeA = zeros(n,m) ; 
        for i = 1 : m
            tildeA(:,i) = y(i) * X(:,i) ;
        end

        A = tildeA' ;
        F = @(w) (-C) * vectOnes' * log(1 ./ ( 1 + exp( - A * w )) ) + 0.5 * norm(w)^2 ;
        gradF = @(w) (-C) * tildeA * ( 1 - 1 ./ ( 1 + exp( - A * w ) ) ) + w ;
        strOk = 1;
        
    elseif strcmp(strFct, 'Beale')
        
        F = @(X) (1.5-X(1)+X(1)*X(2))^2+(2.25-X(1)+X(1)*X(2)^2)^2+(2.625-X(1)+X(1)*X(2)^3)^2 ;
        gradF = @(X) [ (2*(1.5-X(1)+X(1).*X(2))).*(-1+X(2))+(2*(2.25-X(1)+X(1).*X(2).^2)).*(X(2).^2-1)+(2*(2.625-X(1)+X(1).*X(2).^3)).*(X(2).^3-1) , (2*(1.5-X(1)+X(1).*X(2))).*X(1)+(4*(2.25-X(1)+X(1).*X(2).^2)).*X(1).*X(2)+(6*(2.625-X(1)+X(1).*X(2).^3)).*X(1).*X(2).^2] ;
        strOk = 1;
        
    elseif strcmp(strFct, 'Rosenbrock')
        
        F = @(X) 100*(X(1)^2 - X(2))^2 + (X(1)-1)^2 ;
        gradF = @(X) [400*X(1)*(X(1)^2 - X(2)) + 2*(X(1)-1), -200*(X(1)^2 - X(2))];
        strOk = 1;
        
    else
        
        fprintf('Choose another function \n') ; 
    
    end
end

%% Parameters

    eps = input('Tolerance : ') ;
    
    if strcmp(strFct, 'Logistic Regression with L2 regularization') == 1
        u0 = input('Initial cond. as a column vect. of n elements, e.g. ones(n, 1) : ') ;
    else
        u0 = input('Initial cond. as a column vect. of 2 elements : ') ;
    end

%% Optimizers

    fprintf('\n\n Optimizer list ... \n\n') ;
    fprintf('ALV (constant coeff.) \n ALV (increas. seq.) \n ALV (decreas. seq.) \n Polyak \n Modified Polyak \n Laszlo \n Nesterov \n') ;
    
    vectOpt = input('\n\n Vector of optimizers : ') ;       % order is given above (from 1 to 7)
    
    vectOpt = unique(vectOpt) ;
    vectOpt(vectOpt < 1) = [] ;
    vectOpt(vectOpt > 7) = [] ;
    
    
    for k = 1 : length(vectOpt)
       
        if vectOpt(k) == 1
            coeffALV = 1 ;
            stepsizeALV = input('StepSize for ALV with const. coeff : ') ;
            [itALV, objALV, solALV] = ALV(eps,u0, stepsizeALV, F, gradF, coeffALV) ;
        
        elseif vectOpt(k) == 2
            coeffALV = 2 ;
            [itALV, objALV, solALV] = ALV(eps,u0, 1, F, gradF, coeffALV) ;
            
        elseif vectOpt(k) == 3
            coeffALV = 3 ;
            [itALV, objALV, solALV] = ALV(eps,u0, 1, F, gradF, coeffALV) ;
            
        elseif vectOpt(k) == 4
            coeffPolyak = 1 ;
            stepsizePolyak = input('StepSize for Polyak with const. coeff : ') ;
            [itP, objP, solP] = Polyak(eps,u0, stepsizePolyak, F, gradF, coeffPolyak) ;
            
        elseif vectOpt(k) == 5
            coeffPolyak = 2 ;
            stepsizePolyak = input('StepSize for modified Polyak with const. coeff : ') ;
            [itP, objP, solP] = Polyak(eps,u0, stepsizePolyak, F, gradF, coeffPolyak) ;
            
        elseif vectOpt(k) == 6
            stepsizeLaszlo = input('StepSize for Laszlo alg.: ') ;
            [itL, objL, solL] = Laszlo(eps,u0, stepsizeLaszlo, F, gradF) ;
            
        elseif vectOpt(k) == 7
            stepsizeNesterov = input('StepSize for Nesterov alg.: ') ;
            [itN, objN, solN] = Nesterov(eps,u0, stepsizeNesterov, F, gradF) ;
            
        end
        
        
    end
    
 %% Plots
 
    colors = { [1 0.5 0], [1 0 0], [0 0 1], [0 0 0], [0 1 0], [1 1 0], [.61 .51 .74]} ;
    markers = ['s', 'd', 'o', '+', 'x', 'v', '^'] ;
    
 
    if ismember(1, vectOpt) == 1
        semilogy(itALV(2:end), objALV(2:end),'Color', colors{1}, 'Marker', markers(1), 'DisplayName', 'ALV - const. coeff.') ;
        hold on ;
    end
    
    if ismember(2, vectOpt) == 1
        semilogy(itALV(2:end), objALV(2:end), 'Color', colors{2}, 'Marker', markers(2), 'DisplayName', 'ALV - increas. coeff.') ; 
        hold on ;
    end
    
    if ismember(3, vectOpt) == 1
        semilogy(itALV(2:end), objALV(2:end), 'Color', colors{3}, 'Marker', markers(3), 'DisplayName', 'ALV - decreas. coeff.') ; 
        hold on ;
    end
    
    if ismember(4, vectOpt) == 1
        semilogy(itP(2:end), objP(2:end), 'Color', colors{4}, 'Marker', markers(4), 'DisplayName', 'Polyak alg. ') ;
        hold on ;
    end
    
    if ismember(5, vectOpt) == 1
        semilogy(itP(2:end), objP(2:end), 'Color', colors{5}, 'Marker', markers(5), 'DisplayName', 'Modified Polyak alg. ') ; 
        hold on ;
    end
    
    if ismember(6, vectOpt) == 1
        semilogy(itL(2:end), objL(2:end), 'Color', colors{6}, 'Marker', markers(6), 'DisplayName', 'Laszlo alg. ') ; 
        hold on ;
    end
    
    if ismember(7, vectOpt) == 1
        semilogy(itN(2:end), objN(2:end), 'Color', colors{7}, 'Marker', markers(7), 'DisplayName', 'Nesterov alg. ') ; 
        hold on ;
    end

    
    legend() ;
    xlabel('Iterations') ;
    ylabel('$ \| g(x_{n+1}) - g(x_{n}) \| $','interpreter','latex') ;
    title(strFct) ;
    
    