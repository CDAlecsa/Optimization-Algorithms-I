function [itN, ObjN, solN] = Nesterov(eps,u0,h,F,gradF)

    alpha = h ;
    
    ok = 0 ; 
    it = 0 ;
    
    uOld = u0 ;
    uNew = uOld ;
    
    ObjN = [0] ;
    solN = [uOld] ;
    
    while ok == 0 
       
        
        if it == 0 
            v = uOld ;
        else
            beta = 1 * it/(it+3) ;
            v = uOld + beta * (uOld - uOld2) ;      % auxiliary iteration
        end
        
        
        uNew = v - alpha * gradF(v) ;               % new solution 
        

        it = it + 1 ;
         
        solN = [solN uNew] ; 
        ObjN = [ObjN norm(F(uNew)-F(uOld))] ; 
        

        if ObjN(end) <= eps   
           ok = 1 ; 
        end
        
        
        uOld2 = uOld ;
        uOld = uNew ;
        
    end

itN = 1:it+1 ;

end