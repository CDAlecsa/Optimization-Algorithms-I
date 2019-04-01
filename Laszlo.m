function [itL, ObjL, solL] = Laszlo(eps,u0,h, F, gradF)


    alpha = h ;
    
    ok = 0 ; 
    it = 0 ;
    
    uOld = u0 ;
    uNew = uOld ;
    
    ObjL = [0] ;
    solL = [uOld] ;
    
    while ok == 0 
       
        
        if it == 0 
            v = uOld ;
        else
            beta = 0.9 * it/(it+3) ;
            v = uOld + beta * (uOld - uOld2) ;      % auxiliary iteration
        end
        
        
        uNew = v - alpha * gradF(v) ;               % new solution 
        

        it = it + 1 ;
         
        solL = [solL uNew] ; 
        ObjL = [ObjL norm(F(uNew)-F(uOld))] ; 

        

        if ObjL(end) <= eps    
           ok = 1 ; 
        end
        
        
        uOld2 = uOld ;
        uOld = uNew ;
        
    end

itL = 1:it+1 ;

end