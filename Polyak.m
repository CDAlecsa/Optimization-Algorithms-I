function [itP, ObjP, solP] = Polyak(eps, u0, h, F, gradF, coeffPolyak)

    alpha = h ;
    
    ok = 0 ; 
    it = 0 ;
    
    uOld = u0 ;
    uNew = uOld ;
    
    ObjP = [0] ;
    solP = [uOld] ;
    
    while ok == 0 
       
        
        if it == 0 
            v = uOld ;
        else
            
            if coeffPolyak == 1
                beta = it/(it+3) ;
            elseif coeffPolyak == 2
                beta = 0.01 * it/(it+3) ;
            end
            
            v = uOld + beta * (uOld - uOld2) ;      % auxiliary iteration
        end
        
        
        uNew = v - alpha * gradF(uOld) ;               % new solution 
        

        it = it + 1 ;
         
        solP = [solP uNew] ; 
        ObjP = [ObjP norm(F(uNew)-F(uOld))] ; 
        

        if ObjP(end) <= eps   
           ok = 1 ; 
        end
        
        
        uOld2 = uOld ;
        uOld = uNew ;
        
    end

itP = 1:it+1 ;

end