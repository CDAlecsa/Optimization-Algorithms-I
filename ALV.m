function [itALV, ObjALV, solALV] = ALV(eps, u0, h, F, gradF, coeffALV)


    ok = 0 ; 
    it = 0 ;
    
    uOld = u0 ;
    uNew = uOld ;
    
    ObjALV = [0] ;
    solALV = [uOld] ;
    
    while ok == 0 
       
        if coeffALV == 1
            alpha = h ; 
        elseif coeffALV == 2
            alpha = 0.001 * (it+2) / (it+3) ;
        elseif coeffALV == 3
            alpha = 0.001 * (it+4) / (it+3) ;
        end
        
        if it == 0 
            v = uOld ;
        else
            beta = (-0.1) * it/(it+3) ;
            v = uOld + beta * (uOld - uOld2) ;      % auxiliary iteration
        end
        
        
        uNew = v - alpha * gradF(v) ;               % new solution 
        

        it = it + 1 ;
        
        solALV = [solALV uNew] ; 
        ObjALV = [ObjALV norm(F(uNew)-F(uOld))] ; 


        
        if ObjALV(end) <= eps
           ok = 1 ; 
        end
        
        
        uOld2 = uOld ;
        uOld = uNew ;
        
    end

itALV = 1:it+1 ;

end