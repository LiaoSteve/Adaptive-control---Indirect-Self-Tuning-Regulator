function [ p,theta ] = rls_forgetting( p,theta,phi,y ,lambda)    
    k = p *phi'./(lambda*1+phi*p*phi');
    p=p-k*phi*p/lambda;
    theta=theta+k*(y-phi*theta);
end

