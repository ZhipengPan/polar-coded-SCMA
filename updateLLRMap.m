function [L,B] = updateLLRMap(lambda,phi,n,L,B)

    if lambda == 0
        return;
    end
    psi = floor(phi/2);
    if mod(phi,2)==0
        [L,B] = updateLLRMap(lambda-1,psi,n,L,B);
    end
    for omega=0:2^(n-lambda)-1
        if mod(phi,2)==0
            %do sth
            L(phi+omega*2^lambda+1,n+1-lambda) = fFunction(L(psi+2*omega*2^(lambda-1)+1,n+2-lambda),L(psi+(2*omega+1)*2^(lambda-1)+1,n+2-lambda)+B(phi+1+omega*2^lambda+1,n+1-lambda));
        else
            %do sth
            L(phi+omega*2^lambda+1,n+1-lambda) = L(psi+(2*omega+1)*2^(lambda-1)+1,n+2-lambda)+fFunction(L(psi+2*omega*2^(lambda-1)+1,n+2-lambda),B(phi-1+omega*2^lambda+1,n+1-lambda));
        end
    end
end
