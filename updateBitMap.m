function [L,B] = updateBitMap(lambda,phi,n,L,B)

    
    psi = floor(phi/2);
    if mod(phi,2)~=0
        for omega = 0:2^(n-lambda)-1
            B(psi+2*omega*2^(lambda-1)+1,n+2-lambda) = fFunction(B(phi-1+omega*2^(lambda)+1,n+1-lambda),B(phi+omega*2^(lambda)+1,n+1-lambda)+L(psi+(2*omega+1)*2^(lambda-1)+1,n+2-lambda));
            B(psi+(2*omega+1)*2^(lambda-1)+1,n+2-lambda) = B(phi+omega*2^(lambda)+1,n+1-lambda) + fFunction(B(phi-1+omega*2^(lambda)+1,n+1-lambda),L(psi+2*omega*2^(lambda-1)+1,n+2-lambda));
        end
        if mod(psi,2)~=0
            [L,B] = updateBitMap(lambda-1,psi,n,L,B);
        end
    end
end