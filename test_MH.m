function test_MH
    X = randn(81,1);
    [e1,e2] = M_H_algorithm(X, 100000);

end

function [e1,e2] = M_H_algorithm(X, MCsize)

    if isempty(MCsize)
        MCsize = 1000;
    end
    c = 3;
    p = 1;
    N = 10;
    theta1 = 0.3;
    theta2 = 0.3;
    
    D1 = [sparse(-eye(c^2*p*(N - p -1))), sparse(zeros(c^2*p*(N - p -1),c^2*p ))];
    D2 = [sparse(zeros(c^2*p*(N - p -1),c^2*p )), sparse(eye(c^2*p*(N - p -1)))];
    D = D1+D2;
    
    tempR = randn(c^2*p*(N-p));
    R = tempR*tempR';
    e = spdiags(ones(c^2*p*(N-p)),-9:9,c^2*p*(N-p),c^2*p*(N-p));
    R = R.*e;
    
    Q = randn(1,c^2*p*(N-p));
    %%%% proposal distributuion : laplace approximate
    sign1 = sign(X);
   
    DX = D*X;
    sign2 = D'*sign(DX);

    mu1 = theta1*c^2*p*(N-p);
    mu2 = theta2*c^2*p*(N-p-1);

    A1 = sparse((mu1/((sum(abs(X)))^2))*(sign1*sign1'));
    A2 = sparse((mu2/((sum(abs(DX)))^2))*(sign2*sign2'));

    A = 2*R-A1-A2;

    %%%% Monte Carlo
    [L,err] = cholcov(A);
    if err ~= 0
        error(message('stats:mvnrnd:BadCovariance2DSymPos'));
    end
    %%%% problem here to be continue
    randvec = randn(c^2*p*(N-p),MCsize);
    MCsample = L\randvec+repmat(X,1,MCsize);

    %%%% M-H
    p_tilde = @(x) exp(-x'*R*x+Q*x-mu1*log(sum(abs(x)))-mu2*log(sum(abs(D*x))));
    naccept = 0;
    for i = 1:MCsize
        ak = min(1,p_tilde(MCsample(:,i))/p_tilde(X));
        if rand(1) > ak
            MCsample(:,i) = X;
            naccept = naccept + 1;
        end

    end

    e1 = mean(log(sum(abs(MCsample),1)));
    e2 = mean(log(sum(abs(D*MCsample),1)));



    %q_tilde = @(z)


end