function [x_list,time_list] = IDS2SVD(A,b,S,lambda,n,J,t)

scale = vecnorm(A);
A= A ./ scale;

SA = S*A;
SSb = transpose(S)*(S*b);
Ab = transpose(A)*SSb;




[U,svals,V] = svd(SA,'econ'); svals = diag(svals);
svals2 = svals.*svals;

tic
x = V*(((V'*Ab)./svals2));

time = toc;
bonus = x;
x_list = [x ./ scale.'];
time_list = [time];
for i = 1:t-1
    tic
    for j = 1:J
    r1 = (transpose(A)*((b-A*x))-lambda*(x));
    uj = IDS2in(A,r1,V,svals2,lambda,J,i+1);
    x = x+uj;
    end
    x_list = [x_list,x ./ scale.'];
    time = toc;
    time_list = [time_list,time];
end
end


function x = IDS2in(A,Ab,V,svals2,lambda,J,t)


xy = V*(((V'*Ab)./svals2));;
r2 = Ab-(transpose(A)*(A*xy));
xy1 = V*(((V'*r2)./svals2));;
X = [xy,xy1];
ARY=A*X;

[QAX,RAX] = qr(transpose(ARY)*ARY,0);
a = RAX\(transpose(QAX)*(transpose(X)*(Ab)));

x = X*a;

if t>1
    for i = 1:J
    r1 = Ab-(transpose(A)*(A*x)-lambda*(x));
    u =  IDS2in(A,r1,V,svals2,lambda,J,t-1);
    x=x+u;
    end
end

end
