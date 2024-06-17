
function [x_list,time_list] = IDS2(A,b,S,lambda,n,J,t)

SA = S*A;
Ab = transpose(A)*b;

tic


[Q1,R1]=qr(SA,0);
y=transpose(R1)\Ab;
x = R1\y;

time = toc;
bonus = x;
x_list = [x];
time_list = [time];
for i = 1:t-1
    tic
    for j = 1:J
    r1 = (transpose(A)*((b-A*x))-lambda*(x));
    uj = IDS2in(A,r1,R1,lambda,J,i+1);
    x = x+uj;
    end
    x_list = [x_list,x];
    time = toc;
    time_list = [time_list,time];
end
end


function x = IDS2in(A,Ab,R1,lambda,J,t)


xy = R1\(transpose(R1)\(Ab));
r2 = Ab-(transpose(A)*(A*xy));
xy1 = R1\(transpose(R1)\(r2));
X = [xy,xy1];
ARY=A*X;

[QAX,RAX] = qr(transpose(ARY)*ARY,0);
a = RAX\(transpose(QAX)*(transpose(X)*(Ab)));

x = X*a;

if t>1
    for i = 1:J
    r1 = Ab-(transpose(A)*(A*x)-lambda*(x));
    u =  IDS2in(A,r1,R1,lambda,J,t-1);
    x=x+u;
    end
end

end



