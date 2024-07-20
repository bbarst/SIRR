cla
clc;
clear
m = 1000;
n = 30;
s = 100;
lambda = 0;

conds = [-9];
res_sizes = [-1];
method1_result = cell(length(conds), length(res_sizes));
method2_result = cell(length(conds), length(res_sizes));
qr_result = cell(length(conds), length(res_sizes));
k = conds(1);

e_length = 10^(res_sizes(1));
%generate matrix
R = normrnd(0,1,m,n);
[U,~] = qr(R,0);
R = normrnd(0,1,n,n);
[V,~] = qr(R,0);

Sigma = diag(logspace(0,k,n));

A = U*Sigma*transpose(V);
x0 = normrnd(0,1,n,1);
e = normrnd(0,1,m,1);
e = e-U*transpose(U)*e;
e = e/norm(e)*e_length;
b = A*x0 + e;

xstar = x0;
%xstar = (transpose(A)*A+lambda*eye(n))\(transpose(A)*b);
%solve
t1 = 10;
t2 = 5;
J = 2;
J1 = 5;
K = 2;
num = 50;


E1 = [0;0;0];
E2 = [0;0;0];
E3 = [0;0;0];
for i=1:num
x_list1 = inf;

S = normrnd(0,1,s,m)/sqrt(s);
s_list = [];



[x1,x_list1,time_list1] = IDS_solver(A,b,s,K,[J,J1],t2);
[x2,x_list2,~,~] = fossils(A,b,s,[],[],true);
[x_list,time_list] = IDS3(A,b,S,lambda,n,3,50);

x3 = x_list(:,end);
e1t = norm(x1-xstar)/norm(xstar);
e2t = norm(x2-xstar)/norm(xstar);
e3t = norm(x3-xstar)/norm(xstar);
eA1t = norm(A'*(A*x1-b))/norm(A'*b);
eA2t = norm(A'*(A*x2-b))/norm(A'*b);
eA3t = norm(A'*(A*x3-b))/norm(A'*b);
be1t = backward_error_ls(A,b,x1);
be2t = backward_error_ls(A,b,x2);
be3t = backward_error_ls(A,b,x3);

E1(:,end+1) = [e1t;eA1t;be1t];
E2(:,end+1) = [e2t;eA2t;be2t];
E3(:,end+1) = [e3t;eA3t;be3t];

end
ans1=sum(E1,2)/num
ans2=sum(E2,2)/num

plot(log(E1(2,:)),'LineWidth', 3)
hold on
plot(log(E2(2,:)),'LineWidth', 3)
hold on
plot(log(E3(2,:)),'LineWidth', 3)
legend('IDS','fossil','alg2')
name = strcat('kappa:1e',num2str(-k),'   size of error:',num2str(e_length));
xlabel(name)

function [x_list,time_list] = IDS(A,b,s,lambda,n,J,t,x1)
m=size(A,1);
S = sparse_sign_backup(s,m,8);
SA = S*A;
Ab = transpose(A)*b;

tic


[Q1,R1]=qr(SA,0);
y=transpose(R1)\Ab;
x = R1\y;
if x1 ~= 0
    x=x1;
end
time = toc;
bonus = x;
x_list = [x];
time_list = [time];
for i = 1:t-1
    tic
    r1 = A'*(b-A*x);
    uj = IDSin(A,r1,R1,lambda,J,2);
    x = x+uj;
    
    
    x_list = [x_list,x];
    time = toc;
    time_list = [time_list,time];
end
end

function x = IDSin(A,Ab,R1,lambda,J,t)
x = R1\(transpose(R1)\(Ab));
for i = 1:t
    r1 = Ab-(transpose(A)*(A*x)-lambda*(x));
    u =  R1\(transpose(R1)\(r1));
    x=x+u;  
    
end
end

function [x_list,time_list] = IDS3(A,b,S,lambda,n,J,t)

SA = S*A;
SSb = transpose(S)*(S*b);
Ab = transpose(A)*b;


SASA = transpose(SA)*SA;
M = SASA+lambda*eye(n);
[Q,R] = qr(M,0);
tic
x = R\(transpose(Q)*Ab);
[Q1,R1]=qr(SA,0);
y = R1*x;
y=transpose(R1)\Ab;
x = R1\y;
time = toc;
bonus = x;
x_list = [x];
time_list = [time];
for i = 1:t-1
    r1 = A'*(b-A*x);
    uj = IDS2in(A,r1,SASA,Q,R,R1,lambda,J,3);
    x = x+uj;
    x_list = [x_list,x];
    time = toc;
    time_list = [time_list,time];
end
end

function x = IDS2in(A,Ab,SASA,Q,R,R1,lambda,J,t)
xy = R1\(transpose(R1)\(Ab));
r2 = Ab-(transpose(A)*(A*xy));
xy1 = R1\(transpose(R1)\(r2));
X = [xy,xy1];
ARY=A*X;

[QAX,RAX] = qr(transpose(ARY)*ARY,0);
a = RAX\(transpose(QAX)*(transpose(X)*(Ab)));

x = X*a;
global order
order = order+1;
if t>1
    for i = 1:J
    r1 = Ab-(transpose(A)*(A*x)-lambda*(x));
    u =  IDS2in(A,r1,SASA,Q,R,R1,lambda,J,t-1);
    x=x+u;  
    end
end
end