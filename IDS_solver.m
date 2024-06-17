function [x,x_list,time_list] = IDS_solver(A,b,varargin)
%IDS Solve A*x = b in the least-squares sense using IDS
%   Optional parameters (use [] for default value):
%   - d: sketching dimension (default 12*size(A,2))
%   - k: dimension of KryLov space needed in each accleration
%   - J: number of steps of iteration in each refinement
%   - t: number of refinement.
%   - summary: a function of the current iterate to be recorded at each
%     iteration. All summary values will be rows of stats, the second
%     output of this function.
%   - verbose: if true, print at each iteration. (default false)
%   - reproducible: if true, use slower, reproducible implementation of
%     sparse sign embeddings (default false)

    %count time
    tic
    Anum = isnumeric(A);
    if Anum

        m = size(A,1);
        n = size(A,2);
        Afun = @(x,op) mul(A,x,op);
    else
        scale = 1;
        m = size(b,1);
        n = size(A(zeros(size(b,1),0), true),1);
        Afun = A;
    end

    if length(varargin) >= 1 && ~isempty(varargin{1})
        d = varargin{1};
    else
        d = 12*n;
    end
    
    if length(varargin) >= 2 && ~isempty(varargin{2})
        k = varargin{2};
    else
        k = 2;
    end
    
    if length(varargin) >= 3 && ~isempty(varargin{3})
        J = varargin{3};
    else
        J = 2;
    end
    
    if length(varargin) >= 4 && ~isempty(varargin{4})
        t = varargin{4};
    else
        t = 4;
    end
    
%     if length(varargin) >= 5 && ~isempty(varargin{5})
%         summary = varargin{3};
%     else
%         summary = [];
%     end

%     if length(varargin) >= 5 && ~isempty(varargin{5})
%         verbose = varargin{5};
%     else
%         verbose = false;
%     end

    if length(varargin) >= 5 && ~isempty(varargin{5})
        reproducible = varargin{5};
    else
        reproducible = false;
    end
    
    
    
    if ~reproducible && exist('sparsesign','file') == 3
        S = sparsesign(d,m,8);
    else
        warning(['Using slower and slightly incorrect backup ' ...
            'implementation sparse_sign_backup.m. For the better ' ...
            'implementation, build the mex file `sparsesign.c` ' ...
            'using the command `mex sparsesign.c`.']);
        S = sparse_sign_backup(d,m,8);
    end


    stats = [];
    x_list = zeros(n,t);
    time_list = zeros(1,t);
    time=toc;
    if Anum
        SA = full(S*A);
    else
        SA = full(Afun(S',true)');
    end
%     fprintf('sketch time %d\n',toc-time)
    [~,R] = qr(SA,0);
    Ab = Afun(b,true);
    y = transpose(R)\Ab;
    x = R\y;
    time = toc;
    x_list(:,1) = x;
    time_list(1) = time;
    
    for i = 1:t-1
        for j = 1:J
            r = b-Afun(x,false);
            Ar = Afun(r,true);
            uj = IDS2in(A,Ar,R,J,i+1,k,n,Afun);
            x = x+uj;
        end
        x_list(:,i+1) = x;
        time = toc-time;
        time_list(i+1) = time;
        
    end
    
end


function x = IDS2in(A,Ar,R,J,t,K,n,Afun)
%solving Ax=r in y-space y=R*x using KryLov space
X = zeros(n,K);
xy = R\(R'\(Ar));
r2 = Ar-(transpose(A)*(A*xy));
xy1 = R\(transpose(R)\(r2));
X = [xy,xy1];
% X(:,1) = xy;

% for k = 1:K-1
%     r2 = Ar-Afun(Afun(xy,false),true);
%     xy = R\(R'\r2);
%     X(:,k+1) = xy;
% end

AX = A*X;
[QAX,RAX] = qr(AX'*AX,0);
a = RAX\(QAX'*(X'*Ar));
x = X*a;

if t>1
    for i = 1:J
    r1 = Ar-Afun(Afun(x,false),true);
    u =  IDS2in(A,r1,R,J,t-1,K,n,Afun);
    x = x+u;
    end
end

end


function b = mul(A, x, op)
    if op
        b = A'*x;
    else
        b = A*x;
    end
end