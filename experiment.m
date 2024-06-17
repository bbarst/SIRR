cla
clc;
clear
m = 1000;
n = 30;
s = 100;
lambda = 0;

conds = [-12,-8,-4];
res_sizes = [-3,-6];
method1_result = cell(length(conds), length(res_sizes));
method2_result = cell(length(conds), length(res_sizes));
qr_result = cell(length(conds), length(res_sizes));
for idx1 = 1:length(conds)
    k = conds(idx1);
    for idx2 = 1:length(res_sizes)
        e_length = 10^(res_sizes(idx2));
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
        t = 6;
        J = 2;
        num = 1;

        error_1 = zeros(1,t);
        %error_2 = zeros(1,floor(log(t)/log(J)));
        error_2 = zeros(1,t);
        be_list1=[];
        be_list2=[];
        xl1=zeros(n,t);
        %xl2=zeros(n,floor(log(t)/log(J)));
        xl2=zeros(n,t);
        tl1 = zeros(1,t);
        %tl2 = zeros(1,floor(log(t)/log(J)));
        tl2 = zeros(1,t);
        

        for i=1:num
        x_list1 = inf;

        S = normrnd(0,1,s,m)/sqrt(s);
        s_list = [];

        
  
        [~,x_list1,time_list1] = IDS_solver(A,b,s,2,J,t);
        [x_list2,time_list2] = IDS2(A,b,S,lambda,n,J,t);

        xl1 = xl1+x_list1;
        xl2 = xl2+x_list2;

        end
        xl1 = xl1/num;
        xl2 = xl2/num;
        error_1 = vecnorm(xl1-xstar)/norm(xstar);
        error_2 = vecnorm(xl2-xstar)/norm(xstar);
        eA1 = vecnorm(A*(xl1-xstar))/norm(b);
        eA2 = vecnorm(A*(xl2-xstar))/norm(b);
        for j = 1:t
            be_list1 =[be_list1, backward_error_ls(A,b,xl1(:,j))];
        end
        for j = 1:t
            be_list2 =[be_list2, backward_error_ls(A,b,xl2(:,j))];
        end
        method1_result{idx1,idx2}=[error_1;eA1;be_list1];
        method2_result{idx1,idx2}=[error_2;eA2;be_list2];
        
        %[Q,R]=qr(A,0);
        %xqr = R\(Q'*b);
        %xqr = (transpose(A)*A+lambda*eye(n))\(transpose(A)*b);
        xqr = A\b;
        qr_result{idx1,idx2}=[norm(xqr-xstar)/norm(xstar);norm(A*(xqr-xstar))/norm(b);norm(A*(xqr-xstar))/norm(b)];
        

    end
end

%add legend
legend_method = {' IHS',' IDS','Matlab Direct'};
length1 = length(conds);
length2 = length(legend_method);

legend_cond = cell(1,length1);
for i =1:length1
    legend_cond{i}=strcat('kappa=10^',num2str(-conds(i)));
end

legend_list = cell(1,length1*length2);
for i =length1:-1:1
    for j = 1:length2
        legend_list{length2*(3-i)+j} = strcat(legend_cond{i},legend_method{j});
    end
    
end

close all
colors = ['#A2A9DD','#C4A5DE','#F0988C',"#EEB66D","#D94F33","#834026"];
close all
for idx2 = 1:length(res_sizes)
    res_size = res_sizes(idx2);
    for idx1 = length(conds):-1:1
        condA = -conds(idx1);
        condA = 10^condA;
        for j = 1:3
            figure(j+3*(idx2-1));
            %set(gca,'position',[0.1 0.1 0.8 0.92]);
            %axis normal;
            method1_idx = method1_result{idx1,idx2};
            method2_idx = method2_result{idx1,idx2};
            qr_idx = qr_result{idx1,idx2};
            a = semilogy(method1_idx(j,:),'*--', 'LineWidth', 1,...
                'Color', colors(idx1)); hold on
            a.Color(4)=0.5;
            semilogy(method2_idx(j,:),'-', 'LineWidth', 2,...
                'Color', colors(idx1)); hold on
            yline(qr_idx(j),':', 'LineWidth', 2, 'Color', colors(idx1))
            %legend(legend_list)
            if j == 1
                ylabel('Forward Error','FontSize',14);
                legend({'','\kappa = 10^{2}','','','\kappa = 10^{4}','','','\kappa = 10^{8}',''},'Location','best')
            end
            if j == 2
                ylabel('Residual Error','FontSize',14);
            end
            if j == 3
                ylabel('Backward Error','FontSize',14);
            end
            xlabel('Iteration i');

        end
    end
end
%saveas(gca,'final.pdf')




for i=1:2
base = strcat('ihsres/res',num2str(i));
fn = strcat(base,'_fe_saa_ihsids.png');
saveas(figure(3*i-2),fn)
fn = strcat(base,'_re_saa_ihsids.png');
saveas(figure(3*i-1),fn)
fn = strcat(base,'_be_saa_ihsids.png');
saveas(figure(3*i),fn)
end