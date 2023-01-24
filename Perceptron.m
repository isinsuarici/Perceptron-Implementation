
clear all;
close all;
clc;

load("bim488_hw2.mat")
c1x = classes(1,1:100);
c2x = classes(1,101:200);
c1y = classes(2,1:100);
c2y = classes(2,101:200);

plot(c1x, c1y, 'o');
hold
plot(c2x, c2y, '.');

rho=0.01;
w_init1 = [1, 1, -0.5]';
w_init2 = [1, 1, 0.5]';

% input for bias is 1, so I added them to 3rd row
A= ones(1,200);
classes= [classes;A];

[weight_1]=perceptron(classes,class_labels,w_init1,rho);
disp('Final weight vector for w_init1 is:')
disp(weight_1');
x1(1,1)=0;
x1(1,2)=-weight_1(3,1)/weight_1(2,1);
x2(1,1)=-weight_1(3,1)/weight_1(1,1);
x2(1,2)=0;
line(x1,x2,Color="red") % point1 = (0, -b/w2) and point2 = (-b/w1 ,0)

[weight_2]=perceptron(classes,class_labels,w_init2,rho);
disp('Final weight vector for w_init2 is:')
disp(weight_2');
x12(1,1)=0;
x12(1,2)=-weight_2(3,1)/weight_2(2,1); 
x22(1,1)=-weight_2(3,1)/weight_2(1,1); 
x22(1,2)=0; 
line(x12,x22); % point1 = (0, -b/w2) and point2 = (-b/w1 ,0)

function [w]=perceptron(X,y,w_init,rho)

[row,col]=size(X);
max_iter=50000; % for control
w=w_init;        
iter_num=0;         
missclasified=col;     

while(missclasified>0)&&(iter_num<max_iter)
    iter_num=iter_num+1;
    missclasified=0;
    
    sum_symbolpart=zeros(row,1);
    for i=1:col
        if((X(:,i)'*w)<=0 && (y(i)==1) ) % missclasified case-1
            % e(i) = d(i) - y(i)
            % I asssumed that transfer function is +1 for c1 and 0 for c2
            % So, e(i) = 1 - 0 = 1
            % (But if it is +1 for c1 and -1 for c2, e(i) would be 2 for c1
            % and -2 for c2.)
            e(i)=1;
            missclasified=missclasified+1;
            sum_symbolpart=sum_symbolpart+(e(i)*X(:,i));
            
        elseif ((X(:,i)'*w)>0 && (y(i)==2) ) %missclasified case-2
            % e(i) = -1 - 0 = -1
            e(i)=-1;
            missclasified=missclasified+1;
            sum_symbolpart=sum_symbolpart+(e(i)*X(:,i));
        end        
    end
    w=w+rho*sum_symbolpart;
end
end


