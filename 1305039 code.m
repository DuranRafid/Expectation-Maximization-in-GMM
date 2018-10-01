k = 3;

ActualmewArr = {};
ActualsigmaArr = {};
X = [];

for i=1:k
   mew = [-5+(i-1)*5 -5+(i-1)*5];
   ActualmewArr{i} = mew;
   rng(i);
   sigma = randn(2,2);
   sigma = sigma* sigma';
   ActualsigmaArr{i} = sigma;
   rng(i);
   temp = mvnrnd(mew,sigma,500);
   X = [X;temp];
end

disp('Actual ones');
for i=1:k
    disp(ActualmewArr{i});
    disp(ActualsigmaArr{i});
end

plot(X(:,1),X(:,2),'*b');

for i=1:k 
    a = ActualsigmaArr{i}(1,1);
    b = ActualsigmaArr{i}(2,2);
    x0 = ActualmewArr{i}(1);
    y0 = ActualmewArr{i}(2);

    t = -pi:0.1:pi;
    x = x0 + a*cos(t);
    y = y0 + b*sin(t);
    hold on;
    plot(x,y,'b');
end

[n,d] = size(X);

mewArr = {};
sigmaArr = {};
wArr = {};

for i=1:k
   mew = randn(1,2);
   mewArr{i} = mew
   sigma = rand(2,2)*100;
   sigma = sigma* sigma';
   sigmaArr{i} = sigma;
   wArr{i} = rand(1,1);
end


disp('Initially');
for i=1:k
    disp(mewArr{i});
    disp(sigmaArr{i});
end

%calculate log likelihood

log_likelihood = 0;

for j=1:n
    s = 0;
    for i=1:k
        N(i,j) = mvnpdf(X(j,:),mewArr{i},sigmaArr{i});
        s = s + wArr{i}*N(i,j);
    end
    log_likelihood = log_likelihood + log(s);
end

prev_log_likelihood = log_likelihood - 1000;
        

% E step

while(abs(log_likelihood - prev_log_likelihood)>0.0000000001)
    for j=1:n
        sum_w = 0;
        for i=1:k
            P(i,j) = wArr{i}*N(i,j);
            sum_w = sum_w + P(i,j);
        end

        for i=1:k
            P(i,j) = P(i,j)/sum_w;
        end  
    end

%M Step

    for i=1:k
    mewArr{i} = P(i,:)*X / sum(P(i,:)) ;
    end


    for i=1:k
    sumsigma = zeros(2,2);
    for j=1:n
        sumsigma = sumsigma + P(i,j)*(X(j,:) - mewArr{i})'*(X(j,:) - mewArr{i});
    end
    sigmaArr{i} = sumsigma/ sum(P(i,:));
    end


    for i=1:k
        wArr{i} = sum(P(i,:)) / n ;
    end

    prev_log_likelihood = log_likelihood ;
    log_likelihood = 0;

    for j=1:n
    s = 0;
    for i=1:k
        N(i,j) = mvnpdf(X(j,:),mewArr{i},sigmaArr{i});
        s = s + wArr{i}*N(i,j);
    end
    log_likelihood = log_likelihood + log(s);
    end

end

disp('Derived ones');
for i=1:k
    disp(mewArr{i});
    disp(sigmaArr{i});
    disp(wArr{i});
end

for i=1:k 
    a = sigmaArr{i}(1,1);
    b = sigmaArr{i}(2,2);
    x0 = mewArr{i}(1);
    y0 = mewArr{i}(2);

    t = -pi:0.1:pi;
    x = x0 + a*cos(t);
    y = y0 + b*sin(t);
    hold on;
    plot(x,y,'r');
end





        



   


