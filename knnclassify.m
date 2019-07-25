function [classify, Error] = knnclassify(testingset, trainingset, truth_train, truth_test, K)
%Classifier using KNN technique

 n = 1; k = 1;
 [r,~] = size(trainingset);
 [t,~] = size(testingset);
 dist = zeros([t, r]);
 
 %Compute Euclidean distance for testing set
   for n = 1:t
        for k = 1:r
            dist(n,k) = sqrt(sum((trainingset(k,:) - testingset(n,:)).^2));
            k = k+1;
        end
    n = n + 1;
    k = 1;
   end
% 
% Sort Euclidean Distance and assign truth value
 i = 1; n =1;
 assign = zeros(t,4);
 for i = 1:t
     n = 1;
     A = sort(dist(i,:)); %take ith column of distance matrix and sort it by value
     B = A(1:K);          %find first K values
     while n<=K
     [R,C(i,n)] = ismember(B(n), dist(i,:));
     assign(i,n) = truth_train(C(i,n),1);
        if sum(assign(i,:))>= (K/2+0.5)
         assign(i,n+1) = 1;
        end
     n = n+1;
     end
     i = i+1;
 end
 final = assign(:,(K));
 classify = final;
 
 Error_raw = sum(abs(final - truth_test));
 Error = (t - Error_raw)/t;
end

