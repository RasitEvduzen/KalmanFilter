function [K] = Kernel(kernelSelect,X,gamma)
NoD=length(X);
K = zeros(NoD);

switch kernelSelect
    case 'rbf'
        for i=1:NoD
            for j=1:NoD
                K(i,j) = exp(-gamma*(X(i,1)-X(j,1))^2);
            end
        end
end
end