%Convert Locations to Distance matrix
%Assume L is matrix of form NumFrams x 2
function d = Loc2Dist(L)
n = size(L);
d = zeros(n(1)); %Initialise d
for i = 1:n(1)
    for j = 1:n(1)
        d(i,j) = sqrt((L(i,1) - L(j,1))^2 + (L(i,2)-L(j,2))^2); %Compute Euclidean distance.
    end
end
end