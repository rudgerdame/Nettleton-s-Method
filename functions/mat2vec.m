function vec=mat2vec(mat,keep)
% Input: 
% mat - data matrix 
% keep - matrix of ones and zeros specifying which elements are to be kept
% Output:
% vec - column vector
%
% For example:
% X=rand(4,4); vec=mat2vec(X,X>.5);
%
% Modified by pjames@alum.mit.edu, 10/2/15

defval('keep',ones(size(mat)));
keep=logical(keep);
[~,s2]=size(mat);
vec=mat(keep(:,1),1);
for j=2:s2
    vec=vertcat(vec,mat(keep(:,j),j));
end



end