function X = normailize_column(D) 

[ N1 ,~] = size(D) ; 

q = repmat(sum(D.^2).^0.5 , N1,1) ; 
X = D./q ; 