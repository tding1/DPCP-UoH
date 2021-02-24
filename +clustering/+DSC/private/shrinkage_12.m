% Copyright @ Mostafa Rahmani, 2017

function x = shrinkage_12(q,gamma) 

[N1,~] = size(q) ;  
q1 = sum(q.^2).^0.5 ; 
q2 = repmat(q1,N1,1) ; 
G = q - (q./q2)*gamma ; clear q2 ;
q3 = (q1 >= gamma) ; clear q1 ; 
q4 = repmat(q3,N1,1) ; clear q3 ; 

x = G.*q4 ; 




