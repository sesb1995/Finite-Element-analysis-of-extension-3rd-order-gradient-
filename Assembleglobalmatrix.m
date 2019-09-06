function [ ] = Assembleglobalmatrix()

global nnod  ;
global globalmatrix_l globalmatrix_lHELP globalmatrix_r;

globalmatrix_r = zeros(8*nnod,1);
globalmatrix_l = zeros(8*nnod,8*nnod);
for k=1:1:8
    for m=1:1:8
         for i=1:1:nnod%%%swipe the rows of local matrix
            for j=1:1:nnod%%%swipe the columns of local matrix

  globalmatrix_l(i+nnod*(k-1),j+nnod*(m-1)) = globalmatrix_l(i+nnod*(k-1),j+nnod*(m-1)) + globalmatrix_lHELP(i,j,k,m);

            end
         end
       
     end 
end



end