function [ ] = matrixsolver()

global Ex Ey connectivity connectivitynew nnod nrelm H;
global localmatrix_l localmatrix_r globalmatrix_l globalmatrix_r MATERIALPARAMETERS ;


answer = globalmatrix_l\globalmatrix_r ;
MATERIALPARAMETERS(:,H)=answer ;

end