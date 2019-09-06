function [] = iteration()
global  MATERIALPARAMETERS A B H  nrelm connectivitynew profileA profileB nnod;
profileA= zeros(nrelm,4) ;
profileB= zeros(nrelm,4) ;
MATERIALPARAMETERSA=zeros(nnod,1);
MATERIALPARAMETERSB=zeros(nnod,1);
MATERIALPARAMETERSA=MATERIALPARAMETERS( (round(4.01*nnod):5*nnod),H );
MATERIALPARAMETERSB=MATERIALPARAMETERS( (round(5.01*nnod):6*nnod),H );
    for iel=1:nrelm  
        nd=connectivitynew(iel,1:end);         % extract connected node for (iel)-th element
        profileA(iel,:) = ((MATERIALPARAMETERSA(nd))) ;
        profileB(iel,:) = ((MATERIALPARAMETERSB(nd)))  ;
        % extract component value of the node 
    end


A(H+1,:)=  ( profileA(:,1)+ profileA(:,2)+profileA(:,3)+ profileA(:,4) )/4  ; 
B(H+1,:)=  ( profileB(:,1)+ profileB(:,2)+profileB(:,3)+ profileB(:,4) )/4  ;

end