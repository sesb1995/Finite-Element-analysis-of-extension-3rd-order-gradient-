 function [] = postprocessing()
global  nnod nrelm MATERIALPARAMETERS connectivitynew H Ex Ey profileA profileB;
profileA= zeros(nrelm,4) ;
profileB= zeros(nrelm,4) ;
MATERIALPARAMETERSA=MATERIALPARAMETERS( (1:nnod),H-1 );
MATERIALPARAMETERSB=MATERIALPARAMETERS( (round(1.01*nnod):2*nnod),H-1 );    
for iel=1:nrelm  
        nd=connectivitynew(iel,1:end);         % extract connected node for (iel)-th element
        profileA(iel,:) = ((MATERIALPARAMETERSA(nd))) ;
        profileB(iel,:) = ((MATERIALPARAMETERSB(nd)));
        % extract component value of the node 
end  
%%%%%%%%%%%%
 figure(1)
  for i=1:nrelm   %%%length(Edof(:,1))
     plot((Ex(i,[1 3 4 2 1]))+abs(profileA(i,[1 3 4 2 1]))  ,(Ey(i,[1 3 4 2 1]))+abs(profileB(i,[1 3 4 2 1]))  );
     hold on
  end

  hold on
%     for i=1:nrelm     %%%length(Edof(:,1))
%       plot( (Ex(i,[1 3 4 2 1]))  ,(Ey(i,[1 3 4 2 1])) );
%       hold on
%    end

 end




