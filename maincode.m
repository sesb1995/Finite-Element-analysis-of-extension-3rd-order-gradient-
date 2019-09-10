clc ;
clear all ;

%%%%%%%%%%%%%coefficients
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global x y Ex Ey connectivity connectivitynew H C P M z ;
C=150 ; M=5 ; P=150; 
syms x y z;
global nnod nrelm ;

%%%%%%%%%%%organizing elements and nodes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load Elements.txt
load Nodes.txt

nrelm=size(Elements,1);
nnod=size(Nodes,1);
Ex=zeros(nrelm,16);
Ey=zeros(nrelm,16);
connectivity(:,1)=sort(1:nrelm);

for i=1:nrelm
 
    connectivity(i,1:16)=Elements(i,6:21); %%%%creating connectivity matrix 
    Ex(i,:)=[Nodes(connectivity(i,1),2) Nodes(connectivity(i,2),2) Nodes(connectivity(i,3),2) Nodes(connectivity(i,4),2) Nodes(connectivity(i,5),2) Nodes(connectivity(i,6),2) Nodes(connectivity(i,7),2) Nodes(connectivity(i,8),2) Nodes(connectivity(i,9),2) Nodes(connectivity(i,10),2) Nodes(connectivity(i,11),2) Nodes(connectivity(i,12),2) Nodes(connectivity(i,13),2) Nodes(connectivity(i,14),2) Nodes(connectivity(i,15),2) Nodes(connectivity(i,16),2)];
    Ey(i,:)=[Nodes(connectivity(i,1),3) Nodes(connectivity(i,2),3) Nodes(connectivity(i,3),3) Nodes(connectivity(i,4),3) Nodes(connectivity(i,5),3) Nodes(connectivity(i,6),3) Nodes(connectivity(i,7),3) Nodes(connectivity(i,8),3) Nodes(connectivity(i,9),3) Nodes(connectivity(i,10),3) Nodes(connectivity(i,11),3) Nodes(connectivity(i,12),3) Nodes(connectivity(i,13),3) Nodes(connectivity(i,14),3) Nodes(connectivity(i,15),3) Nodes(connectivity(i,16),3)];
    
end

%%%%%%%%%%%%new connectivity
%%%%%%%%%%%%%%%%%%%%%%%%%%new connectivity

connectivitynew=zeros(nrelm,9) ;
connectivitynew(:,1)=connectivity(:,2) ;
connectivitynew(:,2)=connectivity(:,7) ;
connectivitynew(:,3)=connectivity(:,8) ;
connectivitynew(:,4)=connectivity(:,3) ;
connectivitynew(:,5)=connectivity(:,6) ;
connectivitynew(:,6)=connectivity(:,14) ;
connectivitynew(:,7)=connectivity(:,15) ;
connectivitynew(:,8)=connectivity(:,9) ;
connectivitynew(:,9)=connectivity(:,5) ;
connectivitynew(:,10)=connectivity(:,13) ;
connectivitynew(:,11)=connectivity(:,16) ;
connectivitynew(:,12)=connectivity(:,10) ;
connectivitynew(:,13)=connectivity(:,1) ;
connectivitynew(:,14)=connectivity(:,12) ;
connectivitynew(:,15)=connectivity(:,11) ;
connectivitynew(:,16)=connectivity(:,4) ;

for i=1:nrelm
 
    connectivity(i,1:16)=Elements(i,6:21); %%%%creating connectivity matrix 
    Ex(i,:)=[Nodes(connectivitynew(i,1),2) Nodes(connectivitynew(i,2),2) Nodes(connectivitynew(i,3),2) Nodes(connectivitynew(i,4),2) Nodes(connectivitynew(i,5),2) Nodes(connectivitynew(i,6),2) Nodes(connectivitynew(i,7),2) Nodes(connectivitynew(i,8),2) Nodes(connectivitynew(i,9),2) Nodes(connectivitynew(i,10),2) Nodes(connectivitynew(i,11),2) Nodes(connectivitynew(i,12),2) Nodes(connectivitynew(i,13),2) Nodes(connectivitynew(i,14),2) Nodes(connectivitynew(i,15),2) Nodes(connectivitynew(i,16),2)];
    Ey(i,:)=[Nodes(connectivitynew(i,1),3) Nodes(connectivitynew(i,2),3) Nodes(connectivitynew(i,3),3) Nodes(connectivitynew(i,4),3) Nodes(connectivitynew(i,5),3) Nodes(connectivitynew(i,6),3) Nodes(connectivitynew(i,7),3) Nodes(connectivitynew(i,8),3) Nodes(connectivitynew(i,9),3) Nodes(connectivitynew(i,10),3) Nodes(connectivitynew(i,11),3) Nodes(connectivitynew(i,12),3) Nodes(connectivitynew(i,13),3) Nodes(connectivitynew(i,14),3) Nodes(connectivitynew(i,15),3) Nodes(connectivitynew(i,16),3)];
    
end

localmatrix_l = zeros(128,128) ;%%left
localmatrix_r = zeros(128,1) ;%%right
globalmatrix_l = zeros(nnod,nnod) ;%%left
globalmatrix_r = zeros(nnod,1) ;%%right

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%using picard's method
% for H=1:1:30
%%%%%%%%%%%%%%%%%creating local matrix (linear 4*4)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%  for n=1:1:nrelm%%%%swipe the elements
     for k=1:1:2
         for m=1:1:8          
            for i=1:1:16%%%swipe the rows of local matrix      
               for j=1:1:16%%%swipe the columns of local matrix
      
     localmatrix_l(i+16*(k-1),j+16*(m-1))= intquad(PDE8(k,m,i,j)) ;

               end 
           end
         end   
    end
% end
toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Drichlet boundary conditions
%%%%%%%%%%%%%%%%%%%%%%%%%

% for i=1:1:nnod     
%     %%%inlet
%     
%    if(i==7||i==6)
%       globalmatrix_l(i,:)=0 ;
%      globalmatrix_l(i,i)=1 ;
%       globalmatrix_r(i,1)=p0 ;
%    end
%   if(i==9||i==8)
%           globalmatrix_l(i,:)=0 ;  
%           globalmatrix_l(i,i)=1 ;
%           globalmatrix_r(i,1)=p0 ;     
%   end
%   
%   %%%outlet
%     if(i==4||i==3)
%       globalmatrix_l(i,:)=0 ;
%       globalmatrix_l(i,i)=1 ;
%       globalmatrix_r(i,1)=pinfinity ;  
%     end
%    if(i==23||i==22)
%           globalmatrix_l(i,:)=0 ;  
%           globalmatrix_l(i,i)=1 ;
%           globalmatrix_r(i,1)=pinfinity ;
%    end 
% 
% 
%  end
% 
% %%%%%%%%%%%%%%%%%%%%%%%%post processing the pressure 
% global Pressure ;
% global profile ;   
% Pressure = globalmatrix_l\globalmatrix_r ;
% profile = zeros(nrelm,4) ;
%    
%     for iel=1:nrelm  
%         nd=connectivity(iel,1:end);         % extract connected node for (iel)-th element
%         profile(iel,:) = (Pressure(nd)) ;         % extract component value of the node 
%     end
%     figure(1)
%     for i=1:nrelm    
%     fill(Ex(i,[1:end 1]),Ey(i,[1:end 1]),profile(i,[1:end 1]));
%     hold on ;   
%     end
%     colormap 'jet'
%     title('pressure')
%     
%     
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%velocity calculation
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure(3)
%         for i=1:nrelm    
%     fill(Ex(i,[1:end 1]),Ey(i,[1:end 1]),profile3(i,[1:end 1]));
%     hold on ;   
%     end
%     colormap 'jet'
%     title('Vx')
% end