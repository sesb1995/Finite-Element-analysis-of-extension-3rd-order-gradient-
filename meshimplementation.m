function [  ] = meshimplementation()
clc ;
clear all ;

tic
%%%%%%%%%%%%%coefficients

global Ex Ey connectivity connectivitynew nnod nrelm ;

%%%%%%%%%%%organizing elements and nodes

load Elements.txt
load Nodes.txt

nrelm=size(Elements,1);
nnod=size(Nodes,1);
Ex=zeros(nrelm,4);
Ey=zeros(nrelm,4);
connectivity(:,1)=sort(1:nrelm);
% A=zeros(100,1);
% B=zeros(100,1);

for i=1:nrelm
 
    connectivity(i,1:4)=Elements(i,6:9); %%%%creating connectivity matrix 
    Ex(i,:)=[Nodes(connectivity(i,1),2) Nodes(connectivity(i,2),2) Nodes(connectivity(i,3),2) Nodes(connectivity(i,4),2) ];
    Ey(i,:)=[Nodes(connectivity(i,1),3) Nodes(connectivity(i,2),3) Nodes(connectivity(i,3),3) Nodes(connectivity(i,4),3) ];
    
end

%%%%%%%%%%%%new connectivity

connectivitynew=zeros(nrelm,4) ;
connectivitynew(:,3)=connectivity(:,1) ;
connectivitynew(:,1)=connectivity(:,2) ;
connectivitynew(:,2)=connectivity(:,3) ;
connectivitynew(:,4)=connectivity(:,4) ;
for i=1:nrelm
 
    connectivity(i,1:4)=Elements(i,6:9); %%%%creating connectivity matrix 
    Ex(i,:)=[Nodes(connectivitynew(i,1),2) Nodes(connectivitynew(i,2),2) Nodes(connectivitynew(i,3),2) Nodes(connectivitynew(i,4),2) ];
    Ey(i,:)=[Nodes(connectivitynew(i,1),3) Nodes(connectivitynew(i,2),3) Nodes(connectivitynew(i,3),3) Nodes(connectivitynew(i,4),3) ];
    
end

%  figure(1)
%   for i=1:nrelm     %%%length(Edof(:,1))
%      plot(Ex(i,[1 3 4 2 1]) ,Ey(i,[1 3 4 2 1]) );
%      hold on
%   end

  end
