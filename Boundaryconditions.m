function [ ] = Boundaryconditions()

global Ex Ey connectivity connectivitynew nnod nrelm H;
global localmatrix_l localmatrix_r globalmatrix_l globalmatrix_r MATERIALPARAMETERS;

%%%%drichlet boundary conditions
for i=1:1:nnod


if(i==1||i==2||i==5||i==6||i==7||i==8||i==9||i==10||i==11||i==12)%%%%%x1 and x2 are fixed
     globalmatrix_l(i,:)=0 ;
     globalmatrix_l(i,i)=1 ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%fixed x2     (2)
    globalmatrix_l(nnod+i,:)=0 ;
    globalmatrix_l(nnod+i,nnod+i)=1 ;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     globalmatrix_l((3*nnod)+i,:)=0 ;
     globalmatrix_l((3*nnod)+i,(3*nnod)+i)=1 ;
     globalmatrix_r((3*nnod)+i,1)=0;

end

%%%%%%%%%%%%%%%bending on the right side
if(i==3||i==4||i==21||i==22||i==23||i==24||i==25||i==26||i==27||i==28)
%  
  globalmatrix_l((2*nnod)+i,:)=0 ;
  globalmatrix_l((2*nnod)+i,(2*nnod)+i)=1 ;
  globalmatrix_r((2*nnod)+i,1)=-0.1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%
  globalmatrix_l((3*nnod)+i,:)=0 ;
  globalmatrix_l((3*nnod)+i,(3*nnod)+i)=1 ;
  globalmatrix_r((3*nnod)+i,1)=-0.3;

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end

end

%%neuman boundaray conditions ( X1,11=4/15 )
