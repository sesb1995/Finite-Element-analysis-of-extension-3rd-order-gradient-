meshimplementation() ;
%%%%%%%%%%%%%%%%%%%%%%%%%
global  nnod nrelm ; %%%%%%%mesh implementation
global  MATERIALPARAMETERS A B H Ex Ey globalmatrix_l globalmatrix_r profileA profileB;
MATERIALPARAMETERS=zeros(8*nnod,30) ;
A=zeros(30,nrelm);
B=zeros(30,nrelm);
H=1 ;

%%%%%%%%%%%%%%%meshing the profile
 while(1)
  localmatrix() ;
  Assembleglobalmatrix();
  Boundaryconditions();
  matrixsolver()
  iteration()
  H=H+1;
  %%%%%%%%%%%%error analysis
  if(H>2)
  
  max=0 ;
  for i=1:1:nnod
  
  if(abs ( MATERIALPARAMETERS(i,H-1)-MATERIALPARAMETERS(i,H-2) ))
  max= abs ( (MATERIALPARAMETERS(i,H-1)-MATERIALPARAMETERS(i,H-2)) ) ; 
  counter=i ;
  end
  
  end
  %%%%%%%%%%%%%%%%%%

  disp((max/abs( MATERIALPARAMETERS(counter,H-1) ) )) ;
  %%%%%%%%%%%%%%%%%%
  if ( max/ abs( MATERIALPARAMETERS(counter,H-1) )<10^-3 )
      break ;
  end
  
  end
  %%%%%%%%%%%%%%%
 end
postprocessing()