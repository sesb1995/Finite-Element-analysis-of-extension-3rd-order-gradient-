function [  ] = localmatrix()
%%%%%%%coefficients

C=150; P=150;%%(cosserat theory)
%%picard method
%%%%

global Ex Ey connectivitynew nnod nrelm A B H matrix jacobian ;
global localmatrix_l globalmatrix_lHELP;
%%%%basic matrices

localmatrix_l = zeros(4,4) ;%%left
globalmatrix_lHELP = zeros(nnod,nnod,8,8) ;%%right

%%%%%%%%%%%%%%%%%%%%%%%quadrature points

say=zeros(4,1);%shape function
DsayX=zeros(4,1);
DsayY=zeros(4,1);
w=zeros(4,1) ;
p=zeros(4,1) ;
w(1)=.347854 ; w(2)=.652145 ; w(3)=.652145 ; w(4)=.347854 ;
p(1)=-.861136 ; p(2)=-.339981 ; p(3)=.339981 ; p(4)=.861136 ;

%%%%%%%%%%%%%%%%%%creating local matrix (linear 4*4)
  for n=1:1:nrelm                    %%nrelm%%%%swipe the elements
     for k=1:1:8
         for m=1:1:8          
            for i=1:1:4%%%swipe the rows of local matrix      
               for j=1:1:4%%%swipe the columns of local matrix
                 matrix=0;
                 ehsan=0 ;
                   for i1=1:1:4%%quadrature
                    for j1=1:1:4%quadrature
     X=p(i1);Y=p(j1);%assigning each x and y a quadrature point
     
     %%shape functions
     
     say(1) = (1/4)*(1-X)*(1-Y) ;
     say(2) = (1/4)*(1+X)*(1-Y) ;
     say(3) = (1/4)*(1-X)*(1+Y) ;
     say(4) = (1/4)*(1+X)*(1+Y) ;
     
     %%%%derivative of shape function
     
     DsayX(1)=Y/4 - 1/4 ;
     DsayX(2)=1/4 - Y/4 ;
     DsayX(3)=- Y/4 - 1/4 ;
     DsayX(4) =Y/4 + 1/4 ;
     DsayY(1)=X/4 - 1/4 ;
     DsayY(2)=- X/4 - 1/4 ;
     DsayY(3)=1/4 - X/4 ;
     DsayY(4)=X/4 + 1/4 ;
     
     %%%%jacobian & analyzing PDE(transformin integral in master element to real element)
     
     A1 = [DsayX(1) DsayX(2) DsayX(3) DsayX(4); DsayY(1) DsayY(2) DsayY(3) DsayY(4) ] ;                              
     B1 = [Ex(n,1) Ey(n,1) ; Ex(n,2) Ey(n,2) ; Ex(n,3) Ey(n,3) ; Ex(n,4) Ey(n,4) ] ;
     jacobian = A1*B1 ;
     invjacobian=inv(jacobian) ;
     DXi=( ( invjacobian(1,1)*DsayX(i))+( invjacobian(1,2)*DsayY(i)) )     ;
     DXj=( ( invjacobian(1,1)*DsayX(j))+( invjacobian(1,2)*DsayY(j)) )     ;
     DYi=( ( invjacobian(2,1)*DsayX(i))+( invjacobian(2,2)*DsayY(i)) )     ;
     DYj=( ( invjacobian(2,1)*DsayX(j))+( invjacobian(2,2)*DsayY(j)) )     ;
     
%%%8 PDE 
%%%%%%%%%%1st PDE

if(k==1)

if(m==1)    
ehsan=-1*DsayY(i)*DsayY(j)*DYi*DYj*det(jacobian) ;
elseif (m==2) 
  ehsan= ( -( A(H,n)*say(i)*DsayY(j)*DYj)+(B(H,n)*say(i)*DsayX(j)*DXj ) )*det(jacobian) ;
elseif (m==3) 
ehsan= ( (say(i)*say(j))+(C*DsayX(i)*DsayX(j)*DXi*DXj) )*det(jacobian) ;
elseif (m==7) 
ehsan= -P*DsayX(i)*DsayX(j)*DXi*DXj*det(jacobian) ;
else
    ehsan=0 ;
end

end
%%%%%%%%%%%%%%%2nd PDE

if(k==2)
   
if(m==1)
ehsan= ( ( A(H,n)*say(i)*DsayY(j)*DYj )-(B(H,n)*say(i)*DsayX(j)*DXj ) )*det(jacobian) ;
elseif(m==2)
ehsan=-1*DsayY(i)*DsayY(j)*DYi*DYj*det(jacobian)  ; 
elseif(m==4)
ehsan=( (say(i)*say(j))+(C*DsayX(i)*DsayX(j)*DXi*DXj ) )*det(jacobian) ;
elseif(m==8)
ehsan=-P*DsayX(i)*DsayX(j)*DXi*DXj*det(jacobian)  ; 
else
    ehsan=0 ;
end

end
%%%%%%%%%%%%%%3rd PDE

if(k==3)

if(m==1)   
ehsan= DsayX(i)*DsayX(j)*DXi*DXj*det(jacobian)  ;
elseif(m==3)
ehsan=say(i)*say(j)*det(jacobian) ; 
else
    ehsan=0 ;
end

end
%%%%%%%%%%%%4th PDE

if(k==4)

if(m==2)
ehsan=DsayX(i)*DsayX(j)*DXi*DXj*det(jacobian)  ; 
elseif(m==4)
ehsan=say(i)*say(j)*det(jacobian)   ; 
else
    ehsan=0 ;
end

end
%%%%%%%%%%%%5th PDE

if(k==5) 

if(m==3)    
ehsan=DsayX(i)*DsayX(j)*DXi*DXj*det(jacobian) ;
elseif(m==7)
ehsan=say(i)*say(j)*det(jacobian)  ;
else
    ehsan=0 ;
end

end
%%%%%%%%%%%%%%6th PDE

if(k==6)

if(m==4)    
ehsan=DsayX(i)*DsayX(j)*DXi*DXj*det(jacobian)   ; 
elseif(m==8)
ehsan=say(i)*say(j)*det(jacobian)  ;
else
    ehsan=0 ;
end

end
%%%%%%%%%%7th PDE

if(k==7)

if(m==1)    
ehsan=DsayY(i)*DsayY(j)*DYi*DYj*det(jacobian) ; 
elseif(m==3)
ehsan=-say(i)*say(j)*det(jacobian);
elseif(m==5)
ehsan=say(i)*say(j)*det(jacobian)   ;
elseif(m==7)
ehsan=-C*say(i)*say(j)*det(jacobian)  ; 
else
    ehsan=0;
end

end
%%%%%%%%%%8th PDE

if(k==8) 

if(m==2)    
ehsan=DsayY(i)*DsayY(j)*DYi*DYj*det(jacobian) ; 
elseif(m==4)
ehsan=-say(i)*say(j)*det(jacobian)  ;
elseif(m==6)
ehsan=say(i)*say(j)*det(jacobian)  ; 
elseif(m==8)
ehsan=-C*say(i)*say(j)*det(jacobian) ;
else
    ehsan=0;
end

end

%%%%%%%assigning values to local matrix
     matrix= matrix+ ( ehsan*(w(i1)*w(j1)) ) ;

                    end
                   end
 localmatrix_l(i,j)=matrix  ;
 globalmatrix_lHELP(connectivitynew(n,i),connectivitynew(n,j),k,m) = globalmatrix_lHELP(connectivitynew(n,i),connectivitynew(n,j),k,m)+ matrix ; 
                end
            end
         end
     end
  end
  
  
  
  
  
  
  
  
  
  
end












