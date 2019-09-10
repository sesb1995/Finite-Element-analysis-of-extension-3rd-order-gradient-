function [ K ] = PDE3( k,i,j ) 
global x y Ex Ey connectivity connectivitynew H C P M ;
C=150 ; M=5 ; P=150; 
syms x y ;
global nnod nrelm ;
%X1    1
%X2    2  
%Q     3
%R     4
%A     5  
%B     6
%S     7
%T     8
X1_C= derivativeX(i)*derivativeX(j) ;    
X2_C=0 ;    
Q_C=shapefunction(i)*shapefunction(j) ; 
R_C=0  ;
A_C=0  ;
B_C=0  ;
S_C=0  ;  
T_C=0  ;     




end