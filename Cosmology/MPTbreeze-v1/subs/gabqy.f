      SUBROUTINE GABQY(FCT,XL,XU,SUM,TOL,IER)                            
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
      DIMENSION X(10),W(10),S(128),SN(128),L(128),LN(128)   
C  ****  PRINTED OUTPUT OF PARTIAL RESULTS: SET IWR=1.                  
      DATA IWR/0/                                                       
C  ****  COEFFICIENTS FOR GAUSS 20-POINT INTEGRATION.                   
      DATA NP,NP2,NP4/10,20,40/                                         
C  ****  ABSCISAS.                                                      
      DATA X/7.6526521133497334D-02,2.2778585114164508D-01,             
     1       3.7370608871541956D-01,5.1086700195082710D-01,             
     2       6.3605368072651503D-01,7.4633190646015079D-01,             
     3       8.3911697182221882D-01,9.1223442825132591D-01,             
     4       9.6397192727791379D-01,9.9312859918509492D-01/             
C  ****  WEIGHTS.                                                       
      DATA W/1.5275338713072585D-01,1.4917298647260375D-01,             
     1       1.4209610931838205D-01,1.3168863844917663D-01,             
     2       1.1819453196151842D-01,1.0193011981724044D-01,             
     3       8.3276741576704749D-02,6.2672048334109064D-02,             
     4       4.0601429800386941D-02,1.7614007139152118D-02/             
C  ****  CORRECTED TOLERANCE.                                           
      CTOL=DMAX1(TOL,1.0D-13)                                           
      PTOL=0.01D0*CTOL                                                  
      H=XU-XL                                                           
C                                                                       
      IF(IWR.EQ.1) THEN                                                 
      WRITE(6,10)                                                       
   10 FORMAT(///5X,'GAUSS ADAPTIVE-BIPARTITION QUADRATURE')             
      WRITE(6,11) XL,XU,TOL                                             
   11 FORMAT(/5X,'XL = ',1PD15.8,', XU = ',D15.8,', TOL =',             
     1D8.1)                                                             
      ENDIF                                                             
      IER=0                                                             
C  ****  GAUSS INTEGRATION FROM XL TO XU.                               
      A=0.5D0*(XU-XL)                                                   
      B=0.5D0*(XL+XU)                                                   
      C=A*X(1)                                                          
      D=W(1)*(FCT(B+C)+FCT(B-C))                                        
      DO 1 I1=2,NP                                                      
      C=A*X(I1)                                                         
    1 D=D+W(I1)*(FCT(B+C)+FCT(B-C))                                     
      SUM=D*A                                                           
C  ****  ADAPTIVE BIPARTITION SCHEME.                                   
      ICALL=NP2                                                         
      LH=1                                                              
      S(1)=SUM                                                          
      L(1)=1                                                            
    2 HO=H                                                              
      H=0.5D0*H                                                         
      ASUM=SUM                                                          
      LHN=0                                                             
      DO 5 I=1,LH                                                       
      K=L(I)                                                            
      SI=S(I)                                                           
      XA=XL+(K-1)*HO                                                    
      XB=XA+H                                                           
      XC=XA+HO                                                          
      A=0.5D0*(XB-XA)                                                   
      B=0.5D0*(XB+XA)                                                   
      C=A*X(1)                                                          
      D=W(1)*(FCT(B+C)+FCT(B-C))                                        
      DO 3 I2=2,NP                                                      
      C=A*X(I2)                                                         
    3 D=D+W(I2)*(FCT(B+C)+FCT(B-C))                                     
      S1=D*A                                                            
      A=0.5D0*(XC-XB)                                                   
      B=0.5D0*(XC+XB)                                                   
      C=A*X(1)                                                          
      D=W(1)*(FCT(B+C)+FCT(B-C))                                        
      DO 4 I3=2,NP                                                      
      C=A*X(I3)                                                         
    4 D=D+W(I3)*(FCT(B+C)+FCT(B-C))                                     
      S2=D*A                                                            
      ICALL=ICALL+NP4                                                   
      S12=S1+S2                                                         
      SUM=SUM+S12-SI                                                    
      IF(DABS(S12-SI).LT.DMAX1(PTOL*DABS(S12),1.0D-32)) GO TO 5         
      LHN=LHN+2                                                         
      IF(LHN.GT.128.OR.ICALL.GT.9999) GO TO 8                           
      SN(LHN)=S2                                                        
      LN(LHN)=K+K                                                       
      SN(LHN-1)=S1                                                      
      LN(LHN-1)=LN(LHN)-1                                               
    5 CONTINUE                                                          
      ERR=DABS(SUM-ASUM)/DMAX1(DABS(SUM),1.0D-32)                       
      IF(IWR.EQ.1) WRITE(6,12) ICALL,SUM,ERR,LHN                        
   12 FORMAT(5X,'N =',I5,', SUM =',1PD19.12,', ERR =',D8.1,             
     1', LH =',I3)                                                      
      IF(ERR.GT.CTOL.AND.LHN.GT.0) GO TO 6                              
      IF(IWR.EQ.1) WRITE(6,13)                                          
   13 FORMAT(5X,'END OF GAUSS-BIPARTITION PROCEDURE'///)                
      RETURN                                                            
    6 LH=LHN                                                            
      DO 7 I=1,LH                                                       
      S(I)=SN(I)                                                        
    7 L(I)=LN(I)                                                        
      GO TO 2                                                           
C  ****  WARNING (LOW ACCURACY) MESSAGE.                                
    8 WRITE(6,14)                                                       
   14 FORMAT(/5X,'*** LOW ACCURACY IN SUBROUTINE GABQY.')                
      WRITE(6,11) XL,XU,TOL                                             
      WRITE(6,15) SUM,ERR                                               
   15 FORMAT(5X,'SUM =',1PD19.12,', ERR =',D8.1//)                      
      IER=1                                                             
      RETURN                                                            
      END                                                               
