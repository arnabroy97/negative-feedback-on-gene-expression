!!! another gillespie trial for negetive feedback via miRNA-----------------------------------
!!! this is just one gillepie history ------------------------------- 
program trial1
implicit none 

real, dimension(13) :: c
real, dimension(11) :: a
integer ::react_count,mu,nu,seed
real :: m,p,s,com,g_c,T,T2,TINT,TPRINT,summ,R1,R2,R2A0,A0,tao,gamma0,alpha,steady_t
real :: mon,moff,pon,poff,son,soff,comon,comoff
real :: run_t1,run_t2

call cpu_time(run_t1)
!!!! define initial values and reaction constants -------------------------------------------

T2=1050000; TINT=0.5; steady_t = 50000
T=0; mon=0; moff=0; m=0; pon=0; poff=0; p=0
son=0; soff=0; s=0; comon=0; comoff=0; com=0; g_c=0;

gamma0 = 0.0004; alpha = 0.95                           !!! notice here that alpha is real so for alpha = 0.5 write alpha = (1.0/2) 

react_count = 0

 c(1)= 0.0007                                  !!!c(1)=r_m / k_r               (to be varied)
 c(2)= 0.0004                               !!!c(2)=gamma_m / g_r
 c(3)= 0.4                                  !!!c(3)=r_p / k_p
 c(4)= 0.0005                               !!!c(4)=gamma_p /g_p
 c(5)= 0                                 !!!c(5)=k_act              (to be varied)
 c(6)= 0.01                                 !!!c(6)=k_deact
 c(7)= 0.0005                                 !!!c(7)=r_s_0 / k_s_0
 c(8)= 0.5                                      !!!c(8)=r_s  / k_s
 c(9)= 0.0003                                     !!!c(9)=gamma_s / g_s
 c(10)= 0.1                                         !!!c(10)=k_+
 c(11)= 0.1                                    !!!c(11)=k_- 
 c(12)= alpha*gamma0                              !!!c(12)=alpha*gamma_0
 c(13)= (1-alpha)*gamma0                          !!!c(13)=(1-alpha)*gamma_0

!beta = c(13)/c(12)

!!!! assign value of T to TPRINT -------------------------------------------------------------------------------------------------

TPRINT=T

!!! open a data file and write the values--------------------------------------------------------------------------------------------
 
open(2, File= '4th_july_state_values_stoi=1_1_hist_k_r=0.0007_k_a=0_k_d=0.01_k_s_0=0.0005_k_s=0.5_k+=_k-=0.1 &
alpha=0.95_T=[50K,10.5L]_srand_parameter.txt',status='unknown', position='append')

open(3, File= '4th_july_state_values_stoi=1_1_hist_k_r=0.0007_k_a=0_k_d=0.01_k_s_0=0.0005_k_s=0.5_k+=_k-=0.1 &
alpha=0.95_T=[50K,10.5L]_srand_run5.dat',status='unknown', position='append')

write(2,5) c(1),c(2),c(3),c(4),c(5),c(6),c(7),c(8),c(9),c(10),c(11),c(12),c(13),gamma0,alpha,T2,TINT
5 format("##","k_r=",1X,F12.7,1X,"; g_r=",1X,F12.7,1X,"; k_p=",1X,F12.7,1X,"; g_p=",1X,F12.7,1X, &
"k_act=",1X,F15.9,1X,"; k_deact=",1X,F15.9,1X,"; k_s_0=",1X,F12.7,1X,"; k_s=",1X,F12.7,1X, &
"; gamma_s=",1X,F12.7,1X,"k_+ =",1X,F12.7,1X,"; k_- =",1X,F12.7,1X,"; alpha*gamma0 =",1X,F15.9,1X, &
"; (1-alpha)*gamma0 =",1X,F15.9,"; gamma0 =",1X,F12.7,1X,"; alpha =",1X,F12.7,1X,'## T2 = ',1X,F12.4,1X, &
'TINT = ',1X,F12.7,1X, 'Stoichiometry = 1','run = 5')
flush(2)
!write(3,6) c(5),c(6),c(7),c(8),c(9)
!6 format("##","k_act=",1X,F20.15,1X,"; k_deact=",1X,F20.15,1X,"; k_s_0=",1X,F12.9,1X,"; k_s=",1X,F12.9, &
!  1X,"; gamma_s=",1X,F12.9)

!write(3,7) c(10),c(11),c(12),c(13),gamma0,alpha
!7 format("##","k_+ =",1X,F12.9,1X,"; k_- =",1X,F12.9,1X,"; alpha*gamma0 =",1X,F20.15,1X,"; (1-alpha)*gamma0 =",1X,F20.15, &
!"; gamma0 =",1X,F12.9,1X,"; alpha =",1X,F12.9)

!write(3,*) '## T2 = ',T2,'TINT = ',TINT, 'Stoichiometry = 1'

write(2,*)'## (1)TPRINT   (2)mon   (3)moff   (4)m   (5)pon   (6)poff   (7)p   (8)son   (9)soff   (10)s   &
		(11)comon   (12)comoff   (13)com    (14)gon    (15)reaction_count'
flush(2)
!!!! call random number generator -------------------------------------------------------------------

   call SYSTEM_CLOCK(COUNT=seed)          !!! SYSTEM_CLOCK is a inbuild fortran subroutine to call a random number generator 
                                          !!!    and      'count=seed'
                                          !!! gives the seed for the random number generator srand
   call srand(seed)


!!!! find values of a(i) --------------------------------------------------------------------------------------------------
10 a(1) = c(1)                               !!! a(1) = r_m
   a(2) = c(2)*m                             !!! a(2) = m*gamma_m
   a(3) = c(3)*m                             !!! a(3) = m*r_p
   a(4) = c(4)*p                             !!! a(4) = p*gamma_p

   if (g_c.EQ.0) then 
           a(5) = c(5)*p                    !!! a(5) = p*k_act
   else                                     !!! or
           a(5) = c(6)                      !!! a(6) = k_deact
           
   end if 

   if (g_c.EQ.0) then
           a(6) = c(7)                      !!! a(6) = r_s_0
   else                                     !!! or
           a(6) = c(8)                      !!! a(6) = r_s
           
   end if 
   
   a(7) = c(9)*s                            !!! a(7) = s*gamma_s
   a(8) = c(10)*m*s                         !!! a(8) = m*s*k_+
   a(9) = c(11)*com                         !!! a(9) = c*k_-
   a(10) = c(12)*com                        !!! a(10) = c*gamma_c
   a(11) = c(13)*com                        !!! a(11) = c*gamma_cat

 a0 = a(1)+a(2)+a(3)+a(4)+a(5)+a(6)+a(7)+a(8)+a(9)+a(10)+a(11)
 
!!!! calculate R1 & R2 by random number generator ---------------------------------------------------

20  R1 = rand(0)
    R2 = rand(0)
    
    if (R1.LT.1E-30.OR.R2.LT.1E-30) then
       go to 20
    end if

    !R2A0 = R2*a0
   
    !Write(2,*) 'R1=',R1,'R2=',R2
    
!!!! increment T by tao ------------------------------------------------------------------------
   
   tao = alog(1/R1)/a0
   !print*, 'tao=',tao

21 T = T + alog(1/R1)/a0
   !write(3,*) 'increment T=',T,'TPRINT=',TPRINT
   
    
22 if (T.LT.TPRINT) then
   go to 25
   end if

23 if (TPRINT.GE.steady_t) then

	write(3,9) TPRINT,mon,moff,m,pon,poff,p,son,soff,s,comon,comoff,com,g_c
9 	format(F20.10,1X,F30.10,1X,F30.10,1X,F30.10,1X,F30.10,1X,F30.10,1X,F30.10,1X,F30.10,1X,F30.10,1X,F30.10,1X, &
	F30.10,1X,F30.10,1X,F30.10,1X,F30.10)
	flush(3)

   end if

!write(4,9) TPRINT,m,p,s,com,g_c,react_count
!flush(4)
 

!!! increase TPRINT -----------------------------------------------------------------------------------------------------------------

 if (TPRINT.GE.T2) go to 46
  
  TPRINT=TPRINT + TINT
  !write(2,*) 'T=',T,'TPRINT=',TPRINT
  go to 22

!!! calculate mu--------------------------------------------------------------------------------------------------------------------- 
25 R2A0=R2*A0
   summ=0

26 do nu=1,11
      mu=nu
      summ=summ+ a(nu)
      if(summ.GE.R2A0) go to 30
   end do

!!!to find the reaction to happen and calculate step 3 in gillespie algorithm-----------------------------------------------------------

30 go to(31,32,33,34,35,36,37,38,39,40,41), mu

31 m=m+1                                             !!!------------------  gene -> m (r_m)  --------------------------
   react_count= react_count + 1
   if (g_c.EQ.0) then
        moff = m; poff = p; soff = s; comoff = com
        mon = 0; pon = 0; son = 0; comon = 0
   end if
   if (g_c.EQ.1) then
        mon = m; pon = p; son = s; comon = com
        moff = 0; poff = 0; soff = 0; comoff = 0
   end if
   !write(3,99) T,TPRINT,mon,moff,m,pon,poff,p,son,soff,s,comon,comoff,com,g_c,react_count
!99 format('every step->',1X,'T=',F12.3,1X,'Tprint=',F12.3,1X,'mon=',F12.4,1X,'moff=',F12.4,1X,'m=',F12.4,1X, &
          !'pon=',F12.4,1X,'poff=',F12.4,1X,'p=',F12.4,1X,'son=',F12.4,1X,'soff=',F12.4,1X,'s=',F12.4,1X, &
          !'con=',F12.4,1X,'coff=',F12.4,1X,'c=',F12.4,1X,'gon=',F12.4,1X,'react_count=',I10)
   
   go to 45

32 m=m-1                                            !!!------------------  m -> phi (gamma_m)  --------------------------
   react_count= react_count + 1
   if (g_c.EQ.0) then
        moff = m; poff = p; soff = s; comoff = com
        mon = 0; pon = 0; son = 0; comon = 0
   end if
   if (g_c.EQ.1) then
        mon = m; pon = p; son = s; comon = com
        moff = 0; poff = 0; soff = 0; comoff = 0
   end if
   !write(3,99) T,TPRINT,mon,moff,m,pon,poff,p,son,soff,s,comon,comoff,com,g_c,react_count
!99 format (A,1X,F12.3,1X,F12.3,1X,F12.4,1X,F12.4,1X,F12.4,1X,F12.4,1X,F12.4,1X,F12.4,1X,& F12.4,1X,F12.4,1X,F12.4,1X,F12.4,1X,F12.4,1X,F12.4,1X,F12.4,1X,I10)
   go to 45

33 p=p+1                                            !!!------------------ m -> p (r_p)  --------------------------
   react_count= react_count + 1
   if (g_c.EQ.0) then
        moff = m; poff = p; soff = s; comoff = com
        mon = 0; pon = 0; son = 0; comon = 0
   end if
   if (g_c.EQ.1) then
        mon = m; pon = p; son = s; comon = com
        moff = 0; poff = 0; soff = 0; comoff = 0
   end if
  !write(3,99) T,TPRINT,mon,moff,m,pon,poff,p,son,soff,s,comon,comoff,com,g_c,react_count
!99 format (A,1X,F12.3,1X,F12.3,1X,F12.4,1X,F12.4,1X,F12.4,1X,F12.4,1X,F12.4,1X,F12.4,1X,& F12.4,1X,F12.4,1X,F12.4,1X,F12.4,1X,F12.4,1X,F12.4,1X,F12.4,1X,I10)
   go to 45

34 p=p-1                                             !!!------------------  p -> phi (gamma_p) -------------------
   react_count= react_count + 1
   if (g_c.EQ.0) then
        moff = m; poff = p; soff = s; comoff = com
        mon = 0; pon = 0; son = 0; comon = 0
   end if
   if (g_c.EQ.1) then
        mon = m; pon = p; son = s; comon = com
        moff = 0; poff = 0; soff = 0; comoff = 0
   end if
   !write(3,99) T,TPRINT,mon,moff,m,pon,poff,p,son,soff,s,comon,comoff,com,g_c,react_count
!99 format (A,1X,F12.3,1X,F12.3,1X,F12.4,1X,F12.4,1X,F12.4,1X,F12.4,1X,F12.4,1X,F12.4,1X,& F12.4,1X,F12.4,1X,F12.4,1X,F12.4,1X,F12.4,1X,F12.4,1X,F12.4,1X,I10)
   go to 45

35 if (g_c.EQ.0) then
   
   p=p-1
   g_c=g_c+1
   react_count=react_count + 1
   mon = m; pon = p; son = s; comon = com;            !!!-----  miRNA-gene + p -> miRNA-gene-protein complex (k_act) ---------------
   moff = 0; poff = 0; soff = 0; comoff = 0 
   
   !write(3,99) T,TPRINT,mon,moff,m,pon,poff,p,son,soff,s,comon,comoff,com,g_c,react_count
!99 format (A,F12.3,1X,F12.3,1X,F12.4,1X,F12.4,1X,F12.4,1X,F12.4,1X,F12.4,1X,F12.4,1X,& F12.4,1X,F12.4,1X,F12.4,1X,F12.4,1X,F12.4,1X,F12.4,1X,F12.4,1X,I10)
   go to 45

   else
 
   p=p+1                                             !!!-----   miRNA-gene-protein-complex -> miRNA-gene + p (k_deact)---------
   g_c=g_c-1
   react_count = react_count + 1
   moff = m; poff = p; soff = s; comoff = com
   mon = 0; pon = 0; son = 0; comon = 0
   
   !write(3,99) T,TPRINT,mon,moff,m,pon,poff,p,son,soff,s,comon,comoff,com,g_c,react_count
!99 format (A,F12.3,1X,F12.3,1X,F12.4,1X,F12.4,1X,F12.4,1X,F12.4,1X,F12.4,1X,F12.4,1X,& F12.4,1X,F12.4,1X,F12.4,1X,F12.4,1X,F12.4,1X,F12.4,1X,F12.4,1X,I10)
   go to 45

   end if
 
36 if (g_c.EQ.0) then

   s=s+1                                             !!!--------  miRNA-gene -> s (r_s_0)  --------------------------
   react_count=react_count + 1
   moff = m; poff = p; soff = s; comoff = com
   mon = 0; pon = 0; son = 0; comon = 0
   !write(3,99) T,TPRINT,mon,moff,m,pon,poff,p,son,soff,s,comon,comoff,com,g_c,react_count
!99 format (A,F12.3,1X,F12.3,1X,F12.4,1X,F12.4,1X,F12.4,1X,F12.4,1X,F12.4,1X,F12.4,1X,& F12.4,1X,F12.4,1X,F12.4,1X,F12.4,1X,F12.4,1X,F12.4,1X,F12.4,1X,I10)
   go to 45

   else
 
   s=s+1                                             !!!--------  miRNA-gene-protein-complex -> s (r_s) -------
   react_count=react_count + 1
   mon = m; pon = p; son = s; comon = com
   moff = 0; poff = 0; soff = 0; comoff = 0
   !write(3,99) T,TPRINT,mon,moff,m,pon,poff,p,son,soff,s,comon,comoff,com,g_c,react_count
!99 format (A,F12.3,1X,F12.3,1X,F12.4,1X,F12.4,1X,F12.4,1X,F12.4,1X,F12.4,1X,F12.4,1X,& F12.4,1X,F12.4,1X,F12.4,1X,F12.4,1X,F12.4,1X,F12.4,1X,F12.4,1X,I10)
   go to 45 
 
   end if

37 s=s-1                                              !!!--------------  s -> phi (gamma_s)  --------------
   react_count=react_count +1
   if (g_c.EQ.0) then
        moff = m; poff = p; soff = s; comoff = com
        mon = 0; pon = 0; son = 0; comon = 0
   end if
   if (g_c.EQ.1) then
        mon = m; pon = p; son = s; comon = com
        moff = 0; poff = 0; soff = 0; comoff = 0
   end if
   !write(3,99) T,TPRINT,mon,moff,m,pon,poff,p,son,soff,s,comon,comoff,com,g_c,react_count
!99 format (A,F12.3,1X,F12.3,1X,F12.4,1X,F12.4,1X,F12.4,1X,F12.4,1X,F12.4,1X,F12.4,1X,& F12.4,1X,F12.4,1X,F12.4,1X,F12.4,1X,F12.4,1X,F12.4,1X,F12.4,1X,I10)
   go to 45

38 m=m-1
   s=s-1                                              !!!--------------  m + s -> c (k_+)  --------------
   com=com+1                                              
   react_count=react_count +1
   if (g_c.EQ.0) then
        moff = m; poff = p; soff = s; comoff = com
        mon = 0; pon = 0; son = 0; comon = 0
   end if
   if (g_c.EQ.1) then
        mon = m; pon = p; son = s; comon = com
        moff = 0; poff = 0; soff = 0; comoff = 0
   end if
   !write(3,99) T,TPRINT,mon,moff,m,pon,poff,p,son,soff,s,comon,comoff,com,g_c,react_count
!99 format (A,F12.3,1X,F12.3,1X,F12.4,1X,F12.4,1X,F12.4,1X,F12.4,1X,F12.4,1X,F12.4,1X,& F12.4,1X,F12.4,1X,F12.4,1X,F12.4,1X,F12.4,1X,F12.4,1X,F12.4,1X,I10)
   go to 45

39 m=m+1                                              !!!--------------  c -> m + s (k_-)  --------------
   s=s+1
   com=com-1
   react_count=react_count +1
   if (g_c.EQ.0) then
        moff = m; poff = p; soff = s; comoff = com
        mon = 0; pon = 0; son = 0; comon = 0
   end if
   if (g_c.EQ.1) then
        mon = m; pon = p; son = s; comon = com
        moff = 0; poff = 0; soff = 0; comoff = 0
   end if
   !write(3,99) T,TPRINT,mon,moff,m,pon,poff,p,son,soff,s,comon,comoff,com,g_c,react_count
!99 format (A,F12.3,1X,F12.3,1X,F12.4,1X,F12.4,1X,F12.4,1X,F12.4,1X,F12.4,1X,F12.4,1X,& F12.4,1X,F12.4,1X,F12.4,1X,F12.4,1X,F12.4,1X,F12.4,1X,F12.4,1X,I10)
   go to 45

40 com=com-1                                              !!!--------------  c -> phi (gamma_c)  --------------
   react_count=react_count +1
   if (g_c.EQ.0) then
        moff = m; poff = p; soff = s; comoff = com
        mon = 0; pon = 0; son = 0; comon = 0
   end if
   if (g_c.EQ.1) then
        mon = m; pon = p; son = s; comon = com
        moff = 0; poff = 0; soff = 0; comoff = 0
   end if
   !write(3,99) T,TPRINT,mon,moff,m,pon,poff,p,son,soff,s,comon,comoff,com,g_c,react_count
!99 format (A,F12.3,1X,F12.3,1X,F12.4,1X,F12.4,1X,F12.4,1X,F12.4,1X,F12.4,1X,F12.4,1X,& F12.4,1X,F12.4,1X,F12.4,1X,F12.4,1X,F12.4,1X,F12.4,1X,F12.4,1X,I10)
   go to 45

41 com=com-1                                              !!!--------------  c -> s (gamma_cat)  --------------
   s=s+1
   react_count=react_count +1
   if (g_c.EQ.0) then
        moff = m; poff = p; soff = s; comoff = com
        mon = 0; pon = 0; son = 0; comon = 0
   end if
   if (g_c.EQ.1) then
        mon = m; pon = p; son = s; comon = com
        moff = 0; poff = 0; soff = 0; comoff = 0
   end if
   !write(3,99) T,TPRINT,mon,moff,m,pon,poff,p,son,soff,s,comon,comoff,com,g_c,react_count
!99 format (A,F12.3,1X,F12.3,1X,F12.4,1X,F12.4,1X,F12.4,1X,F12.4,1X,F12.4,1X,F12.4,1X,& F12.4,1X,F12.4,1X,F12.4,1X,F12.4,1X,F12.4,1X,F12.4,1X,F12.4,1X,I10)
   go to 45

45 if (T.LT.T2) then 
   go to 10                                  !!!!either this statement or the next statements 40 and 41 ---------------------
   end if

46 call cpu_time(run_t2)
write(2,*) '# ---------------------END OF SIMULATION run = 5 -------------------- time taken =',(run_t2 - run_t1)
write(2,*)

end program trial1


