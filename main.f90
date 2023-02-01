program main
 
!############### ソリトンの崩壊 速度の計算 渦の有無の判定　渦度計算 ##################
implicit none
!------------------ 変数の設定 ------------------------
real(8), parameter::a = 7.2d0*10**3,&!物理量
                    pi=4.d0*atan(1.d0),&
                    pi2=2.d0*pi !by = 1.d0, bx = 1.d0 　外部ポテンシャルなしの一様系
real(8),parameter:: dx= 0.2d0,dy=dx   !微小幅の定義
complex(8), parameter :: uso = (0.d0, 1.d0)!座標の分割に用いる変数
integer(8) :: i, j,i0=30,j0=30
integer(8), parameter :: Nx = 100, Ny = 100
!BOXサイズ
real(8),parameter::x0=0.5d0*dx*dble(Nx),&
                   y0=0.5d0*dy*dble(Ny)   !x0=30 ,y0=30
real(8) :: x(0:Nx+1), y(0:Ny+1) !座標
real(8) :: dt, dt1=0.01d0*dx**2, dt2  ! dt1<<(dx)**2 !時間発展に用いる変数
integer(4)::IR,jikoku(8)
real(8)::fre1,fre2,Randam,no=0.1d0**4
complex(8) :: f(-1:Nx+1,-1:Ny+1), fa(-1:Nx+1,-1:Ny+1), fb(-1:Nx+1,-1:Ny+1) !実時間発展に用いる波動関数
complex(8) :: k(-1:Nx+1,-1:Ny+1), ka(-1:Nx+1,-1:Ny+1), kb(-1:Nx+1,-1:Ny+1) !虚時間発展に用いる波動関数
real(8) :: sum !波動関数の規格化に用いる変数
!ポテンシャルに用いる変数
!real(8) :: Vtrap(0:Nx,0:Ny)
integer::Nite !doループから抜け出すための変巣の定義
real(8) :: iteration
real(8), parameter :: seido = 0.1d0**12.d0 !クランク-ニコルソン法においてdoループを抜け出すために用いる変数
real(8) :: iteration2
real(8), parameter :: seido2 = 0.1d0**8.d0 !虚時間発展から抜け出すために用いる変数
integer(8) :: Nstep, Dstep, Estep
integer(8), parameter :: Fstep = 1*(10)**6 !ステップ数
real(8) :: Etot, Ekin, Etrap, Eint, mu, mua, Etotb,px,py !物理量
complex(8)::uzudoz(1:Nx+1,1:Ny+1) !渦度計算に用いる変数
real(8)::uzudo 
real(8)::histx(-200:200),histy(-200:200),hist(0:400),& !速度分布のために用いる変数
         vx(1:Nx+1,1:Ny+1),vy(1:Nx+1,1:Ny+1),v(1:Nx+1,1:Ny+1), &
         bin=0.005d0,bin1=0.005
integer(8)::m,n,l
real(8)::mark(0:Nx+1,0:Ny+1)=0.d0 !渦の有無に用いる変数
real(8),parameter::vorte(-4:4)=(/0.d0,0.d0,0.d0,1.d0,0.d0,1.d0,0.d0,0.d0,0.d0/)
integer(8)::vorz
real(8)::vortpo(1:Nx+1,1:Ny+1) !渦の個数のカウント
integer::vortnu
character(len=40) filename
!-------------------- 設定終了 ---------------------
!
!---------------- 座標系の設定 --------------------------------
!座標
do i = 0, Nx+1
   x(i) = -x0 + dx*dble(i)
end do
do j = 0, Ny+1
   y(j) = -y0 + dy*dble(j)
end do
!---------------- 設定終了 ---------------------------
!
!---------------- ポテンシャルの設定 -------------------
!外部ポテンシャルなしの系
!do j = 0, Ny
   !do i = 0, Nx
      !  Vtrap(i,j) = bx*( x(i) )**2.d0 + by*( y(j) )**2.d0
   !end do
!end do
!################### 虚時間発展 ################################
!##################### 初期状態の設定　#########################
!
!##################### 初期状態の設定　#########################
!do i=0,Nx+1
  !do j=0,Ny+1
   ! k(i,j)=1.d0
  !enddo
!enddo
   
     
do i=0,Nx+1
   do j=0,Ny+1
      k(i,j)=dtanh(x(i)+x(i0))*dtanh(x(i)-x(i0))
   enddo
enddo
    
!規格化
sum = 0.d0
!
do j = 1, Ny
   do i = 1, Nx
    sum = sum + cdabs(k(i,j))**2.d0
   end do 
end do

sum = sum*dx*dy

do j = 0, Ny+1
   do i = 0, Nx+1
      k(i,j) = k(i,j)/dsqrt(sum)
   end do 
end do
!################# 境界条件　#################################
!周期境界条件
k(Nx+1,:) = k(1,:)
k(:,Ny+1) = k(:,1) 
k(0,:) = k(Nx,:)
k(:,0) = k(:,Ny) 
      
    
!------------------ 虚時間発展開始 -----------------------
Nstep = 0
Estep = 100
Etotb = 0.d0
dt2= 4.d0*(0.1d0**4.d0)
do
   if(i==i0)    then
      do j=1,Ny
         do i=1,Nx
            k(i,j)=0.d0
         enddo
      enddo
      
   else
      mu = 0.d0
      do j = 1, Ny
         do i = 1, Nx
            mu = mu + dconjg(k(i,j))*( -(  (k(i+1,j)+k(i-1,j)-2.0d0*k(i,j))/(dx**2.0d0)&
                    + (k(i,j+1)+k(i,j-1)-2.0d0*k(i,j))/(dy**2.0d0)  ) &
                    + a*(cdabs(k(i,j))**2.0d0)*k(i,j))
         end do
      end do
      mu = mu*dx*dy
      do j = 1, Ny
         do i = 1, Nx
            kb(i,j) = k(i,j) + (0.50d0*dt2/(-1.d0))*( -(  (k(i+1,j)+k(i-1,j)-2.0d0*k(i,j))/(dx**2.0d0)&
                             + (k(i,j+1)+k(i,j-1)-2.0d0*k(i,j))/(dy**2.0d0)  ) &
                             +  a*(cdabs(k(i,j))**2.0d0)*k(i,j)-mu*k(i,j)  )
         end do
      end do
   endif
   
   do  !iteration
      ka(:,:) = k(:,:)
      if(i==i0.or.i==i0+50) then 
         do j=1,Ny
            do i=1,Nx
         
            k(i,j)=0.d0
         
            enddo
         enddo 
         
      
     else
      mua = 0.d0
         do j = 1, Nx
            do i = 1, Ny
               mua = mua + dconjg(ka(i,j))*( -(  (ka(i+1,j)+ka(i-1,j)-2.0d0*ka(i,j))/(dx**2.0d0) &
                   + (ka(i,j+1)+ka(i,j-1)-2.0d0*ka(i,j))/(dy**2.0d0)  ) &
                  +  a*(cdabs(ka(i,j))**2.0d0)*ka(i,j))
            end do
         end do
         mua = mua*dx*dy
         do j = 1, Ny
            do i = 1, Nx
               k(i,j) = kb(i,j) + (0.50d0*dt2/(-1.d0))*( -(  (ka(i+1,j)+ka(i-1,j)-2.0d0*ka(i,j))/(dx**2.0d0) &
                  + (ka(i,j+1)+ka(i,j-1)-2.0d0*ka(i,j))/(dy**2.0d0)  ) &
                  + a*(cdabs(ka(i,j))**2.0d0)*ka(i,j) -mua*ka(i,j) )
            end do
         end do
      endif
      !規格化
      sum = 0.d0
      do j = 1, Ny
         do i = 1, Nx
            sum = sum + cdabs(k(i,j))**2.d0
         end do 
      end do
      sum = sum*dx*dy
      do j = 0, Ny+1
         do i = 0, Nx+1
            k(i,j) = k(i,j)/dsqrt(sum)
         end do 
      end do
     !周期境界条件
        k(Nx+1,:) = k(1,:)
        k(:,Ny+1) = k(:,1) 
        k(0,:) = k(Nx,:)
        k(:,0) = k(:,Ny)
      !doループを抜け出すための条件
      iteration = ddim(cdabs(k(Nx/4,Ny/4)-ka(Nx/4,Ny/4)), seido*cdabs(k(Nx/4,Ny/4)))
      if(iteration > 0.d0) cycle
         exit
      end do                              !iteraiton終了
      Nstep = Nstep + 1
   !............... 物理量の出力...................
   !運動エネルギー
      Ekin = 0.d0
      do j = 1, Ny
         do i = 1, Nx
            Ekin = Ekin + (dconjg(k(i,j)))*( -(  (k(i+1,j)+k(i-1,j)-2*k(i,j))/(dx**2) + (k(i,j+1)+k(i,j-1)-2*k(i,j))/(dy**2)  )  )
         end do
      end do 
      Ekin = Ekin*dx*dy
   !相互作用エネルギー
      Eint = 0.d0
      do j = 1, Ny
         do i = 1, Nx
            Eint = Eint + 0.5d0*a*(cdabs(k(i,j))**4.d0)
         end do
      end do
      Eint = Eint*dx*dy
   !全エネルギー
      Etot = Ekin + Eint
   !出力
      open(40,file="im-total_2Dsoliton-en.txt")
      write(40,*) Nstep, Etot
   !虚時間発展を終わらせるための条件
      iteration2 = ddim( dabs(Etot-Etotb), seido2*dabs(Etotb) )
      Etotb = Etot
      if(iteration2 > 0.d0) cycle
         exit
      enddo                               
      close(40)
   !###################### 虚時間発展終了 ############################
   !################### AVS用のデータの出力 ##########################
   open(55,file='dens.bdat',form='unformatted')
   open(56,file='phas.bdat',form='unformatted')
   !###################### 実時間発展 ################################
   call date_and_time(values=jikoku)
   IR=jikoku(5)+jikoku(6)+jikoku(7)
   !初期状態の設定
   do i=1,Nx
      do j=1,Ny
         call URAND1(1,Randam,IR)
         fre1=Randam
         call URAND1(1,Randam,IR)
         fre2=Randam
         f(i,j)=K(i,j)+no*fre1*cdexp(uso*2.d0*pi*fre2)
      enddo
   enddo
   Nstep = 0
   Estep = 10
   Dstep = 100
   dt = 0.d0
   do !時間発展   
      do j = 1, Ny
         do i = 1, Nx
            fb(i,j) = f(i,j) + (0.50d0*dt/(uso))*( -(  (f(i+1,j)+f(i-1,j)-2.0d0*f(i,j))/(dx**2.0d0) &
                             + (f(i,j+1)+f(i,j-1)-2.0d0*f(i,j))/(dy**2.0d0)  )  &
                             + a*(cdabs(f(i,j))**2.0d0)*f(i,j) )                    !a=1/n n:バルクでの数密度
         end do
      end do
    
      Nite=0   
      do  !iterationを繰り返す
         Nite=Nite+1
         fa(:,:) = f(:,:)
         do j = 1, Ny
            do i = 1, Nx
               f(i,j) = fb(i,j) + (0.50d0*dt/(uso))*( -(  (fa(i+1,j)+fa(i-1,j)-2.0d0*fa(i,j))/(dx**2.0d0) &
                                + (fa(i,j+1)+fa(i,j-1)-2.0d0*fa(i,j))/(dy**2.0d0)  ) &
                                +  a*(cdabs(fa(i,j))**2.0d0)*fa(i,j) )
            end do
         end do
         !周期境界条件
         f(Nx+1,:) = f(1,:)
         f(:,Ny+1) = f(:,1) 
         f(0,:) = f(Nx,:)
         f(:,0) = f(:,Ny)     
   !***************収束判定***************
         if(Nite>6)exit
         cycle
      enddo   !iteation終了
   !########### AVS用のデータの書き込み #####################
   if(mod(Nstep,500)==0 ) then   
      !write(55) ((cdabs(f(i,j))**2,i=1,Nx),j=1,Ny)
      !write(56) ((dimag(cdlog(f(i,j))), i=1,Nx),j=1,Ny)
   
   !output of movie file
      write(filename,"(a,i5.5,a)") "data/dens",int(Nstep),".vtk"
      open(10,file=filename)
      write(10,"('# vtk DataFile Version 3.0')")
      write(10,"('test')")
      write(10,"('ASCII ')")
   
      write(10,"('DATASET STRUCTURED_GRID')")
      write(10,"('DIMENSIONS ',3(1x,i3))") Nx, Ny, 1
   
      write(10,"('POINTS ',i9,' float')") Nx*Ny*1
      do j=0,Ny-1
      do i=0,Nx-1
         write(10,"(3(f9.4,1x))") x(i), y(j), 0.d0
      end do
      end do
   
      write(10,"('POINT_DATA ',i9)") Nx*Ny*1
   
   !date input
      write(10,"('SCALARS dens float')")
      write(10,"('LOOKUP_TABLE default')")
      do j=0,Ny-1
      do i=0,Nx-1
         write(10,*)  cdabs(f(i,j))**2
      end do
      end do
   
      close(10)
   end if
   !............... 物理量の出力...................
   !運動量
   px = 0.d0
   py = 0.d0
   !運動エネルギー
      Ekin = 0.d0
   !トラップエネルギー
      Etrap = 0.d0
   !相互作用エネルギー 
      Eint = 0.d0
   !渦のマーク
      mark(:,:)=0.d0

   if(Estep == 2000) then !############
      do j = 1, Ny
         do i = 1, Nx
            px = px - uso*(dconjg(f(i,j)))*((f(i+1,j)-f(i-1,j))/(2.d0*dx)) 
            py = py - uso*(dconjg(f(i,j)))*((f(i,j+1)-f(i,j-1))/(2.d0*dy)) 
            Ekin = Ekin + (dconjg(f(i,j)))*( -(  (f(i+1,j)+f(i-1,j)-2*f(i,j))/(dx**2) + (f(i,j+1)+f(i,j-1)-2*f(i,j))/(dy**2)  )  )
            !Etrap = Etrap + (Vtrap(i,j))*(  (cdabs(f(i,j)))**2.d0    )
            Eint = Eint + 0.5d0*a*(cdabs(f(i,j))**4.d0)
         end do
      end do
      !運動量
      px = px*dx*dy
      py = py*dx*dy
      !運動エネルギー
      Ekin = Ekin*dx*dy
      !相互作用エネルギー
      Eint = Eint*dx*dy
      !全エネルギー
      Etot = Ekin + Eint 
      !エネルギーの出力
      open(42,file="kin-en_2Dsoliton.txt")
      write(42,*) Nstep, Ekin
      !open(43,file="trap-en_2Dsoliton.txt")    !外部ポテンシャルなしの一様系
      !write(43,*) Nstep, Etrap
      open(45,file="int-en_2Dsoliton.txt")    
      write(45,*) Nstep, Eint
      open(46,file="total-en_2Dsoliton.txt")
      write(46,*) Nstep, Etot
      Estep = 0
   end if !#############
   Estep = Estep + 1
   !運動量の出力
   open(47,file="momentum_2Dsoliton-x.txt")
   write(47,*) Nstep, px
   open(48,file="momentum_2Dsoliton-y.txt")
   write(48,*) Nstep, py
   !############ 渦の有無の判定 #################
   if(mod(Nstep,2000)==0) then !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      do j=1,Ny
         do i=1,Nx
            vorz=idint((dimag(cdlog(f(i+1,j)))-dimag(cdlog(f(i,j))))/pi2) &
               +idint((dimag(cdlog(f(i+1,j+1)))-dimag(cdlog(f(i+1,j))))/pi2) &
               +idint((dimag(cdlog(f(i,j+1)))-dimag(cdlog(f(i+1,j+1))))/pi2) &
               +idint((dimag(cdlog(f(i,j)))-dimag(cdlog(f(i,j+1))))/pi2) 
     
         if(vorz==1.or.vorz==-1) then   !並列化できない
            mark(i,j)=1.d0
            mark(i+1,j)=1.d0
            mark(i,j+1)=1.d0
            mark(i+1,j+1)=1.d0
         end if
           
         !mark(i,j)=vorte(vorz)
         !mark(i+1,j)=vorte(vorz)
         !mark(i,j+1)=vorte(vorz)
         !mark(i+1,j+1)=vorte(vorz)
         end do
      end do
      open(1,file='mark1.bdat',form='unformatted')
      write(1)((mark(i,j),i=1,Nx),j=1,Ny)
   end if
   !###################### 渦の有無の判定計算終了 #################################
   !###################### 渦の個数 ###############################################
   if (mod(Nstep,100)==0) then
      do j=1,Nx
         do i=1,Ny
         !渦の位置
         vortpo(i,j)=  (dimag(cdlog(f(i+1,j)*f(i,j)))/pi2) &
            +(dimag(cdlog(f(i+1,j+1)*f(i+1,j)))/pi2) &
            +(dimag(cdlog(f(i,j+1)*f(i+1,j+1)))/pi2) &
            +(dimag(cdlog(f(i,j)*f(i,j+1)))/pi2) 
            
         !渦のカウント
         vortnu=vortnu+idint(((dimag(cdlog(f(i+1,j)*f(i,j)))/pi2) &
                +(dimag(cdlog(f(i+1,j+1)*f(i+1,j)))/pi2) &
                +(dimag(cdlog(f(i,j+1)*f(i+1,j+1)))/pi2) &
                +(dimag(cdlog(f(i,j)*f(i,j+1)))/pi2) )**2)
         end do   
      end do
   end if
   open(20,file='mark2.bdat',form='unformatted')
   write(20)((vortpo(i,j),i=1,Nx),j=1,Ny)
   open(21,file="vortex2-number.txt")
   write(21,*) Nstep ,vortnu
   
   !###################### 渦度の計算 速度のローテーション#########################
   if(mod(Nstep,10)==0) then !&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      do j=1,Ny
         do i=1,Nx
            uzudoz(i,j)=dimag(dconjg(f(i,j))*f(i+1,j)-dconjg(f(i+1,j))*f(i,j)) &
                              /(cdabs(f(i,j))**2+cdabs(f(i+1,j))**2) &
                        +dimag(dconjg(f(i+1,j))*f(i+1,j+1)-dconjg(f(i+1,j+1))*f(i+1,j)) &
                              /(cdabs(f(i+1,j))**2+cdabs(f(i+1,j+1))**2)&
                        +dimag(dconjg(f(i+1,j+1))*f(i,j+1)-dconjg(f(i,j+1))*f(i+1,j+1)) &
                              /(cdabs(f(i,j+1))**2+cdabs(f(i+1,j+1))**2) &
                        +dimag(dconjg(f(i,j+1))*f(i,j)-dconjg(f(i,j))*f(i,j+1)) &
                              /(cdabs(f(i,j))**2+cdabs(f(i,j+1))**2)
         end do
      end do
      open(2,file='uzudo-distribution.bdat',form='unformatted')
      write(2)((cdabs(uzudoz(i,j))**2.d0,i=1,Nx),j=1,Ny)
   end if
   if (mod(Nstep,2000)==0) then
      open(2,file='uzu-distribution.bdat',form='unformatted')
      write(2)((uzudoz(i,j),i=1,Nx),j=1,Ny)
   end if
   if(mod(Nstep,10)==0) then
      !渦度の大きさの時間発展
      do j=1,Ny
         do i=1,Nx
            uzudo=uzudo+ (cdabs(uzudoz(i,j))**2.d0)
         end do
      end do
      uzudo = uzudo*dx*dy
      open(15,file="uzudo_2Dsoliton.txt")
      write(15,*) Nstep, uzudo
   end if
   !###################### 渦度の計算終了 ######################################
   !###################### 速度分布の計算 #########################################
   !別で計算
   !########################## 速度の計算終了 ####################################################### 
        
   !計算の出力
   if(Nstep==10**4) then
      open(10,file="mark1_2Dsoliton.txt")
      do j = 1, Ny
         do i = 1, Nx
            write(10,*) x(i), y(j), mark(i,j)
         end do 
      end do
   end if
   if(Nstep==1.2*10**5) then
   open(11,file="mark2_2Dsoliton.txt")
      do j = 1, Ny
         do i = 1, Nx
            write(11,*) x(i), y(j), mark(i,j)
         end do 
      end do
   end if
   if(Nstep==1.2*10**5) then
      open(11,file="uzudo_2Dsoliton.txt")
      do j = 1, Ny
         do i = 1, Nx
            write(11,*) x(i), y(j), cdabs(uzudoz(i,j))**2.d0
         end do 
      end do
   end if
   !########### ループの終了条件　##########################
   if(Nstep==Fstep) exit
      Nstep=Nstep+1
      dt=dt1   ! dt1<<(dx)**2
   enddo !時間張って終了
   end program main
   
!*******0-1ransuu*********************************************
SUBROUTINE URAND1(N, X, IR)
!************************************************************************
!* UNIFORM RANDOM NUMBER GENERATOR (MIXED CONGRUENTIAL METHOD)          *
!*     PORTABLE BUT SLOW.  THE PERIOD IS ONLY 1664501.                  *
!* PARAMETERS                                                 ;           *
!*   (1) N      (I) THE NUMBER OF RANDOM NUMBERS TO BE GENERATED        *
!*                  (INPUT)                                             *
!*   (2) X      (D) UNIFORM RANDOM NUMBERS (OUTPUT)                     *
!*   (3) IR     (I) THE INITIAL SEED  (INPUT)                           *
!*                  THE SEED FOR THE NEXT CALL (OUTPUT)                 *
!* COPYRIGHT: Y. OYANAGI, JUNE 30, 1989  V.1                            *
!************************************************************************
!*
DOUBLE PRECISION X(N), INVM
PARAMETER (M = 1664501, LAMBDA = 1229, MU = 351750)
PARAMETER (INVM = 1.0D0 / M)
!*PAREMETER CHECK
IF( N .LE. 0) THEN
   WRITE(6,*) '(SUBR.URAND1) PARAMETER ERROR. N = ', N
   WRITE(6,*) 'RETURN WITH NO FURTHER CALCULATION.'
   RETURN
END IF
IF( IR .LT. 0 .OR. IR .GE. M) THEN
   WRITE(6,*) '(SUBR.URAND1) WARNING. IR = ', IR
END IF
!*MAIN LOOP
DO 10 I = 1, N
   IR = MOD( LAMBDA * IR + MU, M)
   X(I) = IR * INVM
10 CONTINUE
RETURN
END
   