      subroutine get_su2rotation(R, SU2)
      real*8,  parameter :: pi=4.*atan(1.d0)
      complex,intent(out) :: SU2(2,2)
      real*8,intent(in)    :: R(3,3)
      real*8   :: angle, axis(3)
      ! local variables
      complex  :: sig0(2,2), sig1(2,2), sig2(2,2), sig3(2,2)
      integer    :: i
      ! dgeev variables
      integer    :: info, iaxis
      real*8     :: det, R2xR3(3)
      real*8   :: mat(3,3), dvl(3,3), dvr(3,3), wi(3), dwork(12)
      real*8   :: wr(3), arg
      !----
      R2xR3(1)=R(2,2)*R(3,3)-R(3,2)*R(2,3)
      R2xR3(2)=R(3,2)*R(1,3)-R(1,2)*R(3,3)
      R2xR3(3)=R(1,2)*R(2,3)-R(2,2)*R(1,3)
      det = dot_product(R(:,1),R2xR3)
      mat = R*det
      !
      arg=((mat(1,1)+mat(2,2)+mat(3,3)-1.0)*0.5)
      if(arg>1.0)  arg=1.0
      if(arg<-1.0) arg=-1.0
      angle=acos(arg)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! check the rotational angle ! yfxu 08/05/2019
      if(abs(angle)<0.001) then 
          angle=0                                     
      else if(abs(abs(angle)-pi)<0.001) then
          angle=sign(1.d0,angle)*pi            
      else if(abs(abs(angle)-pi/2.0)<0.001) then
          angle=sign(1.d0,angle)*pi/2.0
      else if(abs(abs(angle)-2*pi/3.0)<0.001) then
          angle=sign(1.d0,angle)*2*pi/3.0
      else if(abs(abs(angle)-pi/3.0)<0.001) then
          angle=sign(1.d0,angle)*pi/3.0
      else if(abs(abs(angle)-pi/4.0)<0.001) then
          angle=sign(1.d0,angle)*pi/4.0
      else if(abs(abs(angle)-3*pi/4.0)<0.001) then
          angle=sign(1.d0,angle)*3*pi/4.0
      else if(abs(abs(angle)-pi/6.0)<0.001) then
          angle=sign(1.d0,angle)*pi/6.0
      else if(abs(abs(angle)-5*pi/6.0)<0.001) then
          angle=sign(1.d0,angle)*5*pi/6.0
      else  
      !   write(404,*)"The symmetry operation is not correct!", 180*angle/pi
      !   write(*,*)"The symmetry operation is not correct!", 180*angle/pi
          stop
      endif
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      if(abs(abs(angle) - pi) .lt. 1e-3) then
         ! angle is 180 deg => can't find the axis the
         ! easy way. Diagonalize rotation matrix and 
         ! pick the eigenvector corresponding to 
         ! unity eigenvalue.
          do iaxis=1,3
              axis=0.0
              axis(iaxis)=1.0
              axis=axis+matmul(mat,axis)
              !
              if (maxval(abs(axis))>1e-1) exit
          enddo
          if (maxval(abs(axis))<=1e-1) then
              write(0,*) 'can not find axis'
              write(0,*) axis
              write(0,*) mat
              stop
          endif
          axis=axis/sqrt(dot_product(axis,axis))


      !!!!!!!!!!!!!!!!! The arbitary of C2 axis! yfxu 17/05/2019
      ! use the convertion in msg.txt
          if (abs(axis(1))<1e-2) then
              if (abs(axis(2)) > 1e-2 .and. abs(axis(3))> 1e-2) then
                  if (axis(2)> 0) then
                      axis=-1.0*axis
                  endif
              endif
          endif
          if (abs(axis(2))<1e-2) then
              if (abs(axis(1)) > 1e-2 .and. abs(axis(3))> 1e-2) then
                  if (axis(1)< 0) then
                      axis=-1.0*axis
                  endif
              endif
          endif
          if (abs(axis(3))<1e-2) then
              if (abs(axis(1)) > 1e-2 .and. abs(axis(2))> 1e-2) then
                  if (axis(2)< 0) then
                      axis=-1.0*axis
                  endif
              endif
          endif
      !!!!!!!!!!!!!!!!! yfxu 17/05/2019

      else if(abs(angle) .gt. 1e-3) then
         ! standard case. See Altmann's book
         axis(1)=mat(3,2)-mat(2,3)
         axis(2)=mat(1,3)-mat(3,1)
         axis(3)=mat(2,1)-mat(1,2)
         axis=axis/sin(angle)/2.0
      else if(abs(angle) .lt. 1e-4) then
         axis = 0.0
         axis(1)=1.0
      end if
      
      sig0(:,:) = cmplx( 0.0, 0.0)
      sig1(:,:) = cmplx( 0.0, 0.0)
      sig2(:,:) = cmplx( 0.0, 0.0)
      sig3(:,:) = cmplx( 0.0, 0.0)
      sig0(1,1) = 1.0
      sig0(2,2) = 1.0
      sig1(1,2) = 1.0
      sig1(2,1) = 1.0
      sig2(1,2) = cmplx(0.0,-1.0)
      sig2(2,1) = cmplx(0.0, 1.0)
      sig3(1,1) = 1.0
      sig3(2,2) = -1.0
      SU2=0*SU2
      SU2=SU2+cos(angle/2)*sig0
      SU2=SU2-cmplx(0.0,1.0)*sin(angle/2)*(axis(1)*sig1+axis(2)*sig2)
      SU2=SU2-cmplx(0.0,1.0)*sin(angle/2)*axis(3)*sig3
      end
