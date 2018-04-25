subroutine fbfmmatrix(sd,Us,Ud,V,k,mk,A)
    implicit none
    real(8),intent(in) :: sd,Us,Ud,V
    integer(4),intent(in) :: k,mk
    real(8),intent(out) :: A(mk,mk)
    real(8),parameter :: pi=3.1415926d0
    real(8) :: Q1,ki,kj,q,gi,cgi,cgV,cgU,cg,cg1
    integer(4) :: i,j,h
    Q1=k*2d0*pi/mk
    do i=1,mk
        ki=i*2d0*pi/mk/i
        do j=1,mk
            kj=j*2d0*pi/mk/j
            A(i,j)=-1d0*(Us*conjg(v1(ki*i-Q1))*v1(kj*j-Q1)*conjg(v2(kj*j))*v2(ki*i)&
                        +Ud*u(ki*i-Q1)*u(kj*j-Q1)*u(kj*j)*u(ki*i)&
                        +V*cjl(ki*i-kj*j)*(conjg(v1(ki*i-Q1))*v1(kj*j-Q1)*u(kj*j)*u(ki*i)&
                        +u(ki*i-Q1)*u(kj*j-Q1)*conjg(v2(kj*j))*v2(ki*i)))/mk
            A(j,i)=A(i,j)
        enddo
        gi=0d0
        cgi=0d0
        cgV=0d0
        cgU=0d0
        cg=0d0
        cg1=0d0
        do h=0,mk-1
            q=h*2d0*pi/mk
            gi=gi+Us*v1(q)*conjg(v1(q))*v2(ki*i)*conjg(v2(ki*i))+Ud*u(q)*u(q)*u(ki*i)*u(ki*i)
            cgi=cgi+v1(q)*conjg(v1(q))*u(ki*i)*u(ki*i)+u(q)*u(q)*v2(ki*i)*conjg(v2(ki*i))
            cgV=cgV+v1(q)*conjg(v1(q))*u(ki*i-Q1)*u(ki*i-Q1)
            cgU=cgU+u(q)*u(q)*v1(ki*i-Q1)*conjg(v1(ki*i-Q1))
            cg=cg+cjl(q-ki*i+Q1)*conjg(v1(q))*u(q)*v1(ki*i-Q1)*u(ki*i-Q1)
            cg1=cg1+cjl(q)*conjg(v2(ki*i))*v2(ki*i-q)*u(ki*i-q)*u(ki*i)
        enddo
        gi=gi/mk
        cgi=cgi/mk
        cgV=cgV/mk
        cgU=cgU/mk
        cg=cg/mk
        cg1=cg1/mk
        A(i,i)=gi+2d0*V*cgi-2d0*V*(cgV+cgU)+V*(cg+cg1)&
               -1d0*(Us*v2(ki*i)*conjg(v2(ki*i))*v1(ki*i-Q1)*conjg(v1(ki*i-Q1))&
               +Ud*u(ki*i)**2*u(ki*i-Q1)**2)/mk&
               -2d0*V*(u(ki*i)*u(ki*i)*v1(ki*i-Q1)*conjg(v1(ki*i-Q1))&
                        +v2(ki*i)*conjg(v2(ki*i))*u(ki*i-Q1)**2)/mk
    end do
    contains
        real(8) function cjl(q)
            real(8) :: q
            cjl=2d0*cos(q/2d0)
        end function

        complex(8) function v1(q)
            real(8) :: q
            v1=sd*cos(q/2d0)/abs(cos(q/2d0))/sqrt(sd**2+4d0*cos(q/2d0)**2)
        end function

        complex(8) function v2(q)
            real(8) :: q
            v2=sd*cos(q/2d0)/abs(cos(q/2d0))/sqrt(sd**2+4d0*cos(q/2d0)**2)
        end function

        real(8) function u(q)
            real(8) :: q
            u=2d0*abs(cos(q/2d0))/sqrt(sd**2+4d0*cos(q/2d0)**2)
        end function
end subroutine
