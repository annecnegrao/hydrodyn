program AkanMelinaLPI
!-----------------------------------------------------!
! Anne Caroline Negrao (2015)                         !
!                                                     !
! Modelo Hidrodinamico Unidimensional -               !
! Resolucao das equacoes completas de Saint Venant    !
! atraves do metodo de Preissman adaptado para        !
! escoamentos transcriticos.                          !
!                                                     !
! discretizacao em relacao a A (Akan)                 !
! resolucao do sistema Double Sweep (Melina)          !
! filtro LPI (Fread)                                  !
!                                                     !
! Aplicacao no corrego do Gregorio                    !
!                                                     !
! ATENCAO: Esse script exige ambiente Linux ou Mac!   !
!-----------------------------------------------------!
implicit none

integer:: i,j,k,n,it
integer:: jj,nn,nit,infof,infop,jinter
integer:: problema,m,seg1,seg2
real:: g,long,b,s,nM,dx,dt,teta,tol,Q0
real:: AA,PP,inclinacao
real:: dAdt,dQdx,dQdt,dQ2Adx,Adhdx,ASf,rB,rBN
real,allocatable,dimension(:,:):: h,y,A,dAdh,P,Q,Sf,dSfdQ,dSfdh,Fr,rC,rM
real,allocatable,dimension(:):: x,S0,z,dQ,dh,rterm,norma,yup,yinter,ydown
real,allocatable:: residuos(:),Frmax(:),Frmin(:),vmax(:),vmin(:)
real:: dPdh,dBdQ1,dBdh1,dBdQN,dBdhN,ENS
real,allocatable,dimension(:,:):: W
real,allocatable,dimension(:,:):: dCdQ,dCdh,dMdQ,dMdh
logical:: more
character*200:: arq,arq_evento,dummy,temp,CI,fundo,processo,dir
character*10:: evento
character*10,allocatable,dimension(:):: amd,hms

10 format(i3,4f8.2,2f8.1,f8.2,f8.1,3f14.7,f9.3)
11 format(f6.1,3f10.4)


write(*,*)"Leitura dos dados"

open(1,file="info_canal.txt",action="read")
read(1,*)g            ! aceleracao da gravidade
read(1,*)long         ! comprimento do canal
read(1,*)b            ! largura do fundo do canal
read(1,*)s            ! inclinacao daparede lateral
read(1,*)nM           ! rugosidade de Manning
read(1,*)teta         ! coef. temporal
read(1,*)nit          ! numero maximo de iteracoes
read(1,*)tol          ! tolerancia
read(1,*)infof        ! aproximacao para (1) fundo simples ou (2) fundo com degraus
read(1,*)infop        ! processo de (1) aquecimento ou (2) simulacao
read(1,*)m            ! termo m do filtro LPI
close(1)

! verificando a aproximacao de fundo
if(infof==1)then
    fundo = "fundo_simples"
else if(infof==2)then
    fundo = "fundo_degraus"
else
    write(*,*)"Codigo de aproximacao de fundo errado,"
    write(*,*)"verificar info_canal.txt ..."
    goto 33
end if

! verificando o procedimento
if(infop==1)then
    processo = "aquecimento"
else if(infop==2)then
    processo = "simulacao"
else
    write(*,*)"Codigo de processo errado,"
    write(*,*)"verificar info_canal.txt ..."
    goto 33
end if

write(*,'(5a)')"Fazendo ",trim(processo)," para o ",trim(fundo)," ..."


! Leitura das informacoes de fundo
write(arq,'(2a)')trim(fundo),".txt"
open(1,file=arq,action="read")
read(1,*)
jj = 0
do
    read(1,*,end=99)
    jj = jj+1
end do
99 rewind(1)

allocate(S0(jj),z(jj),x(jj),Frmax(jj),Frmin(jj),vmax(jj),vmin(jj))

read(1,*)
do j = 1,jj
    read(1,*)x(j),S0(j),z(j)
end do
close(1)

! secao intermediaria de monitoramento 
if(fundo=="fundo_simples")then
    jinter = 97
else if(fundo=="fundo_degraus")then
    jinter = 100
end if

! abrindo arquivos de saida
if(processo=="aquecimento")then
    write(arq,'(3a)')"./",trim(fundo),"/arq_eventos_cut.txt"
    open(8,file=arq,action="write")
else if(processo=="simulacao")then
    write(arq,'(3a)')"./",trim(fundo),"/ENS.txt"
    open(8,file=arq,action="write")
    write(8,*)"Coeficiente de Eficiencia de Nash-Sutcliffe p/ profundidade"
end if

! Simulando por evento
if(processo=="aquecimento")then
    open(13,file="arq_eventos.txt",action="read")
else if(processo=="simulacao")then
    write(arq,'(3a)')"./",trim(fundo),"/arq_eventos_cut.txt"
    open(13,file=arq,action="read")
end if
do
    if(processo=="aquecimento")then
        read(13,'(a)',end=97)arq_evento
        read(arq_evento,'(17x,a10)')evento
    else if(processo=="simulacao")then
        read(13,'(a,x,a)',end=97)evento,arq_evento
    end if
    write(*,*)arq_evento

    write(dir,'(5a)')"./",trim(fundo),"/",evento,"/"
    !write(*,*)dir
    
    if(processo=="aquecimento")then
        ! criando diretorios:
        write(arq,'(2a)')"mkdir ",trim(dir)
        call system(arq)
    end if

    ! tamanho do evento:
40  open(1,file=arq_evento,action="read")
    nn = 0
    do
        read(1,*,end=91)
        nn = nn+1
    end do
91  rewind(1)
    nn = nn-1
    write(*,*)"nn = ",nn
    !read(*,*)

    allocate(h(0:nn,jj),y(0:nn,jj),A(0:nn,jj),P(0:nn,jj),Q(0:nn,jj),yup(0:nn),yinter(0:nn),ydown(0:nn))
    allocate(Sf(0:nn,jj),dSfdQ(0:nn,jj),dSfdh(0:nn,jj),dAdh(0:nn,jj),rC(jj,nit),rM(jj,nit),Fr(0:nn,jj))
    allocate(residuos(2*jj),W(jj,2),dQ(jj),dh(jj),rterm(jj),norma(jj))
    allocate(dCdQ(jj,jj),dCdh(jj,jj),dMdQ(jj,jj),dMdh(jj,jj))
    allocate(amd(0:nn),hms(0:nn))

    ! zerando as matrizes:
    h = 0. ; y = 0. ; A = 0. ; P = 0. ; Q = 0.
    Sf = 0. ; dSfdQ = 0. ; dSfdh = 0.
    dCdQ = 0. ; dCdh = 0. ; dMdQ = 0. ; dMdh = 0.

    write(*,*)"Lendo o arquivo do evento..."
    do n = 0,nn
        read(1,*)amd(n),hms(n),yup(n),yinter(n),ydown(n)
    end do
    close(1)

    ! obtendo dt
    read(hms(0),'(6x,i2)')seg1
    read(hms(1),'(6x,i2)')seg2
    dt = real(abs(seg1-seg2))
    write(*,*)"dt = ",dt


    if(processo=="aquecimento")then
        ! calculo da vazao inicial como sendo a normal de montante
        AA = (b+yup(0)*s)*yup(0)
        PP = b+2.*yup(0)*sqrt(1.+s**2)
        Q0 = AA*(AA/PP)**(2./3.)*sqrt(S0(5))/nM

        write(arq,'(5a)')"./",trim(fundo),"/",evento,"/h_y_obs.txt"
        open(1,file=arq,action="write")
        write(1,10)int(x(1)),z(1)+yup(0),yup(0)
        write(1,10)int(x(jinter)),z(jinter)+yinter(0),yinter(0)
        write(1,10)int(x(jj)),z(jj)+ydown(0),ydown(0)
        close(1)

        ! condicao inicial linear
        write(arq,'(5a)')"./",trim(fundo),"/",evento,"/CI_linear.txt"
        open(1,file=arq,action="write")
        write(1,*)"Q0=",Q0
        inclinacao = ((z(1)+yup(0))-(z(jj)+ydown(0)))/(x(jj)-x(1))
        h(0,1) = z(1)+yup(0)
        write(1,11)x(1),S0(1),z(1),h(0,1)
        do j = 2,jj
            h(0,j) = h(0,j-1) - inclinacao*(x(j)-x(j-1))
            write(1,11)x(j),S0(j),z(j),h(0,j)
        end do
        close(1)

    else if(processo=="simulacao")then
        ! condicao inicial aquecida
        write(*,*)"Lendo arquivo de condicoes iniciais..."
        write(arq,'(2a)')trim(dir),"CI_aquecido.txt"
        open(1,file=arq,action="read")
        read(1,*)dummy,Q0
        do j = 1,jj
            read(1,*)dummy,dummy,dummy,h(0,j)
        end do
        close(1)
    end if


    !!!!!!!!!!!!!!!!!!!!!!
    !! Condicao Inicial !!
    !!!!!!!!!!!!!!!!!!!!!!
    y(0,:) = h(0,:)-z
    Q(0,:) = Q0


    !!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! Condicoes de Contorno !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!
    if(processo=="aquecimento")then
        yup = y(0,1)
        ydown = y(0,jj)
        y(:,1) = y(0,1)
        y(:,jj) = y(0,jj)
    else if(processo=="simulação")then
        y(:,1) = yup
        y(:,jj) = ydown
    end if
    h(:,1) = y(:,1)+z(1)
    h(:,jj) = y(:,jj)+z(jj)


    !!!!!!!!!!!!!!!
    !! SIMULACAO !!
    !!!!!!!!!!!!!!!
 
    DO n = 0,nn-1
        write(*,*)"Instante de tempo: n = ",n

        ! Propriedades do canal
        call canal_trapezoidal(g,jj,z,h(n,:),b,s,Q(n,:),nM,A(n,:),dAdh(n,:),P(n,:),&
                         & dPdh,Sf(n,:),dSfdQ(n,:),dSfdh(n,:),Fr(n,:))
!     write(*,'(a)')"  j   Q        z       y       h      A       dA/dh   P        dP/dh    Sf&
!                  &            dSf/dQ        dSf/dh       Fr"
!     do j = 1,jj
!         write(*,10)j,Q(n,j),z(j),y(n,j),h(n,j),A(n,j),dAdh(n,j),P(n,j),dPdh,Sf(n,j),dSfdQ(n,j),dSfdh(n,j),Fr(n,j)
!     end do


        ! Estimativa para o proximo instante (igualando ao instante atual)
        Q(n+1,:) = Q(n,:)
        h(n+1,:) = h(n,:)    
    
        call canal_trapezoidal(g,jj,z,h(n+1,:),b,s,Q(n+1,:),nM,A(n+1,:),dAdh(n+1,:),P(n+1,:),&
                         & dPdh,Sf(n+1,:),dSfdQ(n+1,:),dSfdh(n+1,:),Fr(n+1,:))
                         
                         
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !! ITERACAO - NEWTON GENERALIZADO !!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
        dQ = 0.
        dh = 0.
        do it = 1,nit
            write(*,*)"Iteracao: ",it
	
            ! termo redutor segundo FREAD (filtro LPI)
            do j = 1,jj
                if(Fr(n+1,j)>1.)then
                    rterm(j) = 0.
                else
                    rterm(j) = 1.-Fr(n+1,j)**m
                end if
                !write(*,*)rterm(j)
            end do
            !read(*,*)
            !rterm = 1. ! coso sem filtro
	
	    !!!!!!!!!!!!!!
            !! Residuos !!
	    !!!!!!!!!!!!!!
	
            ! contorno montante
            rB = h(n+1,1) - (z(1)+yup(n+1))
            do j = 1,jj-1
                dx = x(j+1)-x(j)
                ! Aproximações:
                dAdt = ((A(n+1,j+1)+A(n+1,j))-(A(n,j+1)+A(n,j)))/(2.*dt)
                dQdx = teta/dx*(Q(n+1,j+1)-Q(n+1,j)) + (1.-teta)/dx*(Q(n,j+1)-Q(n,j))
                dQdt = rterm(j)*((Q(n+1,j+1)+Q(n+1,j))-(Q(n,j+1)+Q(n,j)))/(2*dt)
                dQ2Adx = rterm(j)*teta/dx*(Q(n+1,j+1)**2/A(n+1,j+1)-Q(n+1,j)**2/A(n+1,j)) &
		   & + rterm(j)*(1.-teta)/dx*(Q(n,j+1)**2/A(n,j+1)-Q(n,j)**2/A(n,j))
                Adhdx = teta*(A(n+1,j+1)+A(n+1,j))/2.*(h(n+1,j+1)-h(n+1,j))/dx &
		  & + (1.-teta)*(A(n,j+1)+A(n,j))/2.*(h(n,j+1)-h(n,j))/dx
                ASf = teta*(A(n+1,j+1)+A(n+1,j))/2.*(Sf(n+1,j+1)+Sf(n+1,j))/2. &
		& +  (1.-teta)*(A(n,j+1)+A(n,j))/2.*(Sf(n,j+1)+Sf(n,j))/2.
                !write(*,*)dAdt,dQdx,dQdt,dQ2Adx,Adhdx,ASf
                rC(j,it) = dAdt + dQdx
                rM(j,it) = dQdt + dQ2Adx + g*Adhdx + g*ASf
                !write(*,*)rC(j,it),rM(j,it)
            end do
            ! contorno jusante
            rBN = h(n+1,jj) - (z(jj)+ydown(n+1))

            residuos(1) = -rB
            j = 0
            do k = 2,2*jj-1
                if(mod(k,2)==0)then !par
                    j = j+1
                    residuos(k) = -rC(j,it)
                else
                    residuos(k) = -rM(j,it)
                end if
            end do
            residuos(2*jj) = -rBN
	
	
            ! Verificando necessidade de mais iteracao
            more = .false.
            do j = 1,2*jj
                if(abs(residuos(j))>tol) more = .true.
                !write(*,*)abs(residuos(j)),tol,more
            end do
            !read(*,*)
            if(.not.more)then
                write(*,*)"Fim da iteração"
                exit
            end if
	
	
	    !!!!!!!!!!!!!!!!!!
            !! Coeficientes !!
	    !!!!!!!!!!!!!!!!!!
	
            ! condicoes de montante
            dBdQ1 = 0.
            dBdh1 = 1.
            do j = 1,jj-1
                dx = x(j+1)-x(j)
                ! Continuidade
                dCdQ(j,j) = -teta/dx
                dCdh(j,j) = 1./(2.*dt)*dAdh(n+1,j)
                dCdQ(j,j+1) = teta/dx
                dCdh(j,j+1) = 1./(2.*dt)*dAdh(n+1,j+1)
                ! Momento
                dMdQ(j,j) = rterm(j)/(2.*dt) - rterm(j)*2.*teta*Q(n+1,j)/(A(n+1,j)*dx) &
                     & + g*teta*(A(n+1,j+1)+A(n+1,j))/4.*dSfdQ(n+1,j)
                dMdh(j,j) = rterm(j)*teta/dx*Q(n+1,j)**2/A(n+1,j)**2*dAdh(n+1,j) - g*teta*(A(n+1,j+1)+A(n+1,j))/(2.*dx) &
                     & + g*teta*(h(n+1,j+1)-h(n+1,j))/(2.*dx)*dAdh(n+1,j) + g*teta*(Sf(n+1,j+1)+Sf(n+1,j))/4.*dAdh(n+1,j) &
                     & + g*teta*(A(n+1,j+1)+A(n+1,j))/4.*dSfdh(n+1,j)
                dMdQ(j,j+1) = rterm(j+1)/(2.*dt) + rterm(j+1)*2.*teta*Q(n+1,j+1)/(A(n+1,j+1)*dx)&
                     & + g*teta*(A(n+1,j+1)+A(n+1,j))/4.*dSfdQ(n+1,j+1)
                dMdh(j,j+1) = -rterm(j+1)*teta/dx*Q(n+1,j+1)**2/A(n+1,j+1)**2*dAdh(n+1,j+1) + g*teta*(A(n+1,j+1)+A(n+1,j))/(2.*dx) &
                     & + g*teta*(h(n+1,j+1)-h(n+1,j))/(2.*dx)*dAdh(n+1,j+1) + g*teta*(Sf(n+1,j+1)+Sf(n+1,j))/4.*dAdh(n+1,j+1) &
                     & + g*teta*(A(n+1,j+1)+A(n+1,j))/4.*dSfdh(n+1,j+1)
            end do
            ! condicoes de jusante
            dBdQN = 0.
            dBdhN = 1.
	  
            !do j = 1,jj-1
            !    write(*,12)j,dCdQ(j,j),dCdh(j,j),dCdQ(j,j+1),dCdh(j,j+1),dMdQ(j,j),dMdh(j,j),dMdQ(j,j+1),dMdh(j,j+1)
            !    12 format(i2,8f10.4)
            !end do
            !read(*,*)
	
	
	    !!!!!!!!!!!!!!!!!!
	    !! DOUBLE SWEEP !!
	    !!!!!!!!!!!!!!!!!!
            call double_sweep(jj,dBdQ1,dBdh1,dBdQN,dBdhN,dCdQ,dCdh,dMdQ,dMdh,residuos,W)
            do j = 1,jj
                dQ(j) = W(j,1)
                dh(j) = W(j,2)
                !write(*,*)dQ(j),dh(j)
            end do


            ! Atualizando valores
            Q(n+1,:) = Q(n+1,:) + dQ
            h(n+1,:) = h(n+1,:) + dh
            !write(*,'(10f8.4)')Q(n+1,:)
            !write(*,'(10f8.4)')h(n+1,:)
            !read(*,*)
	
            call canal_trapezoidal(g,jj,z,h(n+1,:),b,s,Q(n+1,:),nM,A(n+1,:),dAdh(n+1,:),P(n+1,:),&
	                      & dPdh,Sf(n+1,:),dSfdQ(n+1,:),dSfdh(n+1,:),Fr(n+1,:))

        end do

        if(processo=="aquecimento")then
            ! Controle de erro:
	    do j = 1,jj
                if(isnan(h(n,j)))then
                    write(*,'(2a)')"Não foi possivel realizar o aquecimento do evento ",trim(evento)
                    write(*,*)" Os níveis iniciais são muito baixos. "
                    write(*,*)" O evento será cortado e o processo reniciado..."
                    write(arq,'(3a)')"sed -e '1d' ",trim(arq_evento)," > temp"
                    call system(arq)
		    write(arq_evento,'(5a)')"./",trim(fundo),"/eventos_cut/evento_",trim(evento),"_cut.txt"
                    write(*,*)arq_evento
                    write(arq,'(2a)')"mv temp ",trim(arq_evento)
                    call system(arq)
                    deallocate(h,y,A,P,Q,yup,yinter,ydown,Sf,dSfdQ,dSfdh,dAdh,rC,rM,Fr)
                    deallocate(residuos,W,dQ,dh,rterm,norma,dCdQ,dCdh,dMdQ,dMdh,amd,hms)
                    goto 40
	        end if
	    end do
        end if
        
    END DO


    if(processo=="aquecimento")then
        ! Gravando arquivos p/ simulacao
        write(8,'(3a)')evento," ",trim(arq_evento)
    else if(processo=="simulacao")then
        ! Calculando o coeficiente de eficiencia Nash-Sutcliffe p/ nivel
        ! procurando NaN
        do n = 0,nn
            if(isnan(h(n,100)))exit
        end do
        ENS = 1. - sum((yinter(0:n-1)-(h(0:n-1,jinter)-z(jinter)))**2)/sum((yinter(0:n-1)-sum(yinter(0:n-1))/real(n))**2)
        write(8,'(a10,f10.4)')evento,ENS
    end if

    ! Impressao dos Resultados
    if(processo=="aquecimento")then
        write(arq,'(2a)')trim(dir),"CI_aquecido.txt"
        write(*,*)arq
        open(1,file=arq,action="write")
        write(1,*)"Q0=",Q(nn,jj)
        do j =1,jj
            write(1,11)x(j),S0(j),z(j),h(nn,j)
        end do
        close(1)
        write(arq,'(2a)')"cp linha.gnu ",trim(dir)
        call system(arq)

    else if(processo=="simulacao")then
        write(arq,'(2a)')trim(dir),"res_xh.txt"
        open(1,file=arq,action="write")
        write(arq,'(2a)')trim(dir),"res_xQ.txt"
        open(2,file=arq,action="write")
        write(arq,'(2a)')trim(dir),"res_xy.txt"
        open(3,file=arq,action="write")
        write(arq,'(2a)')trim(dir),"res_xFr.txt"
        open(4,file=arq,action="write")
        write(arq,'(2a)')trim(dir),"res_xv.txt"
        open(55,file=arq,action="write")
        do j = 1,jj
            write(1,*)x(j),z(j),(h(n,j),n=0,nn)
            write(2,*)x(j),(Q(n,j),n=0,nn)
            write(3,*)x(j),(h(n,j)-z(j),n=0,nn)
            write(4,*)x(j),(Fr(n,j),n=0,nn)
            write(55,*)x(j),(Q(n,j)/A(n,j),n=0,nn)
        end do
        close(1) ; close(3) ; close(2) ; close(4) ; close(55)
    
        write(arq,'(2a)')trim(dir),"res_tQ.txt"
        open(2,file=arq,action="write")
        write(arq,'(2a)')trim(dir),"res_ty.txt"
        open(3,file=arq,action="write")
        write(arq,'(2a)')trim(dir),"res_tFr.txt"
        open(4,file=arq,action="write")
        write(arq,'(2a)')trim(dir),"res_tv.txt"
        open(55,file=arq,action="write")
        do n = 0,nn
            write(2,*)amd(n)," ",hms(n),(Q(n,j),j=1,jj)
            write(3,*)amd(n)," ",hms(n),(h(n,j)-z(j),j=1,jj)
            write(4,*)amd(n)," ",hms(n),(Fr(n,j),j=1,jj)
            write(55,*)amd(n)," ",hms(n),(Q(n,j)/A(n,j),j=1,jj)
        end do
        close(2) ; close(3) ; close(4) ; close(55)

        write(arq,'(2a)')trim(dir),"curva-chave.txt"
        open(1,file=arq,action="write")
        do n = 0,nn
            write(1,*)h(n,1)-z(1),Q(n,1),h(n,jinter)-z(jinter),Q(n,jinter),h(n,jj)-z(jj),Q(n,jj)
        end do
        close(1)

        ! procurando o valor maximo e minimo de Froude e velocidade para cada secao
        Frmax = -100000. ; Frmin = 100000.
        vmax = -100000. ; vmin = 100000.
        do n = 0,nn
	    do j = 1,jj
	        if(Fr(n,j)>Frmax(j)) Frmax(j) = Fr(n,j)
	        if(Fr(n,j)<Frmin(j)) Frmin(j) = Fr(n,j)
	        if(Q(n,j)/A(n,j)>vmax(j)) vmax(j) = Q(n,j)/A(n,j)
	        if(Q(n,j)/A(n,j)<vmin(j)) vmin(j) = Q(n,j)/A(n,j)
	    end do
        end do
        write(arq,'(2a)')trim(dir),"Fr_maxmin.txt"
        open(1,file=arq,action="write")
        write(arq,'(2a)')trim(dir),"v_maxmin.txt"
        open(2,file=arq,action="write")
        do j = 1,jj
	    write(1,*)j,Frmax(j),Frmin(j)
	    write(2,*)j,vmax(j),vmin(j)
        end do
        close(1) ; close(2)

    end if

    deallocate(h,y,A,P,Q,yup,yinter,ydown,Sf,dSfdQ,dSfdh,dAdh,rC,rM,Fr)
    deallocate(residuos,W,dQ,dh,rterm,norma,dCdQ,dCdh,dMdQ,dMdh,amd,hms)

end do
97 close(13)
close(8)

write(*,*)"Fim!!! =D"

33 continue
end program

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine double_sweep(n,dBdQ1,dBdh1,dBdQN,dBdhN,dCdQ,dCdh,dMdQ,dMdh,res,W)
! Anne
! resolucao do sistema pelo metodo da dupla varredura
implicit none

integer:: n,j,k
real:: res(2*n)
real:: dBdQ1,dBdh1,dBdQN,dBdhN
real,dimension(n,n):: dCdQ,dCdh,dMdQ,dMdh
real,dimension(n,2):: Z,F,W
real,dimension(n,2,2):: D,C,A,E
real:: MM(2,2),iMM(2,2),iD(2,2)

do j = 1,n
    do k = 1,2
	Z(j,k) = res(2*j-2+k)
    end do
end do

do j = 1,n
    if(j==1)then
	D(j,1,1) = dBdQ1
	D(j,1,2) = dBdh1
	D(j,2,1) = dCdQ(j,j)
	D(j,2,2) = dCdh(j,j)
    else if(j==n)then
    	D(j,1,1) = dMdQ(j-1,j)
	D(j,1,2) = dMdh(j-1,j)
	D(j,2,1) = dBdQN
	D(j,2,2) = dBdhN
    else
        D(j,1,1) = dMdQ(j-1,j)
	D(j,1,2) = dMdh(j-1,j)
	D(j,2,1) = dCdQ(j,j)
	D(j,2,2) = dCdh(j,j)
    end if
    if(j==1)then
	A(j,:,:) = 0.
    else
	A(j,1,1) = dMdQ(j-1,j-1)
	A(j,1,2) = dMdh(j-1,j-1)
	A(j,2,1) = 0.
	A(j,2,2) = 0.
    end if
    if(j==n)then
	C(j,:,:) = 0.
    else
	C(j,1,1) = 0.
	C(j,1,2) = 0.
	C(j,2,1) = dCdQ(j,j+1)
	C(j,2,2) = dCdh(j,j+1)
    end if
!     write(*,'(6f10.4)')D(j,1,:),C(j,1,:),A(j,1,:)
!     write(*,'(6f10.4)')D(j,2,:),C(j,2,:),A(j,2,:)
!     read(*,*)
end do

! forward
do j = 1,n
    if(j==1)then
	call matriz_inversa_adj_2x2(D(j,:,:),iD)
! 	write(*,*)D(j,1,:),iD(1,:)
! 	write(*,*)D(j,2,:),iD(2,:)
! 	read(*,*)
	E(j,:,:) = -matmul(iD,C(j,:,:))
	F(j,:) = matmul(iD,Z(j,:))
    else
	MM = D(j,:,:) + matmul(A(j,:,:),E(j-1,:,:))
	call matriz_inversa_adj_2x2(MM,iMM)
! 	write(*,*)MM(1,:),iMM(1,:)
! 	write(*,*)MM(2,:),iMM(2,:)
! 	read(*,*)
	E(j,:,:) = -matmul(iMM,C(j,:,:))
	F(j,:) = matmul(iMM,Z(j,:)-matmul(A(j,:,:),F(j-1,:)))
    end if
!     write(*,*)E(j,1,:),F(j,:)
!     write(*,*)E(j,2,:)
end do
! read(*,*)

! backward
do j = n,1,-1
    if(j==n)then
	W(j,:) = F(j,:)
    else
	W(j,:) = matmul(E(j,:,:),W(j+1,:))+F(j,:)
    end if
    !write(*,*)j,W(j,:)
end do
!read(*,*)


return
end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine matriz_inversa_adj_2x2(A,B)
! Anne 05/01/2015
! Inversão de matriz atraves da matriz adjunta
! considerando APENAS matrizes 2x2
! entrada: A = matriz 2x2 que se deseja obter a inversa
! saída: B = matriz invertida
implicit none

integer:: i,j
real:: det
real,dimension(2,2):: A,B,Acof,Adj


Acof(1,1) = (-1.)**(1+1)*A(2,2)
Acof(1,2) = (-1.)**(1+2)*A(2,1)
Acof(2,1) = (-1.)**(2+1)*A(1,2)
Acof(2,2) = (-1.)**(2+2)*A(1,1)

Adj = Acof
do i = 1,2
    do j = 1,2
	if(i/=j)then
	    Adj(i,j) = Acof(j,i)
	end if
    end do
end do
 
det = Adj(1,1)*Adj(2,2)-Adj(1,2)*Adj(2,1)

B = Adj/det

return

end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine canal_trapezoidal(g,n,z,h,b,s,Q,nM,A,dAdh,P,dPdh,Sf,dSfdQ,dSfdh,Fr)
! Anne
! calculo das propriedades do canal trapezoidal
! s = declividade da parede do canal
implicit none

integer,intent(in):: n
integer:: i
real,intent(in):: g,b,s,nM
real:: ka
real,intent(out):: dPdh
real,dimension(n),intent(in):: z,h,Q
real,dimension(n),intent(out):: A,P,Sf,dSfdQ,dSfdh,dAdh,Fr
real,dimension(n):: y


if(g>10.)then
    ka = 1.49    ! unidades britanicas
else
    ka = 1.      ! SI
end if

y = h-z
A = (b+y*s)*y
dAdh = b+2.*s*(h-z)
P = b+2.*y*sqrt(1.+s**2)
dPdh = 2.*sqrt(1.+s**2)

Sf = (nM**2*P**(4./3.)*Q*abs(Q))/(ka**2*A**(10./3.))
dSfdQ = (2*nM**2*P**(4./3.)*abs(Q))/(ka**2*A**(10./3.))
dSfdh = (nM**2*Q*abs(Q)/(ka**2*A**(20./3.)))* &
      & (4./3.*P**(1./3.)*dPdh*A**(10./3.)-10./3.*A**(7./3.)*dAdh*P**(4./3.))

Fr = Q/A/sqrt(g*A/dAdh)

return
end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
