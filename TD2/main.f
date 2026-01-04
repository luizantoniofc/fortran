
        program TD2
            implicit none

            integer n

            real, dimension(:), allocatable :: vetor_y, vetor_x
            real, dimension(:,:), allocatable :: matriz_a

            print*, "Dimension du Système Linéaire (N):"
            read*, n

            allocate(matriz_a(n,n), vetor_y(n), vetor_x(n))

            call creer_matrice(n, vetor_y, matriz_a)

            call triangularisation(n, matriz_a, vetor_y)
            call resolution(n, vetor_x, vetor_y, matriz_a)
            call ecrite(n, vetor_x, vetor_y, matriz_a)

            print*, "Solução X:", vetor_x

        ! ----------------------------------------------------
        ! eu separo o programa principal das subroutines
        ! ----------------------------------------------------
        contains

            subroutine creer_matrice(n, vetor_y, matriz_a)
                integer, intent(in) :: n !entree
                real, dimension(n), intent(out) :: vetor_y    !sortie
                real, dimension(n,n), intent(out):: matriz_a  !sortie

                integer m, p, h


            ! créer la matrice
                do m = 1, n
                    do p = 1, n
                    print*, "Element de la matrice: A(", m, ",", p, "):"
                    read*, matriz_a(m,p)
                    enddo
                enddo
            ! créer le vecteur y
                print*, "Vecteur y:"
                do h = 1, n
                    print*, "y(", h,")"
                    read*, vetor_y(h)
                enddo

            end subroutine creer_matrice


            subroutine triangularisation(n, matriz_a, vetor_y)
                integer, intent(in) :: n
                real, dimension(n,n), intent(inout) :: matriz_a
                real, dimension(n), intent(inout) :: vetor_y
                integer k, j, l
                real dp, d

                do k = 1, n-1
                    dp = 1.0/matriz_a(k,k)
                    do j = k, n
                        matriz_a(k,j) = matriz_a(k,j)*dp
                    enddo
                    vetor_y(k) = vetor_y(k)*dp

                    do j = k+1, n
                        d = matriz_a(j, k)
                        do l = k, n
                            matriz_a(j,l)= matriz_a(j,l)-d*matriz_a(k,l)
                        enddo
                        vetor_y(j) = vetor_y(j)-(d* vetor_y(k))
                    enddo
                enddo

                dp = 1.0 / matriz_a(n,n)
                matriz_a(n,n) = matriz_a(n,n) * dp
                vetor_y(n) = vetor_y(n) * dp

            end subroutine triangularisation

            ! 2° balayge
            subroutine resolution(n, vetor_x, vetor_y, matriz_a)
                integer, intent(in) :: n
                real, dimension(n), intent(out) :: vetor_x
                real, dimension(n), intent(in) :: vetor_y
                real, dimension(n,n), intent(in) :: matriz_a

                integer i, j, linha
                real soma_termo
                vetor_x(n) = vetor_y(n)/matriz_a(n,n)

                do i = 1, n-1
                    soma_termo = 0.0
                    linha = n-i
                    do j = linha + 1, n
                        soma_termo = soma_termo + (matriz_a(linha, j)*
     +                  vetor_x(j))
                    enddo

                    vetor_x(linha) = (vetor_y(linha) - soma_termo) /
     +               matriz_a(linha, linha)
                enddo
            end subroutine resolution

            subroutine ecrite (n, vetor_x, vetor_y, matriz_a)
                integer, intent(in) :: n
                real, dimension(n, n), intent(in) :: matriz_a
                real, dimension(n), intent(in) :: vetor_x, vetor_y
                integer i, j, unidade
                character(len=50) :: nome_arquivo

                print*, "Saisissez le nom de le fichier pour enregistrer
     +           le résultat:"
                read*, nome_arquivo
                unidade = 10

                open(unit=unidade, file=nome_arquivo, status='unknown')

                write(unidade, *) "Résultats de la méthode de Gauss"
                write(unidade, *) "Dimensions du systme (N): ", n
                write(unidade, *)
                write(unidade, *) '----------------------'

                print*, "Solution de X:"
                do i = 1, n
                    write(unidade, *) "X", i, "=", vetor_x(i)
                enddo
                close(unit=unidade)
                print*, "Données enregistrées avec succès dans le

     +          fichier : ", nome_arquivo
            endsubroutine
        end program TD2
