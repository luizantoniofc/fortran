        program TD1
            implicit none
            integer i, np
            doubleprecision somax, somay, somaxy, somax2, somay2
            doubleprecision x, y, a, b, cc, numerador, denomin

             somax2 = 0
             somax = 0
             somay = 0
             somaxy = 0
             somay2 = 0

            print*, "Nombre de points:"
            read*, np

            do i=1, np
                print*, "x", i, ":"
                read*, x
                somax = somax + x
                somax2 = somax2 + x**2

                print*, "y", i, ":"
                read*, y
                somay = somay + y
                somay2 = somay2 + y**2
                somaxy = somaxy + x*y
            enddo

            a =((np*somaxy)-(somax*somay))/((np*somax2)-(somax*somax))
            b=((somax2*somay)-(somax*somaxy))/(np*somax2-(somax*somax))


            numerador = np * somaxy - somax * somay
            denomin =(np*somax2-somax*somax)*(np*somay2-somay * somay)
            cc = numerador / sqrt(denomin)

            print*, "Somme x^2:", somax2
            print*, "Somme x:", somax
            print*, "Somme y:", somay
            print*, "Somme xy:", somaxy
            print*, "Somme y^2:", somay2

            print*
            print*, "Pente (a):", a
            print*, "Origine (b):", b
            print*, "Coefficient de correlation (r):", cc

        end program

