        program TD4
            implicit none
            integer i
            real R, T_fusion, delta_cp, h_fusion
            real x_max, x_eau, x_ethanol, T_cong

            R = 8.314
            x_max = 0.4

            print*, "Valeur de la Temperature de fusion de l'eau [K]:"
            read*, T_fusion

            print*, "Valeur de Delta Cp [J/mol.K]:"
            read*, delta_cp

            print*, "Valeur de la Enthalpie de fusion de l'eau [J/mol]:"
            read*, h_fusion

            do i = 0, 50 !quebrei o intervalo em 50 pedacinhos
                x_ethanol = DBLE(i) * 0.008d0
                x_eau = 1.0d0 - x_ethanol
                if (x_ethanol<1.0d-6) then ! x = 0, é água pura
                    T_cong = T_fusion
                else
                    call calculer_T_cong(x_eau, T_fusion, h_fusion,
     +              delta_cp, R, T_cong)
                endif
            enddo

            print*, "x ethanol:", x_ethanol
            print*, "x eau:", x_eau
            print*, "Temperature de congelation:", T_cong, "K"



        end program TD4
        !------------------------------
        ! subroutine
        !------------------------------
        subroutine calculer_T_cong(x_eau, T_fusion, h_fusion,
     +      delta_cp, R, T_cong)
            implicit none
            real, intent(in) :: x_eau, T_fusion, h_fusion,
     +      delta_Cp, R
            real, intent (out)  :: T_cong
            real T, f, df, zero
            integer i, i_max

            T = T_fusion - 0.5 !Valeur initiale
            zero = 1.0d-6
            i_max = 100 ! valeur maximal d'itérations

            do i = 1, i_max
            !f est la fonction qu'on doit trouver le zero; df est la derivee, avec Cp cte:
                f = - (h_fusion/R)*(1/T - 1/T_fusion) +
     +          (delta_cp/R)*( (T_fusion/T) - 1.0 + LOG(T/T_fusion) )
     +          - LOG(x_eau)
            !df(T) est la derivee
                df = - (h_fusion/R)*(-1.0d0/(T**2))
     +          + (delta_cp/R)*(-T_fusion/(T**2) + 1.0d0/T)

                if(abs(df) < 1.0d-10) exit
                T_cong = T - f/df

                if (abs(T_cong-T)<zero) return

                T = T_cong
            enddo
        end subroutine calculer_T_cong


