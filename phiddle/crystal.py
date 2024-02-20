import numpy as np

class BaseCrystal():

    def __init__(self, a, b, c, α, β, γ):

        self.a = a
        self.b = b
        self.c = c
        self.α = α
        self.β = β 
        self.γ = γ 


    @property
    def free_param_num(self):
        return 6

    def get_free_params(self):
        return [self.a, self.b, self.c, self.α, self.β, self.γ]

    def set_params(self, a, b, c, α, β, γ):
        self.a = a
        self.b = b
        self.c = c
        self.α = α
        self.β = β
        self.γ = γ



    # Default util functions
    def get_peak(self, h, k, l):
        return (2 * np.pi/self.volume *
               np.sqrt(h**2 * self.b**2 * self.c**2 * np.sin(self.α)**2
                     + k**2 * self.a**2 * self.c**2 * np.sin(self.β)**2
                     + l**2 * self.a**2 * self.b**2 * np.sin(self.γ)**2
                     + 2*h * k * self.a * self.b * self.c**2 * (np.cos(self.α)*np.cos(self.β) - np.cos(self.γ))
                     + 2*k * l * self.a**2 * self.b * self.c * (np.cos(self.β)*np.cos(self.γ) - np.cos(self.α))
                     + 2*h * l * self.a * self.b**2 * self.c * (np.cos(self.α)*np.cos(self.γ) - np.cos(self.β))
                     )
                )

    @property
    def volume(self):
        return (self.a * self.b * self.c
                * np.sqrt( 1+2*np.cos(self.α)*np.cos(self.β)*np.cos(self.γ)
                        - np.cos(self.α)**2 - np.cos(self.β)**2 - np.cos(self.γ)**2 ) )

    def __repr__(self):
        return (f"a: {self.a}, b: {self.b}, c: {self.c}\n"
                f"α: {self.α}, β: {self.β}, γ: {self.γ}")


class Cubic(BaseCrystal):

    def __init__(self, a, b, c, α, β, γ):

        super().__init__(a, b, c, α, β, γ)

    def get_free_params(self):
        return [self.a]

    def set_params(self, a):
        self.a = a
        self.b = a
        self.c = a

    @property
    def free_param_num(self):
        return 1 

    def get_peak(self, h: int, k: int, l: int):
        return 2 * np.pi * np.sqrt(h**2 + k**2 + l**2) / self.a

    @property
    def volume(self):
        return a**3


class Tetragonal(BaseCrystal):

    def __init__(self, a, b, c, α, β, γ):

        super().__init__(a, b, c, α, β, γ)


    def set_params(self, a, c):
        self.a = a
        self.b = a
        self.c = c 

    def get_free_params(self):
        return [self.a, self.c]
    
    @property
    def free_param_num(self):
        return 2

    def get_peak(self, h: int, k: int, l:int):
        return 2 * np.pi * np.sqrt(h**2*self.c**2 + k**2*self.c**2 + l**2*self.a**2) / (self.a*self.c)

    @property
    def volume(self):
        return self.a**2 * self.c


class Orthorhombic(BaseCrystal):

    def __init__(self, a, b, c, α, β, γ):

        super().__init__(a, b, c, α, β, γ)


    def set_params(self, a, b, c):
        self.a = a
        self.b = b 
        self.c = c

    def get_free_params(self):
        return [self.a, self.b, self.c]

    @property
    def free_param_num(self):
        return 3


    def get_peak(self, h: int, k: int, l: int):
        return (2 * np.pi / self.volume *
                np.sqrt(h**2 * self.b**2 * self.c**2
                + k**2 * self.a**2 * self.c**2
                + l**2 * self.a**2 * self.b**2))

    @property
    def volume(self):
        return self.a * self.b * self.c


class Monoclinic(BaseCrystal):

    def __init__(self, a, b, c, α, β, γ):

        super().__init__(a, b, c, α, β, γ)


    def set_params(self, a, b, c, β):
        self.a = a
        self.b = b
        self.c = c
        self.β = β


    def get_free_params(self):
        return [self.a, self.b, self.c, self.β]

    @property
    def free_param_num(self):
        return 4

    def get_peak(self, h: int, k: int, l: int):
        return (2 * np.pi / self.volume 
                * sqrt(h**2 * self.b**2 * self.c**2
                     + k**2 * self.a**2 * self.c**2 * np.sin(self.β)**2
                     + l**2 * self.a**2 * self.b**2
                     - 2*h * l * self.a * self.b**2 * self.c * np.cos(self.β)))

    @property
    def volume(self):
        return self.a * self.b * self.c * np.sin(self.β)


class Hexagonal(BaseCrystal):

    def __init__(self, a, b, c, α, β, γ):

        super().__init__(a, b, c, α, β, γ)


    def set_params(self, a, c,):
        self.a = a
        self.b = a 
        self.c = c

    def get_free_params(self):
        return [self.a, self.c]

    def get_peak(self, h: int, k: int, l: int):
        return ( 2 * np.pi / self.volume *
                 np.sqrt(h**2 * self.b**2 * self.c**2
                       + k**2 * self.a**2 * self.c**2
                       + l**2 * self.a**2 * self.b**2 * np.sin(self.γ)**2
                       + 2*h * k * self.a * self.b * self.c**2 * 
                            (np.cos(self.α)*np.cos(self.β) - np.cos(self.γ))
                         )
                )


    @property
    def free_param_num(self):
        return 2


    @property
    def volume(self):
        return (a * b * c
        * sqrt( 1+2*cos(α)*cos(β)*cos(γ)
                - cos(α)^2 - cos(β)^2 - cos(γ)^2 ) )

class Rhombohedral(BaseCrystal):

    def __init__(self, a, b, c, α, β, γ):

        super().__init__(a, b, c, α, β, γ)


    def set_params(self, a, α):
        self.a = a
        self.b = a
        self.c = a 
        self.α = α
        self.β = α 
        self.γ = α

    def get_free_params(self):
        return [self.a, self.α]

    @property
    def free_param_num(self):
        return 2


    def get_peak(self, h: int, k: int, l: int):
        return (2 * np.pi/self.volume *
                np.sqrt(h**2 * self.a**4 * np.sin(self.α)**2
                      + k**2 * self.a**4 * np.sin(self.α)**2
                      + l**2 * self.a**4 * np.sin(self.α)**2
                      + ((2*self.a**4)*(np.cos(self.α)**2 - np.cos(self.α))
                          * (h * k + k * l + h * l)
                         )
                    ) 
                )

    # Temporarily use the default volume function
    # @property
    # def volume(self):
    #     return (self.a * self.b * self.c
    #             * np.sqrt( 1+2*np.cos(self.α)*np.cos(self.β)*np.cos(self.γ)
    #                     - np.cos(self.α)**2 - np.cos(self.β)**2 - np.cos(self.γ)**2 ) )


class Triclinic(BaseCrystal):

    def __init__(self, a, b, c, α, β, γ):

        super().__init__(a, b, c, α, β, γ)




def get_crystal(a, b, c, α, β, γ):

    if α==β==γ:
        if α==90.:
            if a==b==c:
                return Cubic(a, b, c, α, β, γ) 
            elif a==b and a!=c:
                return Tetragonal(a, b, c, α, β, γ)
            elif a!=b and b!=c:
                return Orthorhombic(a, b, c, α, β, γ) 
        elif a==b==c and α!=90.:
            return Rhombohedral(a, b, c, α, β, γ)
        else:
            return Triclinic(a, b, c, α, β, γ)

    elif α==γ==90. and β!=90.:
        return Monoclinic(a, b, c, α, β, γ)

    elif α==β==90. and γ==120.:
        return Hexagonal(a, b, c, α, β, γ)
    else:
        return Triclinic(a, b, c, α, β, γ)


