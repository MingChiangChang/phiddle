import numpy as np


class BasePeakProfile():

    def __init__(self):
        pass

    # def __call__(self, x):
    #     pass


class Gauss(BasePeakProfile):

    def __init__(self):
        super().__init__()#sigma)
        self.n_param = 0


    def __call__(self, x):
        return np.exp( - x**2)

    @property
    def free_param_num(self):
        return 0


class Lorentz(BasePeakProfile):

    def __init__(self):
        super().__init__()#sigma)
        self.n_param = 0


    def __call__(self, x):
        return  1 / (1 + x**2)

    @property
    def free_param_num(self):
        return 0


class PseudoVoigt(BasePeakProfile):

    def __init__(self, alpha):
        super().__init__()#sigma)
        self.alpha = alpha
        self.n_param = 1 


    def __call__(self, x):
        return  self.alpha * np.exp( - x**2) + (1 - self.alpha) *(1 / (1 + x**2))


    @property
    def free_param_num(self):
        return 1
    # TODO: getter and setter for alpha to constrain it



class FixedPseudoVoigt(BasePeakProfile):

    def __init__(self, alpha):
        super().__init__()# sigma)
        self.alpha = alpha
        self.n_param = 1 


    def __call__(self, x):
        return  self.alpha * np.exp( - x**2) + (1 - self.alpha) *(1 / (1 + x**2))


    @property
    def free_param_num(self):
        return 0
