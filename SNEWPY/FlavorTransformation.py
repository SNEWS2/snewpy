'''
  Copyright (c) 2020 James Kneller

  This file is part of SNEWPY.

  SNEWPY is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  SNEWPY is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
'''

from abc import abstractmethod, ABC

import numpy as np

theta12=np.deg2rad(33.) 
theta13=np.deg2rad(9.) 
theta23=np.deg2rad(45.)

class FlavorTransformation(ABC):

    @abstractmethod
    def p(self):
        pass
    
    @abstractmethod
    def pbar(self):
        pass    


class NoTransformation(FlavorTransformation):
      def p(self):
          return 1.
    
      def pbar(self):
          return 1.


class AdiabaticMSW_NMO(FlavorTransformation):
      def p(self):
          return pow(np.sin(theta13),2.)
    
      def pbar(self):
          return pow(np.cos(theta12)*np.cos(theta13),2.)    



class AdiabaticMSW_IMO(FlavorTransformation):
      def p(self):
          return pow(np.sin(theta12)*np.cos(theta13),2.)
    
      def pbar(self):
          return pow(np.sin(theta13),2.)


class TwoFlavorDecoherence(FlavorTransformation):
      def p(self):
          return 0.5
    
      def pbar(self):
          return 0.5


class ThreeFlavorDecoherence(FlavorTransformation):
      def p(self):
          return 1./3.
    
      def pbar(self):
          return 1./3.
