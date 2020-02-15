# -*- coding: utf-8 -*-
"""Simple supernova oscillation physics.
"""

from abc import abstractmethod, ABC

import numpy as np

theta12 = np.deg2rad(33.) 
theta13 = np.deg2rad(9.) 
theta23 = np.deg2rad(45.)

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
