from snewpy.neutrino import MixingParameters, ThreeFlavorMixingParameters, FourFlavorMixingParameters

class ThreeFlavorTransformation:
    _mixing_params = ThreeFlavorMixingParameters(**MixingParameters())
    @property
    def mixing_params(self):
        return self._mixing_params

    @mixing_params.setter
    def mixing_params(self, val):
        return self._mixing_params.update(**val)

class FourFlavorTransformation(ThreeFlavorTransformation):
    _mixing_params = FourFlavorMixingParameters(**MixingParameters())