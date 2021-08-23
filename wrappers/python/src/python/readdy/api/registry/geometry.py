import readdy._internal.readdybinding.api.geom as _geom
from readdy.api.utils import vec3_of as _v3_of


class Geometries:
    def __init__(self, units):
        self._units = units

    def create_sphere(self, center, radius):
        center = self._units.convert(center, self._units.length_unit)
        radius = self._units.convert(radius, self._units.length_unit)
        return _geom.Sphere(_v3_of(center), radius)

    def create_box(self, v0, v1):
        v0 = self._units.convert(v0, self._units.length_unit)
        v1 = self._units.convert(v1, self._units.length_unit)
        return _geom.Box(_v3_of(v0), _v3_of(v1))

    def create_capsule(self, center, direction, radius, length):
        center = self._units.convert(center, self._units.length_unit)
        direction = self._units.convert(direction, self._units.length_unit)
        radius = self._units.convert(radius, self._units.length_unit)
        length = self._units.convert(length, self._units.length_unit)
        return _geom.Capsule(_v3_of(center), _v3_of(direction), radius, length)
