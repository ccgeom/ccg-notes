import unittest
import pyvista as pv

from ccgeom.manifold import Manifold


class TestDiff(unittest.TestCase):

    def setUp(self):
        self.poly = pv.Polygon()
        self.disc = pv.Disc()
        self.box = pv.Box()
        self.cone = pv.Cone()

    def tearDown(self):
        pass

    def test_poly(self):
        m = Manifold(self.poly)

    def test_disc(self):
        m = Manifold(self.disc)

    def test_box(self):
        m = Manifold(self.box)

    def test_cone(self):
        m = Manifold(self.cone)
