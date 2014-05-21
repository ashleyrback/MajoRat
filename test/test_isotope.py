import ROOT

import rat

import isotope
import physics

import unittest
import os

class TestIsotope(unittest.TestCase):
    def setUp(self):
        self._isotope = isotope.Isotope("Xe136")
    def tearDown(self):
        pass
class TestSNOPlusTe(unittest.TestCase):
    def setUp(self):
        self._te130 = isotope.SNOPlusTe()
        self._converter = physics.ZeroNuConverter("Te130")
        seed = os.getpid()
        self._random = ROOT.TRandom3(seed)
    def test_init(self):
        self.assertEqual(self._te130._av_volume, 903.3)
        self.assertEqual(self._te130._density_lab, 0.862e3)
    def test_get_name(self):
        self.assertEqual(self._te130.get_name(), "Te130")
    def test_mass(self):
        te_mass = 810.475591248 # kg
        self.assertEqual(self._te130.get_mass(), te_mass)
    def test_mass_in_fv(self):
        te_mass = 240.53402185501733 # kg
        self.assertEqual(self._te130.get_mass_in_fv(), te_mass)
    def test_number_nuclei(self):
        self.assertIsNone(self._te130._number_nuclei)
        self.assertEqual(self._te130.get_number_nuclei(), 
                         1.1150580505749745e+27)
        self.assertEqual(self._te130.get_number_nuclei(False), 
                         3.7571705068828865e27)
    def test_decays_from_t_half(self):
        self.assertEqual(self._te130.get_decays_from_t_half(5e+25), 
                         15.457986878334248)
    def test_decays_from_mass(self):
        self.assertEqual(self._te130.get_decays_from_mass(0.0933806512323), 
                         15.45798687833113)
    def test_decays_equality(self):
        t_half = self._random.Uniform(1e24, 10e24)
        decays_from_t_half = self._te130.get_decays_from_t_half(t_half)
        decays_from_mass = self._te130.get_decays_from_mass\
            (self._converter.half_life_to_mass(t_half))
        self.assertAlmostEqual(decays_from_t_half/1.0e27, 
                               decays_from_mass/1.0e27, 12)
    def test_extremes(self):
        effective_mass_max = self._converter.half_life_to_mass(0.0) # when t_half=0
        t_half_min = self._converter.mass_to_half_life(effective_mass_max)
        t_half_max = self._converter.mass_to_half_life(0.0) # when mass=0
        effective_mass_min = self._converter.half_life_to_mass(t_half_max)
        self.assertAlmostEqual(self._te130.get_decays_from_mass\
                                   (effective_mass_max)/1.0e16,
                               1.2947514665383534, 12)
        self.assertAlmostEqual(self._te130.get_decays_from_t_half\
                                   (t_half_min)/1.0e16,
                               1.2947514665383534, 12)
        self.assertAlmostEqual(self._te130.get_decays_from_mass\
                                   (effective_mass_min)/1.0e-5,
                               1.2947514665383534, 12)
        self.assertAlmostEqual(self._te130.get_decays_from_t_half\
                                   (t_half_max)/1.0e-5,
                               1.2947514665383534, 12)
    def tearDown(self):
        pass

if __name__ == "__main__":
    unittest.main()
