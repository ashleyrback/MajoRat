import ROOT

import rat

import physics

import unittest
import os

class TestDoubleBeta(unittest.TestCase):
    def setUp(self):
        self._double_beta = physics.DoubleBeta("Te130")
    def test_init(self):
        self.assertEqual(self._double_beta._phase_space, 1529e-21)
        self.assertEqual(self._double_beta._matrix_element, 3.31)
    def test_get_half_life(self):
        self.assertEqual(self._double_beta.get_half_life(), 
                         5.969480351183813e+16)
    def tearDown(self):
        pass
class TestZeroNuConverter(unittest.TestCase):
    def setUp(self):
        self._converter = physics.ZeroNuConverter("Te130")
        self._converter.set_conversion_factor()
        self._mass = 0.316940621063
        seed = os.getpid()
        self._random = ROOT.TRandom3(seed)
    def test_init(self):
        self.assertEqual(self._converter._two_nu_t_half,
                         5.969480351183813e+16)
        self.assertEqual(self._converter._t_half_min, 5.969480351183813e+10)
        self.assertEqual(self._converter._t_half_max, 5.9694803511838125e+31)
        self.assertEqual(self._converter._phase_space, 14.22e-15)
        self.assertEqual(self._converter._matrix_element, 4.03)
        self.assertEqual(self._converter._coupling_constant, 1.269)
    def test_conversion_factor(self):
        self.assertEqual(self._converter.get_conversion_factor(), 
                         6.603009171798192e+11)
    def test_get_mass_max(self):
        self.assertEqual(self._converter.get_mass_max(), 2702549.374000821)
    def test_get_mass_min(self):
        self.assertEqual(self._converter.get_mass_min(), 8.546211510904833e-05)
    def test_mass_to_half_life(self):
        self.assertEqual(self._converter.mass_to_half_life(self._mass),
                         4.3403823804398253e+24)
    def test_half_life_to_mass(self):
        self.assertEqual(self._converter.half_life_to_mass\
                             (self._converter.mass_to_half_life(self._mass)),
                         self._mass)
    def test_converter_equality(self):
        t_half = self._random.Uniform(1.0e24, 10.0e24)
        t_half_from_mass = self._converter.mass_to_half_life\
            (self._converter.half_life_to_mass(t_half))
        self.assertAlmostEqual(t_half_from_mass/1.0e24,
                               t_half/1.0e24, 12)
    def test_extremes(self):
        effective_mass_max = self._converter.half_life_to_mass(0.0)
        t_half_min = self._converter.mass_to_half_life(effective_mass_max)
        t_half_max = self._converter.mass_to_half_life(0.0)
        effective_mass_min = self._converter.half_life_to_mass(t_half_max)
        self.assertEqual(effective_mass_max, 2702549.374000821)
        self.assertEqual(t_half_min, 59694803511.838104)
        self.assertEqual(t_half_max, 5.969480351183813e+31)
        self.assertEqual(effective_mass_min, 8.546211510904833e-05)
    def tearDown(self):
        pass

if __name__ == "__main__":
    unittest.main()
