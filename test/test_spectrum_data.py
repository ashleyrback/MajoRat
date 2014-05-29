import ROOT

import rat

from spectrum_data import SpectrumData
import physics
import isotope

import unittest

path = "/home/ashley/snoplus/analysis/doubleBeta/MajoRat/test/RAT4-5_1k_decay0_2beta_Te130_0_1_0_3-5.root"

class TestSpectrumData(unittest.TestCase):
    def setUp(self):
        self._spectrum = SpectrumData(path, 5.0e25)
        self._converter = physics.ZeroNuConverter("Te130")
    def test_rat_release(self):
        self.assertEqual(self._spectrum._rat_release, "RAT4.5")
    def test_n_events(self):
        self.assertEqual(self._spectrum._n_events, 1e3)
    def test_generator(self):
        self.assertEqual(self._spectrum._generator.get_generator(), "decay0")
    def test_type(self):
        self.assertEqual(self._spectrum._generator.get_type(), "2beta")
    def test_isotope(self):
        self.assertIsInstance(self._spectrum._generator.get_isotope(),
                              isotope.SNOPlusTe)
        self.assertEqual(self._spectrum._generator.get_isotope().get_name(),
                         "Te130")
    def test_level(self):
        self.assertEqual(self._spectrum._generator.get_level(), 0)
    def test_mode(self):
        self.assertEqual(self._spectrum._generator.get_mode(), 1)
    def test_e_low(self):
        self.assertEqual(self._spectrum._generator.get_e_lo(), 0.0)
    def test_e_high(self):
        self.assertEqual(self._spectrum._generator.get_e_hi(), 3.5)
    def test_spectral_index(self):
        self.assertEqual(self._spectrum._spectral_index, 0)
    def test_label(self):
        self.assertEqual(self._spectrum._label, "0#nu#beta#beta")
    def test_dsreader(self):
        events = rat.dsreader(self._spectrum._path)
        ds, run = events.next()
        self.assertIsInstance(ds, ROOT.RAT.DS.Root)
        self.assertIsInstance(run, ROOT.RAT.DS.Run)
    def test_read_from_file(self):
        entries = 0
        events = 0
        for ds, run in rat.dsreader(self._spectrum._path):
            entries += 1
            events += ds.GetEVCount()
        self.assertEqual(entries, 1000)
        self.assertEqual(events, 1021)
    def test_get_histogram(self):
        self.assertIsNone(self._spectrum._histogram)
        self.assertIsNone(self._spectrum._unscaled_histogram)
        always_recreate=True
        hist = self._spectrum.get_histogram(always_recreate)
        self.assertIsInstance(hist, ROOT.TH1D)
        self.assertEqual(hist.Integral(), 1021)
        self.assertIsInstance(self._spectrum._unscaled_histogram, ROOT.TH1D)
        self.assertEqual(self._spectrum._unscaled_histogram.Integral(), 1021)
    def test_get_number_nuclei(self):
        self.assertEqual(self._spectrum.get_number_nuclei(), 
                         1.1150580505749745e+27)
    def test_set_events_by_t_half(self):
        self._spectrum.set_events_by_t_half()
        self.assertEqual(self._spectrum._events, 15.457986878334248)
        self._spectrum.set_events_by_t_half(0.0)
        self.assertEqual(self._spectrum._events, 1.294751466538353e+16)
        t_half_max = self._converter.mass_to_half_life(0.0)
        self._spectrum.set_events_by_t_half(t_half_max)
        self.assertEqual(self._spectrum._events, 1.2947514665383531e-05)
    def test_set_events_by_mass(self):
        self._spectrum.set_events_by_mass(0.05)
        self.assertEqual(self._spectrum._events, 4.431789725182429)
        self._spectrum.set_events_by_mass(0.0)
        self.assertEqual(self._spectrum._events, 1.2947514665383531e-5)
        effective_mass_max = self._converter.half_life_to_mass(0.0)
        self._spectrum.set_events_by_t_half(effective_mass_max)
        self.assertEqual(self._spectrum._events, 1.294751466538353e16)
    def test_scale_by_t_half(self):
        self._spectrum.scale_by_t_half()
        self.assertEqual(self._spectrum._histogram.Integral(), 
                         15.457986878334248)
        self.assertEqual(self._spectrum._unscaled_histogram.Integral(), 1021)
        self._spectrum.scale_by_t_half(0)
        self.assertAlmostEqual(self._spectrum._histogram.Integral()/1e+16, 
                         1.294751466538353, 12)
        self.assertEqual(self._spectrum._unscaled_histogram.Integral(), 1021)
        self._spectrum.scale_by_t_half()
        self.assertAlmostEqual(self._spectrum._histogram.Integral(), 
                         15.457986878334248, 12)
        self.assertEqual(self._spectrum._unscaled_histogram.Integral(), 1021)
    def test_scale_by_mass(self):
        self._spectrum.scale_by_mass(0.05)
        self.assertEqual(self._spectrum._histogram.Integral(), 
                         4.431789725182429)
        self.assertEqual(self._spectrum._unscaled_histogram.Integral(), 1021)
        self._spectrum.scale_by_mass(0)
        self.assertAlmostEqual(self._spectrum._histogram.Integral()/1e-5,
                               1.2947514665383531, 12)
        self.assertEqual(self._spectrum._unscaled_histogram.Integral(), 1021)
        self._spectrum.scale_by_mass(0.05)
        self.assertAlmostEqual(self._spectrum._histogram.Integral(), 
                         4.431789725182429, 12)
        self.assertEqual(self._spectrum._unscaled_histogram.Integral(), 1021)
    def tearDown(self):
        pass
class TestSpectrumData2Beta(unittest.TestCase):
    def setUp(self):
        self._spectrum = SpectrumData(path, 5.0e25)
    def test_rat_release(self):
        self.assertEqual(self._spectrum._rat_release, "RAT4.5")
    def test_n_events(self):
        self.assertEqual(self._spectrum._n_events, 1e3)
    def test_generator(self):
        self.assertEqual(self._spectrum._generator.get_generator(), "decay0")
    def test_type(self):
        self.assertEqual(self._spectrum._generator.get_type(), "2beta")
    def test_isotope(self):
        self.assertIsInstance(self._spectrum._generator.get_isotope(),
                              isotope.SNOPlusTe)
        self.assertEqual(self._spectrum._generator.get_isotope().get_name(),
                         "Te130")
    def test_level(self):
        self.assertEqual(self._spectrum._generator.get_level(), 0)
    def test_mode(self):
        self.assertEqual(self._spectrum._generator.get_mode(), 1)
    def test_e_low(self):
        self.assertEqual(self._spectrum._generator.get_e_lo(), 0.0)
    def test_e_high(self):
        self.assertEqual(self._spectrum._generator.get_e_hi(), 3.5)
    def test_spectral_index(self):
        self.assertEqual(self._spectrum._spectral_index, 0)
    def test_label(self):
        self.assertEqual(self._spectrum._label, "0#nu#beta#beta")
    def tearDown(self):
        pass
class TestSpectrumDataSolar(unittest.TestCase):
    def setUp(self):
        self._spectrum = SpectrumData("RAT4-5_1k_solar_B8.root")
    def test_rat_release(self):
        self.assertEqual(self._spectrum._rat_release, "RAT4.5")
    def test_n_events(self):
        self.assertEqual(self._spectrum._n_events, 1e3)
    def test_generator(self):
        self.assertEqual(self._spectrum._generator.get_generator(), "solar")
    def test_type(self):
        self.assertEqual(self._spectrum._generator.get_type(), "nue")
    def test_isotope(self):
        self.assertEqual(self._spectrum._generator.get_isotope(), "B8")
    def test_e_low(self):
        self.assertEqual(self._spectrum._generator.get_e_lo(), 0.0)
    def test_e_high(self):
        self.assertEqual(self._spectrum._generator.get_e_hi(), 3.5)
    def test_spectral_index(self):
        self.assertIsNone(self._spectrum._spectral_index)
    def test_label(self):
        self.assertEqual(self._spectrum._label, "B8")
    def tearDown(self):
        pass
class TestSpectrumDataDecayChain(unittest.TestCase):
    def setUp(self):
        self._spectrum = SpectrumData("RAT4-5_1k_decaychain_Rh102.root")
    def test_rat_release(self):
        self.assertEqual(self._spectrum._rat_release, "RAT4.5")
    def test_n_events(self):
        self.assertEqual(self._spectrum._n_events, 1e3)
    def test_generator(self):
        self.assertEqual(self._spectrum._generator.get_generator(), "decaychain")
    def test_type(self):
        self.assertIsNone(self._spectrum._generator.get_type())
    def test_isotope(self):
        self.assertEqual(self._spectrum._generator.get_isotope(), "Rh102")
    def test_e_low(self):
        self.assertEqual(self._spectrum._generator.get_e_lo(), 0.0)
    def test_e_high(self):
        self.assertEqual(self._spectrum._generator.get_e_hi(), 3.5)
    def test_spectral_index(self):
        self.assertIsNone(self._spectrum._spectral_index)
    def test_label(self):
        self.assertEqual(self._spectrum._label, "Rh102")
    def tearDown(self):
        pass
class TestSpectrumDataBackg(unittest.TestCase):
    def setUp(self):
        self._spectrum = SpectrumData("RAT4-5_1k_decay0_K42.root")
    def test_rat_release(self):
        self.assertEqual(self._spectrum._rat_release, "RAT4.5")
    def test_n_events(self):
        self.assertEqual(self._spectrum._n_events, 1e3)
    def test_generator(self):
        self.assertEqual(self._spectrum._generator.get_generator(), "decay0")
    def test_type(self):
        self.assertEqual(self._spectrum._generator.get_type(), "backg")
    def test_isotope(self):
        self.assertEqual(self._spectrum._generator.get_isotope(), "K42")
    def test_e_low(self):
        self.assertEqual(self._spectrum._generator.get_e_lo(), 0.0)
    def test_e_high(self):
        self.assertEqual(self._spectrum._generator.get_e_hi(), 3.5)
    def test_spectral_index(self):
        self.assertIsNone(self._spectrum._spectral_index)
    def test_label(self):
        self.assertEqual(self._spectrum._label, "K42")
    def tearDown(self):
        pass
class TestSpectrumDataNonMajoRatFile(unittest.TestCase):
    def setUp(self):
        self._spectrum = SpectrumData("TeLoadedTe130_0n2b_r10_s0_p1.ntuple.root")
    def test_rat_release(self):
        self.assertIsNone(self._spectrum._rat_release)
    def test_n_events(self):
        self.assertIsNone(self._spectrum._n_events)
    def test_generator(self):
        self.assertIsNone(self._spectrum._generator)
    def test_spectral_index(self):
        self.assertIsNone(self._spectrum._spectral_index)
    def test_label(self):
        self.assertIsNone(self._spectrum._label)
    def test_options(self):
        self.assertIsNotNone(self._spectrum._options)
    def tearDown(self):
        pass

if __name__ == "__main__":
    unittest.main()
