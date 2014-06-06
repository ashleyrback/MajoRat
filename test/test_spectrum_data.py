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
        self.assertEqual(self._spectrum._label, "0#nu#beta#beta-Truth")
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
        hist_label = self._spectrum._label + "-Test"
        self.assertIsNone(self._spectrum._histograms.get(hist_label))
        self.assertIsNone(self._spectrum._unscaled_histograms.get(hist_label))
        always_remake = True
        hist = self._spectrum.get_histogram(hist_label, always_remake)
        self.assertIsInstance(hist, ROOT.TH1D)
        self.assertEqual(hist.Integral(), 1021)
        unscaled_histogram = self._spectrum.get_unscaled_histogram(hist_label)
        self.assertIsInstance(unscaled_histogram, ROOT.TH1D)
        self.assertEqual(unscaled_histogram.Integral(), 1021)
    def test_get_number_nuclei(self):
        self.assertEqual(self._spectrum.get_number_nuclei(), 
                         7.47001779975032e+26)
    def test_set_events_by_t_half(self):
        self._spectrum.set_events_by_t_half()
        self.assertEqual(self._spectrum._events, 10.35564355325908)
        self._spectrum.set_events_by_t_half(0.0)
        self.assertEqual(self._spectrum._events, 8.673823301223735e+15)
        t_half_max = self._converter.mass_to_half_life(0.0)
        self._spectrum.set_events_by_t_half(t_half_max)
        self.assertEqual(self._spectrum._events, 8.673823301223735e-06)
    def test_set_events_by_mass(self):
        self._spectrum.set_events_by_mass(0.05)
        self.assertEqual(self._spectrum._events, 2.9689528822999485)
        self._spectrum.set_events_by_mass(0.0)
        self.assertEqual(self._spectrum._events, 8.673823301223737e-06)
        effective_mass_max = self._converter.half_life_to_mass(0.0)
        self._spectrum.set_events_by_t_half(effective_mass_max)
        self.assertEqual(self._spectrum._events, 8.673823301223735e15)
    def test_scale_by_t_half(self):
        hist_label = self._spectrum._label + "-Test"
        always_remake = True
        histogram = self._spectrum.get_histogram(hist_label, always_remake)
        t_half = "default"
        self._spectrum.scale_by_t_half(t_half, hist_label)
        self.assertEqual(histogram.Integral(), 10.35564355325908)
        unscaled_histogram = self._spectrum.get_unscaled_histogram(hist_label)
        self.assertEqual(unscaled_histogram.Integral(), 1021)
        t_half = 0.0
        self._spectrum.scale_by_t_half(t_half, hist_label)
        self.assertAlmostEqual(histogram.Integral()/1e+16, 
                               0.8673823301223735, 12)
        self.assertEqual(unscaled_histogram.Integral(), 1021)
        t_half = "default"
        self._spectrum.scale_by_t_half(t_half, hist_label)
        self.assertAlmostEqual(histogram.Integral(), 10.35564355325908, 12)
        self.assertEqual(unscaled_histogram.Integral(), 1021)
    def test_scale_by_mass(self):
        hist_label = self._spectrum._label + "-Test"
        always_remake = True
        histogram = self._spectrum.get_histogram(hist_label, always_remake)
        self._spectrum.scale_by_mass(0.05, hist_label)
        self.assertEqual(histogram.Integral(), 2.9689528822999485)
        unscaled_histogram = self._spectrum.get_unscaled_histogram(hist_label)
        self.assertEqual(unscaled_histogram.Integral(), 1021)
        self._spectrum.scale_by_mass(0, hist_label)
        self.assertAlmostEqual(histogram.Integral()/1e-5,
                               0.8673823301223735, 12)
        self.assertEqual(unscaled_histogram.Integral(), 1021)
        self._spectrum.scale_by_mass(0.05, hist_label)
        self.assertAlmostEqual(histogram.Integral(), 2.9689528822999485, 12)
        self.assertEqual(unscaled_histogram.Integral(), 1021)
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
        self.assertEqual(self._spectrum._label, "0#nu#beta#beta-Truth")
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
        self.assertEqual(self._spectrum._label, "B8-Truth")
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
        self.assertEqual(self._spectrum._label, "Rh102-Truth")
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
        self.assertEqual(self._spectrum._label, "K42-Truth")
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
