import ROOT

import rat

from production import Production

import test_spectrum_data

import unittest
import os

path = "/home/ashley/snoplus/analysis/doubleBeta/MajoRat/test/TeLoadedTe130_0n2b_r10_s0_p1.ntuple.root"
path2 = "/home/ashley/snoplus/analysis/doubleBeta/MajoRat/test/TeLoadedTe130_0n2b_r10_s0_p2.ntuple.root"

class TestProduction(unittest.TestCase):
    def setUp(self):
        self._spectrum = Production(path)
        self._spectrum.set_parameters()
    def test_majorat_name(self):
        majorat_name = "RAT4-5_prod_decay0_2beta_Te130_0_1_0-0_3-5"
        self.assertEqual(self._spectrum._majorat_name, majorat_name)
    def test_extension(self):
        self.assertEqual(self._spectrum._ext, ".ntuple.root")
    def testMakeHistogram(self):
        hist_label = self._spectrum._label + "-Test"
        append = False
        always_remake = True
        self._spectrum.make_histogram(hist_label, append, always_remake)
        histogram = self._spectrum.get_histogram(hist_label)
        self.assertEqual(histogram.GetEntries(), 10001)
    def testAppend(self):
        hist_label = self._spectrum._label + "-Test"
        append = False
        always_remake = True
        hist_path = os.environ.get("MAJORAT_TEST")
        self._spectrum.make_histogram(hist_label, append, always_remake)
        self._spectrum.write_histograms(hist_path+"/")
        histogram = self._spectrum.get_histogram(hist_label)
        self.assertEqual(histogram.GetEntries(), 10001)
        append = True
        always_remake = False
        self._spectrum2 = Production(path2)
        self._spectrum2.set_parameters()
        bin_width = "default"
        self._spectrum2.make_histogram(hist_label, append, always_remake, 
                                       bin_width, hist_path+"/")
        histogram2 = self._spectrum2.get_histogram(hist_label)
        self.assertNotEqual(histogram.GetEntries(), histogram2.GetEntries())
        self.assertEqual(histogram2.GetEntries(), 20003)
    def tearDown(self):
        pass
class TestProductionNoSetParameters(test_spectrum_data.TestSpectrumDataNonMajoRatFile):
    def setUp(self):
        self._spectrum = Production(path)
    def tearDown(self):
        pass
class TestProduction2Beta(test_spectrum_data.TestSpectrumData2Beta):
    def setUp(self):
        self._spectrum = Production(path)
        self._spectrum.set_parameters()
    def test_n_events(self):
        pass
    def test_label(self):
        self.assertEqual(self._spectrum._label, "0#nu#beta#beta")
    def tearDown(self):
        pass
class TestProductionSolar(test_spectrum_data.TestSpectrumDataSolar):
    def setUp(self):
        path = "TeLoadedB8_solar_r10_s0_p1.ntuple.root"
        self._spectrum = Production(path)
        self._spectrum.set_parameters()
    def test_n_events(self):
        pass
    def test_label(self):
        self.assertEqual(self._spectrum._label, "B8")
    def tearDown(self):
        pass
class TestProductionDecayChain(test_spectrum_data.TestSpectrumDataDecayChain):
    def setUp(self):
        path = "TeLoadedRh102_r10_s0_p1.ntuple.root"
        self._spectrum = Production(path)
        self._spectrum.set_parameters()
    def test_n_events(self):
        pass
    def test_label(self):
        self.assertEqual(self._spectrum._label, "Rh102")
    def tearDown(self):
        pass
class TestProductionBackg(test_spectrum_data.TestSpectrumDataBackg):
    def setUp(self):
        path = "TeLoadedK42_r10_s0_p1.ntuple.root"
        self._spectrum = Production(path)
        self._spectrum.set_parameters()
    def test_n_events(self):
        pass
    def test_label(self):
        self.assertEqual(self._spectrum._label, "K42")
    def tearDown(self):
        pass
