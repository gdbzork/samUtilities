from unittest import TestCase
import os, logging
from samUtilities.betterHits import BetterHits
from tests import TestHandler

class TestBetterHits(TestCase):

  def setUp(self):
    self.hdlr = TestHandler()
    self.log = logging.getLogger("testLog")
    self.log.setLevel(logging.DEBUG)
    self.log.propagate = False
    self.log.addHandler(self.hdlr)
    self.testPath = os.path.join(os.path.dirname(__file__),"../testData")

  def test_sanity(self):
    bh = BetterHits(log=self.log)
    bh.log.debug("This is a test: %s", os.getcwd())
    bh.log.debug("Who am I?  %s", __file__)
    self.assertTrue(self.hdlr.matches("Who am I\?  .*betterHits_test.py$"))
    self.assertTrue(self.hdlr.matches("This is a test: .*samUtilities$"))

  def test_loadGTF(self):
    bh = BetterHits(log=self.log)
    regions = []
    data = bh.loadGTF(os.path.join(self.testPath,"test_2S_rRNA.gtf"),regions)
    self.assertEqual(3,len(data))
    self.assertEqual(("rDNA",16526,16555),data[0])
    self.assertEqual(("chrX",23291700,23291729),data[1])

  def test_loadHits_genomic(self):
    bh = BetterHits(log=self.log)
    data = bh.loadHits(os.path.join(self.testPath,"test_genomic_mm2.sam"))
    self.assertEqual(8,len(data)) # filters out non-primary hits, so 8 kept
    
  def test_loadHits_unwanted(self):
    bh = BetterHits(log=self.log)
    data = bh.loadHits(os.path.join(self.testPath,"test_unwanted_dm6.sam"))
    self.assertEqual(3,len(data))
    self.assertEqual(27,data["K00252:335:HWMMGBBXX:2:1101:25966:1437"])
    self.assertEqual(24,data["K00252:335:HWMMGBBXX:2:1101:2372:1455"])
    self.assertTrue(self.hdlr.matches("loaded 4 reference sequences \(3 primary\)$"))

  def test_overlap_sanity(self):
    bh = BetterHits(log=self.log)
    with bh.openFile(os.path.join(self.testPath,"test_genomic_mm2.sam")) as fd:
      data = fd.fetch(until_eof=True)
      data.__next__()
      data.__next__()
      rec = data.__next__()
      self.assertTrue(bh.overlaps(rec,[("chr2L",8000,10000)]))

  def test_overlap_wrongContig(self):
    bh = BetterHits(log=self.log)
    with bh.openFile(os.path.join(self.testPath,"test_genomic_mm2.sam")) as fd:
      data = fd.fetch(until_eof=True)
      data.__next__()
      data.__next__()
      rec = data.__next__()
      self.assertFalse(bh.overlaps(rec,[("zork",8000,10000)]))

  def test_overlap_wrongCoords(self):
    bh = BetterHits(log=self.log)
    with bh.openFile(os.path.join(self.testPath,"test_genomic_mm2.sam")) as fd:
      data = fd.fetch(until_eof=True)
      data.__next__()
      data.__next__()
      rec = data.__next__()
      self.assertFalse(bh.overlaps(rec,[("chr2L",6000,8000)]))

  def test_overlap_notFirstCoords(self):
    bh = BetterHits(log=self.log)
    with bh.openFile(os.path.join(self.testPath,"test_genomic_mm2.sam")) as fd:
      data = fd.fetch(until_eof=True)
      data.__next__()
      data.__next__()
      rec = data.__next__()
      self.assertTrue(bh.overlaps(rec,[("thing",0,10),("chr2L",8000,10000)]))

  def test_compare(self):
    bh = BetterHits(log=self.log)
    ref = os.path.join(self.testPath,"test_unwanted_dm6.sam")
    alt = os.path.join(self.testPath,"test_genomic_mm2.sam")
    results = bh.compare(ref,alt,None)
    self.assertEqual(2,len(results))
    self.assertEqual("K00252:335:HWMMGBBXX:2:1101:25966:1437",results[1][0])
    self.assertEqual(28,results[1][1])
    self.assertEqual(27,results[1][2])

  def test_compare_gtf(self):
    bh = BetterHits(log=self.log)
    ref = os.path.join(self.testPath,"test_unwanted_dm6.sam")
    alt = os.path.join(self.testPath,"test_genomic_mm2.sam")
    gtf = os.path.join(self.testPath,"test_2S_rRNA.gtf")
    results = bh.compare(ref,alt,gtf)
    self.assertEqual(1,len(results))
    self.assertEqual("K00252:335:HWMMGBBXX:2:1101:25966:1437",results[0][0])
    self.assertEqual(28,results[0][1])
    self.assertEqual(27,results[0][2])
    
