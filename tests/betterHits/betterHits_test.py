from unittest import TestCase
import os
from samUtilities.betterHits import BetterHits

class TestBetterHits(TestCase):

  def test_sanity(self):
    x = BetterHits()
    x.log.debug("This is a test: %s", os.getcwd())
    x.log.debug("Who am I?  %s", __file__)
    self.assertTrue(x != None)
