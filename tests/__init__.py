import sys
import re
from logging.handlers import BufferingHandler

class TestHandler(BufferingHandler):

  def __init__(self):
    BufferingHandler.__init__(self,0)

  def shouldFlush(self):
    return False

  def emit(self,record):
    self.format(record)
    self.buffer.append(record)

  def logCount(self):
    return len(self.buffer)

  def matches(self,reg):
    pat = re.compile(reg)
    found = False
    for rec in self.buffer:
      if pat.match(rec.message):
        found = True
        break
    return found
