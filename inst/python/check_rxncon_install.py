"""
Checks if rxncon has been installed on the system, exits with an error code if it does not.
"""
import sys

try:
  import rxncon
except ImportError as e:
  print(e)
  sys.exit(1)
