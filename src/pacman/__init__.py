import sys

# HACK: Sets the global encoding to utf-8 for output
sys.stdout.reconfigure(encoding='utf-8')

from . import lib
from . import s00_table
from . import s01_horizons
from . import s02_barycorr
from . import s03_refspectra
from . import s10_direct_images
from . import s20_extract
from . import s21_bin_spectroscopic_lc
from . import s30_run
