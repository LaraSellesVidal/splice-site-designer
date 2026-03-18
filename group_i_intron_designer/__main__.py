"""Allow running the package as ``python -m group_i_intron_designer``."""

import sys

from .cli import main

sys.exit(main())
