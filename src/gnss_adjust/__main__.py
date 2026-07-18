"""Entry point for ``python -m gnss_adjust``."""

import sys

from .cli import main

if __name__ == "__main__":
    sys.exit(main())
