#!/usr/bin/python3
"""Imports code/ from this checkout, not from the copy baked into the image"""

import os, sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
