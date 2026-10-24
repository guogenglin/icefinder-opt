# -*- coding: utf-8 -*-
"""
Created on Fri Oct 23 15:45:35 2026

@author: Genglin Guo
@e-mail: 2019207025@njau.edu.cn
"""

from __future__ import annotations

from pathlib import Path

from .log import ICEfinderError


def meta_prediction(input_file: Path, output_dir: Path, threads: int) -> None:
    '''
    Reject the unsupported legacy metagenome mode before starting tools or writing output files.
    '''
    raise ICEfinderError(
        'The legacy metagenome mode is not supported in icefinder-opt_v1. '
        'Use assembly mode for assembled genomes. '
        'The original legacy implementation remains in src/icefinder-opt/metaICE.py.'
    )
