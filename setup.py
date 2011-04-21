#!/usr/bin/env python

from distutils.core import setup

setup(name='MethylAnalyzer',
      version='0.1.0',
      description='Analyze Methyl-MAPS data',
      author='Yurong Xin',
      author_email='xinyuro@pi.cpmc.columbia.edu',
      url='http://github.com/epigenomics/methylmaps',
      packages=['MethylAnalyzer'],
      scripts = ['scripts/create_cpg_track.py', 'scripts/create_frag_track.py',
                 'scripts/create_microarray_track.py', 'scripts/create_wiggle_track.py',
                 'scripts/filter.py', 'scripts/parse_mates.py', 'scripts/parse_sites.py',
                 'scripts/run_pipeline.py', 'scripts/score.py'],
      package_data={'MethylAnalyzer': ['data/*']},
      requires=['numpy', 'pyfasta', 'pysam'],
      license='GPL',
      )

