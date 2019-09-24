# coding=utf-8

from setuptools import setup, find_packages

LONG_DESCRIPTION = '''Mobile element insertion (MEI) is a major category of structure variations (SVs). 
The rapid development of long read sequencing provides the opportunity to sensitively discover MEIs. 
However, the signals of MEIs implied by noisy long reads are highly complex, due to the repetitiveness 
of mobile elements as well as the serious sequencing errors. Herein, we propose Realignment-based 
Mobile Element insertion detection Tool for Long read (rMETL). rMETL takes advantage of 
its novel chimeric read re-alignment approach to well handle complex MEI signals. 
Benchmarking results on simulated and real datasets demonstrated that rMETL has the ability 
to more sensitivity discover MEIs as well as prevent false positives. 
It is suited to produce high quality MEI callsets in many genomics studies.'''

setup(
    name = "rMETL",
    version = "1.0.4",
    description = "realignment-based Mobile Element insertion detection Tool for Long read",
    author = "Jiang Tao",
    author_email = "tjiang@hit.edu.cn",
    url = "https://github.com/tjiangHIT/rMETL",
    license = "MIT Licence",
    packages = find_packages('src'),
    package_dir = {'': 'src'},
    long_description = LONG_DESCRIPTION,
    zip_safe = False,
    install_requires = ['pysam', 'Biopython', 'Cigar']
)
