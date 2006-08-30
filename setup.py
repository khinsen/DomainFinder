#!/usr/bin/env python

from distutils.core import setup

class Dummy:
    pass
pkginfo = Dummy()
execfile('DomainFinderLib/__pkginfo__.py', pkginfo.__dict__)

setup (name = "DomainFinder",
       version = pkginfo.__version__,
       description = "Domain motion analysis for proteins",
       long_description =
"""DomainFinder is an interactive program for the determination and
characterization of dynamical domains in proteins. Its key features
are

 - computational efficiency: even large proteins can be analyzed
   using a desktop computer in a few minutes

 - ease of use: a state-of-the-art graphical user interface

 - export of results for visualization and further analysis
   (VRML, PDB, and MMTK object format)
""",
       author = "Konrad Hinsen",
       author_email = "hinsen@cnrs-orleans.fr",
       url = "http://dirac.cnrs-orleans.fr/DomainFinder/",
       license = "GPL",
       packages = ['DomainFinderLib'],
       scripts = ['DomainFinder', 'DomainFinderModes', 'TransitionPath'],
       )


