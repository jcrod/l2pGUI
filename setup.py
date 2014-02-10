from distutils.core import setup

setup(
    name='l2pGUI',
    version='0.9',
    author='Jose Rodriguez',
    author_email='josrod@nerc.ac.uk',
    url='www.sgf.rgo.ac.uk',
    packages=['l2pGUI'],
    scripts=['bin/l2pgui_run.py'],
    description='Graphical display client for listen2planes',
    long_description=open('README.txt').read(),
    )
