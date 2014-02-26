from distutils.core import setup

setup(
    name='l2pGUI',
    version='0.9',
    author='jcrod',
    author_email='',
    url='www.sgf.rgo.ac.uk',
    packages=['l2pGUI'],
    scripts=['bin/l2pgui_run.py'],
    data_files=[('l2pGUI', ['conf/l2pGUI.cfg'])],
    description='Graphical display client for listen2planes',
    long_description=open('README.txt').read(),
    )
