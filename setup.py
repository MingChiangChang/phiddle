from setuptools import setup, find_packages

setup(
    name='phiddle',
    version='0.0.1',
    description='Phase labeling gui',
    author='Ming-Chiang Chang',
    author_email='mc2663@cornell.edu',
    python_requires='>=3.6.0',
    url='http://github.com/mingchianghchang/phiddle',
    packages=find_packages(),
    install_requires=['numpy', 'matplotlib', 'julia', 'pyqt6', 'pymatgen', 'xrayutilities', 'astropy']
)
