import versioneer
import setuptools


setuptools.setup(name='xlink_kme_sbml',
      description='Kinetic Simulation Framework for Crosslinks',
      author='Kai-Michael Kammer',
      author_email='kai-michael.kammer@uni-konstanz.de',
      url='https://github.com/stengel-laboratory/xlink-kme-sbml',
      version=versioneer.get_version(),
      cmdclass=versioneer.get_cmdclass(),
      packages=setuptools.find_packages(),
      install_requires=[
          'numpy',
          'pandas',
          'tellurium',
          ''
      ],
      license='MIT',
      python_requires='>=3.10'
      )
