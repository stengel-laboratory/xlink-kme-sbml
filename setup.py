import versioneer
import setuptools


setuptools.setup(name='xlink-kme-sbml',
      description='Kinetic Simulation Framework for Crosslinks',
      author='Kai Kammer',
      author_email='kai-michael.kammer@uni-konstanz.de',
      url='https://git.uni-konstanz.de/kai-michael-kammer/xlink-kme-sbml',
      version=versioneer.get_version(),
      cmdclass=versioneer.get_cmdclass(),
      packages=setuptools.find_packages(),
      install_requires=[
          'numpy',
          'pandas',
          'tellurium'
      ],
      license='MIT',
      python_requires='>=3.6'
      )
