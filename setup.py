from setuptools import find_packages, setup

# read the contents of README file
from os import path
from io import open  # for Python 2 and 3 compatibility

# get __version__ from _version.py
ver_file = path.join('biomni', 'version.py')
with open(ver_file) as f:
    exec(f.read())

this_directory = path.abspath(path.dirname(__file__))


# read the contents of README.md
def readme():
    with open(path.join(this_directory, 'README.md'), encoding='utf-8') as f:
        return f.read()


setup(name='biomni',
      version=__version__,
      license='MIT',
      description='Biomni',
      long_description=readme(),
      long_description_content_type='text/markdown',
      url='https://github.com/snap-stanford/biomni',
      author='Biomni Team',
      author_email='kexinh@cs.stanford.edu',
      packages=find_packages(exclude=['test']),
      zip_safe=False,
      include_package_data=True,
      install_requires=[],
      setup_requires=['setuptools>=38.6.0']
      )