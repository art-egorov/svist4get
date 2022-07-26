from setuptools import setup

with open("README.md", "r") as fh:
      long_description = fh.read()

setup(name='svist4get',
      version='1.3',
      description='A simple visualization tool for genomic tracks from sequencing experiments',
      url='https://github.com/art-egorov/svist4get',
      author='Artyom Egorov',
      author_email='artyom.egorov@hotmail.com',
      license='WTFPL',
      packages=['svist4get'],
      install_requires = ['reportlab', 'biopython','configs','argparse','Pybedtools', 'wand', 'statistics'],
      long_description = long_description,
      long_description_content_type = "text/markdown",
      scripts = ['bin/svist4get', 'bin/svist4get_copier'],
      zip_safe=False, 
      package_data={'svist4get':['svist4get_data/*', 'svist4get_data/fonts/*','svist4get_data/palettes/*','svist4get_data/triplet_codes/*','svist4get_data/help/*']},
      include_package_data=True)
