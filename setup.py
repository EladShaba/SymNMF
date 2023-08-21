from setuptools import Extension, setup

module = Extension("symnmfmod", sources=['symnmf.c','symnmfmodule.c'])
setup(name='symnmfmod',
     version='1.0',
     description='Python wrapper for custom C extension',
     ext_modules=[module])
