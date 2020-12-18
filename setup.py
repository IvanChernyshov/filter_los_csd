import setuptools

with open('README.md', 'r', encoding = 'utf-8') as inpf:
    long_description = inpf.read()

setuptools.setup(
    name = 'filter_los_csd',
    version = '1.4.1',
    author = 'Ivan Yu. Chernyshov',
    author_email = 'ivan.chernyshoff@gmail.com',
    description = 'CML tool filtering line-of-sights contacts',
    long_description = long_description,
    long_description_content_type = 'text/markdown',
    url = 'https://github.com/EPiCs-group/filter_los_csd',
    packages = [],
    scripts = ['bin/filter_los_csd.py'],
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
        'Development Status :: 4 - Beta',
        'Topic :: Scientific/Engineering :: Chemistry'
    ],
    install_requires=[
        'numpy==1.19.3', # TODO: remove restriction after fmod() fix: https://developercommunity.visualstudio.com/content/problem/1207405/fmod-after-an-update-to-windows-2004-is-causing-a.html
        'pandas>=1.0',
        'PyCifRW>=4.0'
    ],
    python_requires = '>=3.6',
)
