Get bolometric estimates of SNe

Note: You need the absolute LATEST version of astroquery to work with SVO

Author: Ian Thornton

Install with pip:

```bash
pip install extrabol
```

Example input:

```python
extrabol ./example/PSc000174_extrabol.dat --verbose -m 1a
```

This example may take a minute the first time you run it as astroquery fetches filter information!

Filters must be specified by their [SVO](http://svo2.cab.inta-csic.es/svo/theory/fps3/) ID.

