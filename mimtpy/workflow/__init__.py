## Dynamic import for modules used in routine workflows, i.e. mimtpyApp.py
## Recommended usage:
##  import mimtpy
##  import mimtpy.workflow

from pathlib import Path
import importlib

# expose the following modules
__all__=[
     'HDFEOS_to_geotiff',
     'plot_geotiff',
     'track_offset',
     'generate_horzvert',
     'H5UNW_to_geotiff',
     'save_geodmod',
]

root_module = Path(__file__).parent.parent.name #mimtpy
for module in __all__:
    importlib.import_module(root_module + '.' + module)
