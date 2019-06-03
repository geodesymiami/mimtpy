# LOGGING
import logging
import os, sys
import importlib

mimt_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(1, mimt_path)
sys.path.insert(1, os.path.join(mimt_path, 'defaults'))
sys.path.insert(1, os.path.join(mimt_path, 'objects'))
sys.path.insert(1, os.path.join(mimt_path, 'utils'))

try:
    os.environ['RSMAS_INSAR']
except KeyError:
    print('Using default MintPy Path: %s' % (insar_path))
    os.environ['RSMAS_INSAR'] = insar_path


logging.basicConfig(filename="example.log",
                            format='%(asctime)s | %(name)-25s | [ %(levelname)s ]'
                                                       ' | %(filename)s:%(lineno)d | %(message)s',
                                                                           level=logging.DEBUG)
ch = logging.StreamHandler()
verbose = False
if verbose:
    ch.setLevel(logging.DEBUG)
else:
    ch.setLevel(logging.ERROR)
warning_logger = logging.getLogger("process_rsmas")
warning_logger.addHandler(ch)
logger = logging.getLogger("process_rsmas." + "__init__")

logger.debug('Starting Logger')
#logger.error('YO WHATS GOOD???')
