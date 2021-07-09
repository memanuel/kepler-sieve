import os
import sys

module_path = os.path.abspath(os.path.join('../src_py'))
if module_path not in sys.path:
    sys.path.append(module_path)