from .pwdftio import *

__all__ = ['PWDFT']

def clean_pycache():
    """Recursively delete all __pycache__ directories and .pyc files in this package."""
    import os
    import shutil
    pkg_dir = os.path.dirname(__file__)
    for root, dirs, files in os.walk(pkg_dir):
        for d in dirs:
            if d == '__pycache__':
                shutil.rmtree(os.path.join(root, d), ignore_errors=True)
        for f in files:
            if f.endswith('.pyc'):
                try:
                    os.remove(os.path.join(root, f))
                except Exception:
                    pass
