__author__ = 'Bas Rustenburg'

import shutil
import os

def before_all(context):
    shutil.rmtree('tmp', ignore_errors=True)
    os.makedirs('tmp')
    context.tmpdir = os.path.abspath('tmp')
