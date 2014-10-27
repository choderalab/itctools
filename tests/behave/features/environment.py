__author__ = 'Bas Rustenburg'

import shutil
import os

def before_all(context):
    """
    Remove the temporary directories that exist and make a clean one.

    :param context:
    :return:
    """
    os.chdir('../..')
    shutil.rmtree('tmp', ignore_errors=True)
    os.makedirs('tmp')
    context.tmpdir = os.path.abspath('tmp')
