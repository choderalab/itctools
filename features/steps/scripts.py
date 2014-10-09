__author__ = 'Bas Rustenburg'

from behave import given, when, then, use_step_matcher

import sys
import os
sys.path.append(os.path.abspath('../../'))

@given(u'that scripts are in the directory "{directory}"')
def step_impl(context, directory):
    """
    Set the directory for running scripts, after making sure it exists.

    :param context:
    :param str directory:
    :return:
    """
    import os
    assert os.path.isdir(directory) is True
    context.scripts = os.path.abspath(directory)

@given(u'the working directory is "{directory}"')
def step_impl(context, directory):
    """
    Check whether a given directory exists.
    If it does not exists, create it in current directory.

    :param context:
    :param str directory:
    :return:
    """
    import os
    if not os.path.isdir(directory) is True:
        os.makedirs(directory)

    context.workdir = os.path.abspath(directory)


@when(u'the script "{scriptname}" is called')
def step_impl(context, scriptname):
    """
    Execute a script on the command line.

    :param context:
    :param str scriptname:
    :param str dirname:
    :return:
    """
    import os
    oldwd = os.getcwd()
    os.chdir(context.workdir)
    print(os.system('%s/%s' % (context.scripts, scriptname)))
    os.chdir(oldwd)

@then(u'a file called "{filename}" is created')
def step_impl(context, filename):
    """
    Check for existence of a file.

    :param context:
    :param str filename:
    :return:
    """
    import os
    assert os.path.isfile('%s/%s' % (context.workdir, filename)) is True

@then(u'"{filename}" is not an empty file')
def step_impl(context, filename):
    """
    Verify that the given file is not empty
    :param context:
    :param str filename:
    :param str format:
    :return:
    """
    import os
    oldwd = os.getcwd()
    os.chdir(context.workdir)
    assert os.stat(filename).st_size > 0
    os.chdir(oldwd)

@then(u'"{filename}" is a .{format} file')
def step_impl(context, filename, format):
    """
    Verify format of given file

    TODO: dummy test
    :param context:
    :param str filename:
    :param str format:
    :return:
    """
    assert False
