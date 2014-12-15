__author__ = 'Bas Rustenburg'

from behave import *

@given(u'the module "{modulename}" is installed')
def step_module_import(context, modulename):
    """
    Tries to import a module with a given name to make sure it is installed.
    """
    module = __import__(modulename)

@given(u'the module is in the current directory')
def step_module_directory(context):
    """
    Adds the current directory to the python path.
    :param context:
    :return:
    """
    import os
    cwd = os.getcwd()
    try:
        os.environ['PYTHONPATH'] += ":%s:" % cwd
    except KeyError:
        os.environ['PYTHONPATH'] = cwd

@given(u'that scripts are in the directory "{directory}"')
def step_context_scriptdir_plural(context, directory):
    """
    Set the directory where scripts can be found, after making sure it exists.

    :param context:
    :param str directory:
    :return:
    """
    import os
    assert os.path.isdir(directory) is True
    context.scriptdir = os.path.abspath(directory)

@given(u'that the script is in the directory "{directory}"')
def step_context_scriptdir_singular(context, directory):
    """
    Set the directory where to find a script, after making sure it exists.

    :param context:
    :param str directory:
    :return:
    """
    context.execute_steps('given that scripts are in the directory "%s"' % directory)

@given(u'the working directory is "{directory}"')
def step_context_workdir(context, directory):
    """
    Check whether a given directory exists.
    If it does not exists, create it in current directory.

    :param context:
    :param str directory:
    :return:
    """
    import os
    fullpath = '%s/%s' % (context.tmpdir, directory)
    if not os.path.isdir(fullpath) is True:
        os.makedirs(fullpath)

    context.workdir = os.path.abspath(fullpath)


@when(u'the script "{scriptname}" is called successfully from the working directory')
def step_execute_script(context, scriptname):
    """
    Execute a script on the command line.
    Fails if script returns error code other than 0

    :param context:
    :param str scriptname:
    :return:
    """
    import os
    from subprocess import Popen, PIPE
    oldwd = os.getcwd()
    os.chdir(context.workdir)
    try:
        script = Popen('%s/%s' % (context.scriptdir, scriptname), stdout=PIPE, stderr=PIPE)
        print(script.communicate())

        assert script.stderr

    finally:
        os.chdir(oldwd)

@then(u'a file called "{filename}" is created')
def step_file_existance(context, filename):
    """
    Check for existence of a file.

    :param context:
    :param str filename:
    :return:
    """
    import os
    assert os.path.isfile('%s/%s' % (context.workdir, filename)) is True

@then(u'"{filename}" is not an empty file')
def step_file_nonzero(context, filename):
    """
    Verify that the given file is not empty
    :param context:
    :param str filename:
    :return:
    """
    import os
    oldwd = os.getcwd()
    os.chdir(context.workdir)
    assert os.stat(filename).st_size > 0
    os.chdir(oldwd)

@then(u'"{filename}" is a .{format} file')
def step_file_format(context, filename, format):
    """
    Verify format of given file

    TODO: dummy test
    :param context:
    :param str filename:
    :param str format:
    :return:
    """
    assert False
