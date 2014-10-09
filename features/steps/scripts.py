__author__ = 'Bas Rustenburg'

from behave import given, when, then, use_step_matcher



@given(u'a directory called "{directory}"')
def step_impl(context, directory):
    """
    Check whether a given directory exists.
    If it does not exists, create it in the current path.

    :param context:
    :param str directory:
    :return:
    """
    pass

@when(u'the {lang} script "{scriptname}" is called')
def step_impl(context, lang, scriptname):
    """
    Execute a script on the command line.
    :param str lang:
    :param context:
    :param str scriptname:
    :return:
    """
    pass

@then(u'a file called "{filename}" is created')
def step_impl(context, filename):
    """
    Check for existence of a file.

    :param context:
    :param str filename:
    :return:
    """
    pass

@then(u'"{filename}" is a .{format} file')
def step_impl(context, filename, format):
    """
    Verify format of given file
    :param context:
    :param str filename:
    :param str format:
    :return:
    """
    pass
