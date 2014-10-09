__author__ = 'Bas Rustenburg'

from behave import *

@given('a directory called "{directory}"')
def step_impl(context, directory):
    """
    Check whether a given directory exists.
    If it does not exists, create it in the current path.

    :param context:
    :param str directory:
    :return:
    """
    pass

@when('the python script "{scriptname}" is called')
def step_impl(context, scriptname):
    """
    Execute a python script on the command line.

    :param context:
    :param str scriptname:
    :return:
    """
    pass

@then('a file called "{filename}" is created')
def step_impl(context, filename):
    """
    Check for existence of a file.

    :param context:
    :param str filename:
    :return:
    """
    pass

@then('"{filename}" is formatted as a? ."{format}" file')
def step_impl(context, format):
    """
    Verify the format of a file.

    :param context:
    :param str format:
    :return:
    """
    pass
