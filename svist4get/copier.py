import shutil
import os


def copy():
    current_directory = os.getcwd()
    internal_directory = os.path.join(os.path.dirname(__file__), 'svist4get_data')
    current_directory = current_directory + '/svist4get_data'

    shutil.copytree(internal_directory, current_directory)
