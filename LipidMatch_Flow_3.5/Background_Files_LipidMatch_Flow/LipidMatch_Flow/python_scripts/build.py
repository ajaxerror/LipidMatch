#!python3.8

import os
import shutil
import subprocess

"""
Starting out this python build script for easy build management of the python packages.
Using python os and subprocess library to try to make the script os-independent.
"""

SHELL = 'cmd.exe'  # cmdline for windows

# define projects by their dirnames and main project filenames
python_projects = [
    {'dir': 'EICgen', 'file': 'EIC_gen.py'}
]

python_scripts_dir = os.path.dirname(os.path.abspath(__file__))

for project in python_projects:
    path = os.path.join(python_scripts_dir, project['dir'])
    file = os.path.join(path, project['file'])
    os.chdir(path)

    print(f"Compiling python project: {project['dir']}")

    # create virtual env
    subprocess.run('python -m venv env', shell=True)

    # activate virtual environment - the command here is shell specific, so
    # you will have to specify your shell at the top and extend the switch
    # to include your shell
    if SHELL == "cmd.exe":
        subprocess.run(f'source {os.path.join(path, "env", "Scripts", "activate.bat")}', shell=True)
    elif SHELL == 'bash' or SHELL == 'zsh':
        subprocess.run(f'{SHELL} -c "source {os.path.join(path, "env", "bin", "activate")}"', shell=True)
    subprocess.run('pip install -r requirements.txt', shell=True)

    # remove old build dirs
    shutil.rmtree(os.path.join(path, 'dist'))
    shutil.rmtree(os.path.join(path, 'build'))

    # compile
    subprocess.run(f'pyinstaller {file}', shell=True)

    # check for build success
    exe_base = project['file'].split('.')[0]
    if os.path.exists(
            os.path.join(path, 'dist', exe_base, f'{exe_base}.exe')
    ) or os.path.exists(
        os.path.join(path, 'dist', exe_base, exe_base)  # if compiled on mac or linux
    ):
        print(f"{project['dir']} compiled successfully.")
    else:
        print(f"{project['dir']} FAILED to compile!")

    # deactivate virtual env
    if SHELL == "cmd.exe":
        subprocess.run('deactivate', shell=True)
    elif SHELL == 'bash' or SHELL == 'zsh':
        subprocess.run(f'{SHELL} -c "deactivate"', shell=True)

    # clean up env dir
    shutil.rmtree(os.path.join(path, 'env'))
