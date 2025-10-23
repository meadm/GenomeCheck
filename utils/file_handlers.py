import os
import shutil
import atexit
from typing import List


def create_temp_directory(directory_path: str = None) -> str:
    """Create temporary directory for uploaded files.

    If directory_path is None, the function will use the environment variable
    SESSION_TEMP_DIR if present, otherwise './temp/'. The function returns the
    path that was created for convenience.
    """
    if directory_path is None:
        directory_path = os.environ.get('SESSION_TEMP_DIR', './temp/')
    try:
        os.makedirs(directory_path, exist_ok=True)
    except OSError as e:
        print(f"Error creating directory: {e}")
    return directory_path


def save_uploaded_files(uploaded_files, directory_path: str = None) -> List[str]:
    """Save uploaded files to temporary directory and return their paths."""
    file_paths = []
    if directory_path is None:
        directory_path = os.environ.get('SESSION_TEMP_DIR', './temp/')
    for f in uploaded_files:
        path = os.path.join(directory_path, f.name)
        with open(path, "wb") as out:
            out.write(f.read())
        file_paths.append(path)
    return file_paths


def cleanup_temp_directory(directory_path: str = None) -> None:
    """Remove temporary directory and all its contents.

    If directory_path is None, the function will use SESSION_TEMP_DIR or
    './temp/' as a fallback.
    """
    if directory_path is None:
        directory_path = os.environ.get('SESSION_TEMP_DIR', './temp/')
    shutil.rmtree(directory_path, ignore_errors=True)


def register_cleanup_on_exit(directory_path: str = None) -> None:
    """Register cleanup function to run when the application exits.

    If directory_path is None, the registered cleanup will use
    SESSION_TEMP_DIR or './temp/' at execution time.
    """
    atexit.register(lambda: cleanup_temp_directory(directory_path))
