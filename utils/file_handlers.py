import os
import shutil
import atexit
from typing import List


def create_temp_directory(directory_path: str = './temp/') -> None:
    """Create temporary directory for uploaded files."""
    try:
        os.makedirs(directory_path, exist_ok=True)
    except OSError as e:
        print(f"Error creating directory: {e}")


def save_uploaded_files(uploaded_files, directory_path: str = './temp/') -> List[str]:
    """Save uploaded files to temporary directory and return their paths."""
    file_paths = []
    for f in uploaded_files:
        path = f"{directory_path}{f.name}"
        with open(path, "wb") as out:
            out.write(f.read())
        file_paths.append(path)
    return file_paths


def cleanup_temp_directory(directory_path: str = './temp/') -> None:
    """Remove temporary directory and all its contents."""
    shutil.rmtree(directory_path, ignore_errors=True)


def register_cleanup_on_exit(directory_path: str = './temp/') -> None:
    """Register cleanup function to run when the application exits."""
    atexit.register(lambda: cleanup_temp_directory(directory_path))
