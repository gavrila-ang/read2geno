import os
import sys

def validate_paths(path_list):
    """
    Validates that all paths in the given list exist on the filesystem.
    
    Parameters:
    - path_list (list of str): A list of file or directory paths to validate.

    If any path does not exist, it prints an error message for each missing path
    and exits the script.
    """
    
    failed_paths = []

    # Check the existence of each path in the list
    for path in path_list:
        if not os.path.exists(path):
            failed_paths.append(path)
            print(f"ERROR: {path} does not exist.")

    # If any paths are missing, exit the script
    if failed_paths:
        print("Exiting pipeline due to missing paths.")
        sys.exit(1)
