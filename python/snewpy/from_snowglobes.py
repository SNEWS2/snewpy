import os

from snewpy.snowglobes import collate as NEW_collate

def collate(Branch, Model_Path, Tarball, detector_input="all", skip_plots=False, return_tables=False, verbose=False, remove_generated_files=True):
    print("[INFO] The `collate` function has been moved to the `snewpy.snowglobes` module.")
    tarball_path = os.path.join(Model_Path, Tarball)
    return NEW_collate(Branch, tarball_path, detector_input, skip_plots, verbose, remove_generated_files)
