import os
from warnings import warn

from snewpy.snowglobes import collate as NEW_collate

def collate(Branch, Model_Path, Tarball, detector_input="all", skip_plots=False, return_tables=False, verbose=False, remove_generated_files=True):
    # Deprecated since SNEWPY v1.1
    warn(f"The 'snewpy.from_snowglobes.collate()' function is deprecated. Use 'snewpy.snowglobes.collate()' instead.", FutureWarning)
    tarball_path = os.path.join(Model_Path, Tarball)
    return NEW_collate(Branch, tarball_path, skip_plots=skip_plots)
