import os
from warnings import warn

from snewpy.snowglobes import simulate

def go(SNOwGLoBESdir, Models_Path, Tarball, detector_input="all", verbose=False):
    # Deprecated since SNEWPY v1.1
    warn(f"The 'run_snowglobes.go()' function is deprecated. Use 'snowglobes.simulate()' instead.", FutureWarning)
    tarball_path = os.path.join(Models_Path, Tarball)
    return simulate(SNOwGLoBESdir, tarball_path, detector_input=detector_input)
