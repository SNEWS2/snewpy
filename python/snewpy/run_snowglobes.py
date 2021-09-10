import os

from snewpy.snowglobes import simulate

def go(SNOwGLoBESdir, Models_Path, Tarball, detector_input="all", verbose=False):
    print("[INFO] The `go` function has been moved to the `snewpy.snowglobes` module.")
    tarball_path = os.path.join(Models_Path, Tarball)
    return simulate(SNOwGLoBESdir, tarball_path, detector_input, verbose)
