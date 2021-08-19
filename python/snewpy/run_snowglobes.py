from snewpy.snowglobes import go as NEW_go

def go(SNOwGLoBESdir, Models_Path, Tarball, detector_input = all, verbose=False):
    print("[INFO] The `go` function has been moved to the `snewpy.snowglobes` module.")
    return NEW_go(SNOwGLoBESdir, Models_Path, Tarball, detector_input, verbose)
