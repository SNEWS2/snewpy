from snewpy.snowglobes import generate_time_series as NEW_generate_time_series
from snewpy.snowglobes import generate_fluence as NEW_generate_fluence


def generate_time_series(model_path, model_file, model_type, transformation_type, transformation_parameters, d, output_filename, ntbins, deltat):
    print("[INFO] The `generate_time_series` function has been moved to the `snewpy.snowglobes` module.")
    return NEW_generate_time_series(model_path, model_file, model_type, transformation_type, transformation_parameters, d, output_filename, ntbins, deltat)


def generate_fluence(model_path, model_file, model_type, transformation_type, d, output_filename, tstart=None, tend=None):
    print("[INFO] The `generate_fluence` function has been moved to the `snewpy.snowglobes` module.")
    return NEW_generate_fluence(model_path, model_file, model_type, transformation_type, d, output_filename, tstart, tend)
